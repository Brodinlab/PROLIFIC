library(lme4)
library(tidyverse)
library(lmerTest) # Add lmerTest package for Type III F-tests
library(ggpubr)
library(ggVennDiagram)
library(rstatix)
library(gghalves)

# load the data
df_all <- read_csv("data/clinical_outcomes_25_bind_eq5d5l.csv") %>%
    mutate(time = factor(Activity.HOPE.CASE.ID, levels = c("baseline", "eot", "followup1", "followup2"))) %>%
    filter(!is.na(treatment_group), !is.na(time), !is.na(value))

df_deltas <- read_csv("data/clinical_outcomes_delta_25_bind_eq5d5l.csv")

get_subset_ids <- function(df_deltas, clinical_outcome, comparison, direction, quantile) {
    df_delta <- df_deltas %>%
        filter(
            clinical_outcome == !!clinical_outcome,
            comparison == !!comparison,
            !is.na(delta)
        )
    pax_ids <- df_delta %>%
        filter(treatment_group == "Paxlovid") %>%
        filter(if (direction == "increase") delta > quantile(delta, quantile) else delta < quantile(delta, quantile)) %>%
        pull(Subject.Id)
    placebo_ids <- df_delta %>%
        filter(treatment_group == "Placebo") %>%
        pull(Subject.Id)
    return(list(pax_ids = pax_ids, placebo_ids = placebo_ids))
}
get_results <- function(clinical_outcome, comparison, step, pos_trend) {
    # Filter data for specific clinical outcome
    df <- df_all %>% filter(clinical_outcome == !!clinical_outcome)
    # Loop through different percentiles
    results <- lapply(seq(0.2, 0.8, step), function(percentile) {
        # Subset the pax group
        pax_ids <- get_subset_ids(df_deltas, clinical_outcome, comparison, pos_trend, percentile)[["pax_ids"]]
        placebo_ids <- get_subset_ids(df_deltas, clinical_outcome, comparison, pos_trend, percentile)[["placebo_ids"]]
        df_subset <- df %>% filter(Subject.Id %in% c(pax_ids, placebo_ids))
        # Fit the model
        model <- lmer(
            value ~ treatment_group + time + (1 | Subject.Id),
            data = df_subset
        )
        # Get results
        anova_res <- anova(model, type = 3)
        coef_res <- summary(model)$coefficients["treatment_groupPlacebo", "Estimate"]
        # Return results
        data.frame(
            percentile = percentile,
            p_value = anova_res["treatment_group", "Pr(>F)"],
            coefficient = coef_res
        )
    }) %>%
        bind_rows() %>%
        mutate(adjusted_p = p.adjust(p_value, method = "BH"))
    return(results)
}

# find the cutoff for those of interest
find_cutoff <- function(res_list, outcome, direction) {
    sig_percentiles <- res_list[[outcome]] %>%
        filter(adjusted_p < 0.05) %>%
        pull(percentile)
    if (length(sig_percentiles) == 0) {
        return(NA)
    }
    if (direction == "increase") {
        return(min(sig_percentiles))
    } else if (direction == "decrease") {
        return(max(sig_percentiles))
    }
}

# linepolot and MEM
create_lineplot <- function(outcome, comparison, direction, percentile) {
    df <- df_all %>%
        filter(clinical_outcome == outcome)
    total_pax_n <- df_deltas %>%
        filter(
            clinical_outcome == !!outcome,
            comparison == !!comparison,
            !is.na(delta)
        ) %>%
        filter(treatment_group == "Paxlovid") %>%
        pull(Subject.Id) %>%
        length()
    # cutoff <- find_cutoff(res_list, outcome, direction)
    pax_ids <- get_subset_ids(df_deltas, outcome, comparison, direction, percentile)[["pax_ids"]]
    placebo_ids <- get_subset_ids(df_deltas, outcome, comparison, direction, percentile)[["placebo_ids"]]
    df_subset <- df %>%
        filter(Subject.Id %in% c(pax_ids, placebo_ids)) %>%
        mutate(n = n_distinct(Subject.Id), .by = treatment_group)
    df_subset %>%
        ggplot(aes(x = time, y = value, color = treatment_group, group = Subject.Id)) +
        geom_line(alpha = 0.2) +
        geom_smooth(aes(group = treatment_group), method = "loess", se = FALSE, linewidth = 2) +
        scale_color_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
        theme_pubr() +
        ggtitle(paste0(
            outcome, "\n",
            "Placebo n=", length(placebo_ids), ", ",
            "Paxlovid n=", length(pax_ids),
            " (out of ", total_pax_n, ")"
        ))
}

res_list <- list(
    get_results("StepsPerDay", "followup1", 0.01, "increase"),
    get_results("RHI", "followup1", 0.01, "increase"),
    get_results("6MWT", "followup1", 0.01, "increase"),
    get_results("6MWT_of_predicted", "followup1", 0.01, "increase"),
    get_results("MOCA", "followup2", 0.01, "increase"),
    get_results("EQ5D5L", "eot", 0.01, "increase"),
    get_results("CAT", "eot", 0.01, "decrease"),
    get_results("Nijmegen", "followup2", 0.01, "decrease"),
    get_results("Depaul", "followup2", 0.01, "decrease"),
    get_results("FSS", "followup2", 0.01, "decrease")
)
names(res_list) <- c("StepsPerDay", "RHI", "6MWT", "6MWT_of_predicted", "MOCA", "EQ5D5L", "CAT", "Nijmegen", "Depaul", "FSS")

# Fig 3a ------------------------------------------------------------------------------

outcomes <- c("EQ5D5L", "Nijmegen", "Depaul", "FSS", "StepsPerDay", "RHI", "6MWT", "CAT", "6MWT_of_predicted", "MOCA")
plots <- lapply(outcomes, function(i) {
    res_list[[i]] %>%
        ggplot(aes(x = percentile, y = adjusted_p)) +
        geom_point() +
        geom_line() +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
        theme_pubr() +
        ggtitle(i)
})

ggarrange(plotlist = plots, ncol = 5, nrow = 2) %>%
    ggexport(
        filename = "figures/Fig3/Fig3a.pdf",
        width = 25, height = 10
    )

# get the ids for significant outcomes and save ---------------------------------------
# Get the patient IDs for each outcome
steps_ids <- get_subset_ids(
    df_deltas, "StepsPerDay", "followup1", "increase",
    find_cutoff(res_list, "StepsPerDay", "increase")
)

fss_ids <- get_subset_ids(
    df_deltas, "FSS", "followup2", "decrease",
    find_cutoff(res_list, "FSS", "decrease")
)

eq5d_ids <- get_subset_ids(
    df_deltas, "EQ5D5L", "eot", "increase",
    find_cutoff(res_list, "EQ5D5L", "increase")
)

# Combine all IDs into a single dataframe
total_steps_ids <- df_deltas %>%
    filter(
        clinical_outcome == "StepsPerDay",
        comparison == "followup1",
        !is.na(delta)
    ) %>%
    filter(treatment_group == "Paxlovid") %>%
    pull(Subject.Id)

total_fss_ids <- df_deltas %>%
    filter(
        clinical_outcome == "FSS",
        comparison == "followup2",
        !is.na(delta)
    ) %>%
    filter(treatment_group == "Paxlovid") %>%
    pull(Subject.Id)

total_eq5d_ids <- df_deltas %>%
    filter(
        clinical_outcome == "EQ5D5L",
        comparison == "eot",
        !is.na(delta)
    ) %>%
    filter(treatment_group == "Paxlovid") %>%
    pull(Subject.Id)

ids_df <- bind_rows(
    data.frame(
        Subject.Id = steps_ids[["pax_ids"]],
        group = "responder",
        source = "StepsPerDay"
    ),
    data.frame(
        Subject.Id = setdiff(total_steps_ids, steps_ids[["pax_ids"]]),
        group = "non_responder",
        source = "StepsPerDay"
    ),
    data.frame(
        Subject.Id = fss_ids[["pax_ids"]],
        group = "responder",
        source = "FSS"
    ),
    data.frame(
        Subject.Id = setdiff(total_fss_ids, fss_ids[["pax_ids"]]),
        group = "non_responder",
        source = "FSS"
    ),
    data.frame(
        Subject.Id = eq5d_ids[["pax_ids"]],
        group = "responder",
        source = "EQ5D5L"
    ),
    data.frame(
        Subject.Id = setdiff(total_eq5d_ids, eq5d_ids[["pax_ids"]]),
        group = "non_responder",
        source = "EQ5D5L"
    )
)

write_csv(ids_df, "data/subset_ids.csv")

# Fig 3b ------------------------------------------------------------------------------
outcome <- "EQ5D5L"

# First calculate all p-values and correct them
df_outcome <- df_all %>%
    filter(
        clinical_outcome == !!outcome
    ) %>%
    left_join(
        ids_df %>% filter(source == !!outcome),
        by = "Subject.Id"
    ) %>%
    mutate(
        group = if_else(treatment_group == "Placebo", "Placebo", group),
        group = factor(group, levels = c("responder", "non_responder", "Placebo"))
    ) %>%
    filter(!is.na(group))

res_p <- df_outcome %>%
    group_by(Activity.HOPE.CASE.ID) %>%
    wilcox_test(value ~ group) %>%
    adjust_pvalue(method = "BH")

# Create plots with adjusted p-values
plots <- lapply(unique(res_p$Activity.HOPE.CASE.ID), function(time_point) {
    df <- df_outcome %>%
        filter(Activity.HOPE.CASE.ID == time_point)

    # Get adjusted p-values for this time point
    time_p_values <- res_p %>%
        filter(Activity.HOPE.CASE.ID == !!time_point) %>%
        add_xy_position(x = "group")

    ggplot(df, aes(x = group, y = value)) +
        geom_half_violin(aes(fill = group, color = group), side = "r", color = NA, nudge = 0.08) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fatten = NULL, fill = NA, color = "black") +
        stat_summary(fun = median, geom = "point", size = 3) +
        geom_half_point(aes(fill = group, color = group), side = "l", alpha = 0.3, range_scale = 0.5, size = 2) +
        theme_pubr() +
        theme(panel.background = element_blank()) +
        scale_fill_manual(values = c("non_responder" = "#888888", "responder" = "#CCB566", "Placebo" = "#888888")) +
        scale_color_manual(values = c("non_responder" = "#888888", "responder" = "#CCB566", "Placebo" = "#888888")) +
        stat_pvalue_manual(
            data = time_p_values,
            label = "p.adj",
            label.scientific = TRUE,
            label.round = 2,
        ) +
        theme(legend.position = "none") +
        ggtitle(sprintf(
            "%s %s\nResponder: %d, Non-responder: %d, Placebo: %d",
            time_point,
            outcome,
            sum(df$group == "responder"),
            sum(df$group == "non_responder"),
            sum(df$group == "Placebo")
        )) +
        theme(legend.position = "none")
})

ggarrange(plotlist = plots, ncol = 4, nrow = 1) %>%
    ggexport(filename = "figures/Fig3/Fig3b.pdf", width = 30, height = 10)

# Fig 3c ------------------------------------------------------------------------------
ggarrange(
    create_lineplot("StepsPerDay", "followup1", "increase", find_cutoff(res_list, "StepsPerDay", "increase")),
    create_lineplot("EQ5D5L", "eot", "increase", find_cutoff(res_list, "EQ5D5L", "increase")),
    create_lineplot("FSS", "followup2", "decrease", find_cutoff(res_list, "FSS", "decrease")),
    ncol = 2, nrow = 2
) %>%
    ggexport(
        filename = "figures/Fig3/Fig3c.pdf",
        width = 15, height = 15
    )

# Fig 3d ------------------------------------------------------------------------------
# Create list for Venn diagram
x <- list(
    StepsPerDay = steps_ids[["pax_ids"]],
    FSS = fss_ids[["pax_ids"]],
    EQ5D5L = eq5d_ids[["pax_ids"]]
)

# Create and export Venn diagram using ggVennDiagram
ggVennDiagram(x,
    category.names = names(x),
    label_alpha = 0,
    edge_size = 1
) +
    scale_fill_gradient(low = "#FFFFFF", high = "#A1665E") +
    theme(legend.position = "none") +
    scale_color_manual(values = "black")

ggsave(
    "figures/Fig3/Fig3d.pdf",
    width = 10,
    height = 10
)

