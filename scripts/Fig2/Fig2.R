library(tidyverse)
library(ggpubr)
library(gghalves)

df_nanostring <- read_csv("data/nanostring_df_long.csv.gz") %>%
    mutate(time = factor(time, levels = c(1, 16, 45)))

# df_nanostring %>%
#     select(Subject.Id, treatment_group) %>%
#     distinct() %>%
#     group_by(treatment_group) %>%
#     tally()

viral_markers <- c("Envelope", "ORF3a", "Nucleocapsid", "ORF7a", "ORF1ab", "SARS-CoV-2_orf1ab_REV", "Spike", "FYN")

# Fig 2a
# Calculate log2 fold change between time 45 and time 1 for each viral marker
df_logfc <- df_nanostring %>%
    filter(gene %in% viral_markers, time %in% c(1, 45)) %>%
    pivot_wider(
        id_cols = c(Subject.Id, gene, treatment_group),
        names_from = time,
        values_from = value
    ) %>%
    mutate(log2fc = log2(`45` / `1`)) %>%
    filter(!is.infinite(log2fc), !is.na(log2fc))
set.seed(42)
df_logfc_sub <- df_logfc %>%
    filter(gene == "ORF1ab")
df_logfc_sub %>%
    ggplot(aes(x = treatment_group, y = log2fc, fill = treatment_group, color = treatment_group)) +
    geom_half_violin(side = "r", color = NA, nudge = 0.08) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fatten = NULL, fill = NA, color = "black") +
    stat_summary(fun = median, geom = "point", size = 3) +
    geom_half_point(side = "l", alpha = 0.3, range_scale = 0.5, size = 2) +
    theme_pubr() +
    theme(panel.background = element_blank()) +
    scale_fill_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
    scale_color_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
    stat_compare_means(aes(group = treatment_group), method = "wilcox.test", label = "p") +
    ggtitle(sprintf(
        "ORF1ab\nPax: %d, Placebo: %d",
        sum(df_logfc_sub$treatment_group == "Paxlovid"),
        sum(df_logfc_sub$treatment_group == "Placebo")
    )) +
    ylab("log2fc 45 vs 1") +
    theme(legend.position = "none")
ggsave("figures/Fig2/Fig2a.pdf", width = 8, height = 8)

# Fig 2b-c
all_pvalues <- lapply(viral_markers, function(viral_marker) {
    df_viral_marker <- df_nanostring %>% filter(gene == viral_marker)
    viral_marker_stratify <- df_viral_marker %>%
        filter(time == 1) %>% # baseline
        mutate(positivity = value > median(value)) %>%
        select(Subject.Id, positivity)
    df_viral_marker_positive <- df_viral_marker %>%
        left_join(viral_marker_stratify, by = "Subject.Id") %>%
        filter(!is.na(positivity), positivity == TRUE)

    # Get p-values for each treatment group and comparison
    df_viral_marker_positive %>%
        group_by(treatment_group) %>%
        summarise(
            p_16vs1 = wilcox.test(
                value[time == 16],
                value[time == 1],
                paired = TRUE
            )$p.value,
            p_45vs1 = wilcox.test(
                value[time == 45],
                value[time == 1],
                paired = TRUE
            )$p.value
        ) %>%
        pivot_longer(
            cols = starts_with("p_"),
            names_to = "comparison",
            values_to = "p.value"
        ) %>%
        mutate(viral_marker = viral_marker)
}) %>%
    bind_rows() %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))

lapply(viral_markers, function(viral_marker) {
    df_viral_marker <- df_nanostring %>% filter(gene == viral_marker)
    viral_marker_stratify <- df_viral_marker %>%
        filter(time == 1) %>% # baseline
        mutate(positivity = value > median(value)) %>%
        select(Subject.Id, positivity)
    df_viral_marker_positive <- df_viral_marker %>%
        left_join(viral_marker_stratify, by = "Subject.Id") %>%
        filter(!is.na(positivity), positivity == TRUE)

    # Get adjusted p-values for this marker
    marker_p_values <- all_pvalues %>%
        filter(viral_marker == !!viral_marker) %>%
        mutate(
            comparison = if_else(comparison == "p_16vs1", "16_1", "45_1"),
            label = sprintf("%.4f", p.adj),
            # Add y position for p-value labels - use max value plus some padding
            y.position = rep(c(max(df_viral_marker_positive$value) * 1.1, max(df_viral_marker_positive$value) * 1.2), 2),
            # Add required group columns for stat_pvalue_manual
            group1 = as.character(1),
            group2 = substr(comparison, 1, 2)
        )

    df_viral_marker_positive %>%
        ggplot(aes(x = time, y = value, color = treatment_group)) +
        geom_point() +
        geom_line(aes(group = Subject.Id), alpha = 0.5) +
        theme_pubr() +
        ylab(paste("Patients with", viral_marker, "above median")) +
        facet_wrap(~treatment_group) +
        stat_pvalue_manual(
            data = marker_p_values,
            label = "label",
        ) +
        scale_color_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
        ggtitle(sprintf(
            "%s\nPlacebo: %d, Paxlovid: %d",
            viral_marker,
            sum((df_viral_marker_positive %>% filter(time == 1))$treatment_group == "Placebo"),
            sum((df_viral_marker_positive %>% filter(time == 1))$treatment_group == "Paxlovid")
        ))
}) %>%
    ggarrange(plotlist = ., ncol = 3, nrow = 3) %>%
    ggexport(
        filename = "figures/Fig2/Fig2b-c.pdf",
        width = 20, height = 20
    )
