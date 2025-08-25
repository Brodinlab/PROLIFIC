library(tidyverse)
library(readxl)
library(ggpubr)
library(gghalves)
library(rstatix)
library(gridExtra)

df_deltas <- read_csv("data/clinical_outcomes_delta_25_bind_eq5d5l.csv") %>%
    distinct()

df_wide <- df_deltas %>%
    mutate(identifier = paste0(clinical_outcome, "_", comparison)) %>%
    select(Subject.Id, treatment_group, identifier, delta) %>%
    pivot_wider(names_from = identifier, values_from = delta) %>%
    filter(!is.na(EQ5D5L_eot))

col_set_1 <- c(
    "RHI_followup1",
    # "6MWT_followup1", "6MWT_followup2", "6MWT_of_predicted_followup1", "6MWT_of_predicted_followup2", # they are largely duplicated
    "6MWT_followup1", "6MWT_followup2",
    "StepsPerDay_followup1",
    "EQ5D5L_eot"
)
col_set_2 <- c(
    "CAT_eot", "CAT_followup1", "CAT_followup2",
    "Nijmegen_followup1", "Nijmegen_followup2",
    "FSS_followup1", "FSS_followup2",
    "Depaul_followup2"
)

df_score <- df_wide %>%
    mutate(
        across(all_of(col_set_1),
            ~ case_when(
                is.na(.) ~ NA_real_,
                . > 0 ~ 1, # increase is improvement
                TRUE ~ 0
            ),
            .names = "{.col}_res"
        )
    ) %>%
    mutate(
        across(all_of(col_set_2),
            ~ case_when(
                is.na(.) ~ NA_real_,
                . < 0 ~ 1, # decrease is improvement
                TRUE ~ 0
            ),
            .names = "{.col}_res"
        )
    ) %>%
    mutate(
        across(all_of("MOCA_followup2"),
            ~ case_when(
                is.na(.) ~ NA_real_,
                . > 1 ~ 1, # increase is improvement, set criteria to more than 1
                TRUE ~ 0
            ),
            .names = "{.col}_res"
        )
    ) %>%
    mutate(
        composite_score = rowMeans(select(., ends_with("_res")), na.rm = TRUE),
        non_na_count = rowSums(!is.na(select(., ends_with("_res"))))
    ) %>%
    mutate(score_group = if_else(composite_score >= 0.75, "high(>=0.75)", "low(<0.75)"))

set.seed(42)
df_score_qc <- df_score %>%
    filter(non_na_count > 6)

g_rainplot <- df_score_qc %>%
    ggplot(aes(x = treatment_group, y = composite_score, fill = treatment_group, color = treatment_group)) +
    geom_half_violin(side = "r", color = NA, nudge = 0.08) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fatten = NULL, fill = NA, color = "black") +
    stat_summary(fun = median, geom = "point", size = 3) +
    geom_half_point(side = "l", alpha = 0.3, range_scale = 0.5, size = 2) +
    theme_pubr() +
    theme(panel.background = element_blank()) +
    scale_fill_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
    scale_color_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
    stat_compare_means() +
    ggtitle(sprintf(
        "Pax: %d, Placebo: %d",
        sum(df_score_qc$treatment_group == "Paxlovid"),
        sum(df_score_qc$treatment_group == "Placebo")
    )) +
    theme(legend.position = "none") +
    geom_text(
        data = df_score_qc %>%
            group_by(treatment_group) %>%
            summarise(median = median(composite_score, na.rm = TRUE)),
        aes(y = median + 0.4, label = sprintf("%.2f", median)),
        vjust = -0.5
    )

tb <- df_score_qc %>%
    group_by(treatment_group, score_group) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = score_group, values_from = n) %>%
    column_to_rownames("treatment_group")

g_rainplot <- g_rainplot +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "red") +
    annotate("text",
        x = 0.52, y = 0.76,
        label = "score = 0.75",
        vjust = -0.5
    ) +
    annotate("text",
        x = 1.5, y = 0.6,
        label = sprintf("Fisher's exact test p = %.4f", fisher_test(tb)$p)
    ) +
    annotation_custom(
        tableGrob(tb, theme = ttheme_minimal())
    )

g_rainplot

ggsave("figures/Fig3/Fig3a.pdf", plot = g_rainplot, width = 8, height = 8)

write_csv(df_score_qc, "data/df_composite_score_0.75.csv")
