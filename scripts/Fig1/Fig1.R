library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(gghalves)

data_impute <- read_csv("data/EQ5D5L_imputed.csv")


# figure 1c
ggplot(data_impute, aes(x = timepoint, y = Form.vasskala, group = Subject.Id, color = treatment_group)) +
    geom_point(alpha = 0.3, size = 3) +
    geom_line(alpha = 0.3) + 
    geom_smooth(aes(group = treatment_group), method = "loess", se = FALSE, linewidth = 2) + # Added trendline per group
    scale_fill_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
    scale_color_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
    scale_x_discrete(labels = c("Baseline", "Day 16", "Day 45", "Day 90")) +
    ylab("EQ-VAS") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/Fig1/Fig1c.pdf", width = 8, height = 8)


# figure 1d
data_delta <- data_impute %>%
    pivot_wider(
        names_from = timepoint,
        values_from = Form.vasskala
    ) %>%
    mutate(
        delta = EQ5D5L_endOfTreatment - EQ5D5L
    ) %>%
    select(Subject.Id, treatment_group, delta)

set.seed(42)
data_delta %>%
    ggplot(aes(x = treatment_group, y = delta, fill = treatment_group, color = treatment_group)) +
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
        "Pax: %d, Placebo: %d",
        sum(data_delta$treatment_group == "Paxlovid"),
        sum(data_delta$treatment_group == "Placebo")
    )) +
    theme(legend.position = "none") +
    geom_text(
        data = data_delta %>%
            group_by(treatment_group) %>%
            summarise(median = median(delta, na.rm = TRUE)),
        aes(y = median + 0.5, label = median),
        vjust = -0.5
    )

ggsave("figures/Fig1/Fig1d.pdf", width = 8, height = 8)
