library(tidyverse)
library(ggpubr)
library(gghalves)
library(rstatix)

df_score_qc <- read_csv("data/df_composite_score_0.75.csv")
df_nanostring <- read_csv("data/nanostring_df_long.csv.gz")

df <- df_nanostring %>%
    filter(time == 1, treatment_group == "Paxlovid") %>%
    mutate(Assay = paste0("nanostring_", gene)) %>%
    select(Subject.Id, Assay, value) %>%
    mutate(Assay = make.names(Assay)) %>%
    left_join(df_score_qc %>% select(Subject.Id, composite_score, score_group), by = "Subject.Id")

viral_marker <- paste0("nanostring_", c("Envelope", "ORF3a", "Nucleocapsid", "ORF7a", "ORF1ab", "SARS.CoV.2_orf1ab_REV", "Spike", "FYN"))

p_values <- lapply(viral_marker, function(feature) {
    df_subset <- df %>% filter(Assay == feature)
    wilcox_test <- wilcox.test(value ~ score_group, data = df_subset)
    data.frame(
        feature = feature,
        p.value = wilcox_test$p.value
    )
}) %>%
    bind_rows() %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))

plots <- viral_marker %>%
    lapply(function(feature) {
        df_subset <- df %>% filter(Assay == feature, !is.na(score_group))
        # Get adjusted p-value for this feature
        p_adj <- p_values %>%
            filter(feature == !!feature) %>%
            pull(p.adj)
        ggplot(df_subset, aes(x = score_group, y = value, fill = score_group, color = score_group, alpha = score_group)) +
            geom_half_violin(side = "r", color = NA, nudge = 0.08) +
            geom_boxplot(width = 0.1, outlier.shape = NA, fatten = NULL, fill = NA, color = "black") +
            stat_summary(fun = median, geom = "point", size = 3) +
            geom_half_point(side = "l", alpha = 0.3, range_scale = 0.5, size = 2) +
            theme_pubr() +
            theme(panel.background = element_blank()) +
            scale_fill_manual(values = c("low(<0.75)" = "#A1665E", "high(>=0.75)" = "#A1665E")) +
            scale_color_manual(values = c("low(<0.75)" = "#A1665E", "high(>=0.75)" = "#A1665E")) +
            scale_alpha_manual(values = c("low(<0.75)" = 0.5, "high(>=0.75)" = 1)) +
            annotate("text",
                x = 1.5, y = max(df_subset$value),
                label = sprintf("FDR adjusted p = %s", format.pval(p_adj, digits = 3))
            ) +
            ggtitle(sprintf(
                "%s\nHigh: %d, Low: %d",
                feature,
                sum(df_subset$score_group == "high(>=0.75)"),
                sum(df_subset$score_group == "low(<0.75)")
            )) +
            ylab(feature) +
            theme(legend.position = "none")
    })

ggarrange(plotlist = plots, ncol = 4, nrow = 4) %>%
    ggexport(filename = "figures/Extended Data Fig/Extended Data Fig3.pdf", width = 20, height = 20)
