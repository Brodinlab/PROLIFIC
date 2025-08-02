library(tidyverse)
library(ggpubr)

df_nanostring <- read_csv("data/nanostring_df_long.csv.gz") %>%
    mutate(time = factor(time, levels = c(1, 16, 45)))

viral_markers <- c("Envelope", "ORF3a", "Nucleocapsid", "ORF7a", "ORF1ab", "SARS-CoV-2_orf1ab_REV", "Spike", "FYN")

# for all patients -------------------------------------------------------------
# First collect all p-values
all_pvalues_all_patients <- lapply(viral_markers, function(viral_marker) {
    df_viral_marker <- df_nanostring %>% filter(gene == viral_marker)
    # Get p-values for each treatment group and comparison
    df_viral_marker %>%
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

# Now create plots using adjusted p-values
lapply(viral_markers, function(viral_marker) {
    df_viral_marker <- df_nanostring %>% filter(gene == viral_marker)
    # Get adjusted p-values for this marker
    marker_p_values <- all_pvalues_all_patients %>%
        filter(viral_marker == !!viral_marker) %>%
        mutate(
            comparison = if_else(comparison == "p_16vs1", "16_1", "45_1"),
            label = sprintf("%.4f", p.adj),
            # Add y position for p-value labels - use max value plus some padding
            y.position = rep(c(max(df_viral_marker$value) * 1.1, max(df_viral_marker$value) * 1.2), 2),
            # Add required group columns for stat_pvalue_manual
            group1 = as.character(1),
            group2 = substr(comparison, 1, 2)
        )

    df_viral_marker %>%
        ggplot(aes(x = time, y = value, color = treatment_group)) +
        geom_point() +
        geom_line(aes(group = Subject.Id), alpha = 0.5) +
        theme_pubr() +
        ylab(paste("All patients -", viral_marker)) +
        facet_wrap(~treatment_group) +
        stat_pvalue_manual(
            data = marker_p_values,
            label = "label",
        ) +
        scale_color_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E"))
}) %>%
    ggarrange(plotlist = ., ncol = 3, nrow = 3) %>%
    ggexport(
        filename = "figures/Extended Data Fig/Extended Data Fig1.pdf",
        width = 20, height = 20
    )
