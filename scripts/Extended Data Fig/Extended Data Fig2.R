library(lme4)
library(tidyverse)
library(lmerTest) # Add lmerTest package for Type III F-tests
library(ggpubr)

# load the data
df_all <- read_csv("data/clinical_outcomes_25_bind_eq5d5l.csv") %>%
    mutate(time = factor(Activity.HOPE.CASE.ID, levels = c("baseline", "eot", "followup1", "followup2"))) %>%
    filter(!is.na(treatment_group), !is.na(time), !is.na(value))
anova_results <- list()

for (outcome in unique(df_all$clinical_outcome)) {
    # Fit the model
    model <- lmer(
        value ~ treatment_group + time + (1 | Subject.Id),
        data = df_all %>% filter(clinical_outcome == outcome)
    )
    # Get Type III F-tests and store results
    anova_results[[outcome]] <- anova(model, type = 3)
}

# Extract p-values for both treatment group and time effects
p_values_treatment <- sapply(anova_results, function(x) x["treatment_group", "Pr(>F)"])
p_values_time <- sapply(anova_results, function(x) x["time", "Pr(>F)"])

# Adjust p-values using BH method
p_values_all <- c(p_values_treatment, p_values_time)
p_adj_all <- p.adjust(p_values_all, method = "BH")
p_adj_treatment <- p_adj_all[1:length(p_values_treatment)]
p_adj_time <- p_adj_all[(length(p_values_treatment) + 1):length(p_values_all)]

# Create plots with both p-values
plot_list <- lapply(unique(df_all$clinical_outcome), function(outcome) {
    df_all %>%
        filter(clinical_outcome == outcome) %>%
        ggplot(aes(x = time, y = value, color = treatment_group, group = Subject.Id)) +
        geom_line(alpha = 0.2) +
        geom_smooth(aes(group = treatment_group), method = "loess", se = FALSE, linewidth = 2) +
        scale_color_manual(values = c("Placebo" = "#888888", "Paxlovid" = "#A1665E")) +
        theme_pubr() +
        ggtitle(outcome) +
        annotate("text",
            x = 1, y = Inf,
            label = paste0(
                "Treatment adj. p = ", format.pval(p_adj_treatment[outcome], digits = 2, eps = 1e-27),
                "\nTime adj. p = ", format.pval(p_adj_time[outcome], digits = 2, eps = 1e-27)
            ),
            vjust = 1.2, hjust = 0
        )
})

ggarrange(plotlist = plot_list, ncol = 3, nrow = 4) %>%
    ggexport(filename = "figures/Extended Data Fig/Extended Data Fig2.pdf", width = 15, height = 20)
