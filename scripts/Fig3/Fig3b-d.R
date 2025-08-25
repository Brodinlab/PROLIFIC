library(tidyverse)
library(ggpubr)
library(gghalves)
library(rstatix)
library(glmnet)

df_score_qc <- read_csv("data/df_composite_score_0.75.csv")
df_nanostring <- read_csv("data/nanostring_df_long.csv.gz")
df_olink <- read_csv("data/olink_df_qc_clean_2.csv.gz")

df <- bind_rows(
    df_nanostring %>%
        filter(time == 1, treatment_group == "Paxlovid") %>%
        mutate(Assay = paste0("nanostring_", gene)) %>%
        select(Subject.Id, Assay, value),
    df_olink %>%
        filter(Visit == "Baseline", treatment_group == "Paxlovid") %>%
        mutate(Assay = paste0("olink_", Assay)) %>%
        select(Assay, NPX, Subject.Id) %>%
        rename(value = NPX)
) %>%
    mutate(Assay = make.names(Assay)) %>%
    left_join(df_score_qc %>% select(Subject.Id, composite_score, score_group), by = "Subject.Id")

df_wide <- df %>%
    select(Subject.Id, Assay, value, composite_score, score_group) %>%
    pivot_wider(names_from = Assay, values_from = value) %>%
    select(-Subject.Id) %>%
    drop_na()

# LASSO -------------------------------------------------------------------

x <- as.matrix(df_wide %>% select(-score_group, -composite_score))
y <- factor(df_wide$score_group, levels = c("low(<0.75)", "high(>=0.75)"))

x <- scale(x)
# Fit LASSO logistic regression with cross-validation
set.seed(42)
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc")

# Fig 4b -------------------------------------------------------------------
pdf(sprintf("figures/Fig3/Fig3b.pdf"))
plot(cvfit)
dev.off()

# Extract coefficients at best lambda
coef_lasso <- coef(cvfit, s = "lambda.min")
selected_features <- rownames(coef_lasso)[which(coef_lasso[, 1] != 0)]
selected_features <- setdiff(selected_features, "(Intercept)")

lasso_features <- data.frame(
    feature = selected_features,
    coef = coef_lasso[selected_features, 1]
) %>%
    mutate(abs_coef = abs(coef)) %>%
    arrange(desc(abs_coef)) %>%
    mutate(feature = fct_reorder(feature, abs_coef))

# Fig 4c -------------------------------------------------------------------
p_values <- lapply(lasso_features$feature, function(feature) {
    df_subset <- df %>% filter(Assay == feature)
    wilcox_test <- wilcox.test(value ~ score_group, data = df_subset)
    data.frame(
        feature = feature,
        p.value = wilcox_test$p.value
    )
}) %>%
    bind_rows() %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))
lasso_features %>%
    mutate(direction = if_else(coef > 0, "high in high(>=0.75)", "high in low(<0.75)")) %>%
    left_join(p_values, by = "feature") %>%
    mutate(significance = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01 ~ "**",
        p.adj < 0.05 ~ "*",
        TRUE ~ NA_character_
    )) %>%
    ggplot(aes(x = feature, y = coef)) +
    geom_segment(aes(xend = feature, y = 0, yend = coef, color = direction, alpha = direction),
        show.legend = FALSE
    ) +
    geom_point(aes(color = direction, fill = direction, alpha = direction),
        size = 4,
    ) +
    geom_text(aes(y = coef + sign(coef) * 0.1, label = significance), vjust = ifelse(lasso_features$coef > 0, -0.5, 1.5)) +
    scale_color_manual(values = c("high in high(>=0.75)" = "#A1665E", "high in low(<0.75)" = "#A1665E")) +
    scale_fill_manual(values = c("high in high(>=0.75)" = "#A1665E", "high in low(<0.75)" = "#A1665E")) +
    scale_alpha_manual(values = c("high in high(>=0.75)" = 1, "high in low(<0.75)" = 0.5)) +
    # coord_flip() +
    labs(x = "Feature", y = "LASSO coefficient") +
    theme_light() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_line(linetype = "dashed"),
        panel.grid.major = element_line(linetype = "dashed"),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 10),
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            # color = ifelse(rev(lasso_features$feature %in% overlapping_features), "red", "black")
        ),
        legend.position = "top"
    )
ggsave("figures/Fig3/Fig3c.pdf", width = 10, height = 6)

# Fig 4d -------------------------------------------------------------------

plots <- p_values %>%
    filter(p.adj < 0.05) %>%
    pull(feature) %>%
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
    ggexport(filename = "figures/Fig3/Fig3d.pdf", width = 20, height = 20)
