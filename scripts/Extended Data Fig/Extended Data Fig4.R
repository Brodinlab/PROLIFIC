library(tidyverse)
library(ggpubr)
library(rstatix)
library(caret)
library(gghalves)
library(glmnet)
source("scripts/function/RandomForest.R")

subset_ids <- read_csv("data/subset_ids.csv")
df_nanostring <- read_csv("data/nanostring_df_long.csv.gz")
df_olink <- read_csv("data/olink_df_qc_clean_2.csv.gz")

outcome <- "FSS" # EQ5D5L, FSS or StepsPerDay

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
    left_join(subset_ids %>% filter(source == outcome) %>% rename(clinical_outcome = source), by = "Subject.Id") %>%
    filter(!is.na(group)) %>%
    mutate(Assay = make.names(Assay))

# Prepare data for modeling
df_wide <- df %>%
    select(Subject.Id, Assay, value, group) %>%
    pivot_wider(names_from = Assay, values_from = value) %>%
    select(-Subject.Id) %>%
    drop_na()
df_wide$group <- factor(df_wide$group, levels = c("non_responder", "responder"))

# LASSO -------------------------------------------------------------------
# LASSO logistic regression for feature selection
# Prepare data for glmnet
x <- as.matrix(df_wide %>% select(-group))
y <- df_wide$group

x <- scale(x)
# Fit LASSO logistic regression with cross-validation
set.seed(42)
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc")

# Plot cross-validated error curve
pdf(sprintf("figures/Extended Data Fig/Extended Data Fig4/Extended Data Fig4a_%s.pdf", outcome))
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

# compare the features selectedby LASSO --------------------------------
# First calculate all p-values and adjust them
p_values <- lapply(lasso_features$feature, function(feature) {
    df_subset <- df %>% filter(Assay == feature)
    wilcox_test <- wilcox.test(value ~ group, data = df_subset)
    data.frame(
        feature = feature,
        p.value = wilcox_test$p.value
    )
}) %>%
    bind_rows() %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))

plots <- p_values %>%
    filter(p.adj < 0.05) %>%
    pull(feature) %>%
    lapply(function(feature) {
        df_subset <- df %>% filter(Assay == feature)
        # Get adjusted p-value for this feature
        p_adj <- p_values %>%
            filter(feature == !!feature) %>%
            pull(p.adj)
        ggplot(df_subset, aes(x = group, y = value, fill = group, color = group)) +
            geom_half_violin(side = "r", color = NA, nudge = 0.08) +
            geom_boxplot(width = 0.1, outlier.shape = NA, fatten = NULL, fill = NA, color = "black") +
            stat_summary(fun = median, geom = "point", size = 3) +
            geom_half_point(side = "l", alpha = 0.3, range_scale = 0.5, size = 2) +
            theme_pubr() +
            theme(panel.background = element_blank()) +
            scale_fill_manual(values = c("non_responder" = "#888888", "responder" = "#A1665E")) +
            scale_color_manual(values = c("non_responder" = "#888888", "responder" = "#A1665E")) +
            annotate("text",
                x = 1.5, y = max(df_subset$value),
                label = sprintf("FDR adjusted p = %s", format.pval(p_adj, digits = 3))
            ) +
            ggtitle(sprintf(
                "%s\nResponder: %d, Non-responder: %d",
                feature,
                sum(df_subset$group == "responder"),
                sum(df_subset$group == "non_responder")
            )) +
            ylab(feature) +
            theme(legend.position = "none")
    })

ggarrange(plotlist = plots, ncol = 4, nrow = 4) %>%
    ggexport(filename = sprintf("figures/Extended Data Fig/Extended Data Fig4/Extended Data Fig4_rainplot_%s.pdf", outcome), width = 20, height = 20)

lasso_features %>%
    mutate(direction = if_else(coef > 0, "high in responder", "high in non-responder")) %>%
    left_join(p_values, by = "feature") %>%
    mutate(significance = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01 ~ "**",
        p.adj < 0.05 ~ "*",
        TRUE ~ NA_character_
    )) %>%
    ggplot(aes(x = feature, y = coef)) +
    geom_segment(aes(xend = feature, y = 0, yend = coef, color = direction),
        alpha = ifelse(lasso_features$coef > 0, 1, 0.5),
        show.legend = FALSE
    ) +
    geom_point(aes(color = direction, fill = direction),
        size = 4,
        alpha = ifelse(lasso_features$coef > 0, 1, 0.5)
    ) +
    geom_text(aes(y = coef + sign(coef) * 0.1, label = significance), vjust = ifelse(lasso_features$coef > 0, -0.5, 1.5)) +
    scale_color_manual(values = c("high in responder" = "#A1665E", "high in non-responder" = "#A1665E")) +
    scale_fill_manual(values = c("high in responder" = "#A1665E", "high in non-responder" = "#A1665E")) +
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
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top"
    )
ggsave(sprintf("figures/Extended Data Fig/Extended Data Fig4/Extended Data Fig4_feature_selection_%s.pdf", outcome), width = 10, height = 6)


# Random forest ------------------------------------------------------------

group <- as.numeric(df_wide$group) - 1 # 0 for non_responder, 1 for responder
res.rf <- rf_fit(df_wide %>% select(-group), group, ntree = 100)

pdf(sprintf("figures/Extended Data Fig/Extended Data Fig4/Extended Data Fig4b_%s.pdf", outcome))
plot(res.rf$roc_curve, print.auc = TRUE)
dev.off()
