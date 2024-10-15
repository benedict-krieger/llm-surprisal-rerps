library(here)
library(ggplot2)
library(dplyr)
library(viridisLite)


aic_df <- read.csv("../results/erp_aic/aic_diffs.csv")

aic_df <- aic_df %>%
  mutate(
    study_labs = factor(
                        study,
                        levels = c("adsbc21", "dbc19", "adbc23"),
                        labels = c("Aurnhammer et al. (2021)",
                                   "Delogu et al. (2019)",
                                   "Aurnhammer et al. (2023)")),

    lme_labs = factor(
                      lme,
                      levels = c("condition",
                                 "leo13b",
                                 "gerpt2large",
                                 "gerpt2"),
                      labels = c("Condition",
                                 "LeoLM",
                                 "GerPT-2 large",
                                 "GerPT-2"))
  )

p <-  ggplot(aic_df, aes(x = time_window, y = norm_aic, fill = lme_labs)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Regression AIC - Null AIC", fill = "LME") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ study_labs, scales = "free_x", space = "free_x") +
  scale_fill_viridis_d()

ggsave("../results/erp_aic/lme_fit_tws.pdf", plot=p, width = 8, height = 6)