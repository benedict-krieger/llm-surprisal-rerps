library(here)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(glue)


aic_df <- read.csv("../results/erp_aic/aic_diffs.csv")

aic_df <- aic_df %>%
  mutate(
    study_labs = factor(
                        study,
                        levels = c("adsbc21",
                                   "dbc19",
                                   "dbc19_corrected",
                                   "adbc23"),
                        labels = c("Aurnhammer et al. (2021)",
                                   "Delogu et al. (2019)",
                                   "Delogu et al. (2019) Re-estimated",
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

dbc19c_df <- aic_df %>% filter(study == "dbc19_corrected", time_window == "P600")

rest_df <- aic_df %>% filter(study != "dbc19_corrected") %>% filter(!(study == "dbc19" & time_window == "P600"))


#print(dbc19c_df)
#print(rest_df)


plot_aics <- function(df,title,width,height)
{

p <-  ggplot(df, aes(x = time_window, y = norm_aic, fill = lme_labs)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Regression AIC - Null AIC", fill = "LME") +
  theme_minimal() +
  facet_grid(. ~ study_labs, scales = "free_x", space = "free_x") +
  scale_fill_viridis_d() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    strip.text = element_text(size = 11)
  )

ggsave(glue("../results/erp_aic/lme_fit_{title}.pdf"), plot=p, width = width, height = height)

}

plot_aics(dbc19c_df,"dbc19c", width = 10, height = 3)
plot_aics(rest_df,"rest", width = 10, height = 3)