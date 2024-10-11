library(here)
library(ggplot2)
library(dplyr)
library(viridisLite)


#aic_df <- read.csv("../results/erp_aic/aic_diffs.csv")

#aic_df <- aic_df %>%
#    mutate(
#        study_labs = factor(study, levels = c("adsbc21", "dbc19", "adbc23"), labels = c("Aurnhammer et al. (2021)", "Delogu et al. (2019)", "Aurnhammer et al. (2023)")),
#        lme_labs = factor(lme, levels = c("condition", "leo13b", "gerpt2large", "gerpt2"), labels = c("Condition", "LeoLM", "GerPT-2 large", "GerPT-2"))
#    )

##############
#### N400 ####
##############

#n4_aic <- aic_df[aic_df$time_window == "N400",]

# Create the bar plot
#n4_p <- ggplot(n4_aic, aes(x = study_labs, y = norm_aic, fill = lme_labs)) +
#  geom_bar(stat = "identity", position = position_dodge()) +
#  scale_x_discrete(limits = c("Aurnhammer et al. (2021)", "Delogu et al. (2019)", "Aurnhammer et al. (2023)")) +   # Custom order for study
#  scale_fill_discrete(limits = c("Condition", "LeoLM", "GerPT-2 large", "GerPT-2")) +  # Custom order for lmes
  #scale_fill_brewer(palette="Set2")+
#  scale_fill_viridis_d()+
#  labs(x = "",y = "Regression AIC - Null AIC", fill = "LME") +
#  theme_minimal() +
#  ggtitle("N400") +
#  theme(plot.title = element_text(size = 16),
#  axis.text.x = element_text(size = 11))

#ggsave("../results/erp_aic/n400_aic_diffs.pdf", plot=n4_p, width = 8, height = 6)

##############
#### P600 ####
##############

#p6_aic <- aic_df[aic_df$time_window == "P600",]

# Create the bar plot
#p6_p <- ggplot(p6_aic, aes(x = study_labs, y = norm_aic, fill = lme_labs)) +
#  geom_bar(stat = "identity", position = position_dodge()) +
#  scale_x_discrete(limits = c("Aurnhammer et al. (2021)", "Delogu et al. (2019)", "Aurnhammer et al. (2023)")) +   # Custom order for study
#  scale_fill_discrete(limits = c("Condition", "LeoLM", "GerPT-2 large", "GerPT-2")) +  # Custom order for lmes
#  #scale_fill_brewer(palette="Set2")+
#  scale_fill_viridis_d()+
#  labs(x = "", y = "Regression AIC - Null AIC", fill = "LME") +
#  theme_minimal() +
#  ggtitle("P600") +
#  theme(plot.title = element_text(size = 16),
#  axis.text.x = element_text(size = 11))

#ggsave("../results/erp_aic/p600_aic_diffs.pdf", plot=p6_p, width = 8, height = 6)


aic_df <- read.csv("../results/erp_aic/aic_diffs.csv")

aic_df <- aic_df %>%
    mutate(
        study_labs = factor(study, levels = c("adsbc21", "dbc19", "adbc23"), labels = c("Aurnhammer et al. (2021)", "Delogu et al. (2019)", "Aurnhammer et al. (2023)")),
        lme_labs = factor(lme, levels = c("condition", "leo13b", "gerpt2large", "gerpt2"), labels = c("Condition", "LeoLM", "GerPT-2 large", "GerPT-2"))
    )

  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Regression AIC - Null AIC", fill = "LME") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ study_labs, scales = "free_x", space = "free_x") +  # Faceting by letter
  scale_fill_viridis_d()

ggsave("../results/erp_aic/test.pdf", plot=test, width = 8, height = 6)