########## BiodiversiTREE ##########

### Analysis of the relationships between earthworm ### 
### and soil SI vales/C:N ratios across soil depth ###
### Author: Gen-Chang Hsu ###

### Libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(modelr)
library(broom.mixed)
library(readxl)
library(lme4)
library(car)
library(sjPlot)

# Data preparation ==============================================================================
### Load the data
Soil_data_by_plot <- readRDS("./Output/Data_clean/Soil_data_by_plot.rds")
Earthworm_data_by_plot <- readRDS("./Output/Data_clean/Earthworm_data_by_plot.rds")


### Convert the soil data into wide format
Soil_data_wide <- Soil_data_by_plot %>%
  pivot_wider(id_cols = Plot, 
              names_from = Depth, 
              values_from = c("Mean_d15N_soil", "Mean_d13C_soil", "Mean_C:N_soil"))


### Join the earthworm and soil data
All_data <- left_join(Earthworm_data_by_plot, 
                      Soil_data_wide,
                      by = "Plot")


# Fit lmer models ==============================================================================
### Fit the lmer models for each species
lmer_out <- All_data %>%
  group_by(Species) %>%
  mutate_at(.vars = vars(starts_with("Mean_")), .funs = scale) %>%
  nest() %>%
  mutate(Mod_fit = map(data, function(x) {
    mods <- list()
    mods[[1]] <-
      lmer(d13C_earthworm ~ `Mean_d13C_soil_0-2` + (1 |
                                                      Plot),
           data = x,
           REML = FALSE)
    mods[[2]] <-
      lmer(d13C_earthworm ~ `Mean_d13C_soil_2-5` + (1 |
                                                      Plot),
           data = x,
           REML = FALSE)
    mods[[3]] <-
      lmer(
        d13C_earthworm ~ `Mean_d13C_soil_5-10` + (1 |
                                                    Plot),
        data = x,
        REML = FALSE
      )
    mods[[4]] <-
      lmer(
        d13C_earthworm ~ `Mean_d13C_soil_10-20` + (1 |
                                                     Plot),
        data = x,
        REML = FALSE
      )
    mods[[5]] <-
      lmer(d13C_earthworm ~ `Mean_d15N_soil_0-2` + (1 |
                                                      Plot),
           data = x,
           REML = FALSE)
    mods[[6]] <-
      lmer(d13C_earthworm ~ `Mean_d15N_soil_2-5` + (1 |
                                                      Plot),
           data = x,
           REML = FALSE)
    mods[[7]] <-
      lmer(
        d13C_earthworm ~ `Mean_d15N_soil_5-10` + (1 |
                                                    Plot),
        data = x,
        REML = FALSE
      )
    mods[[8]] <-
      lmer(
        d13C_earthworm ~ `Mean_d15N_soil_10-20` + (1 |
                                                     Plot),
        data = x,
        REML = FALSE
      )
    mods[[9]] <-
      lmer(d13C_earthworm ~ `Mean_C:N_soil_0-2` + (1 |
                                                     Plot),
           data = x,
           REML = FALSE)
    mods[[10]] <-
      lmer(d13C_earthworm ~ `Mean_C:N_soil_2-5` + (1 |
                                                     Plot),
           data = x,
           REML = FALSE)
    mods[[11]] <-
      lmer(d13C_earthworm ~ `Mean_C:N_soil_5-10` + (1 |
                                                      Plot),
           data = x,
           REML = FALSE)
    mods[[12]] <-
      lmer(
        d13C_earthworm ~ `Mean_C:N_soil_10-20` + (1 |
                                                    Plot),
        data = x,
        REML = FALSE
      )
    return(mods)
  })) %>%
  mutate(Mod_summary = map(Mod_fit, function(X) {
    lapply(X, function(x) {
      coef <- coef(x)$Plot[1, 2]
      lower <- confint(x)[4, 1]
      upper <- confint(x)[4, 2]
      pval <- car::Anova(x)$`Pr(>Chisq)`
      data.frame(
        Beta = coef,
        Lower = lower,
        Upper = upper,
        P_val = pval
      )
    })
  }))

lmer_out <- lmer_out %>%
  mutate(Mod_summary_df = map(Mod_summary, function(x) {
    df <- bind_rows(x) %>%
      mutate(
        Signif = ifelse(P_val < 0.05, "*", ""),
        Var = factor(
          rep(
            c("Mean_d13C_soil", "Mean_d15N_soil", "Mean_C:N_soil"),
            each = 4
          ),
          levels = c("Mean_d13C_soil", "Mean_d15N_soil", "Mean_C:N_soil"),
          ordered = T
        ),
        Depth = factor(
          rep(c("0-2", "2-5", "5-10", "10-20"), 3),
          levels = c("0-2", "2-5", "5-10", "10-20"),
          ordered = T
        )
      )
  })) 


### Plot the standardized coefficients
lmer_out %>% select(Species, Mod_summary_df) %>%
  unnest(cols = Mod_summary_df) %>%
  ggplot(aes(x = fct_rev(Depth), y = Beta, color = Signif)) + 
  geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() + 
  facet_grid(Species ~ Var, 
             labeller = labeller(Species = as_labeller(c(`Allolobophora chlorotica` = "italic(A.~chlorotica)",
                                                         `Aporrectodea caliginosa` = "italic(A.~caliginosa)",
                                                         `Aporrectodea trapezoides` = "italic(A.~trapezoides)",
                                                         `Diplocardia sp.` = "italic(Diplocardia)~sp.",
                                                         `Lumbricus friendi` = "italic(L.~friendi)",
                                                         `Lumbricus rubellus` = "italic(L.~rubellus)"), label_parsed)  
                                 
                                 , Var = as_labeller(
               c(`Mean_d13C_soil` = "delta^{13}*C", 
                 `Mean_d15N_soil` = "delta^{15}*N",
                 `Mean_C:N_soil` = "C:N~ratio"), label_parsed))) + 
  labs(x = "Soil depth (cm)", y = expression(paste("Standardized coefficients (", beta, ")"))) + 
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_rect(fill = "grey90"),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 9),
        legend.position = c(0.5, 1.08),
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.15, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.size = unit(1, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 8),
        legend.text.align = 0,
        legend.box.just = "center",
        legend.justification = c(0.5, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", 
                                         size = 0.5, 
                                         linetype = "solid", 
                                         colour = "transparent")) + 
  scale_color_manual(values = c("black", "red"), guide = F)

ggsave("Output/Figures/Beta_coef.tiff", width = 8, height = 8, dpi = 600)

  
  
