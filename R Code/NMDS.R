########## BiodiversiTREE ##########

### NMDS on the SI values and C:N ratios of earthworms ###
### Author: Gen-Chang Hsu ###

### Libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(modelr)
library(broom)
library(readxl)
library(vegan)
library(ggsci)


# NMDS on the SI values and C:N ratios of earthworms ==============================================================================
### Load the data
Soil_data_by_plot <- readRDS("./Output/Data_clean/Soil_data_by_plot.rds")
Earthworm_data_by_plot <- readRDS("./Output/Data_clean/Earthworm_data_by_plot.rds")


### Background‐corrected earthworm d13C and d15N
Earthworm_data_corrected <- Soil_data_by_plot %>% 
  group_by(Plot) %>%
  summarise_at(vars(starts_with("Mean")), mean, na.rm = T) %>%
  right_join(Earthworm_data_by_plot, by = "Plot") %>% 
  mutate(d15N_earthworm_correct = d15N_earthworm - Mean_d15N_soil,
         d13C_earthworm_correct = d13C_earthworm - Mean_d13C_soil)


### NMDS on SI values and C:N ratios
NMDS_earthworm <- Earthworm_data_corrected %>%
  select(d13C_earthworm_correct,          
         d15N_earthworm_correct,
         `C:N_earthworm`) %>%
  metaMDS(distance = "eu", autotransform = FALSE)     # Euclidean distance used

# Extract the site scores for plotting
NMDS_score <- scores(NMDS_earthworm, display = "sites") %>%
  as.data.frame() %>%
  mutate(Species = Earthworm_data_corrected$Species)

# NMDS plot
ggplot(NMDS_score, aes(x = NMDS1, y = NMDS2, color = Species)) + 
  geom_point(size = 1.5, alpha = 0.7, stroke = 1) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = rep(unit(0.05,"null"), 4),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.225, 0.825),
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.15, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.key.size = unit(1, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 7),
        legend.text.align = 0,
        legend.box.just = "center",
        legend.justification = c(0.5, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent")) +
  stat_ellipse(type = "t", level = 0.5, size = 0.8) +
  scale_color_npg(name = NULL, 
                  label = c(expression(italic("Allolobophora chlorotica")),
                            expression(italic("Aporrectodea caliginosa")),
                            expression(italic("Aporrectodea trapezoides")),
                            expression(paste(italic("Diplocardia"), " sp.")),
                            expression(italic("Lumbricus friendi")),
                            expression(italic("Lumbricus rubellus")))) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-5, 5, 2.5), labels = c("-5.0", "-2.5", "0.0", "2.5", "5.0"))

ggsave("Output/Figures/NMDS.tiff", width = 6, height = 5, dpi = 600)

  










