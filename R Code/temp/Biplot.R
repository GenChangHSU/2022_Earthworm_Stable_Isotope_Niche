## ------------------------------------------------------------------------
## Title: Stable Isotope Biplots 
##
## Author: Gen-Chang Hsu
##
## Date: 2021-02-24
##
## Description: Create earthworm stable isotope biplots for each dataset.
##  
##
## Notes:
##
##
## ------------------------------------------------------------------------


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(ggsci)


# Import files ------------------------------------------------------------
all_data_clean <- readRDS("./Output/Data_clean/all_data_clean.rds")


# ggplot2 theme -----------------------------------------------------------
mytheme <- theme(
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 15, margin = margin(t = 10)),
      axis.title.y = element_text(size = 15, margin = margin(r = 5)),
      plot.title = element_text(hjust = 0.5, size = 18),
      plot.margin = margin(0.2, 0.05, 0.05, 0.05, "null"),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(colour = "transparent"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
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
      legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent"))
  

# Code starts here ---------------------------------------------------------

### Background-adjusted earthworm stable isotope values
adjusted_data <- all_data_clean %>%
  group_by(Dataset) %>%
  
  # Compute the differences between plot-level soil SI values and 
  # site-level (dataset-level) grand mean values (i.e., background references)
  mutate(d13C_soil_grand_mean = mean(d13C_soil_0_5),      
         d15N_soil_grand_mean = mean(d15N_soil_0_5)) %>%
  
  # Adjust the earthworm SI values by shifting back the differences between 
  # plot-level soil SI values and grand mean values 
  mutate(d13C_worm_adjusted = d13C_worm - (d13C_soil_0_5 - d13C_soil_grand_mean),
         d15N_worm_adjusted = d15N_worm - (d15N_soil_0_5 - d15N_soil_grand_mean)) %>%
  ungroup() %>%
  
  # Reorder the levels in column "Dataset" for plotting purpose
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))


### Stable isotope biplots for BiodiversiTREE and BARC datasets
label_df1 <- data.frame(Dataset = c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC"),
                       x = c(-33.5, -33.5, -33.5), 
                       y = c(13.25, 13.25, 13.25),
                       Label = c("(a)", "(b)", "(c)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

adjusted_data %>% 
  filter(Dataset %in% c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC")) %>%
  ggplot(aes(x = d13C_worm_adjusted, y = d15N_worm_adjusted)) +
  geom_point(aes(color = Species), size = 1.2, alpha = 0.7, stroke = 1) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), 
       y = expression(paste(delta^{15}, "N (\u2030)", sep = ""))) +
  
  # 50% data ellipses using a multivariate t-distribution
  stat_ellipse(aes(color = Species), type = "t", level = 0.5, size = 0.6) +
  facet_wrap(~Dataset, scales = "free_y") +
  geom_text(data = label_df1, aes(x = x, y = y, label = Label), size = 5) +
  scale_color_npg(name = NULL, 
                  label = c(expression(italic("Allolobophora chlorotica")),
                            expression(italic("Aporrectodea caliginosa")),
                            expression(italic("Aporrectodea trapezoides")),
                            expression(paste(italic("Diplocardia"), " sp.")),
                            expression(italic("Lumbricus friendi")),
                            expression(italic("Lumbricus rubellus")))) +
  coord_cartesian(xlim = c(-32, -13), ylim = c(0, 12), clip = "off") +
  scale_y_continuous(breaks = seq(1, 11, 2), labels = c("1.0", "3.0", "5.0", "7.0", "9.0", "11.0")) + 
  scale_x_continuous(breaks = seq(-30, -15, 5), labels = c("-30.0", "-25.0", "-20.0", "-15.0")) +
  guides(color = guide_legend(nrow = 2, byrow = T)) +
  mytheme + 
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 1.2)) 
  
ggsave("Output/Figures/Biplot1.tiff", width = 11, height = 4.5, dpi = 600)


### Stable isotope biplots for SERC datasets
label_df2 <- data.frame(Dataset = c("SERC_2011", "SERC_2013"),
                       x = c(-33.5, -33.5), 
                       y = c(13.25, 13.25),
                       Label = c("(a)", "(b)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

adjusted_data %>% 
  filter(Dataset %in% c("SERC_2011", "SERC_2013")) %>%
  ggplot(aes(x = d13C_worm_adjusted, y = d15N_worm_adjusted)) +
  geom_point(aes(color = Species), size = 1.2, alpha = 0.7, stroke = 1) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), 
       y = expression(paste(delta^{15}, "N (\u2030)", sep = ""))) +
  
  # 50% data ellipses using a multivariate t-distribution
  stat_ellipse(aes(color = Species), type = "t", level = 0.5, size = 0.6) +
  facet_wrap(~Dataset, scales = "free_y") +
  geom_text(data = label_df2, aes(x = x, y = y, label = Label), size = 5) +
  scale_color_brewer(name = NULL, palette = "Set1",
                  label = c(expression(italic("Amynthas hilgendorfi")),
                            expression(italic("Aporrectodea caliginosa")),
                            expression(italic("Eisenoides lonnbergi")),
                            expression(italic("Lumbricus rubellus")),
                            expression(italic("Lumbricus terrestris")),
                            expression(italic("Octolasion cyaneum")))) +
  coord_cartesian(xlim = c(-35, -20), ylim = c(-3, 6), clip = "off") +
  scale_y_continuous(breaks = seq(-3, 6, 1.5), labels = c("-3.0", "-1.5", "0", "1.5", "3.0", "4.5", "6.0")) + 
  scale_x_continuous(breaks = seq(-35, -20, 3), labels = c("-35.0", "-32.0", "-29.0", "-26.0", "-23.0", "-20.0")) +
  guides(color = guide_legend(nrow = 2, byrow = T)) +
  mytheme + 
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 1.2)) 

ggsave("Output/Figures/Biplot2.tiff", width = 8, height = 4.5, dpi = 600)




















