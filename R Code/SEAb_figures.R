## ------------------------------------------------------------------------
## Title: Visualization of Species' SEAb and the Percent Overlap
##
## Author: Gen-Chang Hsu
##
## Date: 2021-03-05
##
## Description: 
## 1. Create biplots with SEAbs.
## 2. Visualize species' SEAbs using dot charts with line ranges. 
## 3. Visualize percent overlap in SEAb between species pairs using heatmaps.
##
## Notes:
##
##
## ------------------------------------------------------------------------
set.seed(123)


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(ggsci)


# Import files ------------------------------------------------------------
all_data_clean <- readRDS("./Output/Data_clean/all_data_clean.rds")
SEAb_df <- readRDS("./Output/Data_clean/SEAb.rds")


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


### Stable isotope biplots
# BiodiversiTREE and BARC datasets
SEAb_points_df_1 <- SEAb_df %>% 
  filter(Dataset %in% c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC")) %>%
  select(Dataset, Species, SEAb_points) %>%
  unnest(cols = SEAb_points) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

label_df1 <- data.frame(Dataset = c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC"),
                        x = c(-32.5, -32.5, -32.5), 
                        y = c(13.25, 13.25, 13.25),
                        Label = c("(a)", "(b)", "(c)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

adjusted_data %>% 
  filter(Dataset %in% c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC")) %>%
  ggplot(aes(x = d13C_worm_adjusted, y = d15N_worm_adjusted)) +
  geom_point(aes(color = Species), size = 1.2, alpha = 0.7, stroke = 1) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), 
       y = expression(paste(delta^{15}, "N (\u2030)", sep = ""))) +
  geom_path(data = SEAb_points_df_1, aes(x = x, y = y, color = Species), size = 1) +
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

ggsave("Output/Figures/SEAb_biplot1.tiff", width = 11, height = 4.5, dpi = 600)

# SERC dataset
SEAb_points_df_2 <- SEAb_df %>% 
  filter(Dataset %in% c("SERC_2011", "SERC_2013")) %>%
  select(Dataset, Species, SEAb_points) %>%
  unnest(cols = SEAb_points) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

label_df2 <- data.frame(Dataset = c("SERC_2011", "SERC_2013"),
                        x = c(-35.5, -35.5), 
                        y = c(7, 7),
                        Label = c("(a)", "(b)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

adjusted_data %>% 
  filter(Dataset %in% c("SERC_2011", "SERC_2013")) %>%
  ggplot(aes(x = d13C_worm_adjusted, y = d15N_worm_adjusted)) +
  geom_point(aes(color = Species), size = 1.2, alpha = 0.7, stroke = 1) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), 
       y = expression(paste(delta^{15}, "N (\u2030)", sep = ""))) +
  geom_path(data = SEAb_points_df_2, aes(x = x, y = y, color = Species), size = 1) +
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

ggsave("Output/Figures/SEAb_biplot2.tiff", width = 8, height = 4.5, dpi = 600)


### SEAb dot charts
# BiodiversiTREE and BARC datasets
label_df1 <- data.frame(Dataset = c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC"),
                        x = c(5.85, 5.85, 5.85), 
                        y = c(-1, -1, -1),
                        Label = c("(a)", "(b)", "(c)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

Sp_df1 <- adjusted_data %>% 
  filter(Dataset %in% c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC")) %>%
  .$Species %>% 
  unique() %>% 
  sort()

SEAb_df %>% 
  filter(Dataset %in% c("BiodiversiTREE_1", "BiodiversiTREE_2", "BARC")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T)) %>%
  mutate(Species = factor(Species, levels = Sp_df1, ordered = T)) %>%
  ggplot(aes(x = fct_rev(Species), y = SEAb_mean)) + 
  geom_pointrange(aes(ymin = `SEAb_2.5%`, ymax = `SEAb_97.5%`, color = Species)) +
  coord_flip(xlim = c(1, 5), ylim = c(0, 35), clip = "off") + 
  facet_wrap(~Dataset, scales = "free_x") +
  labs(x = NULL, y = expression(paste(SEA[b], " (", "\u2030"^2, ") ", "(Mean ± 95% CI)"))) +
  geom_text(data = label_df1, aes(x = x, y = y, label = Label), size = 5) +
  scale_x_discrete(drop = F,
                   limit = rev(Sp_df1[Sp_df1 != "Diplocardia sp."]),
                   label = rev(c(expression(italic("Allolobophora chlorotica")),
                                 expression(italic("Aporrectodea caliginosa")),
                                 expression(italic("Aporrectodea trapezoides")),
                                 expression(italic("Lumbricus friendi")),
                                 expression(italic("Lumbricus rubellus"))))) + 
  scale_color_npg(name = NULL, guide = F, drop = F) +
  scale_y_continuous(breaks = seq(0, 35, 5), labels = c("0", "5.0", "10.0", "15.0", "20.0", "25.0", "30.0", "35.0")) + 
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = rep(unit(0.05,"null"), 4),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.y = element_blank())
  
ggsave("Output/Figures/SEAb_dotchart1.tiff", width = 11, height = 4, dpi = 600)

# SERC dataset
label_df2 <- data.frame(Dataset = c("SERC_2011", "SERC_2013"),
                        x = c(5.85, 5.85), 
                        y = c(-0.3, -0.3),
                        Label = c("(a)", "(b)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

Sp_df2 <- adjusted_data %>% 
  filter(Dataset %in% c("SERC_2011", "SERC_2013")) %>%
  .$Species %>% 
  unique() %>% 
  sort()

SEAb_df %>% 
  filter(Dataset %in% c("SERC_2011", "SERC_2013")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T)) %>%
  mutate(Species = factor(Species, levels = Sp_df2, ordered = T)) %>%
  ggplot(aes(x = fct_rev(Species), y = SEAb_mean)) + 
  geom_pointrange(aes(ymin = `SEAb_2.5%`, ymax = `SEAb_97.5%`, color = Species)) +
  coord_flip(xlim = c(1, 5), ylim = c(0, 10), clip = "off") + 
  facet_wrap(~Dataset, scales = "free_x") +
  labs(x = NULL, y = expression(paste(SEA[b], " (", "\u2030"^2, ") ", "(Mean ± 95% CI)"))) +
  geom_text(data = label_df2, aes(x = x, y = y, label = Label), size = 5) +
  scale_x_discrete(drop = F,
                   limit = rev(Sp_df2[Sp_df2 != "Lumbricus terrestris"]),
                   label = rev(c(expression(italic("Amynthas hilgendorfi")),
                                 expression(italic("Aporrectodea caliginosa")),
                                 expression(italic("Eisenoides lonnbergi")),
                                 expression(italic("Lumbricus rubellus")),
                                 expression(italic("Octolasion cyaneum"))))) + 
  scale_color_brewer(name = NULL, guide = F, drop = F, palette = "Set1",
                     label = c(expression(italic("Amynthas hilgendorfi")),
                               expression(italic("Aporrectodea caliginosa")),
                               expression(italic("Eisenoides lonnbergi")),
                               expression(italic("Lumbricus rubellus")),
                               expression(italic("Octolasion cyaneum")))) +
  scale_y_continuous(breaks = seq(0, 10, 2.5), labels = c("0", "2.5", "5.0", "7.5", "10.0")) + 
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = rep(unit(0.05,"null"), 4),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.y = element_blank())

ggsave("Output/Figures/SEAb_dotchart2.tiff", width = 8, height = 4, dpi = 600)


### SEAb percent overlap heatmap
SEAb_df %>% filter(Dataset == "BiodiversiTREE_1") %>%
  select(1:3, SEAb_overlap) %>%
  unnest(cols = SEAb_overlap) %>%
  ggplot() +
  geom_tile(aes(x = Species, y = fct_rev(SEAb_overlap_with), fill = `Percent_overlap_%`)) + 
  scale_x_discrete(position = "top", expand = c(0, 0),
                   label = c(expression(italic("Allolobophora \n   chlorotica")),
                                 expression(italic("Aporrectodea \n   caliginosa")),
                                 expression(italic("Aporrectodea \n  trapezoides")),
                                 expression(italic("Lumbricus \n   friendi")),
                                 expression(italic("Lumbricus \n  rubellus")))) +
  scale_y_discrete(expand = c(0, 0), 
                   label = rev(c(expression(italic("Allolobophora \n   chlorotica")),
                             expression(italic("Aporrectodea \n   caliginosa")),
                             expression(italic("Aporrectodea \n  trapezoides")),
                             expression(italic("Lumbricus \n   friendi")),
                             expression(italic("Lumbricus \n  rubellus"))))) +
  scale_fill_continuous(trans = 'reverse', name = "Percent (%)") + 
  guides(fill = guide_colorbar(reverse = T)) + 
  labs(x = expression(paste(SEA[b], " (", "\u2030"^2, ") of")), 
       y = "Overlapping with") +
  mytheme + 
  theme(axis.title.x.top = element_text(size = 15, margin = margin(b = 15)),
        axis.title.y.left = element_text(size = 17, margin = margin(r = 8)),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black", vjust = 1),
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "null"),
        legend.key.width = unit(0.7, "cm"),
        legend.text = element_text(size = 8))

ggsave("Output/Figures/Overlap.tiff", width = 8, height = 5.5, dpi = 600)
















