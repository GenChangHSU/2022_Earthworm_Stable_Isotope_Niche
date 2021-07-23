## -----------------------------------------------------------------------------
## Title: Visualization of Species' SEAb and Pairwise Percentage SEAb Overlaps
##
## Author: Gen-Chang Hsu
##
## Date: 2021-07-13
##
## Description: 
## 1. Stable isotope scatterplots along with the SEAbs of earthworm species and soil mean and se
## 2. Dot charts of SEAbs. 
## 3. Heatmaps of Pairwise Percentage SEAb overlaps.
## 4. Soil line range plot
## 5. Species big delta plot
## Notes:
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(ggsci)
library(ggpubr)
library(cowplot)


# Import files -----------------------------------------------------------------
all_data_clean <- readRDS("./Output/Data_clean/all_data_clean.rds")
SEAb_df <- readRDS("./Output/Data_clean/SEAb.rds")


# ggplot2 theme ----------------------------------------------------------------
mytheme <- theme(
  axis.text.x = element_text(size = 12, color = "black"),
  axis.text.y = element_text(size = 12, color = "black"),
  axis.title.x = element_text(size = 15, margin = margin(t = 10)),
  axis.title.y = element_text(size = 15, margin = margin(r = 5)),
  plot.title = element_text(hjust = 0.5, size = 25),
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


############################### Code starts here ###############################

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


### 1. Stable isotope scatterplots along with the SEAbs of earthworm species
# Create a fixed color palette for all species across the sites
Species_names1 <- adjusted_data %>% 
  filter(Dataset %in% c("BDTR1", "BDTR2", "BARC")) %>%
  .$Species %>%
  unique() %>% 
  sort()

Species_names2 <- adjusted_data %>% 
  filter(Dataset %in% c("SERC1", "SERC2")) %>%
  .$Species %>%
  unique() %>% 
  sort()

Pallete_colors <- c("#FFD966", 
                    "#843C0C", 
                    "#BF9001", 
                    "#3C5488FF", 
                    "#8FAADC", 
                    "#70AD47", 
                    "#C55A11", 
                    "#984EA3", 
                    "#385723", 
                    "#F4B183")

pal <- set_names(Pallete_colors, union(Species_names1, Species_names2))

# BiodiversiTREE and BARC datasets
SEAb_points_df_1 <- SEAb_df %>% 
  filter(Dataset %in% c("BDTR1", "BDTR2", "BARC")) %>%
  select(Dataset, Species, SEAb_points) %>%
  unnest(cols = SEAb_points) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

label_df1 <- data.frame(Dataset = c("BDTR1", "BDTR2", "BARC"),
                        x = c(-30, -30, -30), 
                        y = c(11.2, 11.2, 11.2),
                        Label = c("(a)", "(b)", "(c)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

adjusted_data %>% 
  filter(Dataset %in% c("BDTR1", "BDTR2", "BARC")) %>%
  ggplot(aes(x = d13C_worm_adjusted, y = d15N_worm_adjusted)) +
  geom_point(aes(color = Species), size = 1, alpha = 0.7, stroke = 0.5) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), 
       y = expression(paste(delta^{15}, "N (\u2030)", sep = ""))) +
  geom_path(data = SEAb_points_df_1, aes(x = x, y = y, color = Species), size = 0.5) +
  facet_wrap(~Dataset, scales = "free", nrow = 2) +
  geom_text(data = label_df1, aes(x = x, y = y, label = Label), size = 5) +
  scale_color_manual(name = NULL, 
                     values = pal,
                     limits = Species_names1,
                     label = c(expression(italic("Allolobophora chlorotica")),
                            expression(italic("Aporrectodea caliginosa")),
                            expression(italic("Aporrectodea trapezoides")),
                            expression(italic("Diplocardia caroliniana")),
                            expression(italic("Lumbricus friendi")),
                            expression(italic("Lumbricus rubellus")))) +
  coord_cartesian(xlim = c(-30, -14), ylim = c(0, 10), clip = "off") +
  scale_y_continuous(breaks = seq(0, 10, 2), labels = c("0", "2.0", "4.0", "6.0", "8.0", "10.0")) + 
  scale_x_continuous(breaks = seq(-30, -15, 5), labels = c("-30.0", "-25.0", "-20.0", "-15.0")) +
  mytheme + 
  theme(legend.direction = "vertical",
        legend.position = c(0.75, 0.225),
        legend.text = element_text(size = 10, margin = margin(t = 0)),
        legend.spacing.y = unit(0, "cm"),
        legend.key.height = unit(0.8, "cm"),
        plot.margin = margin(0.2, 0.05, 0.01, 0.05, "null")) 

ggsave("Output/Figures/SEAb_biplot1.tiff", device = tiff, width = 7, height = 7, dpi = 600)

# SERC dataset
SEAb_points_df_2 <- SEAb_df %>% 
  filter(Dataset %in% c("SERC1", "SERC2")) %>%
  select(Dataset, Species, SEAb_points) %>%
  unnest(cols = SEAb_points) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

label_df2 <- data.frame(Dataset = c("SERC1", "SERC2"),
                        x = c(-32, -32), 
                        y = c(7, 7),
                        Label = c("(a)", "(b)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

adjusted_data %>% 
  filter(Dataset %in% c("SERC1", "SERC2")) %>%
  ggplot(aes(x = d13C_worm_adjusted, y = d15N_worm_adjusted)) +
  geom_point(aes(color = Species), size = 1, alpha = 0.7, stroke = 0.5) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), 
       y = expression(paste(delta^{15}, "N (\u2030)", sep = ""))) +
  geom_path(data = SEAb_points_df_2, aes(x = x, y = y, color = Species), size = 0.5) +
  facet_wrap(~Dataset, scales = "free_y") +
  geom_text(data = label_df2, aes(x = x, y = y, label = Label), size = 5) +
  scale_color_manual(name = NULL, 
                     values = pal,
                     limits = Species_names2,
                     label = c(expression(italic("Aporrectodea caliginosa")),
                               expression(italic("Eisenoides lonnbergi")),
                               expression(italic("Lumbricus rubellus")),
                               expression(italic("Lumbricus terrestris")),
                               expression(italic("Metaphire hilgendorfi")),
                               expression(italic("Octolasion cyaneum")))) +
  coord_cartesian(xlim = c(-32, -22), ylim = c(-3, 6), clip = "off") +
  scale_y_continuous(breaks = seq(-4, 6, 2), labels = c("-4.0", "-2.0", "0", "2.0", "4.0", "6.0")) + 
  scale_x_continuous(breaks = seq(-31, -22, 3), labels = c("-31.0", "-28.0", "-25.0", "-22.0")) +
  guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(shape = 19))) + 
  mytheme + 
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 1.2),
        legend.text = element_text(size = 10, margin = margin(r = 3)),
        legend.spacing.y = unit(0, "cm"),
        plot.margin = margin(0.2, 0.05, 0.01, 0.05, "null")) 

ggsave("Output/Figures/SEAb_biplot2.tiff", device = tiff, width = 7, height = 4, dpi = 600)


### 2. Dot charts of SEAbs
# BiodiversiTREE and BARC datasets
label_df1 <- data.frame(Dataset = c("BDTR1", "BDTR2", "BARC"),
                        x = c(5.9, 5.9, 5.9), 
                        y = c(-0.5, -0.5, -0.5),
                        Label = c("(a)", "(b)", "(c)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

Sp_df1 <- adjusted_data %>% 
  filter(Dataset %in% c("BDTR1", "BDTR2", "BARC")) %>%
  .$Species %>% 
  unique() %>% 
  sort()

SEAb_df %>% 
  filter(Dataset %in% c("BDTR1", "BDTR2", "BARC")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T)) %>%
  mutate(Species = factor(Species, levels = Sp_df1, ordered = T)) %>%
  ggplot(aes(x = fct_rev(Species), y = SEAb_mean)) + 
  geom_pointrange(aes(ymin = `SEAb_hdr_lower`, ymax = `SEAb_hdr_upper`, color = Species)) +
  coord_flip(xlim = c(1, 5), ylim = c(0, 35), clip = "off") + 
  facet_wrap(~Dataset, scales = "free_x") +
  labs(x = NULL, y = expression(paste(SEA[b], " (", "\u2030"^2, ") ", "(Mean ± 95% HDI)"))) +
  geom_text(data = label_df1, aes(x = x, y = y, label = Label), size = 5) +
  scale_x_discrete(drop = F,
                   limit = rev(Sp_df1[Sp_df1 != "Diplocardia sp."]),
                   label = rev(c(expression(italic("Allolobophora chlorotica")),
                                 expression(italic("Aporrectodea caliginosa")),
                                 expression(italic("Aporrectodea trapezoides")),
                                 expression(italic("Lumbricus friendi")),
                                 expression(italic("Lumbricus rubellus"))))) + 
  scale_color_manual(name = NULL, guide = "none", drop = F, values = pal) +
  scale_y_continuous(breaks = seq(0, 35, 5), labels = c("0", "5", "10", "15", "20", "25", "30", "35")) + 
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 22),
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
  
ggsave("Output/Figures/SEAb_dotchart1.tiff", width = 10, height = 3.7, dpi = 600)

# SERC dataset
label_df2 <- data.frame(Dataset = c("SERC1", "SERC2"),
                        x = c(5.9, 5.9), 
                        y = c(-0.1, -0.1),
                        Label = c("(a)", "(b)")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))

Sp_df2 <- adjusted_data %>% 
  filter(Dataset %in% c("SERC1", "SERC2")) %>%
  .$Species %>% 
  unique() %>% 
  sort()

SEAb_df %>% 
  filter(Dataset %in% c("SERC1", "SERC2")) %>%
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T)) %>%
  mutate(Species = factor(Species, levels = Sp_df2, ordered = T)) %>%
  ggplot(aes(x = fct_rev(Species), y = SEAb_mean)) + 
  geom_pointrange(aes(ymin = `SEAb_hdr_lower`, ymax = `SEAb_hdr_upper`, color = Species)) +
  coord_flip(xlim = c(1, 5), ylim = c(0, 8), clip = "off") + 
  facet_wrap(~Dataset, scales = "free_x") +
  labs(x = NULL, y = expression(paste(SEA[b], " (", "\u2030"^2, ") ", "(Mean ± 95% HDI)"))) +
  geom_text(data = label_df2, aes(x = x, y = y, label = Label), size = 5) +
  scale_x_discrete(drop = F,
                   limit = rev(Sp_df2[Sp_df2 != "Lumbricus terrestris"]),
                   label = rev(c(expression(italic("Aporrectodea caliginosa")),
                                 expression(italic("Eisenoides lonnbergi")),
                                 expression(italic("Lumbricus rubellus")),
                                 expression(italic("Metaphire hilgendorfi")),
                                 expression(italic("Octolasion cyaneum"))))) + 
  scale_color_manual(name = NULL, guide = "none", drop = F, values = pal,
                     label = c(expression(italic("Aporrectodea caliginosa")),
                               expression(italic("Eisenoides lonnbergi")),
                               expression(italic("Lumbricus rubellus")),
                               expression(italic("Metaphire hilgendorfi")),
                               expression(italic("Octolasion cyaneum")))) +
  scale_y_continuous(breaks = seq(0, 8, 2), labels = c("0", "2", "4", "6", "8")) + 
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 22),
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

ggsave("Output/Figures/SEAb_dotchart2.tiff", width = 7.1, height = 3.5, dpi = 600)


### 3. Heatmaps of pairwise percentage SEAb overlaps
# Customized palette for the heatmaps
heatmap_pal <- c(rgb(246, 255, 0, 0.3*255, maxColorValue = 255), 
                 rgb(246, 255, 0, 0.9*255, maxColorValue = 255),
                 rgb(250, 200, 0, 1*255, maxColorValue = 255), 
                 rgb(255, 0, 0, 1*255, maxColorValue = 255),
                 rgb(255, 0, 0, 1*255, maxColorValue = 255))

# Extract a common legend for plotting
P_legend <-
  SEAb_df %>% filter(Dataset == "BDTR1") %>%
  select(1:3, SEAb_overlap) %>%
  unnest(cols = SEAb_overlap) %>%
  ggplot() +
  geom_tile(aes(x = Species, y = SEAb_overlap_with, fill = `Percent_overlap_%`)) + 
  scale_fill_gradientn(colors = heatmap_pal,
                       values = c(0, 0.3, 0.6, 0.9, 1.0),
                       name = "Percentage (%)",
                       limits = c(0, 90),
                       breaks = c(5, 25, 45, 65, 85)) +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(legend.direction = "horizontal",
        legend.title = element_text(size = 15, margin = margin(0, 0, 2, 0), hjust = 0.5),
        legend.text = element_text(size = 13, margin = margin(0, 0, 0, 3), hjust = 0.5),
        legend.key.width = unit(1.33, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.justification = c(0.5, 0.5),
        legend.background = element_blank())

legend <- get_legend(P_legend) %>% as_ggplot()

# Plot for BDTR1
P1 <-
  SEAb_df %>% filter(Dataset == "BDTR1") %>%
  select(1:3, SEAb_overlap) %>%
  unnest(cols = SEAb_overlap) %>%
  ggplot() +
  geom_tile(aes(x = Species, y = SEAb_overlap_with, fill = `Percent_overlap_%`), 
            color = "grey90",
            size = 0) + 
  scale_x_discrete(expand = c(0, 0),
                   label = c(expression(italic("Allolobophora \n   chlorotica")),
                             expression(italic("Aporrectodea \n   caliginosa")),
                             expression(italic("Aporrectodea \n  trapezoides")),
                             expression(italic("Lumbricus \n   friendi")),
                             expression(italic("Lumbricus \n  rubellus")))) +
  scale_y_discrete(expand = c(0, 0), 
                   label = c(expression(italic("Allolobophora \n   chlorotica")),
                                 expression(italic("Aporrectodea \n   caliginosa")),
                                 expression(italic("Aporrectodea \n  trapezoides")),
                                 expression(italic("Lumbricus \n   friendi")),
                                 expression(italic("Lumbricus \n  rubellus")))) +
  scale_fill_gradientn(colors = heatmap_pal,
                       values = c(0, 0.3, 0.6, 0.9, 1.0),
                       name = "Percentage (%)",
                       limits = c(0, 90),
                       breaks = c(5, 25, 45, 65, 85),
                       guide = "none") +
  labs(x = "Species A", 
       y = "Species B",
       title = "BDTR1",
       subtitle = "(a)") +
  coord_fixed(ratio = 1) +
  mytheme + 
  theme(axis.title.x = element_text(size = 15, margin = margin(t = -18)),
        axis.title.y.left = element_text(size = 15, margin = margin(r = 6)),
        axis.text.x = element_text(size = 11, color = "black", 
                                   vjust = 1, angle = 30,
                                   margin = margin(23, 0, 0, 0)),
        axis.text.y = element_text(size = 11, color = "black", 
                                   vjust = 0.5, hjust = 1,
                                   margin = margin(0, -3, 0, 0),
                                   angle = 30),
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "null"),
        legend.key.width = unit(0.7, "cm"),
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 18, hjust = 0.5, vjust = -4.5))

# Plot for BDTR2
P2 <-
  SEAb_df %>% filter(Dataset == "BDTR2") %>%
  select(1:3, SEAb_overlap) %>%
  unnest(cols = SEAb_overlap) %>%
  ggplot() +
  geom_tile(aes(x = Species, y = SEAb_overlap_with, fill = `Percent_overlap_%`), 
            color = "grey90",
            size = 0) + 
  scale_x_discrete(expand = c(0, 0),
                   label = c(expression(italic("Allolobophora \n   chlorotica")),
                             expression(italic("Aporrectodea \n   caliginosa")),
                             expression(italic("Aporrectodea \n  trapezoides")),
                             expression(italic("Lumbricus \n   friendi")),
                             expression(italic("Lumbricus \n  rubellus")))) +
  scale_y_discrete(expand = c(0, 0), 
                   label = c(expression(italic("Allolobophora \n   chlorotica")),
                             expression(italic("Aporrectodea \n   caliginosa")),
                             expression(italic("Aporrectodea \n  trapezoides")),
                             expression(italic("Lumbricus \n   friendi")),
                             expression(italic("Lumbricus \n  rubellus")))) +
  scale_fill_gradientn(colors = heatmap_pal,
                       values = c(0, 0.3, 0.6, 0.9, 1.0),
                       name = "Percentage (%)",
                       limits = c(0, 90),
                       breaks = c(5, 25, 45, 65, 85),
                       guide = "none") +
  labs(x = "Species A", 
       y = "Species B",
       title = "BDTR2",
       subtitle = "(b)") +
  coord_fixed(ratio = 1) +
  mytheme + 
  theme(axis.title.x = element_text(size = 15, margin = margin(t = -18)),
        axis.title.y.left = element_text(size = 15, margin = margin(r = 6)),
        axis.text.x = element_text(size = 11, color = "black", 
                                   vjust = 1, angle = 30,
                                   margin = margin(23, 0, 0, 0)),
        axis.text.y = element_text(size = 11, color = "black", 
                                   vjust = 0.5, hjust = 1,
                                   margin = margin(0, -3, 0, 0),
                                   angle = 30),
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "null"),
        legend.key.width = unit(0.7, "cm"),
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 18, hjust = 0.5, vjust = -4.5))

# Plot for BARC
P3 <-
  SEAb_df %>% filter(Dataset == "BARC") %>%
  select(1:3, SEAb_overlap) %>%
  unnest(cols = SEAb_overlap) %>%
  ggplot() +
  geom_tile(aes(x = Species, y = SEAb_overlap_with, fill = `Percent_overlap_%`), 
            color = "grey90",
            size = 0) + 
  scale_x_discrete(expand = c(0, 0),
                   label = c(expression(italic("Allolobophora \n   chlorotica")),
                             expression(italic("Aporrectodea \n   caliginosa")),
                             expression(italic("Aporrectodea \n  trapezoides")),
                             expression(italic("Lumbricus \n   friendi")),
                             expression(italic("Lumbricus \n  rubellus")))) +
  scale_y_discrete(expand = c(0, 0), 
                   label = c(expression(italic("Allolobophora \n   chlorotica")),
                             expression(italic("Aporrectodea \n   caliginosa")),
                             expression(italic("Aporrectodea \n  trapezoides")),
                             expression(italic("Lumbricus \n   friendi")),
                             expression(italic("Lumbricus \n  rubellus")))) +
  scale_fill_gradientn(colors = heatmap_pal,
                       values = c(0, 0.3, 0.6, 0.9, 1.0),
                       name = "Percentage (%)",
                       limits = c(0, 90),
                       breaks = c(5, 25, 45, 65, 85),
                       guide = "none") +
  labs(x = "Species A", 
       y = "Species B",
       title = "BARC",
       subtitle = "(c)") +
  coord_fixed(ratio = 1) +
  mytheme + 
  theme(axis.title.x = element_text(size = 15, margin = margin(t = -18)),
        axis.title.y.left = element_text(size = 15, margin = margin(r = 6)),
        axis.text.x = element_text(size = 11, color = "black", 
                                   vjust = 1, angle = 30,
                                   margin = margin(23, 0, 0, 0)),
        axis.text.y = element_text(size = 11, color = "black", 
                                   vjust = 0.5, hjust = 1,
                                   margin = margin(0, -3, 0, 0),
                                   angle = 30),
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "null"),
        legend.key.width = unit(0.7, "cm"),
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 15),
        plot.title = element_text(size = 18, hjust = 0.5, vjust = -4.5))

# Layout the plots with grey-filled diagonals
ggdraw() + 
  draw_plot(P1 + geom_tile(data = data.frame(x = 1:5, y = 1:5),
                           aes(x = x, y = y), 
                           fill = "grey70"), 
            x = -0.01, y = 0, width = 0.4, height = 1) + 
  draw_plot(P2 + geom_tile(data = data.frame(x = 1:5, y = 1:5),
                           aes(x = x, y = y), 
                           fill = "grey70"), 
            x = 0.38, y = 0, width = 0.4, height = 1) + 
  draw_plot(P3 + geom_tile(data = data.frame(x = 1:2, y = 1:2),
                           aes(x = x, y = y), 
                           fill = "grey70"), 
            x = 0.69, y = 0.42, width = 0.4, height = 0.52) + 
  draw_plot(legend, x = 0.716, y = 0.1675, width = 0.35, height = 0.25)

ggsave("Output/Figures/Overlap_grey.tiff", device = tiff, width = 13, height = 6, dpi = 600)

# Summary tables for the percentage SEAb overlaps
walk(c("BDTR1", "BDTR2"), function(x){
  SEAb_overlap <- SEAb_df %>% 
    filter(Dataset == x) %>%
    select(2, SEAb_overlap) %>%
    unnest(cols = SEAb_overlap) %>%
    mutate(`Percent_overlap_%` = round(`Percent_overlap_%`, 1)) %>%
    pivot_wider(id_cols = Species, 
                names_from = SEAb_overlap_with,
                values_from = `Percent_overlap_%`) %>%
    relocate(`Allolobophora chlorotica`, .after = Species)
  
  write_csv(SEAb_overlap, paste("Output/Data_clean/SEAb_overlap_", x, ".csv"))
})


walk(c("BARC"), function(x){
  SEAb_overlap <- SEAb_df %>% 
    filter(Dataset == x) %>%
    select(2, SEAb_overlap) %>%
    unnest(cols = SEAb_overlap) %>%
    mutate(`Percent_overlap_%` = round(`Percent_overlap_%`, 1)) %>%
    pivot_wider(id_cols = Species, 
                names_from = SEAb_overlap_with,
                values_from = `Percent_overlap_%`)
  
  write_csv(SEAb_overlap, paste("Output/Data_clean/SEAb_overlap_", x, ".csv"))
})







