########## BiodiversiTREE ##########

### Data Preparation ###
### Author: Gen-Chang Hsu ###

### Libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(modelr)
library(broom)
library(readxl)


# Data Wrangling ==============================================================================
### Load the raw datasets
Earthworm_data_raw <- read_xlsx("./data_raw/BiodiversiTREE_data.xlsx", sheet = 1)
Soil_data_raw <- read_xlsx("./data_raw/BiodiversiTREE_data.xlsx", sheet = 2)


### Summarize the soil data by plot
# Separate the two forest plots (Plot 39 and 40)
Soil_data_clean1 <- Soil_data_raw %>% 
  filter(!Plot %in% c(39, 40)) %>%
  select(Plot, Depth, ends_with(c("d15N", "d13C", "C:N"))) %>%
  group_by(Plot, Depth) %>%
  summarise_all(.funs = mean) %>%
  rename("d15N_soil" = "Linear corr d15N", 
         "d13C_soil" = "Linear corr d13C",
         "C:N_soil" = "C:N") %>%
  ungroup() %>%
  mutate(Depth = factor(Depth, levels = c("0-2", "2-5", "5-10", "10-20"), ordered = T)) %>%
  arrange(Plot, Depth)

Soil_data_clean2 <- Soil_data_raw %>% 
  filter(Plot %in% c(39, 40)) %>%
  select(Plot, Depth, ends_with(c("d15N", "d13C", "C:N"))) %>%
  group_by(Plot, Depth) %>%
  summarise_all(.funs = mean) %>%
  rename("d15N_soil" = "Linear corr d15N", 
         "d13C_soil" = "Linear corr d13C",
         "C:N_soil" = "C:N") %>%
  ungroup() %>%
  mutate(Depth = factor(Depth, levels = c("0-2", "2-5", "5-10", "10-20"), ordered = T)) %>%
  arrange(Plot, Depth)

write_rds(Soil_data_clean1, "./Output/Data_clean/Soil_data_clean1.rds")
write_rds(Soil_data_clean2, "./Output/Data_clean/Soil_data_clean2.rds")


### Clean up the earthworm data
# Separate the two forest plots (Plot 39 and 40)
Earthworm_data_clean1 <- Earthworm_data_raw %>% 
  filter(!Plot %in% c(39, 40)) %>%
  select(Plot, Species, ends_with(c("d15N", "d13C", "C:N"))) %>%
  rename("d15N_earthworm" = "Linear corr d15N", 
         "d13C_earthworm" = "Linear corr d13C",
         "C:N_earthworm" = "C:N") %>% 
  ungroup()

Earthworm_data_clean2 <- Earthworm_data_raw %>% 
  filter(Plot %in% c(39, 40)) %>%
  select(Plot, Species, ends_with(c("d15N", "d13C", "C:N"))) %>%
  rename("d15N_earthworm" = "Linear corr d15N", 
         "d13C_earthworm" = "Linear corr d13C",
         "C:N_earthworm" = "C:N") %>% 
  ungroup()

write_rds(Earthworm_data_clean1, "./Output/Data_clean/Earthworm_data_clean1.rds")
write_rds(Earthworm_data_clean2, "./Output/Data_clean/Earthworm_data_clean2.rds")


### Background‐corrected earthworm d13C and d15N
# 1. Use depth-weighted mean (0-2 and 2-5 cm) as plot-level SI values
# 2. Correct for the earthworm SI values by subtracting the difference between 
#    plot-level SI values and grand means
Earthworm_data_corrected1 <- Soil_data_clean1 %>% 
  filter(Depth %in% c("0-2", "2-5")) %>%
  group_by(Plot) %>%
  summarise_at(vars(ends_with("soil")), function(x){
    x[1]*0.4 + x[2]*0.6
  }) %>% 
  rename(Wmean_d15N_soil = d15N_soil, 
         Wmean_d13C_soil = d13C_soil, 
         `Wmean_C:N_soil` = `C:N_soil`) %>%
  mutate(Grand_d15N_soil = mean(Wmean_d15N_soil),
         Grand_d13C_soil = mean(Wmean_d13C_soil)) %>%
  right_join(Earthworm_data_clean1, by = "Plot") %>% 
  mutate(d15N_earthworm_correct = d15N_earthworm - (Wmean_d15N_soil - Grand_d15N_soil),
         d13C_earthworm_correct = d13C_earthworm - (Wmean_d13C_soil - Grand_d13C_soil))

Earthworm_data_corrected2 <- Soil_data_clean2 %>% 
  filter(Depth %in% c("0-2", "2-5")) %>%
  group_by(Plot) %>%
  summarise_at(vars(ends_with("soil")), function(x){
    x[1]*0.4 + x[2]*0.6
  }) %>% 
  rename(Wmean_d15N_soil = d15N_soil, 
         Wmean_d13C_soil = d13C_soil, 
         `Wmean_C:N_soil` = `C:N_soil`) %>%
  mutate(Grand_d15N_soil = mean(Wmean_d15N_soil),
         Grand_d13C_soil = mean(Wmean_d13C_soil)) %>%
  right_join(Earthworm_data_clean2, by = "Plot") %>% 
  mutate(d15N_earthworm_correct = d15N_earthworm - (Wmean_d15N_soil - Grand_d15N_soil),
         d13C_earthworm_correct = d13C_earthworm - (Wmean_d13C_soil - Grand_d13C_soil))

write_rds(Earthworm_data_corrected1, "./Output/Data_clean/Earthworm_data_corrected1.rds")
write_rds(Earthworm_data_corrected2, "./Output/Data_clean/Earthworm_data_corrected2.rds")
















