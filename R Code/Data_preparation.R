## ------------------------------------------------------------------------
## Title: Data Preparation
##
## Author: Gen-Chang Hsu
##
## Date: 2021-02-24
##
## Description: Clean the earthworm and soil stable isotope datasets 
## and organize them into a single tidy file
##
## Notes:
##
##
## ------------------------------------------------------------------------


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(readxl)


# Import files ------------------------------------------------------------
BiodiversiTREE_worm_raw <- read_xlsx("./data_raw/BiodiversiTREE_data.xlsx", sheet = 1)
BiodiversiTREE_soil_raw <- read_xlsx("./data_raw/BiodiversiTREE_data.xlsx", sheet = 2)
SERC_raw <- read_xlsx("./data_raw/SERC_data.xlsx", sheet = 1)
BARC_raw <- read_xlsx("./data_raw/BARC_data.xlsx", sheet = 1)


# Code start here ---------------------------------------------------------

### BiodiversiTREE_1 (excluding Plot 39 and 40)
BiodiversiTREE_soil_clean1 <- BiodiversiTREE_soil_raw %>% 
  filter(!Plot %in% c(39, 40)) %>%
  select(Plot, Depth, ends_with(c("d13C", "d15N"))) %>%
  group_by(Plot, Depth) %>%
  summarise_all(.funs = mean) %>%
  rename("d13C_soil" = "Linear corr d13C",
         "d15N_soil" = "Linear corr d15N") %>%
  ungroup() %>%
  mutate(Depth = factor(Depth, levels = c("0-2", "2-5", "5-10", "10-20"), ordered = T)) %>%
  arrange(Plot, Depth)

BiodiversiTREE_worm_clean1 <- BiodiversiTREE_worm_raw %>% 
  filter(!Plot %in% c(39, 40)) %>%
  select(Plot, Species, ends_with(c("d13C", "d15N"))) %>%
  rename("d13C_worm" = "Linear corr d13C",
         "d15N_worm" = "Linear corr d15N") %>% 
  ungroup()

BiodiversiTREE_clean1 <- BiodiversiTREE_soil_clean1 %>% 
  filter(Depth %in% c("0-2", "2-5")) %>%
  group_by(Plot) %>%
  summarise_at(vars(ends_with("soil")), function(x){
    x[1]*0.4 + x[2]*0.6  # weighted mean by depth
  }) %>%
  right_join(BiodiversiTREE_worm_clean1, by = "Plot") %>%
  mutate(Dataset = "BiodiversiTREE_1") %>%
  mutate(Plot = as.character(Plot)) %>%
  select(Dataset, Plot, Species, d13C_worm, d15N_worm, d13C_soil, d15N_soil)


### BiodiversiTREE_2 (Plot 39 and 40)
BiodiversiTREE_soil_clean2 <- BiodiversiTREE_soil_raw %>% 
  filter(Plot %in% c(39, 40)) %>%
  select(Plot, Depth, ends_with(c("d13C", "d15N"))) %>%
  group_by(Plot, Depth) %>%
  summarise_all(.funs = mean) %>%
  rename("d13C_soil" = "Linear corr d13C",
         "d15N_soil" = "Linear corr d15N") %>%
  ungroup() %>%
  mutate(Depth = factor(Depth, levels = c("0-2", "2-5", "5-10", "10-20"), ordered = T)) %>%
  arrange(Plot, Depth)

BiodiversiTREE_worm_clean2 <- BiodiversiTREE_worm_raw %>% 
  filter(Plot %in% c(39, 40)) %>%
  select(Plot, Species, ends_with(c("d13C", "d15N"))) %>%
  rename("d13C_worm" = "Linear corr d13C",
         "d15N_worm" = "Linear corr d15N") %>% 
  ungroup()

BiodiversiTREE_clean2 <- BiodiversiTREE_soil_clean2 %>% 
  filter(Depth %in% c("0-2", "2-5")) %>%
  group_by(Plot) %>%
  summarise_at(vars(ends_with("soil")), function(x){
    x[1]*0.4 + x[2]*0.6  # weighted mean by depth
  }) %>%
  right_join(BiodiversiTREE_worm_clean2, by = "Plot") %>%
  mutate(Dataset = "BiodiversiTREE_2") %>%
  mutate(Plot = as.character(Plot)) %>%
  select(Dataset, Plot, Species, d13C_worm, d15N_worm, d13C_soil, d15N_soil)


### BARC
BARC_worm_clean <- BARC_raw %>% 
  filter(`sample type` == "earthworm") %>% 
  select(Plot = quadret, 
         Species = species, 
         d13C_worm = d13C, 
         d15N_worm = d15N)

BARC_soil_clean <- BARC_raw %>% 
  filter(`sample type` == "soil" & depth == "0_5") %>% 
  select(Plot = quadret, 
         d13C_soil = d13C, 
         d15N_soil = d15N)

BARC_clean <- BARC_worm_clean %>%
  left_join(BARC_soil_clean, by = "Plot") %>%
  mutate(Dataset = "BARC") %>%
  select(Dataset, Plot, Species, d13C_worm, d15N_worm, d13C_soil, d15N_soil)


### SERC_2011
SERC_2011_worm_clean <- SERC_raw %>% 
  filter(Year == 2011 & `Sample _type` == "earthworm") %>% 
  select(Plot = Quadret, 
         Species = species, 
         d13C_worm = d13C, 
         d15N_worm = d15N)

SERC_2011_soil_clean <- SERC_raw %>% 
  filter(Year == 2011 & `Sample _type` == "soil" & Soil_depth == "0_5") %>% 
  select(Plot = Quadret, 
         d13C_soil = d13C, 
         d15N_soil = d15N) %>% 
  group_by(Plot) %>%
  summarise(d13C_soil = mean(d13C_soil), 
            d15N_soil = mean(d15N_soil))

SERC_2011_clean <- SERC_2011_worm_clean %>%
  left_join(SERC_2011_soil_clean, by = "Plot") %>%
  mutate(Dataset = "SERC_2011") %>%
  select(Dataset, Plot, Species, d13C_worm, d15N_worm, d13C_soil, d15N_soil)


### SERC_2013
SERC_2013_worm_clean <- SERC_raw %>% 
  filter(Year == 2013 & `Sample _type` == "earthworm") %>% 
  select(Plot = Quadret, 
         Species = species, 
         d13C_worm = d13C, 
         d15N_worm = d15N)

SERC_2013_soil_clean <- SERC_raw %>% 
  filter(Year == 2013 & `Sample _type` == "soil" & Soil_depth == "0_5") %>% 
  select(Plot = Quadret, 
         d13C_soil = d13C, 
         d15N_soil = d15N) %>% 
  group_by(Plot) %>%
  summarise(d13C_soil = mean(d13C_soil), 
            d15N_soil = mean(d15N_soil))

SERC_2013_clean <- SERC_2013_worm_clean %>%
  left_join(SERC_2013_soil_clean, by = "Plot") %>%
  mutate(Dataset = "SERC_2013") %>%
  select(Dataset, Plot, Species, d13C_worm, d15N_worm, d13C_soil, d15N_soil)


### Merge all datasets
all_data_clean <- bind_rows(BiodiversiTREE_clean1, 
                            BiodiversiTREE_clean2,
                            BARC_clean,
                            SERC_2011_clean,
                            SERC_2013_clean) %>%
  rename(d13C_soil_0_5 = d13C_soil,
         d15N_soil_0_5 = d15N_soil) %>%
  mutate(Species = plyr::mapvalues(Species, 
                             from = c("caliginosa", 
                                      "trapezoides",
                                      "cyaneum",
                                      "lonnbergi",               
                                      "rubellus",
                                      "terrestris",
                                      "hilgendorfi"), 
                             to = c("Aporrectodea caliginosa", 
                                    "Aporrectodea trapezoides",
                                    "Octolasion cyaneum",
                                    "Eisenoides lonnbergi",
                                    "Lumbricus rubellus",
                                    "Lumbricus terrestris",
                                    "Amynthas hilgendorfi")))

write_rds(all_data_clean, "./Output/Data_clean/all_data_clean.rds")







