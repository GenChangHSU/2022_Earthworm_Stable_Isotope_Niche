## ------------------------------------------------------------------------
## Title: Analysis of Earthworm Isotopic Niches
##
## Author: Gen-Chang Hsu
##
## Date: 2021-02-26
##
## Description: 
## 1. Compute the total niche areas, corrected standard ellipse areas (SEAc), 
##    SEAs, and 95% ellipse areas for each species
## 2. Calculate the percent overlap in SEAs and 95% ellipse areas between 
##    species pairs. 
## 3. Test the differences in species isotopic niches using PERMANOVA and 
##    PERMADISP.
##
## Notes:
##
##
## ------------------------------------------------------------------------


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(hdrcde)
library(vegan)
library(SIBER)


# Import files ------------------------------------------------------------
all_data_clean <- readRDS("./Output/Data_clean/all_data_clean.rds")


# Code start here ---------------------------------------------------------

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


### Summary statistics of species isotopic niches by dataset

# A vector of the dataset names
dataset_list <- unique(as.character(adjusted_data$Dataset))

Isoniche_list <- lapply(dataset_list, function(dataset){
    
  set.seed(123)  # For PERMANOVA and PERMADISP
  
  # Create SIBER object
  SIBER_data <- adjusted_data %>% 
    filter(Dataset == dataset) %>%
    group_by(Species) %>%
    mutate(n.ind = n()) %>%
    filter(n.ind > 3) %>%
    select(iso1 = d13C_worm_adjusted, 
           iso2 = d15N_worm_adjusted,
           group = Species) %>%
    arrange(group) %>% 
    ungroup() %>%
    mutate(group = plyr::mapvalues(group, 
                                   unique(group), 
                                   1:length(unique(group)))) %>%
    mutate(community = 1) %>% 
    as.data.frame() %>% 
    createSiberObject()
  
  # Get the species names
  Sp_names <- adjusted_data %>% 
    filter(Dataset == dataset) %>%
    group_by(Species) %>%
    mutate(n.ind = n()) %>%
    filter(n.ind > 3) %>%
    .$Species %>% 
    unique() %>% 
    sort()
  
  # Total isoniche area and SEAc
  total_area <- groupMetricsML(SIBER_data) %>% 
    as.data.frame() %>%
    `names<-`(Sp_names) %>%
    t() %>% 
    as.data.frame() %>%
    rownames_to_column(var = "Species") %>%
    select(Species, TA, SEAc)
  
  # SEA overlap between species pairs
  SEA_overlap <- data.frame(
    Sp1 = rep(Sp_names, each = length(Sp_names)),
    Sp2 = rep(Sp_names),
    Ellipse1 = paste(1, rep(1:length(Sp_names), each = length(Sp_names)), sep = "."),
    Ellipse2 = paste(1, rep(1:length(Sp_names)), sep = ".")) %>%
    bind_cols(., map2(.x = .$Ellipse1, .y = .$Ellipse2, function(x, y) {
      maxLikOverlap(x, y, SIBER_data, p.interval = NULL, n = 100)
    }) %>% bind_rows()) %>%
    select(-Ellipse1, -Ellipse2, -area.2) %>%
    rename(SEA = area.1,
           Area = overlap) %>%
    mutate(Percent_overlap = Area/SEA) %>%
    filter(Sp1 != Sp2) %>%
    rename(Species = Sp1) %>%
    group_by(Species) %>%
    rename(SEA_overlap_with = Sp2) %>%
    nest(SEA_overlap = c(SEA_overlap_with, Area, Percent_overlap))
  
  # 95% ellipse area overlap between species pairs
  `95%_EA_overlap` <- data.frame(
    Sp1 = rep(Sp_names, each = length(Sp_names)),
    Sp2 = rep(Sp_names),
    Ellipse1 = paste(1, rep(1:length(Sp_names), each = length(Sp_names)), sep = "."),
    Ellipse2 = paste(1, rep(1:length(Sp_names)), sep = ".")) %>%
    bind_cols(., map2(.x = .$Ellipse1, .y = .$Ellipse2, function(x, y) {
      maxLikOverlap(x, y, SIBER_data, p.interval = 0.95, n = 100)
    }) %>% bind_rows()) %>%
    select(-Ellipse1, -Ellipse2, -area.2) %>%
    rename(`95%_EA` = area.1,
           Area = overlap) %>%
    mutate(Percent_overlap = Area/`95%_EA`) %>%
    filter(Sp1 != Sp2) %>%
    rename(Species = Sp1) %>%
    group_by(Species) %>%
    rename(`95%_EA_overlap_with` = Sp2) %>%
    nest(`95%_EA_overlap` = c(`95%_EA_overlap_with`, Area, Percent_overlap))
  
  # PERMANOVA and PERMADISP tests for multivariate differences in species isoniches
  PERM_tests <- data.frame(
    Sp1 = rep(Sp_names, each = length(Sp_names)),
    Sp2 = rep(Sp_names)) %>%
    select(Sp1, Sp2) %>%
    mutate(PERMANOVA_pval = map2(Sp1, Sp2, function(x, y){
      if (x == y){
        return(NA)
      } else {
        data = filter(filter(adjusted_data, Dataset == dataset), Species %in% c(x, y))
        PERMANOVA_out <- adonis(d15N_worm_adjusted + d13C_worm_adjusted ~ Species,
                                data = data,
                                method = "eu", 
                                permutations = 999)
        return(PERMANOVA_out$aov.tab$`Pr(>F)`[1])
      }
    })) %>%
    mutate(PERMADISP_pval = map2(Sp1, Sp2, function(x, y){
      if (x == y){
        return(NA)
      } else {
        data <- filter(filter(adjusted_data, Dataset == dataset), Species %in% c(x, y))
        Dist <- dist(data[, c("d15N_worm_adjusted", "d13C_worm_adjusted")])  
        PERMADISP_out <- betadisper(Dist, group = data$Species)
        PERMADISP_out_perm <- permutest(PERMADISP_out, pairwise = TRUE, permutations = 999)
        return(PERMADISP_out_perm$tab$`Pr(>F)`[1])
      }
    })) %>%
    filter(Sp1 != Sp2) %>%
    rename(Species = Sp1) %>%
    group_by(Species) %>%
    rename(Test_with = Sp2) %>%
    nest(PERM_tests = c(Test_with, PERMANOVA_pval, PERMADISP_pval))
  
  # Combine all the results
  niche_sum_stat <- total_area %>%
    left_join( SEA_overlap, by = "Species") %>%
    left_join(`95%_EA_overlap`, by = "Species") %>%
    left_join(PERM_tests, by = "Species")
  
  return(niche_sum_stat)
  
}) %>% 
  `names<-`(dataset_list)


### Convert the list into a dataframe
Isoniche_df <- Isoniche_list %>% 
  bind_rows(.id = "Dataset")

write_rds(Isoniche_df, "./Output/Data_clean/Isoniche.rds")




