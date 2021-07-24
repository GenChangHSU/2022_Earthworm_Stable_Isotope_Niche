## -----------------------------------------------------------------------------
## Title: Pairwise Comparisons of Earthworm Species Total Isotopic Niches
##
## Author: Gen-Chang Hsu
##
## Date: 2021-07-13
##
## Description: 
## Section 1. Test the pairwise differences in species total isotopic niches using 
##            PERMANOVA and PERMDISP.
## 
## Notes:
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(hdrcde)
library(vegan)
library(SIBER)


# Import files -----------------------------------------------------------------
all_data_clean <- readRDS("./Output/Data_clean/all_data_clean.rds")


############################### Code starts here ###############################

# Section 1 ---------------------------------------------------------------
### Background-adjusted earthworm stable isotope values
adjusted_data <- all_data_clean %>%
  group_by(Dataset) %>%
  
  # Compute the differences between plot-level soil SI values and 
  # site-level (dataset-level) grand mean values (i.e., background references)
  mutate(d13C_soil_grand_mean = mean(unique(d13C_soil_0_5)),      
         d15N_soil_grand_mean = mean(unique(d15N_soil_0_5))) %>%
  
  # Adjust the earthworm SI values by shifting back the differences between 
  # plot-level soil SI values and grand mean values 
  mutate(d13C_worm_adjusted = d13C_worm - (d13C_soil_0_5 - d13C_soil_grand_mean),
         d15N_worm_adjusted = d15N_worm - (d15N_soil_0_5 - d15N_soil_grand_mean)) %>%
  ungroup() %>%
  
  # Reorder the levels in column "Dataset" for plotting purpose
  mutate(Dataset = factor(Dataset, levels = unique(Dataset), ordered = T))


### Summary statistics of species isotopic niches
# A vector of the dataset names
dataset_list <- unique(as.character(adjusted_data$Dataset))

Isoniche_test_list <- lapply(dataset_list, function(dataset){
    
  # Get the species names
  Sp_names <- adjusted_data %>% 
    filter(Dataset == dataset) %>%
    group_by(Species) %>%
    mutate(n.ind = n()) %>%
    filter(n.ind > 3) %>%
    .$Species %>% 
    unique() %>% 
    sort()
  
  # PERMANOVA and PERMDISP tests for multivariate differences in species isoniches
  PERM_tests <- combn(Sp_names, 2) %>%
    t() %>% 
    `colnames<-`(c("Sp1", "Sp2")) %>%
    as_tibble() %>%
    mutate(PERMANOVA = map2_chr(Sp1, Sp2, function(x, y){
      
        data <- filter(filter(adjusted_data, Dataset == dataset), Species %in% c(x, y))
        PERMANOVA_out <- adonis(d15N_worm_adjusted + d13C_worm_adjusted ~ Species,
                                data = data,
                                method = "eu", 
                                permutations = 999)
        
        PERMANOVA_F <- ifelse(PERMANOVA_out$aov.tab$F.Model[1] > 0.01, 
                              round(PERMANOVA_out$aov.tab$F.Model[1], 2),
                              round(PERMANOVA_out$aov.tab$F.Model[1], 3))
        PERMANOVA_P <- ifelse(PERMANOVA_out$aov.tab$`Pr(>F)`[1] > 0.01, 
                              round(PERMANOVA_out$aov.tab$`Pr(>F)`[1], 2),
                              round(PERMANOVA_out$aov.tab$`Pr(>F)`[1], 3))
        
        return(paste0("F = ", PERMANOVA_F, ", P = ", PERMANOVA_P))
    })) %>%
    mutate(PERMDISP = map2_chr(Sp1, Sp2, function(x, y){
      
        data <- filter(filter(adjusted_data, Dataset == dataset), Species %in% c(x, y))
        Dist <- dist(data[, c("d15N_worm_adjusted", "d13C_worm_adjusted")])  
        PERMDISP_out <- betadisper(Dist, group = data$Species)
        PERMDISP_out_perm <- permutest(PERMDISP_out, pairwise = TRUE, permutations = 999)
        
        PERMDISP_F <- ifelse(PERMDISP_out_perm$tab$`F`[1] > 0.01, 
                              round(PERMDISP_out_perm$tab$`F`[1], 2),
                              round(PERMDISP_out_perm$tab$`F`[1], 3))
        PERMDISP_P <- ifelse(PERMDISP_out_perm$tab$`Pr(>F)`[1] > 0.01, 
                              round(PERMDISP_out_perm$tab$`Pr(>F)`[1], 2),
                              round(PERMDISP_out_perm$tab$`Pr(>F)`[1], 3))
        
        return(paste0("F = ", PERMDISP_F, ", P = ", PERMDISP_P))
    })) 
  
  return(PERM_tests)
  
}) %>% 
  `names<-`(dataset_list)


### Convert the list into a dataframe
Isoniche_test_df <- Isoniche_test_list %>% 
  bind_rows(.id = "Dataset") 

write_rds(Isoniche_test_df, "./Output/Data_clean/Isoniche_test_df.rds")
write_csv(Isoniche_test_df, "./Output/Data_clean/Isoniche_test_df.csv")


