## ------------------------------------------------------------------------
## Title: Compare soil stable isotope signatures of plot blocks in the 
##        BiodiversiTREE research area
##
## Author: Gen-Chang Hsu
##
## Date: 2021-06-29
##
## Description: 
## 1. Permutation test comparing the pairwise difference in d13C, d15N, and both
##    ratios of the soil (0-2 and 2-5 cm soil depth) between plot blocks in the
##    BiodiversiTREE research area
## 
## Notes:
## 
## ------------------------------------------------------------------------
set.seed(123)


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(readxl)


# Import files ------------------------------------------------------------
BiodiversiTREE_soil_raw <- read_xlsx("./Data_raw/BiodiversiTREE_data.xlsx", sheet = 2)


# Code starts here ---------------------------------------------------------

### Soil data
BDTR_soil_df <- BiodiversiTREE_soil_raw %>%
  select(Plot, Depth, ends_with(c("d13C", "d15N"))) %>%
  rename("d13C_soil" = "Linear corr d13C",
         "d15N_soil" = "Linear corr d15N") %>%
  mutate(Depth = factor(Depth, levels = c("0-2", "2-5", "5-10", "10-20"), ordered = T)) %>%
  arrange(Plot, Depth) %>%
  filter(Plot %in% c(16, 17, 28, 30, 37, 68, 70, 71, 72, 39, 40)) %>%
  filter(Depth %in% c("0-2", "2-5")) %>%
  mutate(Block = case_when(Plot %in% c(16, 17) ~ "Block1",
                           Plot %in% c(28, 30, 37) ~ "Block2",
                           Plot %in% c(68, 70, 71, 72) ~ "Block3",
                           Plot %in% c(39, 40) ~ "Block4"))

### Permutation
Perm_out_list <- map(c("0-2", "2-5"), function(DEPTH){
  Perm_out <- combn(unique(BDTR_soil_df$Block), m = 2) %>%
    t() %>%
    as_tibble() %>%
    transmute(Block1 = V1, Block2 = V2) %>%
    mutate(
      d13C_perm = map2_dbl(Block1, Block2, function(x, y){
        
        B1 <- filter(BDTR_soil_df, Depth == DEPTH & Block == x)$`d13C_soil`
        B2 <- filter(BDTR_soil_df, Depth == DEPTH & Block == y)$`d13C_soil`
        
        dif_obs <- mean(B1) - mean(B2)
        
        dif_perm <- sapply(1:999, function(x){
          samp <- sample(c(B1, B2), replace = F)
          B1_perm <- samp[length(B1)]
          B2_perm <- samp[(length(B1) + 1):(length(B1) + length(B2))]
          dif_perm <- mean(B1_perm) - mean(B2_perm)
        })
        
        P_val <- sum(abs(dif_obs) < abs(dif_perm))/(length(dif_perm) + 1)
      }),
      
      d15N_perm = map2_dbl(Block1, Block2, function(x, y){
        
        B1 <- filter(BDTR_soil_df, Depth == DEPTH & Block == x)$`d15N_soil`
        B2 <- filter(BDTR_soil_df, Depth == DEPTH & Block == y)$`d15N_soil`
        
        dif_obs <- mean(B1) - mean(B2)
        
        dif_perm <- sapply(1:999, function(x){
          samp <- sample(c(B1, B2), replace = F)
          B1_perm <- samp[length(B1)]
          B2_perm <- samp[(length(B1) + 1):(length(B1) + length(B2))]
          dif_perm <- mean(B1_perm) - mean(B2_perm)
        })
        
        P_val <- sum(abs(dif_obs) < abs(dif_perm))/(length(dif_perm) + 1)
      }))
}) %>% `names<-`(c("0-2", "2-5"))

### Write the file
Perm_out_df <- Perm_out_list %>%
  bind_rows(.id = "Depth") %>%
  write_csv(., "./Output/Data_clean/Soil_perm.csv")


