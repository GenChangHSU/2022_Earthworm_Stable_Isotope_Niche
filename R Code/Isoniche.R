## ------------------------------------------------------------------------
## Title: Analysis of Earthworm Isotopic Niches
##
## Author: Gen-Chang Hsu
##
## Date: 2021-02-26
##
## Description: Compute corrected standard ellipse area (SEAc) for 
## each species and percent SEAc overlap between species pairs. 
## Test the difference in species isotopic niches using PERMANOVA and 
## PERMADISP.
##
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





