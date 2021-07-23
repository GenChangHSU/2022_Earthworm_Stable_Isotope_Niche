## -----------------------------------------------------------------------------
## Title: Bayesian Standard Ellipses Areas and Overlap between Species
##
## Author: Gen-Chang Hsu
##
## Date: 2021-07-13
##
## Description: 
## Section 1. Compute SEAb for each species and pairwise percentage SEAb overlaps.                
##  
## Notes:
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(SIBER)
library(ellipse)
library(coda)
library(ggsci)
library(ggmcmc)


# Import files -----------------------------------------------------------------
all_data_clean <- readRDS("./Output/Data_clean/all_data_clean.rds")


############################### Code starts here ###############################

# Section 1 ---------------------------------------------------------------
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


### Bayesian SEA for each species
dataset_list <- unique(as.character(adjusted_data$Dataset))

SEAb_list <- lapply(dataset_list, function(dataset) {

  # Create SIBER object
  SIBER_data <- adjusted_data %>%
    filter(Dataset == dataset) %>%
    group_by(Species) %>%
    mutate(n.ind = n()) %>%
    filter(n.ind > 3) %>%
    select(
      iso1 = d13C_worm_adjusted,
      iso2 = d15N_worm_adjusted,
      group = Species
    ) %>%
    arrange(group) %>%
    ungroup() %>%
    mutate(group = plyr::mapvalues(
      group,
      unique(group),
      1:length(unique(group))
    )) %>%
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

  # Bayesian model setup
  parms <- list()
  parms$n.iter <- 30000  # Number of iterations
  parms$n.burnin <- 3000  # Burn-in
  parms$n.thin <- 10  # Thinning
  parms$n.chains <- 3  # Number of MCMC chains
  parms$save.output <- TRUE  # Save jags output for diagnostics
  parms$save.dir <- paste0(getwd(), "/", "Output/Jags")  # Folder to save jags output
  
  # Prior (the inverse Wishart prior)
  priors <- list()
  priors$R <- 1*diag(2)
  priors$k <- 2
  priors$tau.mu <- 0.001

  # Fit the model
  ellipses.posterior <- siberMVN(SIBER_data, parms, priors)
  SEA.B <- siberEllipses(ellipses.posterior) %>%
    as.data.frame() %>%
    `names<-`(Sp_names)
  
  # Model diagnostics
  all.files <- dir(parms$save.dir, full.names = TRUE)
  model.files <- all.files[grep("jags_output", all.files)]
  map(1:length(Sp_names), function(a){
    sp <- Sp_names[a]
    path <- paste0("./Output/Diagnostics/", dataset, "_", sp, ".pdf")
    load(model.files[a]) 
    S <- ggs(output)
    ggmcmc(S, plot = c("traceplot", "geweke"), param_page = 2, file = path)
  })
  
  # Compute standard ellipses data points from the posterior draws (for plotting)
  SEAb_points <- map(ellipses.posterior, function(x) {
    sigma <- matrix(apply(x[, 1:4], 2, mean), 2, 2) # Covariance matrix of the ellipses
    mu <- apply(x[, 5:6], 2, mean)  # Center of the ellipses
    sea <- ellipse(sigma, centre = mu, level = pchisq(1, 2), npoints = 100)  # Standard ellipses: pchisq(1, 2)
    return(sea)
  }) %>%
    `names<-`(Sp_names)
  
  SEAb_points <- data.frame(Species = Sp_names) %>%
    mutate(SEAb_points = map(SEAb_points, function(w){as.data.frame(w)}))
  
  # Summary statistics for the posterior draws
  SEAb_sum_stats <- SEA.B %>%
    pivot_longer(cols = everything(), names_to = "Species", values_to = "SEAb") %>%
    group_by(Species) %>%
    summarise(
      SEAb_mean = mean(SEAb),
      SEAb_sd = sd(SEAb),
      SEAb_mode = quantile(SEAb, 0.5),
      `SEAb_hdr_lower` = HPDinterval(as.mcmc(SEAb), 0.95)[1],
      `SEAb_hdr_upper` = HPDinterval(as.mcmc(SEAb), 0.95)[2]
    )

  # SEAb overlap
  SEAb_overlap <- data.frame(
    Sp1 = rep(Sp_names, each = length(Sp_names)),
    Sp2 = rep(Sp_names),
    Ellipse1 = paste(1, rep(1:length(Sp_names), each = length(Sp_names)), sep = "."),
    Ellipse2 = paste(1, rep(1:length(Sp_names)), sep = ".")
  ) %>%
    bind_cols(., map2(.x = .$Ellipse1, .y = .$Ellipse2, function(x, y) {
      bayesianOverlap(x, y,
        ellipses.posterior,
        draws = 100,  # Use first 100 posterior draws to calculate overlap
        p.interval = NULL,
        n = 100  # 100 points to form the ellipses
      ) %>%
        apply(., MARGIN = 2, FUN = function(z) mean(z)) %>%
        t() %>%
        as.data.frame() %>%
        `names<-`(c("area.1", "area.2", "overlap"))
    }) %>% bind_rows()) %>%
    select(-Ellipse1, -Ellipse2, -area.2) %>%
    rename(
      SEAb = area.1,
      Area = overlap
    ) %>%
    mutate(`Percent_overlap_%` = Area / SEAb * 100) %>%
    filter(Sp1 != Sp2) %>%
    rename(Species = Sp1) %>%
    group_by(Species) %>%
    rename(SEAb_overlap_with = Sp2) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    nest(SEAb_overlap = c(SEAb_overlap_with, Area, `Percent_overlap_%`)) %>%
    select(-SEAb)

  # Combine all the results
  SEAb_all <- SEAb_sum_stats %>%
    left_join(SEAb_overlap, by = "Species") %>%
    left_join(SEAb_points, by = "Species")

  return(SEAb_all)
  
}) %>% 
  `names<-`(dataset_list)


### Convert the list into a dataframe
SEAb_df <- SEAb_list %>% 
  bind_rows(.id = "Dataset") 

write_rds(SEAb_df, "./Output/Data_clean/SEAb.rds")
write_csv(select(SEAb_df, -SEAb_overlap, -SEAb_points), "./Output/Data_clean/SEAb.csv")


