
# Race x Lifestyle Differential Vulnerability & Exposure Project
# Sensitivity Analysis 2 - SES not included as covariate


# load libraries
library(tidyverse)  # data management
library(timereg)    # additive survival models
library(VGAM)       # multinomial regression, needed for causal mediation
library(MASS)       # needed for causal mediation functions
library(foreach)    # to bootstrap


# Personal computer; specify locations 
data   <- "C:/.../Processed data/"            # Location of data
model  <- "C:/.../Sensitivity/"  # Location of model output
output <- "C:/.../Output/"

# SCC; ; specify locations 
# setwd("/external/mgmt3/imaging/scratch/Imhpr/kpuka/nhis/")
# data    <- "Data/"
# model  <- "Model/"


# Load data
nhis        <- readRDS (paste0(data, "nhis18_85.rds"))
nhis_male   <- filter(nhis, female==0)
nhis_female <- filter(nhis, female==1)
source("5_CMed sensitivity 2 function.R")


# Set up parallel processing ************************************************************************************************

# Linux 
# library(doMC)
# foreach::getDoParWorkers()      # Identify # of cores that will be used
# registerDoMC(5)                 # Specify number of cores to use  
# foreach::getDoParWorkers()      # Identify # of cores that will be used


# Windows 
# library(parallel)   
# library(doParallel) 
# foreach::getDoParWorkers()                # Identify # of cores that will be used
# cl <- makeCluster(4, outfile = "log.txt") # Specify number of cores to use  
# registerDoParallel(cl)                    # Specify number of cores to use  
# foreach::getDoParWorkers()                # Identify # of cores that will be used


# WOMEN: Bootstrap Causal Mediation *****************************************************************************************

# Analysis
# set.seed(1235)
# CMed_boot_w <- bootstrap_CMed_noEDU(nhis_female, reps=1000, prop=0.20)  # Run analysis using bootstrap
# saveRDS(CMed_boot_w, file.path(model, "CMed_boot_noEDU_w.rds"))        # Save bootstrap results
CMed_boot_w <- readRDS(file.path(model, "CMed_boot_noEDU_w.rds"))        # load bootstrap results

# Results 
CMed_women <- as.data.frame(do.call(cbind, CMed_boot_w)) %>% 
  format_CMed_noEDU()                                                # Compute CI and format results 
CMed_women                                                           # print results 



# MEN: Bootstrap Causal Mediation ********************************************************************************************

# Analysis
# set.seed(1235)
# CMed_boot_m <- bootstrap_CMed_noEDU(nhis_male, reps=1000, prop=0.20)  # Run analysis using bootstrap
# saveRDS(CMed_boot_m, file.path(model, "CMed_boot_noEDU_m.rds"))       # Save bootstrap results
CMed_boot_m <- readRDS(file.path(model, "CMed_boot_noEDU_m.rds"))       # load bootstrap results


# Results 
CMed_men <- as.data.frame(do.call(cbind, CMed_boot_m)) %>%
  format_CMed_noEDU()                                                    # Compute CI and format results 
CMed_men                                                           # print results 


# COMBINE Results
colnames(CMed_men)   <- paste0("men_", colnames(CMed_men))
colnames(CMed_women) <- paste0("women_", colnames(CMed_women))
CMed_table <- cbind(CMed_men, CMed_women) %>% rename(race = men_race, term = men_term)
CMed_table <- CMed_table[c(1:4, 7:8)]
CMed_table
write.csv(CMed_table, file=paste0(output, "Causal Mediation Results noEDU.csv")) # save results

