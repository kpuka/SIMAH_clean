
# Race x Lifestyle Differential Vulnerability & Exposure Project
## Causal Mediation Sensitivity 3
# Analysis for separate age subgroups

# LOAD DATA AND SET FILE LOCATIONS

# load libraries
library(tidyverse)  # data management
library(timereg)    # additive survival models
library(VGAM)       # multinomial regression, needed for causal mediation
library(MASS)       # needed for causal mediation functions
library(foreach)    # to bootstrap


# Specify file locations 
data   <- "C:/.../Processed data/"  # Location of data
model  <- "C:/.../Sensitivity/"     # Location of model output
output <- "C:/.../Output/"

# Load data
nhis        <- readRDS (paste0(data, "nhis18_85.rds"))
nhis_male   <- filter(nhis, female==0)
nhis_female <- filter(nhis, female==1)
source("7_CMed sensitivity 3 function.R")


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


# WOMEN: Bootstrap Causal Mediation ----------------------------------------------------------------------------------------

# Ages 18-59 **************************************************************************************************************
# Analysis 
set.seed(1235)
CMed_boot_w_18_59 <- bootstrap_CMed_ageStrat(nhis_female, reps=1000, prop=0.20, start_time=18, end_time=59.99) # Run analysis
saveRDS(CMed_boot_w_18_59, file.path(model, "CMed_boot_w_18_59.rds"))        # Save bootstrap results
#CMed_boot_w_18_59 <- readRDS(file.path(model, "CMed_boot_w_18_59.rds"))      # load bootstrap results

CMed_women_18_59 <- as.data.frame(do.call(cbind, CMed_boot_w_18_59)) %>% format_CMed_ageStrat () # Compute CI and format results 
CMed_women_18_59                                                                                 # print results 


# Ages 60-69 **************************************************************************************************************
# Analysis 
set.seed(1235)
CMed_boot_w_60_69 <- bootstrap_CMed_ageStrat(nhis_female, reps=1000, prop=0.20, start_time=60, end_time=69.99) # Run analysis
saveRDS(CMed_boot_w_60_69, file.path(model, "CMed_boot_w_60_69.rds"))        # Save bootstrap results
#CMed_boot_w_60_69 <- readRDS(file.path(model, "CMed_boot_w_60_69.rds"))     # load bootstrap results

CMed_women_60_69 <- as.data.frame(do.call(cbind, CMed_boot_w_60_69)) %>% format_CMed_ageStrat () # Compute CI and format results 
CMed_women_60_69                                                                                 # print results 

    
# Ages 70-85 **************************************************************************************************************
# Analysis 
set.seed(1235)
CMed_boot_w_70_85 <- bootstrap_CMed_ageStrat(nhis_female, reps=1000, prop=0.20, start_time=70, end_time=85) # Run analysis
saveRDS(CMed_boot_w_70_85, file.path(model, "CMed_boot_w_70_85.rds"))        # Save bootstrap results
#CMed_boot_w_70_85 <- readRDS(file.path(model, "CMed_boot_w_70_85.rds"))     # load bootstrap results

CMed_women_70_85 <- as.data.frame(do.call(cbind, CMed_boot_w_70_85)) %>% format_CMed_ageStrat () # Compute CI and format results 
CMed_women_70_85                                                                                 # print results 


# MEN: Bootstrap Causal Mediation ----------------------------------------------------------------------------------------

# Ages 18-59 **************************************************************************************************************
# Analysis 
set.seed(1235)
CMed_boot_m_18_59 <- bootstrap_CMed_ageStrat(nhis_male, reps=1000, prop=0.20, start_time=18, end_time=59.99) # Run analysis
saveRDS(CMed_boot_m_18_59, file.path(model, "CMed_boot_m_18_59.rds"))        # Save bootstrap results
#CMed_boot_m_18_59 <- readRDS(file.path(model, "CMed_boot_m_18_59.rds"))     # load bootstrap results

CMed_men_18_59 <- as.data.frame(do.call(cbind, CMed_boot_m_18_59)) %>% format_CMed_ageStrat () # Compute CI and format results 
CMed_men_18_59                                                                                 # print results 


# Ages 60-69 **************************************************************************************************************
# Analysis 
set.seed(1235)
CMed_boot_m_60_69 <- bootstrap_CMed_ageStrat(nhis_male, reps=1000, prop=0.20, start_time=60, end_time=69.99) # Run analysis
saveRDS(CMed_boot_m_60_69, file.path(model, "CMed_boot_m_60_69.rds"))         # Save bootstrap results
#CMed_boot_m_60_69 <- readRDS(file.path(model, "CMed_boot_m_60_69.rds"))      # load bootstrap results

CMed_men_60_69 <- as.data.frame(do.call(cbind, CMed_boot_m_60_69)) %>% format_CMed_ageStrat () # Compute CI and format results 
CMed_men_60_69                                                                                 # print results 


# Ages 70-85 **************************************************************************************************************
# Analysis 
set.seed(1235)
CMed_boot_m_70_85 <- bootstrap_CMed_ageStrat(nhis_male, reps=1000, prop=0.20, start_time=70, end_time=85) # Run analysis
saveRDS(CMed_boot_m_70_85, file.path(model, "CMed_boot_m_70_85.rds"))        # Save bootstrap results
#CMed_boot_m_70_85 <- readRDS(file.path(model, "CMed_boot_m_70_85.rds"))     # load bootstrap results

CMed_men_70_85 <- as.data.frame(do.call(cbind, CMed_boot_m_70_85)) %>% format_CMed_ageStrat () # Compute CI and format results 
CMed_men_70_85                                                                                 # print results 




# COMBINE
colnames(CMed_men_18_59) <- paste0("men_18_59_", colnames(CMed_men_18_59))
colnames(CMed_men_60_69) <- paste0("men_60_69_", colnames(CMed_men_60_69))
colnames(CMed_men_70_85) <- paste0("men_70_85_", colnames(CMed_men_70_85))

colnames(CMed_women_18_59) <- paste0("women_18_59_", colnames(CMed_women_18_59))
colnames(CMed_women_60_69) <- paste0("women_60_69_", colnames(CMed_women_60_69))
colnames(CMed_women_70_85) <- paste0("women_70_85_", colnames(CMed_women_70_85))

CMed_table <- cbind(CMed_women_18_59, CMed_women_60_69, CMed_women_70_85,
                    CMed_men_18_59, CMed_men_60_69, CMed_men_70_85)

CMed_table<-cbind(CMed_women_18_59, CMed_women_60_69) %>% 
  rename(race = women_18_59_race, term = women_18_59_term) %>% 
  dplyr::select(race, term, contains("deaths")) %>%

write.csv(CMed_table, file=paste0(output, "Causal Mediation Results age stratified.csv")) # save results


