
# Race x Lifestyle Differential Vulnerability & Exposure Project
# Sensitivity Analysis - One Mediator at a time


# load libraries
library(tidyverse)  # data management
library(timereg)    # additive survival models
library(VGAM)       # multinomial regression, needed for causal mediation
library(MASS)       # needed for causal mediation functions

# Specify file locations 
data   <- "C:/.../Processed data/"            # Location of data
output <- "C:/.../Output/"



# Load data
nhis        <- readRDS (paste0(data, "nhis18_85.rds"))
nhis_male   <- filter(nhis, female==0)
nhis_female <- filter(nhis, female==1)
source("5_CMed sensitivity 1 function.R")


# MEN: Causal mediation with one mediator: 
alcohol_men <- CMed_oneVar (nhis_male, alcohol5v2, smoking4, bmi_cat, phy_act3, 3, c(1,4,35), c(2,5,39), c(3,6,43)) %>% mutate (mediator = "Alcohol") %>% relocate (mediator)
smoking_men <- CMed_oneVar (nhis_male, smoking4, alcohol5v2, bmi_cat, phy_act3, 1, c(1,4,36), c(2,5,40), c(3,6,44)) %>% mutate (mediator = "Smoking") %>% relocate (mediator)
bmi_men     <- CMed_oneVar (nhis_male, bmi_cat, alcohol5v2, smoking4, phy_act3, 2, c(1,4,36), c(2,5,40), c(3,6,44)) %>% mutate (mediator = "BMI") %>% relocate (mediator)
phy_men     <- CMed_oneVar (nhis_male, phy_act3, alcohol5v2, smoking4, bmi_cat, 3, c(1,4,37), c(2,5,41), c(3,6,45)) %>% mutate (mediator = "Physical") %>% relocate (mediator)

CMed_men <- rbind(alcohol_men, smoking_men, bmi_men, phy_men)
view(CMed_men)


# WOMEN: Causal mediation with one mediator: 
alcohol_women <- CMed_oneVar (nhis_female, alcohol5v2, smoking4, bmi_cat, phy_act3, 3, c(1,4,35), c(2,5,39), c(3,6,43)) %>% mutate (mediator = "Alcohol") %>% relocate (mediator)
smoking_women <- CMed_oneVar (nhis_female, smoking4, alcohol5v2, bmi_cat, phy_act3, 1, c(1,4,36), c(2,5,40), c(3,6,44)) %>% mutate (mediator = "Smoking") %>% relocate (mediator)
bmi_women     <- CMed_oneVar (nhis_female, bmi_cat, alcohol5v2, smoking4, phy_act3, 2, c(1,4,36), c(2,5,40), c(3,6,44)) %>% mutate (mediator = "BMI") %>% relocate (mediator)
phy_women     <- CMed_oneVar (nhis_female, phy_act3, alcohol5v2, smoking4, bmi_cat, 3, c(1,4,37), c(2,5,41), c(3,6,45)) %>% mutate (mediator = "Physical") %>% relocate (mediator)

CMed_women <- rbind(alcohol_women, smoking_women, bmi_women, phy_women)
view(CMed_women)


# Combined
colnames(CMed_men)   <- paste0("men_", colnames(CMed_men))
colnames(CMed_women) <- paste0("women_", colnames(CMed_women))

CMed_table <- cbind(CMed_men, CMed_women) 
CMed_table <- CMed_table[c(1,3:5,9:10)] %>% 
  rename(mediator = men_mediator, label = men_label)

view(CMed_table)
write.csv(CMed_table, file=paste0(output, "Table_e2 Causal Mediation, one mediator.csv")) # save results

