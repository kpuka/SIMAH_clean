
# SES x Lifestyle Differential Vulnerability & Exposure Project
# Objective 2: Causal Mediation File 


# LOAD DATA AND SET FILE LOCATIONS

# load libraries
library(tidyverse)  # data management
library(timereg)    # additive survival models
library(VGAM)       # multinomial regression, needed for causal mediation
library(MASS)       # needed for causal mediation functions
library(knitr)      # formatted table



# Specify the data and output file locations

data    <- "C:/.../nhis/Data/"
output  <- "C:/.../nhis/SES x lifestyle/Output/"
source("Function - Format Results.R")
source("Function - CausalMed Results.R")


# Load data
nhis        <- readRDS (file.path(data, "nhis.rds"))
nhis_male   <- readRDS (file.path(data, "nhis_male.rds"))
nhis_female <- readRDS (file.path(data, "nhis_female.rds"))



# OBJECTIVE 2: Causal Mediation

# The causal mediation analyses involves four steps:
# 1) Fit separate multinomial logistic regressions with each mediator (M1, M2, M3, and M4) as the outcome.
# 2) Create copies of the dataset to account for all possible combinations of the exposure and mediators (3^4*= 81); the dataset was expanded from 229,994 to 18,629,514 (women) and 185,770 to 15,047,370 (men) pseudobservations. 
# 3) Using the expanded dataset, calculate weights for each mediator using the predicted probabilities from Step 1. 
# 4) Fit a marginal structural model using Aalen additive hazards with the weight as weight and the id as a cluster level; this ensures thatrobust standard errors are calculated. The model with robust variance and resampling (robust=TRUE) was not used because of computation limiations.  

# For more details and theoretical justification/description see:
# Lange et al. 2014 https//doi.org/10.1093/aje/kwt270
# Lange et al. 2012 https//doi.org/10.1093/aje/kwr525
# Lange et al. 2011 https//doi.org/10.1097/EDE.0b013e31821c680c


# Data Preparation ----------------------------------------------------------------------------------------
causal_mediation_prep(nhis_female) %>% saveRDS(paste0(output, "expandedData_fem.rds"))
causal_mediation_prep(nhis_male)   %>% saveRDS(paste0(output, "expandedData_male.rds"))


# Run Analyses, WOMEN ----------------------------------------------------------------------------------------------------------------

# Load data
expandedData <-readRDS(file.path(output, "expandedData_fem.rds"))

hist(newMyData$weightM)

# Run model
CMed_f <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) + 
                                                        const(A.edu) * const(edu_M2.smk) +
                                                        const(A.edu) * const(edu_M3.bmi) +
                                                        const(A.edu) * const(edu_M4.phy) +
                                                        const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                          data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
                
saveRDS(CMed_f, file.path(output, "CMed_f.rds"))  # Save model results
     
# Load model and view results
CMed_model <-readRDS(file.path(output, "CMed_f.rds"))  # load model (if needed)
Low_Education <- c(1,3,5,7,9,29,33,37,41)              # List the coefficients of interest for 'low education'
format_CMed (CMed_model, Low_Education) %>% kable()    # print formatted results



 
     

# Run Analyses, MEN ----------------------------------------------------------------------------------------------------------------

# Load data
expandedData <- readRDS(file.path(output, "expandedData_male.rds"))


# Run Model
CMed_m <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) +
                                                        const(A.edu) * const(edu_M2.smk) +
                                                        const(A.edu) * const(edu_M3.bmi) +
                                                        const(A.edu) * const(edu_M4.phy) +
                                                        const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                          data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)

saveRDS(CMed_m, file.path(output, "CMed_m.rds"))       # Save model results


# Load model and view results
CMed_model <-readRDS(file.path(output, "CMed_m.rds"))  # load model (if needed)
Low_Education <- c(1,3,5,7,9,29,33,37,41)              # List the coefficients of interest for 'low education'
format_CMed (CMed_model, Low_Education) %>% kable()    # print formatted results

