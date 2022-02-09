
# SES x Lifestyle Differential Vulnerability & Exposure Project
# Objective 1: Aalen additive hazard models file


# LOAD DATA AND SET FILE LOCATIONS

# load libraries
library(tidyverse)  # data management
library(skimr)      # descriptive statistics
library(survival)   # surivval analyses
library(survminer)  # surivval analyses
library(timereg)    # additive survival models
memory.limit(size=1e+13)


# Specify the data and output file locations
data    <- "C:/.../nhis/Data/"
output  <- "C:/.../nhis/SES x lifestyle/Output/"
source("Function - Format Results.R")

    
# Load data
nhis        <- readRDS (file.path(data, "nhis.rds"))
nhis_male   <- readRDS (file.path(data, "nhis_male.rds"))
nhis_female <- readRDS (file.path(data, "nhis_female.rds"))


# OBJECTIVE 1: Joint Effects, Hazard Models - Stratified by Sex

# The effect estimates from the model can be directly interpreted as the number of additional events (deaths) per 1 person year at risk
# Two different versions of the model were ran identify the interaction effect (model with the interaction term) and the joint effect (model with interacting variable)

## Create function specifying the Aalen models --------------------------------------------------------------------------------------

interaction_model <- function(data, lifestyle){
  data <- data %>%  mutate (lifestyle = {{lifestyle}})
  model <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(edu.factor)*const(lifestyle) + 
                  const(married.factor) + ethnicity.factor + const(factor(srvy_yr)),  data = data)
 return(model)  
}

jointeffect_model <- function(data, lifestyle_edu){
  data <- data %>%  mutate (lifestyle_edu = {{lifestyle_edu}})
  model <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(lifestyle_edu) + 
                  const(married.factor) + ethnicity.factor + const(factor(srvy_yr)),  data = data)
  return(model)
}



# First, Run each model and save results -------------------------------------------------------------------------------------------

# Alcohol Use 
interaction_model(nhis_female, alcohol5v2.factor) %>% saveRDS(paste0(output, "aalen_alc_f.rds"))
jointeffect_model(nhis_female, edu.alc)           %>% saveRDS(paste0(output, "aalen_alc_f2.rds"))

interaction_model(nhis_male, alcohol5v2.factor) %>% saveRDS(paste0(output, "aalen_alc_m.rds"))
jointeffect_model(nhis_male, edu.alc)           %>% saveRDS(paste0(output, "aalen_alc_m2.rds"))

  
# Smoking 
interaction_model(nhis_female, smoking4.factor) %>% saveRDS(paste0(output, "aalen_smk_f.rds"))
jointeffect_model(nhis_female, edu.smk)         %>% saveRDS(paste0(output, "aalen_smk_f2.rds"))

interaction_model(nhis_male, smoking4.factor) %>% saveRDS(paste0(output, "aalen_smk_m.rds"))
jointeffect_model(nhis_male, edu.smk)         %>% saveRDS(paste0(output, "aalen_smk_m2.rds"))

      
# BMI
interaction_model(nhis_female, bmi_cat.factor) %>% saveRDS(paste0(output, "aalen_bmi_f.rds"))
jointeffect_model(nhis_female, edu.bmi)        %>% saveRDS(paste0(output, "aalen_bmi_f2.rds"))

interaction_model(nhis_male, bmi_cat.factor) %>% saveRDS(paste0(output, "aalen_bmi_m.rds"))
jointeffect_model(nhis_male, edu.bmi)        %>% saveRDS(paste0(output, "aalen_bmi_m2.rds"))


# Physical Activity
interaction_model(nhis_female, phy_act3.factor) %>% saveRDS(paste0(output, "aalen_phy_f.rds"))
jointeffect_model(nhis_female, edu.phy)         %>% saveRDS(paste0(output, "aalen_phy_f2.rds"))

interaction_model(nhis_male, phy_act3.factor) %>% saveRDS(paste0(output, "aalen_phy_m.rds"))
jointeffect_model(nhis_male, edu.phy)         %>% saveRDS(paste0(output, "aalen_phy_m2.rds"))


      


# Second, Load and view model results -------------------------------------------------------------------------------------------

# A "aalen_10000py" function was created to extract the coefficients and multiple them by 10,000 to get estimates per 10,000 person years
# The results from the function to extract results (as written below) pertain to the effect of:
    # "Low SES & 'best' lifestyle", 
    # "High SES & 'poor' lifestyle", 
    # "Low SES & 'poor' behavior", 
    # "Interaction effect"
  
      
# Alcohol, Women       
aalen  <- readRDS(paste0(output, "aalen_alc_f.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_alc_f2.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31) # print results of interest

# Alcohol, Men 
aalen  <- readRDS(paste0(output, "aalen_alc_m.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_alc_m2.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31)



#Smoking, Women 
aalen  <- readRDS(paste0(output, "aalen_smk_f.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_smk_f2.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)

#Smoking, Men 
aalen  <- readRDS(paste0(output, "aalen_smk_m.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_smk_m2.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)



#BMI, Women 
aalen  <- readRDS(paste0(output, "aalen_bmi_f.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_bmi_f2.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)

#BMI, Men 
aalen  <- readRDS(paste0(output, "aalen_bmi_m.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_bmi_m2.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)



#Physical activity, Women 
aalen  <- readRDS(paste0(output, "aalen_phy_f.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_phy_f2.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 3); aalen_10000py(aalen2, 4); aalen_10000py(aalen, 23)

#Physical activity, Mn 
aalen  <- readRDS(paste0(output, "aalen_phy_m.rds"))        
aalen2 <- readRDS(paste0(output, "aalen_phy_m2.rds"))        
aalen_10000py(aalen, 1); aalen_10000py(aalen, 3); aalen_10000py(aalen2, 4); aalen_10000py(aalen, 23)




