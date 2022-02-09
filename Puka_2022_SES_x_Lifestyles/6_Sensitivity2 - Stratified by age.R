
# SES x Lifestyle Differential Vulnerability & Exposure Project
# Sensitivity Analyses 2: Stratified by age


# LOAD DATA AND SET FILE LOCATIONS

# load libraries
library(tidyverse)  # data management
library(skimr)      # descriptive statistics
library(gmodels)    # CrossTable command
library(tableone)   # create table one
library(survival)   # surivval analyses
library(survminer)  # surivval analyses
library(timereg)    # additive survival models
library(survey)     # for survey weighted cox model
library(VGAM)       # multinomial regression, needed for causal mediation
library(MASS)       # needed for causal mediation functions


# Specify the data and output file locations
data    <- "C:/.../nhis/Data/"
output  <- "C:/.../nhis/SES x lifestyle/Output/"
source("Function - Format Results.R")
source("Function - CausalMed Results.R")


# Load data
nhis        <- readRDS (file.path(data, "nhis.rds"))
nhis_male   <- readRDS (file.path(data, "nhis_male.rds"))
nhis_female <- readRDS (file.path(data, "nhis_female.rds"))



# Aalen Models Stratified by age (Objective 1) --------------------------------------------------------------------------------------------------

## Create function specifying the Aalen models ---------------------------------------------------------------------------------------------------

interaction_model3 <- function(data, lifestyle, start_time, end_time){
  data <- data %>%  mutate (lifestyle = {{lifestyle}})
  model <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(edu.factor)*const(lifestyle) + 
                  const(married.factor) + ethnicity.factor + const(factor(srvy_yr)),  data = data, start.time=start_time, max.time=end_time)
  return(model)  
}

jointeffect_model3 <- function(data, lifestyle_edu, start_time, end_time){
  data <- data %>%  mutate (lifestyle_edu = {{lifestyle_edu}})
  model <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(lifestyle_edu) + 
                  const(married.factor) + ethnicity.factor + const(factor(srvy_yr)),  data = data, start.time=start_time, max.time=end_time)
  return(model)
}



## First, Run all models and save results -----------------------------------------------------------------------------------------------------------

# Alcohol Use, Women  
interaction_model3(nhis_female, alcohol5v2.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_alc_f_25_60.rds"))
jointeffect_model3(nhis_female, edu.alc, 25, 59.999)           %>% saveRDS(paste0(output, "aalen_alc_f2_25_60.rds"))

interaction_model3(nhis_female, alcohol5v2.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_alc_f_60_70.rds"))
jointeffect_model3(nhis_female, edu.alc, 60, 69.999)           %>% saveRDS(paste0(output, "aalen_alc_f2_60_70.rds"))

interaction_model3(nhis_female, alcohol5v2.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_alc_f_70_85.rds"))
jointeffect_model3(nhis_female, edu.alc, 70, 84.999)           %>% saveRDS(paste0(output, "aalen_alc_f2_70_85.rds"))

# Alcohol Use, Men  
interaction_model3(nhis_male, alcohol5v2.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_alc_m_25_60.rds"))
jointeffect_model3(nhis_male, edu.alc, 25, 59.999)           %>% saveRDS(paste0(output, "aalen_alc_m2_25_60.rds"))

interaction_model3(nhis_male, alcohol5v2.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_alc_m_60_70.rds"))
jointeffect_model3(nhis_male, edu.alc, 60, 69.999)           %>% saveRDS(paste0(output, "aalen_alc_m2_60_70.rds"))

interaction_model3(nhis_male, alcohol5v2.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_alc_m_70_85.rds"))
jointeffect_model3(nhis_male, edu.alc, 70, 84.999)           %>% saveRDS(paste0(output, "aalen_alc_m2_70_85.rds"))



# Smoking, Women  
interaction_model3(nhis_female, smoking4.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_smk_f_25_60.rds"))
jointeffect_model3(nhis_female, edu.smk, 25, 59.999)         %>% saveRDS(paste0(output, "aalen_smk_f2_25_60.rds"))

interaction_model3(nhis_female, smoking4.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_smk_f_60_70.rds"))
jointeffect_model3(nhis_female, edu.smk, 60, 69.999)         %>% saveRDS(paste0(output, "aalen_smk_f2_60_70.rds"))

interaction_model3(nhis_female, smoking4.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_smk_f_70_85.rds"))
jointeffect_model3(nhis_female, edu.smk, 70, 84.999)         %>% saveRDS(paste0(output, "aalen_smk_f2_70_85.rds"))

# Smoking, Men  
interaction_model3(nhis_male, smoking4.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_smk_m_25_60.rds"))
jointeffect_model3(nhis_male, edu.smk, 25, 59.999)         %>% saveRDS(paste0(output, "aalen_smk_m2_25_60.rds"))

interaction_model3(nhis_male, smoking4.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_smk_m_60_70.rds"))
jointeffect_model3(nhis_male, edu.smk, 60, 69.999)         %>% saveRDS(paste0(output, "aalen_smk_m2_60_70.rds"))

interaction_model3(nhis_male, smoking4.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_smk_m_70_85.rds"))
jointeffect_model3(nhis_male, edu.smk, 70, 84.999)         %>% saveRDS(paste0(output, "aalen_smk_m2_70_85.rds"))



# BMI, Women  
interaction_model3(nhis_female, bmi_cat.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_bmi_f_25_60.rds"))
jointeffect_model3(nhis_female, edu.bmi, 25, 59.999)        %>% saveRDS(paste0(output, "aalen_bmi_f2_25_60.rds"))

interaction_model3(nhis_female, bmi_cat.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_bmi_f_60_70.rds"))
jointeffect_model3(nhis_female, edu.bmi, 60, 69.999)        %>% saveRDS(paste0(output, "aalen_bmi_f2_60_70.rds"))

interaction_model3(nhis_female, bmi_cat.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_bmi_f_70_85.rds"))
jointeffect_model3(nhis_female, edu.bmi, 70, 84.999)        %>% saveRDS(paste0(output, "aalen_bmi_f2_70_85.rds"))

# BMI, Men  
interaction_model3(nhis_male, bmi_cat.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_bmi_m_25_60.rds"))
jointeffect_model3(nhis_male, edu.bmi, 25, 59.999)        %>% saveRDS(paste0(output, "aalen_bmi_m2_25_60.rds"))

interaction_model3(nhis_male, bmi_cat.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_bmi_m_60_70.rds"))
jointeffect_model3(nhis_male, edu.bmi, 60, 69.999)        %>% saveRDS(paste0(output, "aalen_bmi_m2_60_70.rds"))

interaction_model3(nhis_male, bmi_cat.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_bmi_m_70_85.rds"))
jointeffect_model3(nhis_male, edu.bmi, 70, 84.999)        %>% saveRDS(paste0(output, "aalen_bmi_m2_70_85.rds"))



# Physical Activity, Women  
interaction_model3(nhis_female, phy_act3.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_phy_f_25_60.rds"))
jointeffect_model3(nhis_female, edu.phy, 25, 59.999)         %>% saveRDS(paste0(output, "aalen_phy_f2_25_60.rds"))

interaction_model3(nhis_female, phy_act3.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_phy_f_60_70.rds"))
jointeffect_model3(nhis_female, edu.phy, 60, 69.999)         %>% saveRDS(paste0(output, "aalen_phy_f2_60_70.rds"))

interaction_model3(nhis_female, phy_act3.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_phy_f_70_85.rds"))
jointeffect_model3(nhis_female, edu.phy, 70, 84.999)         %>% saveRDS(paste0(output, "aalen_phy_f2_70_85.rds"))

# Physical Activity, Men  
interaction_model3(nhis_male, phy_act3.factor, 25, 59.999) %>% saveRDS(paste0(output, "aalen_phy_m_25_60.rds"))
jointeffect_model3(nhis_male, edu.phy, 25, 59.999)         %>% saveRDS(paste0(output, "aalen_phy_m2_25_60.rds"))

interaction_model3(nhis_male, phy_act3.factor, 60, 69.999) %>% saveRDS(paste0(output, "aalen_phy_m_60_70.rds"))
jointeffect_model3(nhis_male, edu.phy, 60, 69.999)         %>% saveRDS(paste0(output, "aalen_phy_m2_60_70.rds"))

interaction_model3(nhis_male, phy_act3.factor, 70, 84.999) %>% saveRDS(paste0(output, "aalen_phy_m_70_85.rds"))
jointeffect_model3(nhis_male, edu.phy, 70, 84.999)         %>% saveRDS(paste0(output, "aalen_phy_m2_70_85.rds"))





## Second: Load and view model results --------------------------------------------------------------------------------------------

# Alcohol, Women
aalen <-readRDS(paste0(output, "aalen_alc_f_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_alc_f2_25_60.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31)

      aalen <-readRDS(paste0(output, "aalen_alc_f_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_alc_f2_60_70.rds"))       
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31)   
      
      aalen <-readRDS(paste0(output, "aalen_alc_f_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_alc_f2_70_85.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31)   


 
# Smoking, Women
aalen <-readRDS(paste0(output, "aalen_smk_f_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_smk_f2_25_60.rds"))                           
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   

      aalen <-readRDS(paste0(output, "aalen_smk_f_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_smk_f2_60_70.rds"))   
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      
      aalen <-readRDS(paste0(output, "aalen_smk_f_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_smk_f2_70_85.rds"))    
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   


# BMI, Women
aalen <-readRDS(paste0(output, "aalen_bmi_f_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_bmi_f2_25_60.rds"))                              
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   

      aalen <-readRDS(paste0(output, "aalen_bmi_f_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_bmi_f2_60_70.rds"))                                
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   

      aalen <-readRDS(paste0(output, "aalen_bmi_f_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_bmi_f2_70_85.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      

# Physical Activity, Women
aalen <-readRDS(paste0(output, "aalen_phy_f_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_phy_f2_25_60.rds"))     
aalen_10000py(aalen, 1); aalen_10000py(aalen, 3);aalen_10000py(aalen2, 4);  aalen_10000py(aalen, 23)   
      
      aalen <-readRDS(paste0(output, "aalen_phy_f_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_phy_f2_60_70.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 3);aalen_10000py(aalen2, 4);  aalen_10000py(aalen, 23)   
      
      aalen <-readRDS(paste0(output, "aalen_phy_f_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_phy_f2_70_85.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 3);aalen_10000py(aalen2, 4);  aalen_10000py(aalen, 23)   




# Alcohol, Men
aalen <-readRDS(paste0(output, "aalen_alc_m_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_alc_m2_25_60.rds"))
aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31) # print results of interest
      
      aalen <-readRDS(paste0(output, "aalen_alc_m_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_alc_m2_60_70.rds"))       
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31)   
      
      aalen <-readRDS(paste0(output, "aalen_alc_m_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_alc_m2_70_85.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 6); aalen_10000py(aalen2, 13); aalen_10000py(aalen, 31)   
      
      
      
# Smoking, Men
aalen <-readRDS(paste0(output, "aalen_smk_m_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_smk_m2_25_60.rds"))                           
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      
      aalen <-readRDS(paste0(output, "aalen_smk_m_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_smk_m2_60_70.rds"))   
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      
      aalen <-readRDS(paste0(output, "aalen_smk_m_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_smk_m2_70_85.rds"))    
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      
      
# BMI, Men
aalen <-readRDS(paste0(output, "aalen_bmi_m_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_bmi_m2_25_60.rds"))                              
aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      
      aalen <-readRDS(paste0(output, "aalen_bmi_m_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_bmi_m2_60_70.rds"))                                
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      
      aalen <-readRDS(paste0(output, "aalen_bmi_m_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_bmi_m2_70_85.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 5); aalen_10000py(aalen2, 10); aalen_10000py(aalen, 28)   
      
      
# Physical Activity, Men
aalen <-readRDS(paste0(output, "aalen_phy_m_25_60.rds"))         
aalen2 <-readRDS(paste0(output, "aalen_phy_m2_25_60.rds"))     
aalen_10000py(aalen, 1); aalen_10000py(aalen, 3);aalen_10000py(aalen2, 4);  aalen_10000py(aalen, 23)   
      
      aalen <-readRDS(paste0(output, "aalen_phy_m_60_70.rds"))                                                      
      aalen2 <-readRDS(paste0(output, "aalen_phy_m2_60_70.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 3);aalen_10000py(aalen2, 4);  aalen_10000py(aalen, 23)   
      
      aalen <-readRDS(paste0(output, "aalen_phy_m_70_85.rds"))        
      aalen2 <-readRDS(paste0(output, "aalen_phy_m2_70_85.rds")) 
      aalen_10000py(aalen, 1); aalen_10000py(aalen, 3);aalen_10000py(aalen2, 4);  aalen_10000py(aalen, 23)   
      
      




      
# WOMEN - Causal Mediation Analysis Stratified by age (Objective 2) ----------------------------------------------------------------------------------

# Load data (used the same data as that generated from the main analysis)
expandedData <-readRDS(file.path(output, "expandedData_fem.rds"))

# AGES 25 - 59 Years - Run Causal Mediation Model
CMed_f_25_59 <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) + 
                                                              const(A.edu) * const(edu_M2.smk) +
                                                              const(A.edu) * const(edu_M3.bmi) +
                                                              const(A.edu) * const(edu_M4.phy) +
                                                              const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                            start.time=25, max.time=59.999,
                            data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
          saveRDS(CMed_f_25_59, file.path(output, "CMed_f_25_59.rds"))       # Save model results
          CMed_f_model_25_59 <-readRDS(file.path(output, "CMed_f_25_59.rds"))  # Load model results
          
          
          
          # Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
          summary(CMed_f_model_25_59)   
          getTE_NotRobust(CMed_f_model_25_59, c(1,3,5,7,9,29,33,37,41))  
          getIE_NotRobust(CMed_f_model_25_59, c(3,5,7,9,29,33,37,41))    
          getTE_IE_NotRobust(CMed_f_model_25_59, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41)) 



# AGES 60 - 69 Years - Run Causal Mediation Model
CMed_f_60_69 <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) + 
                                                              const(A.edu) * const(edu_M2.smk) +
                                                              const(A.edu) * const(edu_M3.bmi) +
                                                              const(A.edu) * const(edu_M4.phy) +
                                                              const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                        start.time=60, max.time=69.999,
                        data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
          saveRDS(CMed_f_60_69, file.path(output, "CMed_f_60_69.rds"))       # Save model results
          CMed_f_model_60_69 <-readRDS(file.path(output, "CMed_f_60_69.rds"))  # Load model results
          
          
          # Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
          summary(CMed_f_model_25_59)   
          getTE_NotRobust(CMed_f_model_25_59, c(1,3,5,7,9,29,33,37,41))  
          getIE_NotRobust(CMed_f_model_25_59, c(3,5,7,9,29,33,37,41))    
          getTE_IE_NotRobust(CMed_f_model_25_59, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41)) 




# AGES 70 - 85 Years - Run Causal Mediation Model
CMed_f_70_85 <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) + 
                                                              const(A.edu) * const(edu_M2.smk) +
                                                              const(A.edu) * const(edu_M3.bmi) +
                                                              const(A.edu) * const(edu_M4.phy) +
                                                              const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                        start.time=70, max.time=84.999,
                        data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
          saveRDS(CMed_f_70_85, file.path(output, "CMed_f_70_85.rds"))       # Save model results
          CMed_f_model_70_85 <-readRDS(file.path(output, "CMed_f_70_85.rds"))  # Load model results
          
          
          
          # Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
          summary(CMed_f_model_25_59)   
          getTE_NotRobust(CMed_f_model_25_59, c(1,3,5,7,9,29,33,37,41))  
          getIE_NotRobust(CMed_f_model_25_59, c(3,5,7,9,29,33,37,41))    
          getTE_IE_NotRobust(CMed_f_model_25_59, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41)) 






# MEN - Causal Mediation Analysis Stratified by age (Objective 2) ----------------------------------------------------------------------------------
          
# Load data (used the same data as that generated from the main analysis)
expandedData <-readRDS(file.path(output, "expandedData_male.rds"))

          
          
# AGES 25 - 59 Years - Run Causal Mediation Model   
CMed_m_25_59 <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) + 
                                                              const(A.edu) * const(edu_M2.smk) +
                                                              const(A.edu) * const(edu_M3.bmi) +
                                                              const(A.edu) * const(edu_M4.phy) +
                                                              const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                        start.time=25, max.time=59.999,
                        data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
          saveRDS(CMed_m_25_59, file.path(output, "CMed_m_25_59.rds"))       # Save model results
          CMed_m_model_25_59 <-readRDS(file.path(output, "CMed_m_25_59.rds"))  # Load model results
          
          
          
          # Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
          summary(CMed_m_model_25_59)  
          getTE_NotRobust(CMed_m_model_25_59, c(1,3,5,7,9,29,33,37,41))  
          getIE_NotRobust(CMed_m_model_25_59, c(3,5,7,9,29,33,37,41))    
          getTE_IE_NotRobust(CMed_m_model_25_59, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41)) 



          

# AGES 60 - 69 Years - Run Causal Mediation Model
CMed_m_60_69 <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) + 
                                                              const(A.edu) * const(edu_M2.smk) +
                                                              const(A.edu) * const(edu_M3.bmi) +
                                                              const(A.edu) * const(edu_M4.phy) +
                                                              const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                          start.time=60, max.time=69.999,
                          data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
          saveRDS(CMed_m_60_69, file.path(output, "CMed_m_60_69.rds"))       # Save model results
          CMed_m_model_60_69 <-readRDS(file.path(output, "CMed_m_60_69.rds"))  # Load model results
          
          
          # Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
          summary(CMed_m_model_25_59)  
          getTE_NotRobust(CMed_m_model_25_59, c(1,3,5,7,9,29,33,37,41))  
          getIE_NotRobust(CMed_m_model_25_59, c(3,5,7,9,29,33,37,41))    
          getTE_IE_NotRobust(CMed_m_model_25_59, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41)) 





# AGES 70 - 85 Years - Run Causal Mediation Model
CMed_m_70_85 <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(A.edu) * const(edu_M1.alc) + 
                                                              const(A.edu) * const(edu_M2.smk) +
                                                              const(A.edu) * const(edu_M3.bmi) +
                                                              const(A.edu) * const(edu_M4.phy) +
                                                              const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                        start.time=70, max.time=84.999,
                        data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
          saveRDS(CMed_m_70_85, file.path(output, "CMed_m_70_85.rds"))       # Save model results
          CMed_m_model_70_85 <-readRDS(file.path(output, "CMed_m_70_85.rds"))  # Load model results
          
          
          
          # Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
          summary(CMed_m_model_25_59)  
          getTE_NotRobust(CMed_m_model_25_59, c(1,3,5,7,9,29,33,37,41))  
          getIE_NotRobust(CMed_m_model_25_59, c(3,5,7,9,29,33,37,41))    
          getTE_IE_NotRobust(CMed_m_model_25_59, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41)) 




