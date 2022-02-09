
# NHIS Data Management File

library(haven)      # Read SAS file
library(tidyverse)  # data management
library(janitor)    # clean variable names
library(skimr)      # descriptive statistics



# Specify the data file location
data_orig <- "C:/.../Original data/"
data_new  <- "C:/.../Processed data/"

# Import data form SAS and edit/re-categorize variables 

nhis_all <- read_sas (paste0(data_orig, "nhis_mort_clean.sas7bdat")) %>%
  zap_formats() %>% zap_label() %>% clean_names() %>%   # removes labels/formats form SAS and clean names
    mutate (# Recode and create variables
          bl_age = age,                  # Baseline age
          end_age = age + yrs_followup,  # Age at death/censor
          income4 = recode(income, `1`=1, `2`=2, `3`=2, `4`=3, `5`=4),              # merge the two middle categories
          alcohol4 = recode(alcohol5, `1`=1, `2`=2, `3`=3, `4`=4, `5`=4),           # merge high risk and very high risk group
          alcohol5v2 = recode(alcohol6, `1`=1, `2`=2, `3`=3, `4`=4, `5`=5, `6`=5),  # separate former/never drinkers and merge high risk and very high risk group
          smoking4 = recode(smoking, `0`=1, `1`=2, `2`=3, `3`=4),                   # recoded such that '1' is the first/smallest category
          allcause_death = mortstat,
      
          # LABEL VARIABLES 
          # Exposures
          edu.factor = factor(edu, levels=c(3,1,2), labels = c("Bachelors", "Highschool", "Some college")),
          income.factor = factor(income, levels=c(4,1,2,3,5), labels = c("High income", "Poor","Near poor","Middle income", "Missing")),
          income4.factor = factor(income4, levels=c(3,1,2,4), labels = c( "High", "Low","Medium", "Missing")),
  
          # Mediator 1 - Alcohol
          alcohol5.factor   = factor(alcohol5, levels=c(1,2,3,4,5),   labels = c("Abstinence", "Low risk","Medium risk","High risk","Very high risk")),
          alcohol6.factor   = factor(alcohol6, levels=c(1,2,3,4,5,6), labels = c("Never Drinker", "Former Drinker", "Low risk","Medium risk","High risk","Very high risk")),
          alcohol4.factor   = factor(alcohol4, levels=c(1,2,3,4),     labels = c("Abstinence", "Low risk","Medium risk","High risk")),
          alcohol5v2.factor = factor(alcohol5v2, levels=c(3,1,2,4,5), labels = c("Low risk","Never Drinker", "Former Drinker", "Medium risk","High risk")), 
          drink_hist.factor = factor(drink_hist, levels=c(0,1,2),     labels = c("Never Drinker", "Former Drinker", "Current Drinker")),
          hed.factor        = factor(hed, levels=c(1,2,3,4),          labels = c("No HED", "HED <1/month", "HED >1/month, <1/week", "HED >=1/week")),
          
          # Mediator 2 - BMI
          bmi_cat.factor = factor(bmi_cat, levels=c(2,1,3,4), labels = c("Healthy weight","Underweight", "Overweight", "Obese")),
                      
          # Mediator 3 - Smoking
          smoking4.factor = factor(smoking4, levels=c(1,2,3,4), labels = c("Never smoker", "Former smoker", "Current some day smoker", "Current everyday smoker")),
                
          # Mediator 4 - Physicial Activity
          phy_act3.factor =factor(phy_act3,levels=c(3,1,2), labels = c("Active", "Sedentary", "Somewhat active")),
  
          # Covariates
          ethnicity.factor = factor(ethnicity, levels=c(1,2,3,4), labels = c("Non-Hispanic White", "Non-Hispanic Black","Hispanic", "Other")),
          ethnicity_detail = factor(ethnicity_detail, levels=c(1,2,3,4,5,6), labels = c("Non-Hispanic White", "Non-Hispanic Black","Hispanic", "AI/AN", "API", "Other, including multiple")),
          female.factor = factor(female, levels=c(0,1), labels = c("Male", "Female")),
          married.factor = factor(married, levels=c(0,1), labels = c("Not married/living togeter", "Married/cohabitating")),
          employed.factor = factor(employed, levels=c(0,1), labels = c("Not employed", "Paid employment, student or retired")),
          diabet.factor = factor(diabet, levels=c(0,1,2), labels = c("No","Borderline","Yes")),
          
          # Outcome        
          allcause_death.factor = factor(allcause_death, levels=c(0,1), labels = c("Alive","Deceased")))
                
         str(nhis_all)
        # Check recoding                                     
        # count(nhis_all, alcohol5, alcohol4)
        # count(nhis_all, alcohol6, alcohol5_2)
        
        
# Create an 'interaction' variable, combining the SES and Health behavior variables
# For main analyses
nhis_all$edu.alc <- interaction(nhis_all$edu.factor, nhis_all$alcohol5v2.factor)
nhis_all$edu.smk <- interaction(nhis_all$edu.factor, nhis_all$smoking4.factor)
nhis_all$edu.bmi <- interaction(nhis_all$edu.factor, nhis_all$bmi_cat.factor)
nhis_all$edu.phy <- interaction(nhis_all$edu.factor, nhis_all$phy_act3.factor)

# For sensitivity analyses
nhis_all$inc.alc <- interaction(nhis_all$income4.factor, nhis_all$alcohol5v2.factor)
nhis_all$inc.smk <- interaction(nhis_all$income4.factor, nhis_all$smoking4.factor)
nhis_all$inc.bmi <- interaction(nhis_all$income4.factor, nhis_all$bmi_cat.factor)
nhis_all$inc.phy <- interaction(nhis_all$income4.factor, nhis_all$phy_act3.factor)

nhis_all$edu.hed <- interaction(nhis_all$edu.factor, nhis_all$hed.factor)
        
        
      
       
# Create subset of data with relevant participants        
# Remove those outside our age range
nhis_all2 <- filter (nhis_all, age>=25 & age <85) %>% 
  mutate(lost = factor(ifelse(complete.cases(yrs_followup, allcause_death, alcohol5v2, bmi_cat, smoking4, phy_act3, edu, age, female, married, ethnicity), "Complete", "Lost")))

nhis <- filter(nhis_all2, lost=="Complete")

# Create database specific to males or females
nhis_female <- filter(nhis, female==1)
nhis_male <- filter(nhis, female==0)
        

# Descriptives of how many were lost to follow-up
nrow(nhis_all)                        # n completed NHIS
nrow(nhis_all) - nrow(nhis_all2)      # Aged <25 or >85
nrow(nhis_all2)                       # n eligible 
count(nhis_all, allcause_death)       # No mortality data
nrow(nhis_all2) - nrow(nhis)          # missing any covariates
nrow(nhis_all2) - nrow(nhis) - 23474  # of those with mortality data, n missing covariates
nrow(nhis)                            # final sample size



# Save copy of final datasets  
saveRDS(nhis_all2,    paste0(data_new, "nhis_all.rds"))      # NHIS data with all participants
saveRDS(nhis,        paste0(data_new, "nhis.rds"))          # NHIS data to be analyzed
saveRDS(nhis_male,   paste0(data_new, "nhis_male.rds"))     # NHIS data to be analyzed (males only)
saveRDS(nhis_female, paste0(data_new, "nhis_female.rds"))   # NHIS data to be analyzed (females only)

