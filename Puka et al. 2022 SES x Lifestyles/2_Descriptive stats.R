
# SES x Lifestyle Differential Vulnerability & Exposure Project
# Descriptive Statistics

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
library(biostat3)   # survRate command


# Specify the data and output file locations
data    <- "C:/.../Processed data/"   # Location of data
output  <- "C:/.../Output/"           # Location of figures/tables


# Load data
nhis_all    <- readRDS (paste0(data, "nhis_all.rds"))
nhis        <- readRDS (paste0(data, "nhis.rds"))
nhis_male   <- readRDS (paste0(data, "nhis_male.rds"))
nhis_female <- readRDS (paste0(data, "nhis_female.rds"))




# DESCRIPTIVE STATISTICS 

# Table 1: Participant characteristics - STRATIFIED BY SEX
tab1 <-CreateTableOne(vars= c("age", "yrs_followup","allcause_death.factor", 
                              "alcohol5v2.factor","smoking4.factor",  "bmi_cat.factor", "phy_act3.factor",
                               "ethnicity.factor", "income4.factor",  "married.factor"),
                      factorVars = c("allcause_death.factor",
                                      "alcohol5v2.factor", "smoking4.factor", "bmi_cat.factor", "phy_act3.factor",
                                     "ethnicity.factor", "income4.factor",  "married.factor"), 
                      strata= c("edu.factor", "female.factor"), addOverall = TRUE, data=nhis)
  table1_v1 <- print(tab1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE)  # Shows sample size and %
  table1_v2 <- print(tab1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE, format="p") # shows % only
  write.csv(table1_v1, file = file.path(output, "Table1 Demographics_V1.csv"))
  write.csv(table1_v2, file = file.path(output, "Table1 Demographics_V2.csv"))
  kableone(table1_v1)

  
# Person years and death rate: tstop = person years; event = # events; rate = events per person year 
survRate(Surv(yrs_followup, allcause_death)~1, data=nhis)          # overall 
survRate(Surv(yrs_followup, allcause_death)~edu, data=nhis)        # for each SES category 
survRate(Surv(yrs_followup, allcause_death)~female+edu, data=nhis) # for each SES * Sex category 
 

# Composition of "Other" race/ethnicity group
nhis %>% filter(ethnicity==4) %>%
  count(ethnicity_detail) %>%
  mutate(percent = n/sum(n))


# eFigure 1: Survival plot 
ggsurvplot_facet(fit = survfit(Surv(bl_age, end_age, allcause_death) ~ edu, data = nhis), 
  data=nhis, facet.by="female.factor", censor = FALSE,xlim = c(25, 100), 
  conf.int = TRUE, 
  legend.labs = c("Low education", "Medium education", "High education"),
  xlab = "Age (years)", 
  ylab = "Overall survival probability") 

      # Age Medium Survival:
      survfit(Surv(bl_age, end_age, allcause_death) ~ edu, data = nhis)
      survfit(Surv(bl_age, end_age, allcause_death) ~ edu, data = nhis_female)
      survfit(Surv(bl_age, end_age, allcause_death) ~ edu, data = nhis_male)


      

# eTable 1: Participant characteristics - STRATIFIED BY SEX
nhis_all_male <- filter(nhis_all, female==0)
nhis_all_female <- filter(nhis_all, female==1)
      
all_vars <- c("age", "edu.factor", "alcohol5v2.factor","smoking4.factor",  "bmi_cat.factor", "phy_act3.factor", "ethnicity.factor",  "married.factor")
factor_vars <- c("edu.factor", "alcohol5v2.factor", "smoking4.factor", "bmi_cat.factor", "phy_act3.factor", "ethnicity.factor",  "married.factor")

tab_e1 <-CreateTableOne(all_vars, factor_vars, strata="lost", data=nhis_all_male)
    table_e1 <- print(tab_e1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE, smd=TRUE, format="p") # shows % only
    write.csv(table_e1, file = file.path(output, "eTable 1 (men) Attrition.csv"))
    kableone(table_e1)
      
    
tab_e1 <-CreateTableOne(all_vars, factor_vars, strata="lost", data=nhis_all_female)
    table_e1 <- print(tab_e1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE, smd=TRUE, format="p") # shows % only
    write.csv(table_e1, file = file.path(output, "eTable 1 (women) Attrition.csv"))
    kableone(table_e1)
    
      