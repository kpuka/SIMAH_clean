
# SES x Lifestyle Differential Vulnerability & Exposure Project
# Objective 2: Causal Mediation File 


# LOAD DATA AND SET FILE LOCATIONS

# load libraries
library(tidyverse)  # data management
library(timereg)    # additive survival models
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


# Create function for data preparation, to repeat analyses for men and women ------------------------------------------
causal_mediation_prep <-function(data) {

      # Data Preparation -----------------------------------------------------------------------------------------------------------------------------
      
      ### Step 0: Select data to use *****************************************************************************************************************
      # **********************************************************************************************************************************************
      mydata <- data %>%
        mutate(A.edu = edu,
          M1.alc = alcohol5v2,
          M2.smk = smoking4,
          M3.bmi = bmi_cat,
          M4.phy = phy_act3) %>%
        dplyr::select(A.edu, M1.alc, M2.smk, M3.bmi, M4.phy, allcause_death, bl_age, end_age, married, ethnicity, srvy_yr)
      
      # specifies the reference category
      mydata$A.edu <- factor(mydata$A.edu, levels=c(1,2,3), labels = c("Low", "Med", "High"))
      mydata$A.edu <- relevel(mydata$A.edu, ref = "High")
      
      # NOTE: For technical reasons, the mediators should be coded as integers starting with 1
      
      
      
      ### Step 1: Fit a model for each mediator, , conditioning on exposure and all confounders *******************************************************
      # ***********************************************************************************************************************************************
      
      # Fit model for each mediator, conditioning on exposure (education) and all confounders
      
      mydata$ATemp <- mydata$A.edu # first, create and use a copy of the exposure variable (for technical reasons related to R)
      fitM1 <- vglm(M1.alc ~ ATemp + bl_age + married + factor(ethnicity) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 3))
      fitM2 <- vglm(M2.smk ~ ATemp + bl_age + married + factor(ethnicity) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 1))
      fitM3 <- vglm(M3.bmi ~ ATemp + bl_age + married + factor(ethnicity) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 2))
      fitM4 <- vglm(M4.phy ~ ATemp + bl_age + married + factor(ethnicity) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 3))
      
      
      
      
      ### Step 2: Construct copies of ID and exposure *************************************************************************************************
      # ***********************************************************************************************************************************************
      
      #Create ID Variable
      mydata$ID <- 1:nrow(mydata) # construct id variable
      
      # Create counterfactual version of exposure (education); repeated 4 times because there are 4 mediators
      levelsOfEDU <- unique(mydata$A.edu)
      myData1 <- mydata
      myData2 <- mydata
      myData3 <- mydata
      myData1$edu_M1.alc <- levelsOfEDU[1]
      myData2$edu_M1.alc <- levelsOfEDU[2]
      myData3$edu_M1.alc <- levelsOfEDU[3]
      tempMyData <- rbind(myData1, myData2, myData3)
      
      myData1 <- tempMyData
      myData2 <- tempMyData
      myData3 <- tempMyData
      myData1$edu_M2.smk <- levelsOfEDU[1]
      myData2$edu_M2.smk <- levelsOfEDU[2]
      myData3$edu_M2.smk <- levelsOfEDU[3]
      tempMyData <- rbind(myData1, myData2, myData3)
      
      myData1 <- tempMyData
      myData2 <- tempMyData
      myData3 <- tempMyData
      myData1$edu_M3.bmi <- levelsOfEDU[1]
      myData2$edu_M3.bmi <- levelsOfEDU[2]
      myData3$edu_M3.bmi <- levelsOfEDU[3]
      tempMyData <- rbind(myData1, myData2, myData3)
      
      myData1 <- tempMyData
      myData2 <- tempMyData
      myData3 <- tempMyData
      myData1$edu_M4.phy <- levelsOfEDU[1]
      myData2$edu_M4.phy <- levelsOfEDU[2]
      myData3$edu_M4.phy <- levelsOfEDU[3]
      newMyData <- rbind(myData1, myData2, myData3)
      
      
      
      
      ### Step 3: Construct weights  *********************************************************************************************************************
      # **************************************************************************************************************************************************
      
      # M1: alcohol
      newMyData$ATemp <- newMyData$A.edu
      tempDir1 <- as.matrix(predict(fitM1,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M1.alc)]
      
      newMyData$ATemp <- newMyData$edu_M1.alc
      tempIndir1 <- as.matrix(predict(fitM1,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M1.alc)]
      
      newMyData$weight1 <- tempIndir1/tempDir1
      
      
      #M2: Smoking
      newMyData$ATemp <- newMyData$A.edu
      tempDir2 <- as.matrix(predict(fitM2,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M2.smk)]
      
      newMyData$ATemp <- newMyData$edu_M2.smk
      tempIndir2 <- as.matrix(predict(fitM2,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M2.smk)]
      
      newMyData$weight2 <- tempIndir2/tempDir2
      
      
      #M3: BMI
      newMyData$ATemp <- newMyData$A.edu
      tempDir3 <- as.matrix(predict(fitM3,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M3.bmi)]
      
      newMyData$ATemp <- newMyData$edu_M3.bmi
      tempIndir3 <- as.matrix(predict(fitM3,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M3.bmi)]
      
      newMyData$weight3 <- tempIndir3/tempDir3
      
      
      #M4: Physical activity
      newMyData$ATemp <- newMyData$A.edu
      tempDir4 <- as.matrix(predict(fitM4,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M4.phy)]
      
      newMyData$ATemp <- newMyData$edu_M4.phy
      tempIndir4 <- as.matrix(predict(fitM4,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M4.phy)]
      
      newMyData$weight4 <- tempIndir4/tempDir4
      
      # Final weight
      newMyData$weightM <- newMyData$weight1 * newMyData$weight2 * newMyData$weight3 * newMyData$weight4
      
      return(newMyData)
}
    

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
                
# Save model results
saveRDS(CMed_f, file.path(output, "CMed_f.rds"))       


# Load model results
CMed_model <-readRDS(file.path(output, "CMed_f.rds"))  

# The 'summary(model)' command will produce the direct effect, indirect effects, and mediated interaction effects. 
# Functions were used to extract the other details. 
              
# Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
summary(CMed_model)   
getTE_NotRobust(CMed_model, c(1,3,5,7,9,29,33,37,41))                             # Function to get the total effect and proportion mediated
getIE_NotRobust(CMed_model, c(3,5,7,9,29,33,37,41))                               # Function to get the total combined indirect effect  
getTE_IE_NotRobust(CMed_model, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41))  # Function to get the proportion mediated of the combined indirect effect
      

 
     

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


# Save model results
saveRDS(CMed_m, file.path(output, "CMed_m.rds"))       

# Load model results
CMed_model <-readRDS(file.path(output, "CMed_m.rds"))  


# The 'summary(model)' command will produce the direct effect, indirect effects, and mediated interaction effects. 
# Functions were used to extract the other details. 

# Get final results. NOTE: THE NUMBERS BELOW MAY HAVE TO BE CHANGED IF A DIFFERENT MODEL IS USED
summary(CMed_model) 
getTE_NotRobust(CMed_model, c(1,3,5,7,9,29,33,37,41))                              # Function to get the total effect and proportion mediated
getIE_NotRobust(CMed_model, c(3,5,7,9,29,33,37,41))                                # Function to get the total combined indirect effect  
getTE_IE_NotRobust(CMed_model, c(1,3,5,7,9,29,33,37,41), c(3,5,7,9,29,33,37,41))   # Function to get the proportion mediated of the combined indirect effect


