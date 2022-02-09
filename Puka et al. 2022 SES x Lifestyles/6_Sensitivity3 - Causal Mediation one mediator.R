# SES x Lifestyle Differential Vulnerability & Exposure Project
# Sensitivity Analyses 3: Causal Mediation using one mediator at a time 


# LOAD DATA AND SET FILE LOCATIONS

# load libraries
library(tidyverse)  # data management
library(timereg)    # additive survival models
library(VGAM)       # multinomial regression, needed for causal mediation
library(MASS)       # needed for causal mediation functions


#Personal computer
data    <- "C:/.../nhis/Data/"
output  <- "C:/.../nhis/SES x lifestyle/Output/"
source("Function - Format Results.R")
source("Function - CausalMed Results.R")
    

# Load data
nhis        <- readRDS (file.path(data, "nhis.rds"))
nhis_male   <- readRDS (file.path(data, "nhis_male.rds"))
nhis_female <- readRDS (file.path(data, "nhis_female.rds"))



# Create function to Run Causal Mediation with a single mediator -----------------------------------------------------------------------

causal_med <- function(data, mediator, cov1, cov2, cov3, M_ref_level, loc_dir, loc_indir, loc_med_int){

    ### Step 0: Select data to use
    mydata <- data %>%
      mutate (Mediator = {{mediator}},
              Cov1 = {{cov1}}, 
              Cov2 = {{cov2}},
              Cov3 = {{cov3}}) %>%
      dplyr::select(edu, Mediator, Cov1, Cov2, Cov3, allcause_death, bl_age, end_age, married, ethnicity, srvy_yr)
    
      # specifies the reference category
      mydata$edu <- factor(mydata$edu, levels=c(1,2,3), labels = c("Low", "Med", "High"))
      mydata$edu <- relevel(mydata$edu, ref = "High")
      
      # NOTE: For technical reasons, the mediators should be coded as integers starting with 1


    ### Step 1: Fit a model for the mediator, conditioning on exposure (education) and all confounders
    mydata$ATemp <- mydata$edu # first, create and use a copy of the exposure variable (for technical reasons related to R)
    fit_Med <- vglm(Mediator ~ ATemp + bl_age + married + factor(ethnicity) + factor(srvy_yr) + factor(Cov1) + factor(Cov2) + factor(Cov3), data = mydata, family=multinomial(refLevel = M_ref_level))


    
    ### Step 2: Construct copies of ID and exposure
    mydata$ID <- 1:nrow(mydata) # construct id variable
    
    # Create counterfactual version of exposure (education)
    levelsOfEDU <- unique(mydata$edu)
    myData1 <- mydata
    myData2 <- mydata
    myData3 <- mydata
    myData1$edu_M <- levelsOfEDU[1]
    myData2$edu_M <- levelsOfEDU[2]
    myData3$edu_M <- levelsOfEDU[3]
    newMyData <- rbind(myData1, myData2, myData3)



    ### Step 3: Construct weights
    newMyData$ATemp <- newMyData$edu
    tempDir1 <- as.matrix(predict(fit_Med,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$Mediator)]
    
    newMyData$ATemp <- newMyData$edu_M
    tempIndir1 <- as.matrix(predict(fit_Med,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$Mediator)]
    
    newMyData$weightM <- tempIndir1/tempDir1
    
    

    #### Step 4: Run Causal Mediation Analysis
    CMed_single <- aalen(Surv(bl_age, end_age, allcause_death) ~ const(edu) * const(edu_M) + 
                                                              const(factor(Cov1)) +  const(factor(Cov2)) + const(factor(Cov3)) +  
                                                              const(married) + factor(ethnicity) + const(factor(srvy_yr)),
                                data=newMyData, weights=newMyData$weightM, clusters=newMyData$ID, robust=0)  
                      

    # The 'summary(model)' command will produce the direct effect, indirect effects, and mediated interaction effects. 
    # Functions were used to extract the other details. 
                    
    # Get final results. 
    # NOTE: The numbers below pertain to the location of the direct effect, indirect effect, and mediated interaction effect
    
    summary <- summary(CMed_single)
    TE      <- getTE_NotRobust(CMed_single, c(loc_dir,loc_indir,loc_med_int))           # Function to get the total effect and proportion mediated
    IE      <- getIE_NotRobust(CMed_single, c(loc_indir,loc_med_int))        # Function to get the total combined indirect effect
    TE_IE   <- getTE_IE_NotRobust(CMed_single, c(loc_dir,loc_indir,loc_med_int), c(loc_indir,loc_med_int))  # Function to get the proportion mediated of the combined indirect effect

    results <- list(summary, TE, IE, TE_IE)
    
    return (results)
    
}



# Run Causal Mediation analyses with single mediator -------------------------------------------------------------------------------------------
## Alcohol use
causal_med(nhis_male, alcohol5v2, smoking4, bmi_cat, phy_act3, 3, 1, 3, 31)
causal_med(nhis_female, alcohol5v2, smoking4, bmi_cat, phy_act3, 3, 1, 3, 31)

## Smoking
causal_med(nhis_male, smoking4, alcohol5v2, bmi_cat, phy_act3, 1, 1, 3, 32)
causal_med(nhis_female, smoking4, alcohol5v2, bmi_cat, phy_act3, 1, 1, 3, 32)

## BMI
causal_med(nhis_male, bmi_cat, alcohol5v2, smoking4, phy_act3, 2, 1, 3, 32)
causal_med(nhis_female, bmi_cat, alcohol5v2, smoking4, phy_act3, 2, 1, 3, 32)

## Physical activity
causal_med(nhis_male, phy_act3, alcohol5v2, smoking4, bmi_cat, 3, 1, 3, 33)
causal_med(nhis_female, phy_act3, alcohol5v2, smoking4, bmi_cat, 3, 1, 3, 33)
