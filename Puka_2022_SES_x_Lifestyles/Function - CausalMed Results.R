
# SES x Lifestyle Differential Vulnerability & Exposure Project
## Functions for Causal Mediation 


# Data preparation for Causal Mediation -------------------------------------------------------------
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


# Extract Data from Causal Mediation ----------------------------------------------------------------

# Direct, indirect and mediated interactive effects and standard errors are derived directly from the summary() command below
# Total effect is obtained by the sum of the three separate effects
# Confidence intervals for total effects and mediated proportions are computed using the code below:

# Function to get the (Not-Robust) total effect and proportion mediated
getTE_NotRobust <- function(CMed_model, v){
  TE <- sum(CMed_model$gamma[v])
  mu <- CMed_model$gamma[v]
  Omega <- CMed_model$var.gamma[v,v]    # To obtain non-robust estimates
  temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
  temp_TE <- apply(temp,1,sum)
  med_prop <- c(mu/TE,1)
  med_prop_CI <- rbind(t(apply(temp/temp_TE, 2, quantile, c(0.025, 0.975))), c(1,1))
  output <- cbind(c(mu,TE), c(apply(temp,2,sd),sd(temp_TE)), med_prop, med_prop_CI)
  colnames(output) <- c("Est.", "SE", "med_prop", "lowerCI", "UpperCI")
  rownames(output) <- c(rownames(CMed_model$gamma)[v],"TE")
  return(output)}

        # # Function to get the (Robust) total effect and proportion mediated
        getTE_Robust <- function(CMed_model, v){
          TE <- sum(CMed_model$gamma[v])
          mu <- CMed_model$gamma[v]
          Omega <- CMed_model$robvar.gamma[v,v] # To obtain robust estimates
          temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
          temp_TE <- apply(temp,1,sum)
          med_prop <- c(mu/TE,1)
          med_prop_CI <- rbind(t(apply(temp/temp_TE, 2, quantile, c(0.025, 0.975))), c(1,1))
          output <- cbind(c(mu,TE), c(apply(temp,2,sd),sd(temp_TE)), med_prop, med_prop_CI)
          colnames(output) <- c("Est.", "SE", "med_prop", "lowerCI", "UpperCI")
          rownames(output) <- c(rownames(CMed_model$gamma)[v],"TE")
          return(output)}


# Function to get the (Not-Robust) total combined indirect effect  
getIE_NotRobust <- function(CMed_model, v){
  IE <- sum(CMed_model$gamma[v])
  mu <- CMed_model$gamma[v]
  Omega <- CMed_model$var.gamma[v,v] # To obtain non-robust estimates
  require(MASS)
  temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
  temp_IE <- apply(temp,1,sum)
  med_prop <- c(mu/IE,1)
  med_prop_CI <- rbind(t(apply(temp/temp_IE, 2, quantile, c(0.025, 0.975))), c(1,1))
  output <- cbind(c(mu,IE), c(apply(temp,2,sd),sd(temp_IE)), med_prop, med_prop_CI)
  colnames(output) <- c("Est.", "SE", "med_prop", "lowerCI", "UpperCI")
  rownames(output) <- c(rownames(CMed_model$gamma)[v],"IE")
  return(output)}

        # Function to get the (Robust) total combined indirect effect  
        getIE_Robust <- function(CMed_model, v){
          IE <- sum(CMed_model$gamma[v])
          mu <- CMed_model$gamma[v]
          Omega <- CMed_model$robvar.gamma[v,v] # To obtain robust estimates
          require(MASS)
          temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
          temp_IE <- apply(temp,1,sum)
          med_prop <- c(mu/IE,1)
          med_prop_CI <- rbind(t(apply(temp/temp_IE, 2, quantile, c(0.025, 0.975))), c(1,1))
          output <- cbind(c(mu,IE), c(apply(temp,2,sd),sd(temp_IE)), med_prop, med_prop_CI)
          colnames(output) <- c("Est.", "SE", "med_prop", "lowerCI", "UpperCI")
          rownames(output) <- c(rownames(CMed_model$gamma)[v],"IE")
          return(output)}



# Function to get the (Not Robust)  proportion mediated of the combined indirect effect
getTE_IE_NotRobust <- function(CMed_model, v, z){
  #total effect
  TE <- sum(CMed_model$gamma[v])
  mu <- CMed_model$gamma[v]
  # Omega <- CMed_model$robvar.gamma[v,v]  # To obtain robust estimates
  Omega <- CMed_model$var.gamma[v,v]     # To obtain non-robust estimates
  require(MASS)
  temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
  temp_TE <- apply(temp,1,sum)
  IE <- sum(CMed_model$gamma[z])
  muIE <- CMed_model$gamma[z]
  #OmegaIE <- CMed_model$robvar.gamma[z,z] # To obtain robust estimates
  OmegaIE <- CMed_model$var.gamma[z,z]    # To obtain non-robust estimates
  require(MASS)
  tempIE <- mvrnorm(n=10^4, mu=muIE, Sigma=OmegaIE)
  temp_IE <- apply(tempIE,1,sum)
  med_prop <- c(IE/TE,1)
  med_prop_CI <- (temp_IE/temp_TE)
  output <- cbind(IE, med_prop, quantile)
  quantile <- quantile(med_prop_CI, c(0.025, 0.975))
  output <- cbind(IE, med_prop, quantile)
  return(output)}

      
      # Function to get the (Not Robust)  proportion mediated of the combined indirect effect
      getTE_IE_Robust <- function(CMed_model, v, z){
        #total effect
        TE <- sum(CMed_model$gamma[v])
        mu <- CMed_model$gamma[v]
        Omega <- CMed_model$robvar.gamma[v,v]  # To obtain robust estimates
        require(MASS)
        temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
        temp_TE <- apply(temp,1,sum)
        IE <- sum(CMed_model$gamma[z])
        muIE <- CMed_model$gamma[z]
        OmegaIE <- CMed_model$robvar.gamma[z,z] # To obtain robust estimates
        require(MASS)
        tempIE <- mvrnorm(n=10^4, mu=muIE, Sigma=OmegaIE)
        temp_IE <- apply(tempIE,1,sum)
        med_prop <- c(IE/TE,1)
        med_prop_CI <- (temp_IE/temp_TE)
        output <- cbind(IE, med_prop, quantile)
        quantile <- quantile(med_prop_CI, c(0.025, 0.975))
        output <- cbind(IE, med_prop, quantile)
        return(output)}

      
      
# Format Data from Causal Mediation ----------------------------------------------------------------

format_CMed <- function (model, coef_list) {
  
  group <- enexpr(coef_list)
  
  # All_coef model coefficients   
  All_coef <- coef(model) %>%
    as.data.frame() %>% 
    rownames_to_column(var = "variable") %>% 
    mutate (deaths_10000py = round(Coef. * 10000, 1),
      lower = round(`lower2.5%` * 10000, 1),
      upper = round(`upper97.5%` * 10000, 1),
      deaths_10000py_CI = paste0(deaths_10000py, " (", lower, ", ", upper, ")")) %>% 
    dplyr::select (variable, deaths_10000py_CI)
  
  coef <- slice(All_coef, coef_list)
  
  
  # Function to get the total effect and proportion mediated
  TE_prop <- getTE_NotRobust(model, coef_list) %>% 
    as.data.frame() %>% rownames_to_column(var = "variable") 
  
  TE <- filter(TE_prop, variable =="TE") %>% 
    mutate (lower = round((`Est.` - (1.96 * SE))*10000,1), 
      upper = round((`Est.` + (1.96 * SE))*10000,1), 
      deaths_10000py = round(`Est.` * 10000,1),
      deaths_10000py_CI = paste0(deaths_10000py, " (", lower, ", ", upper, ")")) %>%
    dplyr::select(variable, deaths_10000py_CI)
  
  prop <- filter(TE_prop, variable !="TE") %>%
    mutate(prop = round(med_prop * 100,0),
      lower = round(lowerCI * 100,0),
      upper = round(UpperCI * 100,0),
      prop_CI = paste0(prop, " (", lower, ", ", upper, ")")) %>% 
    dplyr::select(variable, prop_CI) 
  
  
  # Function to get the total combined indirect effect  (coef_list excluding 1st item)
  IE <- getIE_NotRobust(model, coef_list[-1]) %>% 
    as.data.frame() %>% rownames_to_column(var = "variable") %>% 
    filter (variable == "IE") %>% 
    mutate (lower = round((`Est.` - (1.96 * SE))*10000,1), 
      upper = round((`Est.` + (1.96 * SE))*10000,1), 
      deaths_10000py = round(`Est.` * 10000,1),
      deaths_10000py_CI = paste0(deaths_10000py, " (", lower, ", ", upper, ")")) %>%
    dplyr::select(variable, deaths_10000py_CI) 
  
  
  # Function to get the proportion mediated of the combined indirect effect
  IE_prop <- getTE_IE_NotRobust(model, coef_list, coef_list[-1]) %>% 
    as.data.frame() %>% rownames_to_column(var = "variable") %>%
    pivot_wider(names_from="variable", values_from=c("IE", "med_prop", "quantile")) %>%
    mutate (IE_prop = round(`med_prop_2.5%` * 100, 0),
      lower = round(`quantile_2.5%` *100, 0), 
      upper = round(`quantile_97.5%` *100, 0),
      prop_CI = paste0(IE_prop, " (", lower, ", ", upper, ")"),
      variable = "IE") %>% 
    dplyr::select(variable, prop_CI) 
  
  
  one <- full_join(coef, prop, by="variable") %>%
    mutate(label = c(paste0("02 Direct effect of ", group),
      "04 Alcohol use: differential exposure",
      "06 Smoking: differential exposure",
      "08 BMI: differential exposure",
      "10 Physical activity: differential exposure", 
      "05 Alcohol use: differential vulnerability ",
      "07 Smoking: differential vulnerability ",
      "09 BMI: differential vulnerability ",
      "11 Physical activity: differential vulnerability ")) %>% 
    dplyr::select(label, deaths_10000py_CI, prop_CI) 
  
  two <- TE %>%
    mutate (label = ifelse(variable == "TE", paste0("01 Total effect of ", group), NA),
      prop_CI = "100") %>% 
    dplyr::select(label, deaths_10000py_CI, prop_CI)
  
  three <- full_join(IE, IE_prop, by="variable") %>%
    mutate (label = ifelse(variable == "IE", paste0("03 Indirect effect of ", group), NA)) %>% 
    dplyr::select(label, deaths_10000py_CI, prop_CI) 
  
  
  final <- rbind(one, two, three) %>%
    arrange(label) 
  
  return(final)
}
