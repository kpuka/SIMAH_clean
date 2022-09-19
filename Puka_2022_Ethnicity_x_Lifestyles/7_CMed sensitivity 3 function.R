
# Race x Lifestyle Differential Vulnerability & Exposure Project
## Functions for Causal Mediation Sensitivity 3
# Analysis for separate age subgroups

# Main functions ------------------------------------------------------------------------------------------------------------------
# Bootstrap Causal Mediation Analysis
bootstrap_CMed_ageStrat <- function(data, reps, prop, start_time, end_time) {
  
  # Set up loop to re-run analyses multiple times
  foreach(i = 1:reps, .combine="c", 
          # need to specify the packages and functions that will be used
          .packages = c("tidyverse", "timereg", "VGAM", "MASS"), 
          .export = c("CMed_prep", "getTE", "getIE", "getTE_IE", "format_CMed_bootstrap")) %dopar% {    
    
    cat(paste0("Iteration: ", i, "\n")) # progress indicator
    
    # Prepare data for analysis (function defined below) 
    expandedData <- CMed_prep(data, prop) 
    
    # Run model
    model <- aalen(Surv(bl_age, end_age, allcause_death) ~  const(A.race) * const(race_M1.alc) + 
                                                            const(A.race) * const(race_M2.smk) +
                                                            const(A.race) * const(race_M3.bmi) +
                                                            const(A.race) * const(race_M4.phy) +
                                                            const(married) + const(factor(edu)) + const(factor(srvy_yr)),
                                     start.time=start_time, max.time=end_time,
                                     data=expandedData, weights=expandedData$weightM, clusters=expandedData$ID, robust=0)  
    
    cat(paste0("     Step 4 Complete (data analysis)", "\n")) # progress indicator
    
    
    # Extract coefficients of interest
      # First, list the coefficients of interest for each race/ethnicity
      Black    <- c(1,4,7,10,13,36,45,54,63)
      Hispanic <- c(2,5,8,11,14,40,49,58,67)
      Other    <- c(3,6,9,12,15,44,53,62,71)
      
            # Review model coefficients and selected coefficients to ensure they align
            if (i==1) {
              model_coefficients <- coef(model) %>%
                row.names() %>% 
                as.data.frame() %>%
                rownames_to_column()
              
              print("model coefficients"); print(model_coefficients)
              print("Black"); print(Black)
              print("Hispanic"); print(Hispanic)
              print("Other"); print(Other)
            }
      
      # Second, get and format results for each level of exposure (function defined below)
      results_black <- format_CMed_bootstrap (model, Black)  
      results_hisp  <- format_CMed_bootstrap (model, Hispanic)
      results_other <- format_CMed_bootstrap (model, Other)
    
      final <- rbind (results_black, results_hisp, results_other)
  } 
}


# Format results from the bookstrapping
format_CMed_ageStrat <- function(bootstrap_results){
  
  # Compute CI and format results 
  mean <- rowMeans(bootstrap_results) # get mean estimate
  
  ci <- apply(bootstrap_results, 1, quantile, probs=c(0.025, 0.975)) %>% t() #get 95% CI
  
  results <- cbind(mean, ci) %>% as.data.frame() %>% 
    
    # add labels
    mutate (term = rep(c( "Total effect of race (ref=White)", 
      "Direct effect of raca (ref=White)", 
      "Indirect effect of race (ref=White)", 
      "     Alcohol use: differential exposure", 
      "     Alcohol use: differential vulnerability ", 
      "     Smoking: differential exposure", 
      "     Smoking: differential vulnerability ", 
      "     BMI: differential exposure", 
      "     BMI: differential vulnerability ", 
      "     Physical activity: differential exposure", 
      "     Physical activity: differential vulnerability"), 6),
      race = rep(c("Black", "Hispanic", "Other"), each=22),
      type = rep(c("deaths", "prop"), 3, each=11)) %>% 
    
    # Separate the 'additional deaths' and 'proportion' estimates 
    pivot_wider (names_from="type", values_from=c("mean", "2.5%", "97.5%")) %>%
    
    # reformat 
    mutate (deaths = round(mean_deaths*10000,1),
      deaths_lower = round(`2.5%_deaths`*10000,1),
      deaths_upper = round(`97.5%_deaths`*10000,1),
      prop = round(mean_prop*100,0),
      prop_lower = round(`2.5%_prop`*100,0),
      prop_upper = round(`97.5%_prop`*100,0),
      deaths_10000py_ci = paste0(deaths, " (", deaths_lower, ", ", deaths_upper, ")"),
      prop_ci = paste0(prop, " (", prop_lower, ", ", prop_upper, ")"))%>%
    dplyr::select(race, term, deaths_10000py_ci, prop_ci)
}



# Supporting functions ------------------------------------------------------------------------------------------------------------------
# Functions below have been utilized in the two functions above 

# Data preparation for Causal Mediation
CMed_prep <-function(data, prop = 1) {
  
  ### Step 0: Select data to use *****************************************************************************************************************
  # **********************************************************************************************************************************************
  mydata <- data %>%
    mutate(
      A.race = factor(ethnicity, levels=c(1,2,3,4), labels = c("White", "Black", "Hispanic", "Other")),
      M1.alc = alcohol5v2,
      M2.smk = smoking4,
      M3.bmi = bmi_cat,
      M4.phy = phy_act3) %>%
    dplyr::select(A.race, M1.alc, M2.smk, M3.bmi, M4.phy, allcause_death, bl_age, end_age, married, edu, srvy_yr) %>%
    sample_frac(prop, replace=FALSE) # sample a fraction of the data
  
  cat(paste0("     Step 0 Complete (select data)", "\n")) # progress indicator
  cat(paste0("     Original sample size: ", nrow(mydata)), "\n") 
  
  # NOTE: For technical reasons, the mediators should be coded as integers starting with 1
  
  
  ### Step 1: Fit a model for each mediator, , conditioning on exposure and all confounders *******************************************************
  # ***********************************************************************************************************************************************
  
  # Fit model for each mediator, conditioning on exposure (race) and all confounders
  
  mydata$ATemp <- mydata$A.race # first, create and use a copy of the exposure variable (for technical reasons related to R)
  fitM1 <- vglm(M1.alc ~ ATemp + bl_age + married + factor(edu) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 3))
  fitM2 <- vglm(M2.smk ~ ATemp + bl_age + married + factor(edu) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 1))
  fitM3 <- vglm(M3.bmi ~ ATemp + bl_age + married + factor(edu) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 2))
  fitM4 <- vglm(M4.phy ~ ATemp + bl_age + married + factor(edu) + factor(srvy_yr), data = mydata, family=multinomial(refLevel = 3))
  
  cat(paste0("     Step 1 Complete (Fit model for each mediator)", "\n")) # progress indicator
 
  ### Step 2: Construct copies of ID and exposure *************************************************************************************************
  # ***********************************************************************************************************************************************
  
  #Create ID Variable
  mydata$ID <- 1:nrow(mydata) # construct id variable
  
  # Create counterfactual version of exposure (race); repeated 4 times because there are 4 mediators
  levelsOfRACE <- unique(mydata$A.race)
  myData1 <- mydata
  myData2 <- mydata
  myData3 <- mydata
  myData4 <- mydata
  myData1$race_M1.alc <- levelsOfRACE[1]
  myData2$race_M1.alc <- levelsOfRACE[2]
  myData3$race_M1.alc <- levelsOfRACE[3]
  myData4$race_M1.alc <- levelsOfRACE[4]
  tempMyData <- rbind(myData1, myData2, myData3, myData4)
  
  myData1 <- tempMyData
  myData2 <- tempMyData
  myData3 <- tempMyData
  myData4 <- tempMyData
  myData1$race_M2.smk <- levelsOfRACE[1]
  myData2$race_M2.smk <- levelsOfRACE[2]
  myData3$race_M2.smk <- levelsOfRACE[3]
  myData4$race_M2.smk <- levelsOfRACE[4]
  tempMyData <- rbind(myData1, myData2, myData3, myData4)
  
  myData1 <- tempMyData
  myData2 <- tempMyData
  myData3 <- tempMyData
  myData4 <- tempMyData
  myData1$race_M3.bmi <- levelsOfRACE[1]
  myData2$race_M3.bmi <- levelsOfRACE[2]
  myData3$race_M3.bmi <- levelsOfRACE[3]
  myData4$race_M3.bmi <- levelsOfRACE[4]
  tempMyData <- rbind(myData1, myData2, myData3, myData4)
  
  myData1 <- tempMyData
  myData2 <- tempMyData
  myData3 <- tempMyData
  myData4 <- tempMyData
  myData1$race_M4.phy <- levelsOfRACE[1]
  myData2$race_M4.phy <- levelsOfRACE[2]
  myData3$race_M4.phy <- levelsOfRACE[3]
  myData4$race_M4.phy <- levelsOfRACE[4]
  newMyData <- rbind(myData1, myData2, myData3, myData4)
  
  cat(paste0("     Step 2 Complete (expand data)", "\n")) # progress indicator
  cat(paste0("     Expanded sample size: ", nrow(newMyData)), "\n") 
  
  ### Step 3: Construct weights  *********************************************************************************************************************
  # **************************************************************************************************************************************************
  
  # M1: alcohol
  newMyData$ATemp <- newMyData$A.race
  tempDir1 <- as.matrix(predict(fitM1,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M1.alc)]
  
  newMyData$ATemp <- newMyData$race_M1.alc
  tempIndir1 <- as.matrix(predict(fitM1,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M1.alc)]
  
  newMyData$weight1 <- tempIndir1/tempDir1
  cat(paste0("     Step 3.1 Complete (weghts for alcohol)", "\n")) # progress indicator
  
  
  #M2: Smoking
  newMyData$ATemp <- newMyData$A.race
  tempDir2 <- as.matrix(predict(fitM2,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M2.smk)]
  
  newMyData$ATemp <- newMyData$race_M2.smk
  tempIndir2 <- as.matrix(predict(fitM2,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M2.smk)]
  
  newMyData$weight2 <- tempIndir2/tempDir2
  cat(paste0("     Step 3.2 Complete (weghts for smoking)", "\n")) # progress indicator
  
  
  #M3: BMI
  newMyData$ATemp <- newMyData$A.race
  tempDir3 <- as.matrix(predict(fitM3,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M3.bmi)]
  
  newMyData$ATemp <- newMyData$race_M3.bmi
  tempIndir3 <- as.matrix(predict(fitM3,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M3.bmi)]
  
  newMyData$weight3 <- tempIndir3/tempDir3
  cat(paste0("     Step 3.3 Complete (weghts for BMI)", "\n")) # progress indicator
  
  
  #M4: Physical activity
  newMyData$ATemp <- newMyData$A.race
  tempDir4 <- as.matrix(predict(fitM4,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M4.phy)]
  
  newMyData$ATemp <- newMyData$race_M4.phy
  tempIndir4 <- as.matrix(predict(fitM4,type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$M4.phy)]
  
  newMyData$weight4 <- tempIndir4/tempDir4
  cat(paste0("     Step 3.4 Complete (weghts for physical activity)", "\n")) # progress indicator
  
  
  # Final weight
  newMyData$weightM <- newMyData$weight1 * newMyData$weight2 * newMyData$weight3 * newMyData$weight4
  
  newMyData <- newMyData %>%
    dplyr::select(ID, bl_age, end_age, allcause_death, A.race, race_M1.alc, race_M2.smk, race_M3.bmi, race_M4.phy, married, edu, srvy_yr, weightM)
  
  
  return(newMyData)
  cat(paste0("     Step 3 Complete (weghts combined)", "\n")) # progress indicator
  
}

# Function to get the total effect and proportion mediated
getTE <- function(CMed_model, v){
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

# Function to get the total combined indirect effect  
getIE <- function(CMed_model, v){
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

# Function to get the proportion mediated of the combined indirect effect
getTE_IE <- function(CMed_model, v, z){
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

# Format results to export bootstrapped results
format_CMed_bootstrap <- function (model, coef_list) {
        
        group <- enexpr(coef_list)
        
        # All_coef model coefficients   
        All_coef <- coef(model) %>%
          as.data.frame() %>% 
          rownames_to_column(var = "term") %>% 
          mutate (estimate = Coef.) %>% 
          dplyr::select (term, estimate) 
        
        coef <- slice(All_coef, coef_list)
        
        # Function to get the total effect and proportion mediated
        TE_prop <- getTE(model, coef_list) %>% 
          as.data.frame() %>% rownames_to_column(var = "term")
        
        TE <- filter(TE_prop, term =="TE") %>% 
          mutate (estimate = Est.) %>% 
          dplyr::select (term, estimate) 
        
        prop <- filter(TE_prop, term !="TE") %>%
          dplyr::select(term, med_prop) 
        
        
        # Function to get the total combined indirect effect  (coef_list excluding 1st item)
        IE <- getIE(model, coef_list[-1]) %>% 
          as.data.frame() %>% rownames_to_column(var = "term") %>% 
          filter (term == "IE") %>% 
          mutate (estimate = Est.) %>% 
          dplyr::select (term, estimate) 
        
        
        # Function to get the proportion mediated of the combined indirect effect
        IE_prop <- getTE_IE(model, coef_list, coef_list[-1]) %>% 
          as.data.frame() %>% rownames_to_column(var = "term") %>%
          mutate (term = ifelse(term=="2.5%", "IE", term)) %>%
          filter (term=="IE") %>%
          dplyr::select(term, med_prop)
        
        
        one <- full_join(coef, prop, by="term") %>%
          mutate(label = c(paste0("02 Direct effect of ", group, " (ref=White)"),
            "04 Alcohol use: differential exposure",
            "06 Smoking: differential exposure",
            "08 BMI: differential exposure",
            "10 Physical activity: differential exposure", 
            "05 Alcohol use: differential vulnerability",
            "07 Smoking: differential vulnerability",
            "09 BMI: differential vulnerability",
            "11 Physical activity: differential vulnerability")) %>% 
          dplyr::select(label, estimate, med_prop)
        
        two <- TE %>%
          mutate (label = ifelse(term == "TE", paste0("01 Total effect of ", group, " (ref=White)"), NA),
            med_prop = "1") %>% 
          dplyr::select(label, estimate, med_prop)
        
        three <- full_join(IE, IE_prop, by="term") %>% 
          mutate (label = ifelse(term == "IE", paste0("03 Indirect effect of ", group, " (ref=White)"), NA)) %>% 
          dplyr::select(label, estimate, med_prop) 
        
        
        final <- rbind(one, two, three) %>%
          mutate (med_prop = as.numeric(med_prop)) %>% 
          pivot_longer(cols = c("estimate", "med_prop"), names_to = "type", values_to = "estimate") %>%
          arrange(type, label) %>%
          dplyr::select(estimate)
        
        return(final)
}
      

      
      
