
# Race x Lifestyle Differential Vulnerability & Exposure Project
## Sensitivity Analysis (one mediator at a time)

# Causal Mediation Analysis
CMed_oneVar <-function(data, mediator, cov1, cov2, cov3, M_ref_level, black_coef, hispanic_coef, other_coef) {

  
  ### Step 0: Select data to use *****************************************************************************************************************
  # **********************************************************************************************************************************************
  mydata <- data %>%
    mutate(
      A.race = factor(ethnicity, levels=c(1,2,3,4), labels = c("White", "Black", "Hispanic", "Other")),
      Mediator = {{mediator}},
      Cov1 = {{cov1}}, 
      Cov2 = {{cov2}},
      Cov3 = {{cov3}}) %>%
    dplyr::select(A.race, Mediator, Cov1, Cov2, Cov3, allcause_death, bl_age, end_age, married, edu, srvy_yr) 
  
  
  cat(paste0("     Step 0 Complete (select data)", "\n")) # progress indicator
  cat(paste0("     Original sample size: ", nrow(mydata)), "\n") 
  
  # NOTE: For technical reasons, the mediators should be coded as integers starting with 1
  
  
  ### Step 1: Fit a model for each mediator, , conditioning on exposure and all confounders *******************************************************
  # ***********************************************************************************************************************************************
  
  # Fit model for each mediator, conditioning on exposure (race) and all confounders
  
  mydata$ATemp <- mydata$A.race # first, create and use a copy of the exposure variable (for technical reasons related to R)
  fit_Med <- vglm(Mediator ~ ATemp + bl_age + married + factor(edu) + factor(srvy_yr) +
                             factor(Cov1) + factor(Cov2) + factor(Cov3),
                             data = mydata, family=multinomial(refLevel = M_ref_level))
  
  cat(paste0("     Step 1 Complete (Fit model for mediator)", "\n")) # progress indicator
  
  ### Step 2: Construct copies of ID and exposure *************************************************************************************************
  # ***********************************************************************************************************************************************
  
  #Create ID Variable
  mydata$ID <- 1:nrow(mydata) # construct id variable
  
  # Create counterfactual version of exposure (race)
  levelsOfRACE <- unique(mydata$A.race)
  myData1 <- mydata
  myData2 <- mydata
  myData3 <- mydata
  myData4 <- mydata
  myData1$race_M <- levelsOfRACE[1]
  myData2$race_M <- levelsOfRACE[2]
  myData3$race_M <- levelsOfRACE[3]
  myData4$race_M <- levelsOfRACE[4]
  newMyData <- rbind(myData1, myData2, myData3, myData4)
  
  
  cat(paste0("     Step 2 Complete (expand data)", "\n")) # progress indicator
  cat(paste0("     Expanded sample size: ", nrow(newMyData)), "\n") 
  
  
  ### Step 3: Construct weights  *********************************************************************************************************************
  # **************************************************************************************************************************************************
  
  newMyData$ATemp <- newMyData$A.race
  tempDir1 <- as.matrix(predict(fit_Med, type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$Mediator)]
  
  newMyData$ATemp <- newMyData$race_M
  tempIndir1 <- as.matrix(predict(fit_Med, type = "response", newdata=newMyData))[cbind(1:nrow(newMyData),newMyData$Mediator)]
  
  newMyData$weightM <- tempIndir1/tempDir1
  cat(paste0("     Step 3 Complete (weghts for mediator)", "\n")) # progress indicator
  
  
  # Keep relavent variables 
  newMyData <- newMyData %>%
    dplyr::select(ID, bl_age, end_age, allcause_death, A.race, race_M, married, edu, Cov1, Cov2, Cov3, srvy_yr, weightM)
  
  
  ### Step 4: Run Causal Mediation Analysis  *********************************************************************************************************
  # **************************************************************************************************************************************************
  
  # Run model
  model <- aalen(Surv(bl_age, end_age, allcause_death) ~  const(A.race) * const(race_M) + 
                                                          const(factor(Cov1)) +  const(factor(Cov2)) + const(factor(Cov3)) + 
                                                          const(married) + const(factor(edu)) + const(factor(srvy_yr)),
                data=newMyData, weights=newMyData$weightM, clusters=newMyData$ID, robust=0)  
  
  cat(paste0("     Step 4 Complete (data analysis)", "\n")) # progress indicator

  
  
  # List the coefficients of interest for each race/ethnicity
  Black    <- black_coef
  Hispanic <- hispanic_coef
  Other    <- other_coef
  
  # Get and format results (Function defined below)
  results_black <- format_CMed (model, Black)  
  results_hisp  <- format_CMed (model, Hispanic)
  results_other <- format_CMed (model, Other)
  
  # Combine
  combined <- rbind (results_black, results_hisp, results_other)
  return(combined)
}


# Format results to export bootstrapped results
format_CMed <- function (model, coef_list) {
  
  group <- enexpr(coef_list)
  
  # All_coef model coefficients   
  All_coef <- coef(model) %>%
    as.data.frame() %>% 
    rownames_to_column(var = "term") %>% 
    mutate (deaths_10000py = round(Coef. * 10000, 1),
      lower = round(`lower2.5%` * 10000, 1),
      upper = round(`upper97.5%` * 10000, 1),
      deaths_10000py_CI = paste0(deaths_10000py, " (", lower, ", ", upper, ")")) %>% 
    dplyr::select (term, deaths_10000py_CI)
  
  coef <- slice(All_coef, coef_list)
  
  
  # Function to get the total effect and proportion mediated
  TE_prop <- getTE(model, coef_list) %>% 
    as.data.frame() %>% rownames_to_column(var = "term")
  
  TE <- filter(TE_prop, term =="TE") %>% 
    mutate (lower = round((`Est.` - (1.96 * SE))*10000,1), 
      upper = round((`Est.` + (1.96 * SE))*10000,1), 
      deaths_10000py = round(`Est.` * 10000,1),
      deaths_10000py_CI = paste0(deaths_10000py, " (", lower, ", ", upper, ")")) %>%
    dplyr::select(term, deaths_10000py_CI)
  
  prop <- filter(TE_prop, term !="TE") %>%
    mutate(prop = round(med_prop * 100,0),
      lower = round(lowerCI * 100,0),
      upper = round(UpperCI * 100,0),
      prop_CI = paste0(prop, " (", lower, ", ", upper, ")")) %>% 
    dplyr::select(term, prop_CI) 
  
  
  
  # Function to get the total combined indirect effect  (coef_list excluding 1st item)
  IE <- getIE(model, coef_list[-1]) %>% 
    as.data.frame() %>% rownames_to_column(var = "term") %>% 
    filter (term == "IE") %>% 
    mutate (lower = round((`Est.` - (1.96 * SE))*10000,1), 
      upper = round((`Est.` + (1.96 * SE))*10000,1), 
      deaths_10000py = round(`Est.` * 10000,1),
      deaths_10000py_CI = paste0(deaths_10000py, " (", lower, ", ", upper, ")")) %>%
    dplyr::select(term, deaths_10000py_CI) 
  
  
  
  # Function to get the proportion mediated of the combined indirect effect
  IE_prop <- getTE_IE(model, coef_list, coef_list[-1]) %>% 
    as.data.frame() %>% rownames_to_column(var = "term") %>%
    pivot_wider(names_from="term", values_from=c("IE", "med_prop", "quantile")) %>%
    mutate (IE_prop = round(`med_prop_2.5%` * 100, 0),
      lower = round(`quantile_2.5%` *100, 0), 
      upper = round(`quantile_97.5%` *100, 0),
      prop_CI = paste0(IE_prop, " (", lower, ", ", upper, ")"),
      term = "IE") %>% 
    dplyr::select(term, prop_CI) 
  
  
  
  # Edit and combine results
  one <- full_join(coef, prop, by="term") %>% 
    mutate(label = c(paste0("02 Direct effect of ", group, " (ref=White)"),
      "04 Differential exposure",
      "05 Differential vulnerability")) %>% 
    dplyr::select(term, label, deaths_10000py_CI, prop_CI)
  
  two <- TE %>%
    mutate (label = ifelse(term == "TE", paste0("01 Total effect of ", group, " (ref=White)"), NA),
      prop_CI = "1") %>% 
    dplyr::select(term, label, deaths_10000py_CI, prop_CI)
  
  three <- full_join(IE, IE_prop, by="term") %>% 
    mutate (label = ifelse(term == "IE", paste0("03 Indirect effect of ", group, " (ref=White)"), NA)) %>% 
    dplyr::select(term, label, deaths_10000py_CI, prop_CI) 
  
  
  final <- rbind(one, two, three) %>% 
    arrange(label)
  
  return(final)
}

# Function to get the (Not-Robust) total effect and proportion mediated
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

# Function to get the (Not-Robust) total combined indirect effect  
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

# Function to get the (Not Robust)  proportion mediated of the combined indirect effect
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

