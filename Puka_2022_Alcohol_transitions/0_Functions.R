
# SIMAH - NESARC Alcohol Transitions
# Functions used

# Table 2 - Extract Transition Probabilities (TP) and correct CI to original sample size 
predicted_TP <- function(model, year) {
  
  table <- data.frame(print(pmatrix.msm(model, t=year, ci="norm"))) %>%
    mutate(From = row.names(.)) %>%
    pivot_longer(cols = -From, names_to = "To") %>%
    separate(value, into=c("Estimate","Lower","Upper", NA), sep="\\(|\\,|\\)", convert=TRUE) %>%  # separated based on "(" "," and  ")"  convert=TRUE names variables numeric
    mutate(
      SE = (Upper - Lower) / 3.92,
      SD = SE * sqrt(nrow(nesarc_expanded)),  
      newSE = SD / sqrt(nrow(nesarc)), 
      newLower = round((Estimate - (newSE * 1.96))*100, digits=1),
      newUpper = round((Estimate + (newSE * 1.96))*100, digits=1),
      Estimate = round(Estimate*100, digits=1), 
      EstimateCI = paste0(Estimate, " (", newLower, ", ", newUpper, ")")) %>%
    select (From, To, EstimateCI) %>%
    pivot_wider(names_from = "To", values_from = "EstimateCI")
  
  return(table)
}


# Table 3 - Extract HR results, rearrange, and correct CI to original sample size 
HR_table <- function (model){
  
  table <- data.frame(hazard.msm(model))%>%
    mutate(transition = row.names(.)) %>%
    pivot_longer(cols=-transition) %>%
    extract(name, into=c("Variable","Type"), regex="(.*)\\.(.*)") %>%   # Separate the string (name) into the variable and type of estimate (HR, Upper, Lower), separate at last occuring period .
    pivot_wider(names_from="Type", values_from = "value") %>%
    mutate(
      SE=(log(U) - log(L))/3.92, 
      SD = SE * sqrt(nrow(nesarc_expanded)),  
      newSE = SD / sqrt(nrow(nesarc)),
      newLower = round(exp((log(HR)) - (newSE * 1.96)), digits=2), 
      newUpper = round(exp((log(HR)) + (newSE * 1.96)), digits=2), 
      HR = round(HR, digits=2),
      EstimateCI = paste0(HR, " (", newLower, ", ", newUpper, ")")) %>%
    select(transition, Variable, EstimateCI) %>%
    pivot_wider(names_from = "transition", values_from = EstimateCI) 

  return (table) 
}



# Transition probabilities, all covariates, multiple years
predicted_TP_covs <- function(model, year) {

  age_cat <- unique(nesarc_expanded$age7)
  sex <- unique(nesarc_expanded$female_w1)
  race <- unique(nesarc_expanded$race_w1)
  edu <- unique(nesarc_expanded$edu3)
  
  probs <- list()
  for (i in age_cat){
    for (j in sex){
      for (l in race){
        for (k in edu){
          for (m in 1:year){

          # extract the probabilities
          probs[[paste(i,j,l,k)]] <- data.frame(print(pmatrix.msm(model, t=m, ci="norm",
                                                covariates = list(j, i, k, l)))) %>%  # order is same as that used in the model: female, age, edu, race 
            # modify the output presentation
            mutate(From = row.names(.)) %>%
            pivot_longer(cols = -From, names_to = "To") %>%
            separate(value, into=c("Probability","Lower","Upper", NA), sep="\\(|\\,|\\)", convert=TRUE) %>%  # separated based on "(" "," and  ")"  convert=TRUE names variables numeric
            mutate(SE = (Upper - Lower) / 3.92,
                   SD = SE * sqrt(nrow(nesarc_expanded)),  
                   newSE = SD / sqrt(nrow(nesarc)), 
                   lowerCI = Probability - (newSE * 1.96),
                   upperCI = Probability + (newSE * 1.96),
                   age_cat = i,
                   sex = j, 
                   race = l,
                   edu = k,
                   year = m) %>%
            select (year, age_cat, sex, edu, race, From, To, Probability, lowerCI, upperCI)
          }
        }
      }
    }
  }
  probs <- do.call(rbind, probs)
  row.names(probs) <- NULL  # remove row names
  return(probs)
}



# Transition probabilities by age, multiple years
predicted_TP_age <- function(model, year) {
  
  age_cat <- unique(nesarc_expanded$age7)
  
  probs <- list()
  for (i in 1:year) {
    for (j in age_cat){
      # extract the probabilities
      probs[[paste(i,j)]] <- data.frame(print(pmatrix.msm(model, t=i, ci="norm", covariates = list(age7 = j)))) %>%  
        
        # modify the output presentation
        mutate(From = row.names(.)) %>%
        pivot_longer(cols = -From, names_to = "To") %>%
        separate(value, into=c("Probability","Lower","Upper", NA), sep="\\(|\\,|\\)", convert=TRUE) %>%  # separated based on "(" "," and  ")"  convert=TRUE names variables numeric
        mutate(SE = (Upper - Lower) / 3.92,
              SD = SE * sqrt(nrow(nesarc_expanded)),  
              newSE = SD / sqrt(nrow(nesarc)), 
              lowerCI = Probability - (newSE * 1.96),
              upperCI = Probability + (newSE * 1.96),
              age_cat = j,
              year = i) %>%
        select (year, age_cat, From, To, Probability, lowerCI, upperCI)
    }
  }
  probs <- do.call(rbind, probs)
  row.names(probs) <- NULL  # remove row names
  return(probs)
}


# Transition probabilities by edu, multiple years
predicted_TP_edu <- function(model, year) {
  
  edu <- unique(nesarc_expanded$edu3)
  
  probs <- list()
  for (i in 1:year) {
    for (j in edu){
      # extract the probabilities
      probs[[paste(i,j)]] <- data.frame(print(pmatrix.msm(model, t=i, ci="norm", covariates = list(edu3 = j)))) %>%  
        
        # modify the output presentation
        mutate(From = row.names(.)) %>%
        pivot_longer(cols = -From, names_to = "To") %>%
        separate(value, into=c("Probability","Lower","Upper", NA), sep="\\(|\\,|\\)", convert=TRUE) %>%  # separated based on "(" "," and  ")"  convert=TRUE names variables numeric
        mutate( SE = (Upper - Lower) / 3.92,
                SD = SE * sqrt(nrow(nesarc_expanded)),  
                newSE = SD / sqrt(nrow(nesarc)), 
                lowerCI = Probability - (newSE * 1.96),
                upperCI = Probability + (newSE * 1.96),
                edu = j,
                year = i) %>%
        select (year, edu, From, To, Probability, lowerCI, upperCI)
    }
  }
  probs <- do.call(rbind, probs)
  row.names(probs) <- NULL  # remove row names
  return(probs)
}




# Transition probabilities by sex, multiple years
predicted_TP_sex <- function(model, year) {
  
  sex <- unique(nesarc_expanded$female_w1)
  
  probs <- list()
  for (i in 1:year) {
    for (j in sex){
      # extract the probabilities
      probs[[paste(i,j)]] <- data.frame(print(pmatrix.msm(model, t=i, ci="norm", covariates = list(female_w1 = j)))) %>%  
        
        # modify the output presentation
        mutate(From = row.names(.)) %>%
        pivot_longer(cols = -From, names_to = "To") %>%
        separate(value, into=c("Probability","Lower","Upper", NA), sep="\\(|\\,|\\)", convert=TRUE) %>%  # separated based on "(" "," and  ")"  convert=TRUE names variables numeric
        mutate( SE = (Upper - Lower) / 3.92,
                SD = SE * sqrt(nrow(nesarc_expanded)),  
                newSE = SD / sqrt(nrow(nesarc)), 
                lowerCI = Probability - (newSE * 1.96),
                upperCI = Probability + (newSE * 1.96),
                sex = j,
                year = i) %>%
        select (year, sex, From, To, Probability, lowerCI, upperCI)
    }
  }
  probs <- do.call(rbind, probs)
  row.names(probs) <- NULL  # remove row names
  return(probs)
}




# Transition probabilities by race, multiple years
predicted_TP_race <- function(model, year) {
  
  race <- unique(nesarc_expanded$race_w1)
  
  probs <- list()
  for (i in 1:year) {
    for (j in race){
      # extract the probabilities
      probs[[paste(i,j)]] <- data.frame(print(pmatrix.msm(model, t=i, ci="norm", covariates = list(race_w1 = j)))) %>%  
        
        # modify the output presentation
        mutate(From = row.names(.)) %>%
        pivot_longer(cols = -From, names_to = "To") %>%
        separate(value, into=c("Probability","Lower","Upper", NA), sep="\\(|\\,|\\)", convert=TRUE) %>%  # separated based on "(" "," and  ")"  convert=TRUE names variables numeric
        mutate(SE = (Upper - Lower) / 3.92,
               SD = SE * sqrt(nrow(nesarc_expanded)),  
               newSE = SD / sqrt(nrow(nesarc)), 
               lowerCI = Probability - (newSE * 1.96),
               upperCI = Probability + (newSE * 1.96),
               race = j,
               year = i) %>%
        select (year, race, From, To, Probability, lowerCI, upperCI)
    }
  }
  probs <- do.call(rbind, probs)
  row.names(probs) <- NULL  # remove row names
  return(probs)
}




# Transition probabilities, no covariates (average level of covariates), multiple years
predicted_TP_noCovs <- function(model, year) {
  probs <- list()
  for (i in 1:year) {
      # extract the probabilities
      probs[[paste(i)]] <- data.frame(print(pmatrix.msm(model, t=i, ci="norm"))) %>%  
        # modify the output presentation
        mutate(From = row.names(.)) %>%
        pivot_longer(cols = -From, names_to = "To") %>%
        separate(value, into=c("Probability","Lower","Upper", NA), sep="\\(|\\,|\\)", convert=TRUE) %>%  # separated based on "(" "," and  ")"  convert=TRUE names variables numeric
        mutate( SE = (Upper - Lower) / 3.92,
                SD = SE * sqrt(nrow(nesarc_expanded)),  
                newSE = SD / sqrt(nrow(nesarc)), 
                lowerCI = Probability - (newSE * 1.96),
                upperCI = Probability + (newSE * 1.96),
                year = i) %>%
        select (year, From, To, Probability, lowerCI, upperCI)
    }
  probs <- do.call(rbind, probs)
  row.names(probs) <- NULL  # remove row names
  return(probs)
}













# Function to apply transition probabilities

set.seed(123)  # To get consistent results

transition_alc4 <- function(data, transitions){
  selected <- unique(data$cat)
  rates <- transitions %>% filter(cat == selected)
  data$predicted_cat <- ifelse(data$prob<=rates$cumsum[1], "Non-drinker",
    ifelse(data$prob<=rates$cumsum[2] & data$prob>rates$cumsum[1], "Category I",
      ifelse(data$prob<=rates$cumsum[3] & data$prob>rates$cumsum[2],"Category II",
        ifelse(data$prob<=rates$cumsum[4] & data$prob>rates$cumsum[3],"Category III",NA))))
  return(data)
}


transition_hed <- function(data, transitions){
  selected <- unique(data$cat)
  rates <- transitions %>% filter(cat == selected)
  data$predicted_cat <- ifelse(data$prob <= rates$cumsum[1], "Non-drinker",
    ifelse(data$prob <= rates$cumsum[2] & data$prob > rates$cumsum[1], "Drinker, no HED",
      ifelse(data$prob <= rates$cumsum[3] & data$prob > rates$cumsum[2],"Occasional HED",
        ifelse(data$prob <= rates$cumsum[4] & data$prob > rates$cumsum[3],"Monthly HED",
          ifelse(data$prob <= rates$cumsum[5] & data$prob > rates$cumsum[4],"Weekly HED",NA)))))
  return(data)
}



# Simulate Transitions
simulate_pop <- function(initial_pop, aTP, apply_transitions, years) {
  pop <- initial_pop
  for (i in 1:years) {
    pop <- pop %>% 
      mutate( year= i,
              cat = paste(sex, age_cat, edu, race, predicted_cat, sep="_"),
              prob = runif(nrow(.))) %>%  # generate random prob
      group_by(cat) %>%
        do(apply_transitions(., aTP)) %>% # use 'do( )' to run the function defined earlier
      ungroup() %>% 
      select (-cat, -prob) 
  }
  return(pop)
}
    



# Compare observed vs predicted
compare_pct <- function(predicted_pop, observed_pop, strata="All") {

  # Data prep
  predicted_pop <- predicted_pop %>% mutate(strata = {{strata}})
  observed_pop  <- observed_pop %>% mutate(strata = {{strata}})
  
  # Analyses
  predicted_count <- predicted_pop %>%
    group_by(strata) %>% 
    count(predicted_cat) %>%                                          
    rename(predicted = n, Category = predicted_cat) %>%
    mutate(pred_pct = round(predicted / sum(predicted) * 100,1)) 
  
  observed_count <- observed_pop %>% 
    group_by(strata) %>% 
    count(observed_cat) %>%                                           
    rename(observed = n, Category = observed_cat) %>%
    mutate(obs_pct = round(observed / sum(observed) * 100,1))
  
  comparison <- full_join (observed_count, predicted_count, by=c("Category", "strata")) %>% 
    mutate(diff_pct = round((pred_pct-obs_pct), 2))%>% 
    select (strata, Category, obs_pct, pred_pct, diff_pct) 
  
  return(comparison)
}
  
  

