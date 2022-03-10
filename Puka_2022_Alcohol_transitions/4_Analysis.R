
# SIMAH - NESARC Alcohol Transitions
# Data Analysis

library(tidyverse)  # data management
library(skimr)      # descriptive statistics
library(janitor)    # data management
library(msm)        # model transition probabilities
library(tableone)   # create descriptives table
library(knitr)      # create descriptives table
# options(scipen=999) # prevent the use of scientific notation
memory.limit(size=1e+13)

 
# Specify the data and output file locations
data    <- "C:/.../Processed data/"  # Location of data
output  <- "C:/.../nesarc/Output/"   # Location of tables and figures 
models  <- "C:/.../Models/"          # Location of saved MSM models


# Load data / functions
nesarc_all       <- readRDS(paste0(data, "nesarc_all.rds")) # Contains those with missing data 
nesarc           <- readRDS(paste0(data, "nesarc_clean.rds")) 
nesarc_expanded  <- readRDS(paste0(data, "nesarc_clean_expanded.rds")) 
nesarc3_expanded <- readRDS(paste0(data, "nesarc3_clean_expanded.rds")) 
source("0_Functions.R")


# Load MSM Models (created in 'Model Selection' file)
alc4_crude.msm <- readRDS(paste0(models, "msm3b_crude.RDS"))
alc4.msm       <- readRDS(paste0(models, "msm3b.RDS"))
hed_crude.msm  <- readRDS(paste0(models, "hed_crude.msm.RDS"))
hed.msm        <- readRDS(paste0(models, "hed.msm.RDS"))


# 1) Descriptives ------------------------------------------------------------------------------------------

# years follow-up
nesarc %>%
  filter(wave==2) %>%
  select (years) %>% 
  skim()


# Variables of interest
variables <- c("female", "age", "age7", "race.factor", "edu3", "alc4.factor", "hed.factor")
factor_vars <- c("female", "age7", "race.factor", "edu3", "alc4.factor", "hed.factor")

# Descriptives of expanded data
tab_exp <- CreateTableOne(vars= variables, factorVars = factor_vars, strata="wave.factor", data=nesarc_expanded)
table1_exp <- print(tab_exp, noSpaces = TRUE, catDigits = 0, contDigits = 1, pDigits = 2, printToggle = FALSE, test=FALSE, format="p")  # Shows only % 
write.csv(table1_exp, file=paste0(output,"Table 1 - Descriptives of expanded data.csv")) 
kable(table1_exp)


# Descriptives at basleine and follow-up of included participants 
tab1 <- CreateTableOne(vars= variables, factorVars = factor_vars, strata="wave.factor", data=nesarc)
table1 <- print(tab1, noSpaces = TRUE, catDigits = 0, contDigits = 1, pDigits = 2, printToggle = FALSE, test=FALSE)   
write.csv(table1, file=paste0(output,"Table S1 - Descriptives of original data.csv"))  # export to excel, to copy/paste into manuscript
kable(table1)


# Descriptives for attrition 
nesarc1_all <- filter(nesarc_all, wave==1) # select baseline data
tab1_attr <- CreateTableOne(vars= variables, factorVars = factor_vars, strata="lost.factor", data=nesarc1_all)
table1_attr <- print(tab1_attr, noSpaces = TRUE, catDigits = 0, contDigits = 1, pDigits = 2, printToggle = FALSE, test=FALSE, smd=TRUE)   
write.csv(table1_attr, file=paste0(output,"Table S2 - Attrition Descriptives.csv"))  # export to excel, to copy/paste into manuscript
kable(table1_attr)                             # view in R; R Markdown friendly version


# Descriptives NESARC I and III (expanded data)
tab1_exp3 <- CreateTableOne(vars= variables, factorVars = factor_vars, data=nesarc3_expanded)
table1_exp3 <- print(tab1_exp3, noSpaces = TRUE, catDigits = 0, contDigits = 1, pDigits = 2, printToggle = FALSE, test=FALSE, format="p")  # Shows only % 
write.csv(table1_exp3, file=paste0(output,"Table S6 - NESARC III Descriptives of expanded data.csv")) 
kable(table1_exp3)



# 2) ALCOHOL CONSUMPTION  ---------------------------------------------------------------------------------------
#   2.1) Average annual TP (Table 2a)   ---------------------------------------------------------------------------

    # The MSM Model are loaded above, and were created in Model Selection file
    
    # 2.1.1) Preliminary view of model results 
    alc4.msm                              # model summary
    pmatrix.msm(alc4.msm, t=1, ci="norm") # Transition probabilities at year = t; covariates set to their mean value
    hazard.msm(alc4.msm)                  # Hazard ratios for transition
    
    
    
    # 2.1.2) Function to extract Annual Transition Probabilities (aTP) and correct CI to original sample size 
    # Final, adjusted model
    predicted_TP(model=alc4.msm, year=1) %>% 
      write_csv(paste0(output, "Table 2a - AlcUse Annual TP.csv")) %>% 
      kable()
    
    
    # Crude (unadjusted) model
    predicted_TP(model=alc4_crude.msm, year=1) %>% 
      write_csv(paste0(output, "Table S5a - AlcUse Unadjusted Annual TP.csv")) %>%
      kable()




#   2.2) Generate TP and baseline population ---------------------------------------------------------------

# Function to extract annual TP (Supplement 1)
aTP_alc4 <- predicted_TP_covs (alc4.msm, 1) %>%
  mutate(From = recode(From,"State 1" = "Non-drinker",  # Rename states
                            "State 2" = "Category I",
                            "State 3" = "Category II",
                            "State 4" = "Category III"),
         To = recode(To, "State.1" = "Non-drinker",
                        "State.2" = "Category I",
                        "State.3" = "Category II",
                        "State.4" = "Category III")) %>% 
  write_csv(paste0(output, "Supplement 1 - AlcUse Annual Transition Probabilities.csv")) %>%    # Save a copy of the TP
  mutate(cat = paste(sex, age_cat, edu, race, From, sep="_")) %>% 
  select(cat, To, Probability) %>% 
  group_by(cat) %>% 
  mutate(cumsum = cumsum(Probability)) %>%
  ungroup() 



# Load and set up the initial population (based on NESARC wave 1)
AlcUse_basepop <- nesarc_expanded %>%
  select(idnum, wave, age, age7, female_w1, race_w1, edu3, alc4.factor) %>%  
  rename( sex = female_w1, 
          race = race_w1) %>%
  pivot_wider(names_from="wave", values_from=c("alc4.factor", "age", "age7", "edu3")) %>%  
  mutate (predicted_cat = alc4.factor_1,
          AlcUse_1 = alc4.factor_1,
          AlcUse_2 = alc4.factor_2,
          age = age_1,
          age7 = age7_1,
          edu = edu3_1) %>%
  select(idnum, AlcUse_1, AlcUse_2, sex, age, age7, edu, race, predicted_cat)


# Load and set up the NESARC 3 population
AlcUse_observed_pop <- nesarc3_expanded %>%  # starts with observed data
  mutate(observed_cat = alc4.factor,
         sex = female.factor,
         edu = edu3.factor,
         race = race.factor) %>% 
  select(idnum, observed_cat, sex, age, age7, edu, race)



#   2.3) Simulate and compare to NESARC 3 (Table S7) ----------------------------------------------------------

# Function to simulate population at 11 years follow-up
AlcUse_predicted_pop <- simulate_pop(AlcUse_basepop, aTP_alc4, transition_alc4, 11)


# Function to compare proportions
compare_pct(AlcUse_predicted_pop, AlcUse_observed_pop) 

all  <- compare_pct(AlcUse_predicted_pop, AlcUse_observed_pop) 
sex  <- compare_pct(AlcUse_predicted_pop, AlcUse_observed_pop, sex)
edu  <- compare_pct(AlcUse_predicted_pop, AlcUse_observed_pop, edu)
race <- compare_pct(AlcUse_predicted_pop, AlcUse_observed_pop, race)
age  <- compare_pct(AlcUse_predicted_pop, AlcUse_observed_pop, age7)
    
rbind(all, sex, edu, race, age) %>% write_csv(paste0(output, "Table S7.csv"))



#   2.4) Plot TP over time  --------------------------------------------------------------------

# TP OVERALL ***********************************************************************************
TP_noCovs <- predicted_TP_noCovs(alc4.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker",  "State 2" = "Category I", "State 3" = "Category II", "State 4" = "Category III"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Category I", "State.3" = "Category II", "State.4" = "Category III"))

TP_noCovs_init <- TP_noCovs %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_noCovs_init, TP_noCovs) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker"), 
    To = fct_relevel(To, "Non-drinker")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure 1a - TP of AlcUse over time.tiff"), dpi=600, width=12, height = 4)



# View TP for non-drinkers 
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Non-drinker" & To == "Non-drinker" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0),
          inverse = 100 - Probability) 


# View TP for Category I drinkers 
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Category I" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0)) 


# View TP for Category II drinkers 
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Category II" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0)) 

# View TP for Category III drinkers 
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Category III" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0)) 


# TP by AGE ***********************************************************************************
TP_by_age <- predicted_TP_age(alc4.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker",  "State 2" = "Category I", "State 3" = "Category II", "State 4" = "Category III"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Category I", "State.3" = "Category II", "State.4" = "Category III"),
    age = paste(age_cat, "years"))

TP_by_age_init <- TP_by_age %>% 
  filter(year==1) %>%
  mutate(year = 0,
        Probability = ifelse(From == To, 1, 0),
        lowerCI = ifelse(From == To, 1, 0),
        upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_age_init, TP_by_age) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker"), 
    To = fct_relevel(To, "Non-drinker")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(age~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
  ggsave(paste0(output, "Figure S2 - AlcUse TP over time, stratified by age.tiff"), dpi=600, width=8.5, height = 11, units="in")

  
  

  
# TP by EDUCATION ***********************************************************************************
TP_by_edu <- predicted_TP_edu(alc4.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker",  "State 2" = "Category I", "State 3" = "Category II", "State 4" = "Category III"),
    To =  recode(To, "State.1" = "Non-drinker", "State.2" = "Category I", "State.3" = "Category II", "State.4" = "Category III"),
    edu = recode(edu, "Low" = "Low education", "Med" = "Medium education", "High" = "High education"))

TP_by_edu_init <- TP_by_edu %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_edu_init, TP_by_edu) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker"), 
    To = fct_relevel(To, "Non-drinker"),
    edu = fct_relevel(edu, "Low education", "Medium education", "High education")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(edu~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure S3 - AlcUse TP over time, stratified by edu.tiff"), dpi=600, width=8.5, height = 6, units="in")

  


# TP by RACE ***********************************************************************************
TP_by_race <- predicted_TP_race(alc4.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker",  "State 2" = "Category I", "State 3" = "Category II", "State 4" = "Category III"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Category I", "State.3" = "Category II", "State.4" = "Category III"))

TP_by_race_init <- TP_by_race %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_race_init, TP_by_race) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker"), 
    To = fct_relevel(To, "Non-drinker"),
    race = fct_relevel(race, "White, non-Hispanic")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(race~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure S4 - AlcUse TP over time, stratified by race.tiff"), dpi=600, width=8.5, height = 9, units="in")




# TP by SEX ***********************************************************************************
TP_by_sex <- predicted_TP_sex(alc4.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker",  "State 2" = "Category I", "State 3" = "Category II", "State 4" = "Category III"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Category I", "State.3" = "Category II", "State.4" = "Category III"))

TP_by_sex_init <- TP_by_sex %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_sex_init, TP_by_sex) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker"), 
    To = fct_relevel(To, "Non-drinker")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(sex~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure S5 - AlcUse TP over time, stratified by sex.tiff"), dpi=600, width=8.5, height = 4, units="in")

  
#   2.5) Hazard ratios  (Table 4) -------------------------------------------------------------------------------
# Function to extract HR results, rearrange, and correct CI to original sample size 
HR_table(alc4.msm)  %>%
  rename( "Non-drinker->Category I"   = "State 1 - State 2",
          "Category I->Category II"   = "State 2 - State 3",
          "Category II->Category III" = "State 3 - State 4",
          "Category III->Category II" = "State 4 - State 3",
          "Category II->Category I"   = "State 3 - State 2",
          "Category I->Non-drinker"   = "State 2 - State 1") %>%
  select(Variable, "Non-drinker->Category I", "Category I->Category II", "Category II->Category III", 
                    "Category III->Category II", "Category II->Category I", "Category I->Non-drinker") %>% 
  write_csv(paste0(output, "Table 4 - AlcUse Hazard Ratios.csv")) %>% # save results for paper
  kableone()


#   2.6) Transition frequencies (Tables S3a, S4a) --------------------------------------------------------------

statetable.msm(alc4, idnum, data=nesarc) %>% as.data.frame.array() %>% 
  write_csv(file=paste0(output, "Table S3a - AlcUse Transition Frequecies.csv")) %>%
  kable()

statetable.msm(alc4, idnum, data=nesarc_expanded) %>% as.data.frame.array() %>% 
  write_csv(file=paste0(output, "Table S4a - AlcUse Transitio Frequencies, expanded data.csv")) %>%
  kable()


# 3) HEAVY EPISODIC DRINKING  ----------------------------------------------------------------------------------
#   3.1) Average annual TP (Table 2b)   ---------------------------------------------------------------------------

    # The MSM Model are loaded above, and were created in Model Selection file
    
    # 3.1.1) Preliminary view of model results 
    hed.msm                              # model summary
    pmatrix.msm(hed.msm, t=1, ci="norm") # Transition probabilities at year = t; covariates set to their mean value
    hazard.msm(hed.msm)                  # Hazard ratios for transition


    # 3.1.2) Function to extract Annual Transition Probabilities (aTP) and correct CI to original sample size 
    # Final, adjusted model
    predicted_TP(model=hed.msm, year=1) %>%
      write_csv(paste0(output, "Table 2b - HED Annual TP.csv")) %>% 
      kable()

    # Crude (unadjusted) model
    predicted_TP(model=hed_crude.msm, year=1) %>% 
      write_csv(paste0(output, "Table S5b - HED Unadjusted Annual TP.csv")) %>% 
      kable()


#   3.2) Generate TP and baseline population ---------------------------------------------------------------
    
# Function to extract and format the annual TP
aTP_hed <- predicted_TP_covs (hed.msm, 1) %>%
  mutate(From = recode(From,"State 1" = "Non-drinker",  # Rename states
                            "State 2" = "Drinker, no HED",
                            "State 3" = "Occasional HED",
                            "State 4" = "Monthly HED",
                            "State 5" = "Weekly HED"),
          To = recode(To, "State.1" = "Non-drinker",
                          "State.2" = "Drinker, no HED",
                          "State.3" = "Occasional HED",
                          "State.4" = "Monthly HED",
                          "State.5" = "Weekly HED")) %>%
  write_csv(paste0(output, "Supplement 2 - HED Annual Transition Probabilities.csv")) %>%   # Save TP
  mutate(cat = paste(sex, age_cat, edu, race, From, sep="_")) %>% 
  select(cat, To, Probability) %>% 
  group_by(cat) %>% 
  mutate(cumsum = cumsum(Probability)) %>%
  ungroup()

# Load and set up the initial population (based on NESARC wave 1)
hed_basepop <- nesarc_expanded %>%
  select(idnum, wave, age, age7, female_w1, race_w1, edu3, hed.factor) %>%
  rename(sex = female_w1, 
    race = race_w1) %>%
  pivot_wider(names_from="wave", values_from=c("hed.factor", "age", "age7", "edu3")) %>%
  mutate (predicted_cat = hed.factor_1,
    hed_1 = hed.factor_1,
    hed_2 = hed.factor_2,
    age = age_1,
    age7 = age7_1,
    edu = edu3_1) %>%
  select(idnum, hed_1, hed_2, sex, age, age7, edu, race, predicted_cat)

# 2.2.4) Load and set up the NESARC 3 population
hed_observed_pop <- nesarc3_expanded %>%  # starts with observed data
  mutate( observed_cat = hed.factor,
    sex = female.factor,
    edu = edu3.factor,
    race = race.factor) %>% 
  select(idnum, observed_cat, sex, age, age7, edu, race)


#   3.3) Simulate and compare to NESARC 3 (Table S8) ----------------------------------------------------------

# Function to simulate population at 11 years follow-up
hed_predicted_pop <- simulate_pop(hed_basepop, aTP_hed, transition_hed, 11)


# Function to compare proportions
compare_pct(hed_predicted_pop, hed_observed_pop) %>% kable()

all  <- compare_pct(hed_predicted_pop, hed_observed_pop)
sex  <- compare_pct(hed_predicted_pop, hed_observed_pop, sex)
edu  <- compare_pct(hed_predicted_pop, hed_observed_pop, edu)
race <- compare_pct(hed_predicted_pop, hed_observed_pop, race)
age  <- compare_pct(hed_predicted_pop, hed_observed_pop, age7)

rbind(all, sex, edu, race, age) %>% write_csv(paste0(output, "Table S8.csv"))


#   3.4) Plot TP over time  --------------------------------------------------------------------

# TP OVERALL ***********************************************************************************
TP_noCovs <- predicted_TP_noCovs(hed.msm, 10) %>% 
  mutate( From = recode(From,"State 1" = "Non-drinker", "State 2" = "Drinker, no HED", "State 3" = "Occasional HED", "State 4" = "Monthly HED", "State 5" = "Weekly HED"),
          To = recode(To, "State.1" = "Non-drinker", "State.2" = "Drinker, no HED", "State.3" = "Occasional HED", "State.4" = "Monthly HED", "State.5" = "Weekly HED"))

TP_noCovs_init <- TP_noCovs %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_noCovs_init, TP_noCovs) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker", "Initial state: Drinker, no HED", "Initial state: Occasional HED", "Initial state: Monthly HED", "Initial state: Weekly HED"), 
    To = fct_relevel(To, "Non-drinker", "Drinker, no HED", "Occasional HED", "Monthly HED", "Weekly HED")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure 1b - TP of HED over time.tiff"), dpi=600, width=12, height = 3.5)



# View TP for non-drinkers 
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Non-drinker" & To %in% c("Non-drinker", "Drinker, no HED") & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 1),
    inverse = 100 - Probability) 

# View TP for Drinker, no HED 
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Drinker, no HED" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0)) 

# View TP for Occasional HED
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Occasional HED" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0)) 

# View TP for Monthly HED
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Monthly HED" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0)) 

# View TP for Weekly HED
TP_noCovs %>% 
  select(-lowerCI, - upperCI) %>% 
  filter (From == "Weekly HED" & year %in% c(1,5)) %>%
  mutate (Probability = round((Probability *100), 0)) 



# TP by AGE ***********************************************************************************
TP_by_age <- predicted_TP_age(hed.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker", "State 2" = "Drinker, no HED", "State 3" = "Occasional HED", "State 4" = "Monthly HED", "State 5" = "Weekly HED"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Drinker, no HED", "State.3" = "Occasional HED", "State.4" = "Monthly HED", "State.5" = "Weekly HED"),
    age = paste(age_cat, "years"))

TP_by_age_init <- TP_by_age %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_age_init, TP_by_age) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker", "Initial state: Drinker, no HED", "Initial state: Occasional HED", "Initial state: Monthly HED", "Initial state: Weekly HED"), 
    To = fct_relevel(To, "Non-drinker", "Drinker, no HED", "Occasional HED", "Monthly HED", "Weekly HED")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(age~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure S6 - HED TP over time, stratified by age.tiff"), dpi=600, width=8.5, height = 11, units="in")





# TP by EDUCATION ***********************************************************************************
TP_by_edu <- predicted_TP_edu(hed.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker", "State 2" = "Drinker, no HED", "State 3" = "Occasional HED", "State 4" = "Monthly HED", "State 5" = "Weekly HED"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Drinker, no HED", "State.3" = "Occasional HED", "State.4" = "Monthly HED", "State.5" = "Weekly HED"),
    edu = recode(edu, "Low" = "Low education", "Med" = "Medium education", "High" = "High education"))

TP_by_edu_init <- TP_by_edu %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_edu_init, TP_by_edu) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker", "Initial state: Drinker, no HED", "Initial state: Occasional HED", "Initial state: Monthly HED", "Initial state: Weekly HED"), 
    To = fct_relevel(To, "Non-drinker", "Drinker, no HED", "Occasional HED", "Monthly HED", "Weekly HED"), 
    edu = fct_relevel(edu, "Low education", "Medium education", "High education")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(edu~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure S7 - HED TP over time, stratified by edu.tiff"), dpi=600, width=8.5, height = 6, units="in")




# TP by RACE ***********************************************************************************
TP_by_race <- predicted_TP_race(hed.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker", "State 2" = "Drinker, no HED", "State 3" = "Occasional HED", "State 4" = "Monthly HED", "State 5" = "Weekly HED"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Drinker, no HED", "State.3" = "Occasional HED", "State.4" = "Monthly HED", "State.5" = "Weekly HED"))
    
TP_by_race_init <- TP_by_race %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_race_init, TP_by_race) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker", "Initial state: Drinker, no HED", "Initial state: Occasional HED", "Initial state: Monthly HED", "Initial state: Weekly HED"), 
    To = fct_relevel(To, "Non-drinker", "Drinker, no HED", "Occasional HED", "Monthly HED", "Weekly HED"),
    race = fct_relevel(race, "White, non-Hispanic")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(race~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure S8 - HED TP over time, stratified by race.tiff"), dpi=600, width=8.5, height = 9, units="in")




# TP by SEX ***********************************************************************************
TP_by_sex <- predicted_TP_sex(hed.msm, 10) %>% 
  mutate( 
    From = recode(From,"State 1" = "Non-drinker", "State 2" = "Drinker, no HED", "State 3" = "Occasional HED", "State 4" = "Monthly HED", "State 5" = "Weekly HED"),
    To = recode(To, "State.1" = "Non-drinker", "State.2" = "Drinker, no HED", "State.3" = "Occasional HED", "State.4" = "Monthly HED", "State.5" = "Weekly HED"))
    
TP_by_sex_init <- TP_by_sex %>% 
  filter(year==1) %>%
  mutate(year = 0,
    Probability = ifelse(From == To, 1, 0),
    lowerCI = ifelse(From == To, 1, 0),
    upperCI = ifelse(From == To, 1, 0))

rbind (TP_by_sex_init, TP_by_sex) %>%
  mutate(
    From = paste("Initial state:", From),
    From = fct_relevel(From, "Initial state: Non-drinker", "Initial state: Drinker, no HED", "Initial state: Occasional HED", "Initial state: Monthly HED", "Initial state: Weekly HED"), 
    To = fct_relevel(To, "Non-drinker", "Drinker, no HED", "Occasional HED", "Monthly HED", "Weekly HED")) %>% 
  ggplot(aes(x=year, y=Probability, group=To, color=To, fill=To)) + 
  geom_line(aes(color=To), size=0.75) + geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha=0.2, linetype=0) +
  facet_grid(sex~From) + 
  labs(x = "Time (years)", y="Transition Probability", color="State at Follow-up:", fill="State at Follow-up:") +
  theme(legend.position = "top",
    panel.grid.major=element_line(color="grey90"), 
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA)) + 
  scale_y_continuous(breaks=seq(0, 1, by= .2)) + 
  scale_x_continuous(breaks=seq(0, 10, by= 2))
ggsave(paste0(output, "Figure S9 - HED TP over time, stratified by sex.tiff"), dpi=600, width=8.5, height = 4, units="in")




#   3.5) Hazard ratios  (Table 5) -------------------------------------------------------------------------------

# Function to extract HR results, rearrange, and correct CI to original sample size 
HR_table(hed.msm)  %>%
  rename( "Non-drinker->Drinker, no HED"    = "State 1 - State 2",
          "Drinker, no HED->Non-drinker"    = "State 2 - State 1",
          "Drinker, no HED->Occasional HED" = "State 2 - State 3",
          "Occasional HED->Drinker, no HED" = "State 3 - State 2",
          "Occasional HED->Monthly HED"     = "State 3 - State 4",
          "Monthly HED->Occasional HED"     = "State 4 - State 3",
          "Monthly HED->Weekly HED"         = "State 4 - State 5",
          "Weekly HED->Monthly HED"         = "State 5 - State 4") %>%
  select(Variable,  "Non-drinker->Drinker, no HED", "Drinker, no HED->Occasional HED", "Occasional HED->Monthly HED", "Monthly HED->Weekly HED",
                    "Weekly HED->Monthly HED", "Monthly HED->Occasional HED", "Occasional HED->Drinker, no HED", "Drinker, no HED->Non-drinker") %>%
  write_csv(paste0(output, "Table 5 - HED Hazard Ratios.csv")) %>%  # save results for paper
  kable()




#   3.6) Transition frequencies (Tables S3b, S4b) --------------------------------------------------------------

statetable.msm(hed, idnum, data=nesarc) %>% as.data.frame.array() %>% 
  write_csv(file=paste0(output, "Table S3b - HED Transition Frequecies.csv")) %>%
  kable()

statetable.msm(hed, idnum, data=nesarc_expanded) %>% as.data.frame.array() %>% 
  write_csv(file=paste0(output, "Table S4b - HED Transitio Frequencies, expanded data.csv")) %>%
  kable()



