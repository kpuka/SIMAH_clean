
# SIMAH - NESARC Alcohol Transitions
# Data Extraction (and minor edits)

library(haven)      # Read STATA and SAS data
library(tidyverse)  # data management
library(skimr)      # descriptive statistics
library(survey)     # to accomodate survey weights


# Specify the data and output file locations
data_orig  <- "C:/.../Original data/NESARC I and II/"
data_orig3 <- "C:/.../Original data/NESARC III/"
data_new   <- "C:/.../Processed data/"


# Load data 
nesarc1_orig <- read_dta(paste0(data_orig, "NESARCWave1.dta")) %>% zap_formats() %>% zap_label() %>% zap_labels()
nesarc2_orig <- read_dta(paste0(data_orig, "NESARCWave2.dta")) %>% zap_formats() %>% zap_label() %>% zap_labels()
nesarc3_orig <- read_sas (paste0(data_orig3, "publicfinal_102015.sas7bdat")) %>% zap_formats() %>% zap_label() %>% zap_labels()


# Some variables names are different compared to original codebook; create a file with the available variable names to double check
  # write_csv(as.data.frame(names(nesarc1_orig)), paste0(data_orig, "Wave1_variables.csv"))
  # write_csv(as.data.frame(names(nesarc2_orig)), paste0(data_orig, "Wave2_variables.csv"))


# WAVE 1 -------------------------------------------------------------------------------------------------------------  

nesarc1 <- nesarc1_orig %>%
  
  # Select the variables of interest 
  select(idnum, psu, stratum, weight, CDAY, CMON, CYEAR, 
    sex, olds1q1c, olds1q1d3, olds1q1d5, MARITAL, S1Q6A, olds1q11b, CONSUMER,
    S2AQ4A, S2AQ4B, s2aq4cr, S2AQ4D, S2AQ4E, S2AQ4F, S2AQ4G, coolecf, 
    S2AQ5A, S2AQ5B, s2aq5cr, S2AQ5D, S2AQ5E, S2AQ5F, S2AQ5G, beerecf, 
    S2AQ6A, S2AQ6B, s2aq6cr, S2AQ6D, S2AQ6E, S2AQ6F, S2AQ6G, wineecf, 
    S2AQ7A, S2AQ7B, s2aq7cr, S2AQ7D, S2AQ7E, S2AQ7F, S2AQ7G, liqrecf, S2AQ8E, S2AQ9) %>%
  
  # recode / create new variables
  mutate (
    wave = 1,
    age_diff = 0, # age difference from wave 1
    age = NA,    # placeholder variable - age as NA since it will be extracted from Wave 2
    nw1age = NA, # Placeholder variable - matches  age1 from wave 2
    female = recode(sex, `1` = 0, `2` = 1),
    race = case_when(olds1q1d5==1 & olds1q1d3 ==2 & olds1q1c==2  ~ 1, # white, non-hispanic
                     olds1q1d3==1 & olds1q1c==2 ~ 2, # black, non-hispanic
                     olds1q1c==1 ~ 3,                # Hispanic
                     TRUE ~ 4)) %>%                  # Other, non-hispanic  
  
   # rename variables to align with Wave 2
  rename(
    marital_stat = MARITAL,
    edu = S1Q6A,
    fam_income = olds1q11b,
    drinking_stat = CONSUMER,
    drank5plus_freq = S2AQ8E,
    drank4plus_freq = S2AQ9)
  
    
  # check
  # count(nesarc1, olds1q1d5, olds1q1c, olds1q1d3, race) 
  # count(nesarc1, sex, female)
  
  # Remove extra variables  
  nesarc1 <- select(nesarc1, -olds1q1d5, -olds1q1d3, -olds1q1c, -sex)
  
# Wave 2  -------------------------------------------------------------------------------------------------------------  

nesarc2 <- nesarc2_orig %>%
    
  # Select the variables of interest 
  select(idnum, w2psu, w2stratum, w2weight, w2cday, w2cmon, w2cyear, 
    nw1age, w2AGE, w2SEX, w2ethrace, W2MARITAL, w2s1q15ar, w2s1q19br, w2CONSUMER, 
    w2S2AQ5A, w2S2AQ5B, w2s2aq5cr, w2S2AQ5D, w2S2AQ5E, w2S2AQ5F, w2S2AQ5G, w2coolecf, 
    w2S2AQ6A, w2S2AQ6B, w2s2aq6cr, w2S2AQ6D, w2S2AQ6E, w2S2AQ6F, w2S2AQ6G, w2beerecf, 
    w2S2AQ7A, w2S2AQ7B, w2s2aq7cr, w2S2AQ7D, w2S2AQ7E, w2S2AQ7F, w2S2AQ7G, W2WINEECF, 
    w2S2AQ8A, w2S2AQ8B, w2s2aq8cr, w2S2AQ8D, w2S2AQ8E, w2S2AQ8F, w2S2AQ8G, w2liqrecf, w2S2AQ4F, w2S2AQ4E) %>%
    
  # recode / create new variables
  mutate (
    wave = 2,
    age_diff = w2AGE - nw1age, # follow-up years
    female = recode(w2SEX, `1` = 0, `2` = 1),
    race = recode(w2ethrace, `1`=1, `2`=2, `3`=4, `4`=4, `5`=3)) %>%
  
  # rename variables to align with Wave 1
  rename(
    psu = w2psu, stratum = w2stratum, weight = w2weight, CDAY = w2cday, CMON = w2cmon, CYEAR = w2cyear, 
    age = w2AGE, marital_stat = W2MARITAL, edu = w2s1q15ar,
    fam_income = w2s1q19br, drinking_stat = w2CONSUMER, 
    S2AQ4A = w2S2AQ5A,  S2AQ4B = w2S2AQ5B,  s2aq4cr = w2s2aq5cr,  S2AQ4D = w2S2AQ5D,  S2AQ4E = w2S2AQ5E,  S2AQ4F = w2S2AQ5F,  S2AQ4G = w2S2AQ5G,  coolecf = w2coolecf,  
    S2AQ5A = w2S2AQ6A,  S2AQ5B = w2S2AQ6B,  s2aq5cr = w2s2aq6cr,  S2AQ5D = w2S2AQ6D,  S2AQ5E = w2S2AQ6E,  S2AQ5F = w2S2AQ6F,  S2AQ5G = w2S2AQ6G,  beerecf = w2beerecf,  
    S2AQ6A = w2S2AQ7A,  S2AQ6B = w2S2AQ7B,  s2aq6cr = w2s2aq7cr,  S2AQ6D = w2S2AQ7D,  S2AQ6E = w2S2AQ7E,  S2AQ6F = w2S2AQ7F,  S2AQ6G = w2S2AQ7G,  wineecf = W2WINEECF,  
    S2AQ7A = w2S2AQ8A,  S2AQ7B = w2S2AQ8B,  s2aq7cr = w2s2aq8cr,  S2AQ7D = w2S2AQ8D,  S2AQ7E = w2S2AQ8E,  S2AQ7F = w2S2AQ8F,  S2AQ7G = w2S2AQ8G,  liqrecf = w2liqrecf,  
    S2AQ8E = w2S2AQ4F,  S2AQ9 = w2S2AQ4E, drank5plus_freq = w2S2AQ4F, drank4plus_freq = w2S2AQ4E)

  # check
  # count(nesarc2, w2SEX, female) 
  # count(nesarc2, w2ethrace, race) 
 
  
  # remove extra variables
  nesarc2 <- select(nesarc2, -w2ethrace, -w2SEX)
  

# Merge and dave Wave 1 & 2 Data -----------------------------------------------------------------------------------------------------------
nesarc_raw <- rbind(nesarc1, nesarc2) %>%
    # Fill in age1 data from the nw1age variable from wave2
    arrange (idnum, wave) %>%  
    group_by(idnum) %>%     
      mutate (age_w1 = lead (nw1age, n=1)) %>% 
    ungroup() %>%
    mutate (age = if_else(wave==1, age_w1, age)) %>%
    select(-nw1age, -age_w1)
  

# Save data Wave 1 & 2
saveRDS(nesarc_raw, paste0(data_new, "nesarc_raw.rds"))

  
  
# WAVE 3 ------------------------------------------------------------------------------------------------------------------------

nesarc3 <- nesarc3_orig %>%
  # specify variables to keep
  select( idnum, AUDWEIGHT, NAGE, NSEX, nethrace, NMARITAL, NEDUC, 
          n1q19br, nconsumer, n2aq4h, n2aq4f,
          n2aq5a, n2aq5b, n2aq5cr, n2aq5d, n2aq5e, n2aq5f, n2aq5g, ncoolecf,
          n2aq6a, n2aq6b, n2aq6cr, n2aq6d, n2aq6e, n2aq6f, n2aq6g, nbeerecf, 
          n2aq7a, n2aq7b, n2aq7cr, n2aq7d, n2aq7e, n2aq7f, n2aq7g, nwineecf, 
          n2aq8a, n2aq8b, n2aq8cr, n2aq8d, n2aq8e, n2aq8f, n2aq8g, nliqrecf) 
  
  
# Save data Wave 3
saveRDS(nesarc3, paste0(data_new, "nesarc3_raw.rds"))
