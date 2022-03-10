
# SIMAH - NESARC Alcohol Transitions
# Data Management

library(tidyverse)       # data management
library(skimr)           # descriptive statistics
library(survey)          # to accomodate survey weights
library(lubridate)       # to work with dates/calculate follow-up time
library(splitstackshape) # To replicate data based on sampling weight


# Specify the data file locations
data <- "C:/.../Processed data/"


# NESARC WAVE 1 and 2 --------------------------------------------------------------------------------------------------------
# Edit data - recode and recategorize variables
nesarc <- readRDS(paste0(data, "nesarc_raw.rds")) %>%

  # Re-code variables
  mutate(
    married = recode(marital_stat, `1`= 1, `2`=1, `3`=0, `4`=0, `5`=0, `6`=0),
    edu3 = case_when(edu %in% c(1,2,3,4,5,6,7,8,9) ~ 1,
                     edu %in% c(10,11) ~ 2,
                     edu %in% c(12,13,14) ~ 3),
    income3 = case_when(fam_income %in% c(1,2,3,4,5,6) ~ 1, # split into evenly sized groups based on Wave 1 data
                        fam_income %in% c(7,8,9,10,11) ~ 2,
                        fam_income %in% c(12,13,14,15,16,17,18,19,20,21) ~ 3),
    age3 = cut(age, breaks = c(-Inf, 30, 50, Inf), labels = c("18-29", "30-49", "50+"), right=FALSE),
    age7 = factor(case_when(age < 21 ~ "18-20",
                            age >= 21 & age <26 ~ "21-25",
                            age >= 26 & age <30 ~ "26-29",
                            age >= 30 & age <40 ~ "30-39",
                            age >= 40 & age <50 ~ "40-49",
                            age >= 50 & age <65 ~ "50-64",
                            age >= 65 ~ "65+")),
    age7 = relevel(age7, ref="65+"), 

    # Calculate follow-up time and baseline version of age, sex, and edu
    intv_date = make_date(year=CYEAR, month=CMON, day=CDAY)) %>%
  
  arrange(idnum, wave) %>%
  group_by(idnum) %>%
      mutate (prev_intv_date = lag(intv_date, 1),
              female_w1 = lag(female, 1),
              age_w1 = lag(age, 1),
              race_w1 = lag(race, 1),
              drinking_stat_wave1 = lag(drinking_stat, 1),
              weight_wave2 = lead(weight, 1)) %>%
  ungroup() %>%
  mutate (years = as.numeric(prev_intv_date %--% intv_date /dyears(1)),
          years = if_else(wave==1, 0, years),
          female_w1 = if_else(wave==1,female, female_w1),
          age_w1 = if_else(wave==1,age, age_w1),
          race_w1 = if_else(wave==1, race, race_w1),
          drinking_stat_wave1 = if_else(wave==1, drinking_stat, drinking_stat_wave1),
          weight_wave2 = if_else(wave==2, weight, weight_wave2),

    # Calculating alcohol intake:
        # Recode "99" (missing) to NA
        across(c(S2AQ4B, s2aq4cr, S2AQ4D, S2AQ4E, S2AQ4F, S2AQ4G,
                 S2AQ5B, s2aq5cr, S2AQ5D, S2AQ5E, S2AQ5F, S2AQ5G,
                 S2AQ6B, s2aq6cr, S2AQ6D, S2AQ6E, S2AQ6F, S2AQ6G,
                 S2AQ7B, s2aq7cr, S2AQ7D, S2AQ7E, S2AQ7F, S2AQ7G),
          ~ recode(.x, `99`=NA_real_)),
    
        # Recode frequencies to number drinking days (variable: avg # of drinking days)
        across(c(S2AQ4B, S2AQ4F, S2AQ4G, S2AQ5B, S2AQ5F, S2AQ5G, S2AQ6B, S2AQ6F, S2AQ6G, S2AQ7B, S2AQ7F, S2AQ7G), 
          ~ recode(.x, `1`=365, `2`=273, `3`=182, `4`=104, `5`=52, `6`=30, `7`=12, `8`=9, `9`=4.5, `10`=1.5, `11`=0), .names = "{.col}_days"),
        

        # Calculate total drinks per year (as per the manual)
        coolers_yearly = ifelse(S2AQ4E <= 5, 
                                  (S2AQ4D * (S2AQ4B_days - S2AQ4F_days)) + (S2AQ4E * S2AQ4F_days), # quantity if largest drinks <=5
                                  (S2AQ4D * (S2AQ4B_days - S2AQ4G_days)) + ((S2AQ4G_days - S2AQ4F_days)*(exp((log(pmax(5, S2AQ4D)) + log(S2AQ4E - 1)) /2))) + (S2AQ4E*S2AQ4F_days)), # quantity if largest drinks >5
        
        beers_yearly = ifelse(S2AQ5E <= 5, 
                                  (S2AQ5D * (S2AQ5B_days - S2AQ5F_days)) + (S2AQ5E * S2AQ5F_days),
                                  (S2AQ5D * (S2AQ5B_days - S2AQ5G_days)) + ((S2AQ5G_days - S2AQ5F_days)*(exp((log(pmax(5, S2AQ5D)) + log(S2AQ5E - 1)) /2))) + (S2AQ5E*S2AQ5F_days)),
        
        wine_yearly = ifelse(S2AQ6E <= 5, 
                                  (S2AQ6D * (S2AQ6B_days - S2AQ6F_days)) + (S2AQ6E * S2AQ6F_days), 
                                  (S2AQ6D * (S2AQ6B_days - S2AQ6G_days)) + ((S2AQ6G_days - S2AQ6F_days)*(exp((log(pmax(5, S2AQ6D)) + log(S2AQ6E - 1)) /2))) + (S2AQ6E*S2AQ6F_days)),
        
        liquor_yearly = ifelse(S2AQ7E <= 5, 
                              (S2AQ7D * (S2AQ7B_days - S2AQ7F_days)) + (S2AQ7E * S2AQ7F_days), 
                              (S2AQ7D * (S2AQ7B_days - S2AQ7G_days)) + ((S2AQ7G_days - S2AQ7F_days)*(exp((log(pmax(5, S2AQ7D)) + log(S2AQ7E - 1)) /2))) + (S2AQ7E*S2AQ7F_days)), 
                        
        # Calculate daily ethnanol intake 
        coolers_daily_oz = (coolers_yearly * (s2aq4cr * coolecf))/365,
        beers_daily_oz = (beers_yearly * (s2aq5cr * beerecf))/365,
        wine_daily_oz = (wine_yearly * (s2aq6cr * wineecf))/365,
        liquor_daily_oz = (liquor_yearly * (s2aq7cr * liqrecf))/365) %>%

    mutate(
        #calculate total alcohol ounces per day 
        alc_daily_oz = rowSums(select(., coolers_daily_oz, beers_daily_oz, wine_daily_oz, liquor_daily_oz), na.rm = TRUE), # those with only NAs get a value of 0
        
        # Code as NA those who reported drinking a cooler, beer, wine, or liquor but their alc_daily_oz = 0
        alc_daily_oz = if_else(alc_daily_oz==0 & (S2AQ4A==1 | S2AQ5A==1 | S2AQ6A==1 | S2AQ7A==1), NA_real_, alc_daily_oz),
        
        # Code as NA those who with 'unknown' in at least one of cooler, beer, wine, or liquor and 'No' or 'unknown' to all the others 
        alc_daily_oz = if_else(alc_daily_oz==0 & (S2AQ4A==9 & S2AQ5A%in%c(2,9) & S2AQ6A%in%c(2,9) & S2AQ7A%in%c(2,9)), NA_real_, alc_daily_oz),
        alc_daily_oz = if_else(alc_daily_oz==0 & (S2AQ4A%in%c(2,9) & S2AQ5A==9 & S2AQ6A%in%c(2,9) & S2AQ7A%in%c(2,9)), NA_real_, alc_daily_oz),
        alc_daily_oz = if_else(alc_daily_oz==0 & (S2AQ4A%in%c(2,9) & S2AQ5A%in%c(2,9) & S2AQ6A==9 & S2AQ7A%in%c(2,9)), NA_real_, alc_daily_oz),
        alc_daily_oz = if_else(alc_daily_oz==0 & (S2AQ4A%in%c(2,9) & S2AQ5A%in%c(2,9) & S2AQ6A%in%c(2,9) & S2AQ7A==9), NA_real_, alc_daily_oz),
        
        # Impute 0 ounces daily for non-drinkers
        alc_daily_oz = if_else(drinking_stat%in%c(2,3), 0, alc_daily_oz), 
      
      # Recategorize daily alcohol use using other units
      alc_daily_g = alc_daily_oz * 28.3495,   # Convert daily ounces to grams  
      alc_daily_drinks = alc_daily_oz  / 0.60, # Coverty to daily # of drinks, assuming 0.60oz per drink, as per NESARC guidelines
      

      # Categorize alcohol use as per NESARC guidelines
      alc4_nesarc = case_when(
        # Men:
        female==0 & alc_daily_oz==0 ~ 1,                         # non-drinkers
        female==0 & alc_daily_oz>0 & alc_daily_oz <= 0.257 ~ 2,  # light drinker, i.e., 3 or fewer drinks per week
        female==0 & alc_daily_oz>0.257 & alc_daily_oz <=1.2 ~ 3, # moderate drinker, i.e., 3 to 14 drinks per week
        female==0 & alc_daily_oz>1.2 ~ 4,                        # heavier drinker, i.e., more than 2 drinks
          
        # Women:
        female==1 & alc_daily_oz==0 ~ 1,                         # non-drinkers
        female==1 & alc_daily_oz>0 & alc_daily_oz <= 0.257 ~ 2,  # light drinker, i.e., 3 or fewer drinks per week
        female==1 & alc_daily_oz>0.257 & alc_daily_oz <=0.6 ~ 3, # moderate drinker, i.e., 3 to 7 drinks per week
        female==1 & alc_daily_oz>0.6 ~ 4),                       # heavier drinker, i.e., more than 1 drink
        
      
      # Categorize alcohol use as per SIMAH protocol
      alc6 = case_when(
        # Men & Women:
        wave==1 & drinking_stat==3 ~ 1,                          # lifetime abstinence
        wave==2 & drinking_stat_wave1==3 & drinking_stat==3 ~ 1, # lifetime abstinence
        
        wave==2 & drinking_stat_wave1!=3 & drinking_stat==2 ~ 2, # former drinker
        wave==2 & drinking_stat_wave1!=3 & drinking_stat==3 ~ 2, # former drinker
        
        drinking_stat==2 ~ 2,                                # former drinker 
        drinking_stat==1 & alc_daily_g==0 ~ 2,               # former drinker (indicated they drink, but had 0 grams of alcohol)

        # Men
        female==0 & alc_daily_g >0 & alc_daily_g <=40 ~ 3,   # low risk
        female==0 & alc_daily_g >40 & alc_daily_g <=60 ~ 4,  # medium risk
        female==0 & alc_daily_g >60 & alc_daily_g <=100 ~ 5, # high risk
        female==0 & alc_daily_g >100 ~ 6,                    # very high risk
        
        # Women:
        female==1 & alc_daily_g >0 & alc_daily_g <=20 ~ 3,   # low risk
        female==1 & alc_daily_g >20 & alc_daily_g <=40 ~ 4,  # medium risk
        female==1 & alc_daily_g >40 & alc_daily_g <=60 ~ 5,  # high risk
        female==1 & alc_daily_g >60 ~ 6),                    # very high risk
      
      alc5 = recode(alc6, `6`=5),       # merge high-risk and very-high-risk categories
      alc4 = recode(alc6, `1`=1, `2`=1, `3`=2, `4`=3, `5`=4, `6`=4), # merge absteiners/former and merge high/very-high categories
      alc5_v2 = recode(alc6, `1`=1, `2`=1, `3`=2, `4`=3, `5`=4, `6`=5), # merge absteiners/former categories
         
    # Calculate heavy episodic drinking (HED)
    hed = case_when(
                    # Men:
                    female==0 & drank5plus_freq %in% c(1,2,3,4,5) ~ 5, # HED >= 1/week
                    female==0 & drank5plus_freq %in% c(6,7) ~ 4,       # HED >= 1/month but < 1/week
                    female==0 & drank5plus_freq %in% c(8,9,10) ~ 3,    # HED < 1 / month
                    female==0 & drank5plus_freq %in% c(11) ~ 2,        # No HED in last year
                    
                    # Women:
                    female==1 & drank4plus_freq %in% c(1,2,3,4,5) ~ 5, # HED >= 1/week
                    female==1 & drank4plus_freq %in% c(6,7) ~ 4,       # HED >= 1/month but < 1/week
                    female==1 & drank4plus_freq %in% c(8,9,10) ~ 3,    # HED < 1 / month
                    female==1 & drank4plus_freq %in% c(11) ~ 2),       # No HED in last year
    
    hed = ifelse(alc5 %in% c(1,2), 1, hed)) #Non-drinker


  # Check  
  # select(nesarc, idnum, wave, intv_date, prev_intv_date, years, age_diff) %>% view()
  # select(nesarc, idnum, wave, female, female_w1) %>% view()
  # select(nesarc, idnum, wave, race, race_w1) %>% view()
  # count(nesarc, marital_stat, married)
  # count(nesarc, edu, edu3)
  # view(count(nesarc, fam_income, income3))
  # count(nesarc, S2AQ4B, S2AQ4B_days)
  # count(nesarc, S2AQ5B, S2AQ5B_days)
  # select(nesarc, coolers_yearly,beers_yearly, wine_yearly, liquor_yearly) %>% skim()
  #     select(nesarc, idnum, wave, drinking_stat, S2AQ4A, S2AQ4D, S2AQ4B_days, S2AQ4F_days, S2AQ4E, S2AQ4G_days, coolers_yearly) %>% write.csv(paste0(data_orig, "check_coolers.csv", na=""))
  # select(nesarc, coolers_daily_oz, beers_daily_oz, wine_daily_oz, liquor_daily_oz, alc_daily_oz, alc_daily_g) %>% skim()
  #     select(nesarc, idnum, wave, drinking_stat, alc_daily_g, alc_daily_oz, coolers_daily_oz, beers_daily_oz, wine_daily_oz, liquor_daily_oz, 
  #       S2AQ4A, S2AQ4B, s2aq4cr, S2AQ4D, S2AQ4E, S2AQ4F, S2AQ4G, coolecf, 
  #       S2AQ5A, S2AQ5B, s2aq5cr, S2AQ5D, S2AQ5E, S2AQ5F, S2AQ5G, beerecf, 
  #       S2AQ6A, S2AQ6B, s2aq6cr, S2AQ6D, S2AQ6E, S2AQ6F, S2AQ6G, wineecf, 
  #       S2AQ7A, S2AQ7B, s2aq7cr, S2AQ7D, S2AQ7E, S2AQ7F, S2AQ7G, liqrecf) %>% write.csv(paste0(data_orig, "check_alc.csv", na = ""))
  # count(nesarc, alc4_nesarc)
  # count(nesarc, alc6, alc5, alc4)
  #check if anyone who's previously drank is coded as a lifetime abstainer in wave 2
        # nesarc %>%
        #   select(idnum, wave, alc5, alc_daily_g) %>%
        #   pivot_wider(names_from="wave", values_from = c("alc5", "alc_daily_g")) %>%
        #   filter(alc5_2==1 & alc5_1%in%c(2,3,4,5)) %>%
        #   view()
        # 
        # nesarc %>% 
        #   filter(is.na(alc6)) %>%
        #   select (idnum, wave, alc6, drinking_stat, drinking_stat_wave1, alc_daily_g) %>%
        #   view()
    # filter(nesarc, wave==1) %>% count(hed)



# Select variables of interest and identify those lost to follow-up or with incomplete data
nesarc <- nesarc %>%
  
  # order data by ID, then by wave (timepoint)
  arrange(idnum, wave) %>% 
  
  # Identify the variables to keep
  select(idnum, wave, psu, stratum, weight, weight_wave2,  years, age, age_w1, age3, age7, female, female_w1, race, race_w1, 
    married, edu3, income3, alc_daily_oz, alc_daily_g, alc_daily_drinks, alc4_nesarc, alc6, alc5, alc5_v2, alc4, hed) %>% 
  
  # remove those with missing data 
  group_by(idnum) %>%
    mutate(lost = ifelse(n()>1, 0, 1),                    # Identify those with data at one time point (8,440 observations removed; n=69,306)
           lost = ifelse(any(is.na(alc5)), 1, lost),      # Identify those with missing alcohol data in either time point (424 observations removed; n=68,882)
           lost = ifelse(any(is.na(hed)), 1, lost)) %>%   # Identify those with missing HED in either time point (166 observations removed; n=68,716)
    ungroup() 

  

# label values (first listed category is the reference)

nesarc <- nesarc %>%
  mutate (
    alc6.factor        = factor(alc6, levels=c(1,2,3,4,5,6), labels=c("Lifetime abstainer", "Former drinker", "Category I", "Category II", "Category III", "Category IV")),
    alc5.factor        = factor(alc5, levels=c(1,2,3,4,5), labels=c("Abstainer", "Former", "Category I", "Category II", "Category III")),
    alc5_v2.factor     = factor(alc5_v2, levels=c(1,2,3,4,5), labels=c("Non-drinker", "Category I", "Category II", "Category III", "Category IV")),
    alc4.factor        = factor(alc4, levels=c(1,2,3,4), labels=c("Non-drinker", "Category I", "Category II", "Category III")),
    alc4_nesarc.factor = factor(alc4_nesarc, levels=c(1,2,3, 4), labels=c("Non-drinkers", "Light drinker", "Moderate drinker", "Heavy drinker")),
    hed.factor         = factor(hed, levels=c(1,2,3,4,5), labels=c("Non-drinker", "Drinker, no HED", "Occasional HED", "Monthly HED", "Weekly HED")),
    female.factor      = factor(female, levels=c(0,1), labels=c("Men", "Women")),
    female_w1          = factor(female_w1, levels=c(0,1), labels=c("Men", "Women")), # Sex at wave 1
    race.factor        = factor(race, levels=c(1,2,3,4), labels=c("White, non-Hispanic", "Black, non-Hispanic", "Hispanic", "Other, non-Hispanic")),
    race_w1            = factor(race_w1, levels=c(1,2,3,4), labels=c("White, non-Hispanic", "Black, non-Hispanic", "Hispanic", "Other, non-Hispanic")), # race at wave 1
    married.factor     = factor(married, levels=c(1,0), labels=c("Married/cohab.", "Single")),
    edu3               = factor(edu3, levels=c(3,1,2), labels=c("High", "Low", "Med")),
    income3.factor     = factor(income3, levels=c(3,1,2), labels=c("High", "Low", "Med")),
    wave.factor        = factor(wave, levels=c(1,2), labels=c("Wave 1", "Wave 2")),
    lost.factor        = factor(lost, levels=c(0,1), labels=c("Completed Follow-up", "Lost to Follow-up")))


    
# Create dataframe with complete data at both time points (cc=complete case)
nesarc_cc <- nesarc %>% 
  filter(lost==0)  # Remove those with missing data
      
      
      
# Replicate data (and create unique ID variable) to adjust for sampling weight
nesarc_cc_expanded <- nesarc_cc %>%
  mutate (new_weight = weight_wave2 / 100) %>%  # because original weight variable ranged from 455 to 73,192
  expandRows(., "new_weight") %>%  # replicates data
  
  # Generate unique ID
  group_by(idnum, wave) %>%
    mutate(iter = sprintf("%04d", 1:n()),  # the sprintf("%04d", X) command is used to add leading 0s to make it a variable with 4 digits
           idnum = as.numeric(paste0(idnum, iter))) %>%
  ungroup() %>%
  arrange(idnum, wave) # order data by ID then wave (needed for the MSM model)

# check
# filter(nesarc_cc, wave==1) %>% count(hed)


  
# Save data
saveRDS(nesarc, paste0(data, "nesarc_all.rds"))
saveRDS(nesarc_cc, paste0(data, "nesarc_clean.rds"))
saveRDS(nesarc_cc_expanded, paste0(data, "nesarc_clean_expanded.rds"))

# NESARC WAVE 3 --------------------------------------------------------------------------------------------------------
# Edit data - recode and recategorize variables
nesarc3 <- readRDS(paste0(data, "nesarc3_raw.rds")) %>% 
  
  # Re-code /rename variables
  rename(age = NAGE,
        drinking_stat = nconsumer, 
        drank5plus_freq = n2aq4h,
        drank4plus_freq = n2aq4f,
        weight = AUDWEIGHT) %>% 
  
  mutate(
    female = recode(NSEX, `1` = 0, `2` = 1),
    race = recode(nethrace, `1`=1, `2`=2, `3`=4, `4`=4, `5`=3),
    edu3 = case_when(NEDUC %in% c(1,2,3,4,5,6,7,8,9) ~ 1,
                     NEDUC %in% c(10,11) ~ 2,
                     NEDUC %in% c(12,13,14) ~ 3),
    age3 = cut(age, breaks = c(-Inf, 30, 50, Inf), labels = c("18-29", "30-49", "50+"), right=FALSE),
    age7 = factor(case_when(age < 21 ~ "18-20",
                            age >= 21 & age <26 ~ "21-25",
                            age >= 26 & age <30 ~ "26-29",
                            age >= 30 & age <40 ~ "30-39",
                            age >= 40 & age <50 ~ "40-49",
                            age >= 50 & age <65 ~ "50-64",
                            age >= 65 ~ "65+")),
    married = recode(NMARITAL, `1`= 1, `2`=1, `3`=0, `4`=0, `5`=0, `6`=0),
    
    # Calculating alcohol intake:
    # Recode "99" (missing) to NA
    across( c(n2aq5b, n2aq5cr, n2aq5d, n2aq5e, n2aq5f, n2aq5g,
              n2aq6b, n2aq6cr, n2aq6d, n2aq6e, n2aq6f, n2aq6g,
              n2aq7b, n2aq7cr, n2aq7d, n2aq7e, n2aq7f, n2aq7g,
              n2aq8b, n2aq8cr, n2aq8d, n2aq8e, n2aq8f, n2aq8g), 
        ~ recode(.x, `99`=NA_real_)),
    
    # Recode frequencies to number drinking days (variable: avg # of drinking days)
    across(c(n2aq5b, n2aq5f, n2aq5g,  n2aq6b, n2aq6f, n2aq6g,  n2aq7b, n2aq7f, n2aq7g, n2aq8b, n2aq8f, n2aq8g), 
        ~ recode(.x, `1`=365, `2`=273, `3`=182, `4`=104, `5`=52, `6`=30, `7`=12, `8`=9, `9`=4.5, `10`=1.5, `11`=0), .names = "{.col}_days"),
    
    
    # Calculate total drinks per year (as per the manual)
    coolers_yearly = ifelse(n2aq5e <= 5, 
      (n2aq5d * (n2aq5b_days - n2aq5f_days)) + (n2aq5e * n2aq5f_days), # quantity if largest drinks <=5
      (n2aq5d * (n2aq5b_days - n2aq5g_days)) + ((n2aq5g_days - n2aq5f_days)*(exp((log(pmax(5, n2aq5d)) + log(n2aq5e - 1)) /2))) + (n2aq5e*n2aq5f_days)), # quantity if largest drinks >5
    
    beers_yearly = ifelse(n2aq6e <= 5, 
      (n2aq6d * (n2aq6b_days - n2aq6f_days)) + (n2aq6e * n2aq6f_days), 
      (n2aq6d * (n2aq6b_days - n2aq6g_days)) + ((n2aq6g_days - n2aq6f_days)*(exp((log(pmax(5, n2aq6d)) + log(n2aq6e - 1)) /2))) + (n2aq6e*n2aq6f_days)),
    
    wine_yearly = ifelse(n2aq7e <= 5, 
      (n2aq7d * (n2aq7b_days - n2aq7f_days)) + (n2aq7e * n2aq7f_days), 
      (n2aq7d * (n2aq7b_days - n2aq7g_days)) + ((n2aq7g_days - n2aq7f_days)*(exp((log(pmax(5, n2aq7d)) + log(n2aq7e - 1)) /2))) + (n2aq7e*n2aq7f_days)),
    
    liquor_yearly = ifelse(n2aq8e <= 5, 
      (n2aq8d * (n2aq8b_days - n2aq8f_days)) + (n2aq8e * n2aq8f_days), 
      (n2aq8d * (n2aq8b_days - n2aq8g_days)) + ((n2aq8g_days - n2aq8f_days)*(exp((log(pmax(5, n2aq8d)) + log(n2aq8e - 1)) /2))) + (n2aq8e*n2aq8f_days)),
    
    # Calculate daily ethnanol intake 
    coolers_daily_oz = (coolers_yearly * (n2aq5cr * ncoolecf))/365,
    beers_daily_oz   = (beers_yearly   * (n2aq6cr * nbeerecf))/365,
    wine_daily_oz    = (wine_yearly    * (n2aq7cr * nwineecf))/365,
    liquor_daily_oz  = (liquor_yearly  * (n2aq8cr * nliqrecf))/365) %>%
     
    
  mutate(
    #calculate total alcohol ounces per day 
    alc_daily_oz = rowSums(select(., coolers_daily_oz, beers_daily_oz, wine_daily_oz, liquor_daily_oz), na.rm = TRUE), # those with only NAs get a value of 0
    
    # Code as NA those who reported drinking a cooler, beer, wine, or liquor but their alc_daily_oz = 0
    alc_daily_oz = if_else(alc_daily_oz==0 & (n2aq5a==1 | n2aq6a==1 | n2aq7a==1 | n2aq8a==1), NA_real_, alc_daily_oz),
    
    # Code as NA those who with 'unknown' in at least one of cooler, beer, wine, or liquor and 'No' or 'unknown' to all the others 
    alc_daily_oz = if_else(alc_daily_oz==0 & (n2aq5a==9        & n2aq6a%in%c(2,9) & n2aq7a%in%c(2,9) & n2aq8a%in%c(2,9)), NA_real_, alc_daily_oz),
    alc_daily_oz = if_else(alc_daily_oz==0 & (n2aq5a%in%c(2,9) & n2aq6a==9        & n2aq7a%in%c(2,9) & n2aq8a%in%c(2,9)), NA_real_, alc_daily_oz),
    alc_daily_oz = if_else(alc_daily_oz==0 & (n2aq5a%in%c(2,9) & n2aq6a%in%c(2,9) & n2aq7a==9        & n2aq8a%in%c(2,9)), NA_real_, alc_daily_oz),
    alc_daily_oz = if_else(alc_daily_oz==0 & (n2aq5a%in%c(2,9) & n2aq6a%in%c(2,9) & n2aq7a%in%c(2,9) & n2aq8a==9),        NA_real_, alc_daily_oz),
    
    # Impute 0 ounces daily for non-drinkers
    alc_daily_oz = if_else(drinking_stat%in%c(2,3), 0, alc_daily_oz), 
    
    # Recategorize daily alcohol use using other units
    alc_daily_g = alc_daily_oz * 28.3495,   # Convert daily ounces to grams  
    alc_daily_drinks = alc_daily_oz  / 0.60, # Coverty to daily # of drinks, assuming 0.60oz per drink, as per NESARC guidelines
    
    
     # Categorize alcohol use as per SIMAH protocol
    alc6 = case_when(
      # Men & Women:
      drinking_stat==3 ~ 1,                          # lifetime abstinence
      drinking_stat==2 ~ 2,                          # former drinker
      drinking_stat==1 & alc_daily_g==0 ~ 2,         # former drinker (indicated they drink, but had 0 grams of alcohol)
      
      # Men
      female==0 & alc_daily_g >0 & alc_daily_g <=40 ~ 3,   # low risk
      female==0 & alc_daily_g >40 & alc_daily_g <=60 ~ 4,  # medium risk
      female==0 & alc_daily_g >60 & alc_daily_g <=100 ~ 5, # high risk
      female==0 & alc_daily_g >100 ~ 6,                    # very high risk
      
      # Women:
      female==1 & alc_daily_g >0 & alc_daily_g <=20 ~ 3,   # low risk
      female==1 & alc_daily_g >20 & alc_daily_g <=40 ~ 4,  # medium risk
      female==1 & alc_daily_g >40 & alc_daily_g <=60 ~ 5,  # high risk
      female==1 & alc_daily_g >60 ~ 6),                    # very high risk
    
    alc5 = recode(alc6, `6`=5),       # merge high-risk and very-high-risk categories
    alc4 = recode(alc6, `1`=1, `2`=1, `3`=2, `4`=3, `5`=4, `6`=4), # merge absteiners/former and merge high/very-high categories
    alc5_v2 = recode(alc6, `1`=1, `2`=1, `3`=2, `4`=3, `5`=4, `6`=5), # merge absteiners/former categories
    
    
    # Calculate heavy episodic drinking (HED)
    hed = case_when(
      # Men:
      female==0 & drank5plus_freq %in% c(1,2,3,4,5) ~ 5, # HED >= 1/week
      female==0 & drank5plus_freq %in% c(6,7) ~ 4,       # HED >= 1/month but < 1/week
      female==0 & drank5plus_freq %in% c(8,9,10) ~ 3,    # HED < 1 / month
      female==0 & drank5plus_freq %in% c(11) ~ 2,        # No HED in last year
      
      # Women:
      female==1 & drank4plus_freq %in% c(1,2,3,4,5) ~ 5, # HED >= 1/week
      female==1 & drank4plus_freq %in% c(6,7) ~ 4,       # HED >= 1/month but < 1/week
      female==1 & drank4plus_freq %in% c(8,9,10) ~ 3,    # HED < 1 / month
      female==1 & drank4plus_freq %in% c(11) ~ 2),       # No HED in last year
    
    hed = ifelse(alc5 %in% c(1,2), 1, hed)) #Non-drinker




# Select variables of interest and identify those lost to follow-up or with incomplete data
nesarc3 <- nesarc3 %>%
  
  # Identify the variables to keep
  select(idnum, weight, age, age3, age7, female, race, married, edu3, 
         alc_daily_oz, alc_daily_g, alc_daily_drinks, alc6, alc5, alc5_v2, alc4, hed) %>%
  
  # remove those with missing data 
  filter(!is.na(alc5)) %>%   # 88 observations removed; n=36,221
  filter (!is.na(hed))       # 83 observations removed; n=36,138


# label values (first listed category is the reference)
nesarc3 <- nesarc3 %>%
  mutate (
    alc6.factor    = factor(alc6, levels=c(1,2,3,4,5,6), labels=c("Lifetime abstainer", "Former drinker", "Category I", "Category II", "Category III", "Category IV")),
    alc5.factor    = factor(alc5, levels=c(1,2,3,4,5), labels=c("Abstainer", "Former", "Category I", "Category II", "Category III")),
    alc5_v2.factor = factor(alc5_v2, levels=c(1,2,3,4,5), labels=c("Non-drinker", "Category I", "Category II", "Category III", "Category IV")),
    alc4.factor    = factor(alc4, levels=c(1,2,3,4), labels=c("Non-drinker", "Category I", "Category II", "Category III")),
    hed.factor     = factor(hed, levels=c(1,2,3,4,5), labels=c("Non-drinker", "Drinker, no HED", "Occasional HED", "Monthly HED", "Weekly HED")),
    female.factor  = factor(female, levels=c(0,1), labels=c("Men", "Women")),
    race.factor    = factor(race, levels=c(1,2,3,4), labels=c("White, non-Hispanic", "Black, non-Hispanic", "Hispanic", "Other, non-Hispanic")),
    married.factor = factor(married, levels=c(1,0), labels=c("Married/cohab.", "Single")),
    edu3.factor    = factor(edu3, levels=c(3,1,2), labels=c("High", "Low", "Med")))


nesarc3 %>% select(weight) %>% skim()


# Replicate data (and create unique ID variable) to adjust for sampling weight
nesarc3_expanded <- nesarc3 %>%
  mutate (new_weight = weight / 100) %>%  # because original weight variable ranged from 586 to 49,403
  expandRows(., "new_weight") %>%         # replicates data
  
  # Generate unique ID
  group_by(idnum) %>%
  mutate(iter = sprintf("%04d", 1:n()),  # the sprintf("%04d", X) command is used to add leading 0s to make it a variable with 4 digits
    idnum = as.numeric(paste0(idnum, iter))) %>%
  ungroup() %>%
  arrange(idnum) # order data by ID 

# check
# filter(nesarc_cc, wave==1) %>% count(hed)


# Save data
saveRDS(nesarc3, paste0(data, "nesarc3_clean.rds"))
saveRDS(nesarc3_expanded, paste0(data, "nesarc3_clean_expanded.rds"))
