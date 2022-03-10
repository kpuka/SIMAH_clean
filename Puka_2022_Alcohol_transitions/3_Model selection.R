
# SIMAH - NESARC Alcohol Transitions
# Model selection

# Load packages and data
library(tidyverse)  # data management
library(msm)        # model transition probabilities
memory.limit(size=1e+13)

# Specify the data and output file locations
data    <- "C:/.../Processed data/"  # Location of data
models  <- "C:/.../Models/"          # Location of saved MSM models


# Load data / functions
nesarc_all      <- readRDS(paste0(data, "nesarc_all.rds")) # Contains those with missing data 
nesarc          <- readRDS(paste0(data, "nesarc_clean.rds")) 
nesarc_expanded <- readRDS(paste0(data, "nesarc_clean_expanded.rds")) 


# Data management -------------------------------------------------------------------------------------------------

# Data when age is continuous (ages less than 90 (le90); remove those aged >90 since the exact age is unknown)
nesarc_le90 <- nesarc %>%
  filter (age<90) %>% group_by(idnum) %>% filter(n() > 1) %>% ungroup() %>% 
  mutate (age_sq = age ^ 2, 
    age.c = (age - mean(age)) / sd(age),
    age_sq.c = (age_sq - mean(age_sq)) / sd(age_sq)) 

nesarc_expanded_le90 <- nesarc_expanded %>%
  filter (age<90) %>% group_by(idnum) %>% filter(n() > 1) %>% ungroup() %>% 
  mutate (age_sq = age ^ 2, 
    age.c = (age - mean(age)) / sd(age),
    age_sq.c = (age_sq - mean(age_sq)) / sd(age_sq)) 


# MSM 1: AlcUse (6levels) ---------------------------------------------------------------------------------------------

# Specify allowed transition
# Only allow adjacent transitions, except for transitions back to lifetime abstainers, and abstainer->former drinker
Q <- rbind ( c(0,    0,     0.25, 0,    0,    0),
             c(0,    0,     0.25, 0,    0,    0),
             c(0,    0.25,  0,    0.25, 0,    0),
             c(0,    0,     0.25, 0,    0.25, 0),
             c(0,    0,     0,    0.25, 0,    0.25), 
             c(0,    0,     0,    0,    0.25, 0))

# Specify initial values
Q_ageCat  <- crudeinits.msm(alc6 ~ years, idnum, data=nesarc_expanded, qmatrix=Q)       # When age is categorical
Q_ageCont <- crudeinits.msm(alc6 ~ years, idnum, data=nesarc_expanded_ageC, qmatrix=Q)  # When age is continuous 


# MSM 1: All ages **************************************************************************************
# MSM 1A: Age (3 categories) (all ages)
msm1a <- msm ( alc6 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q_ageCat, 
              center=FALSE, control = list(trace=1, maxit=1000, fnscale = 3000000),
              covariates = ~ female_w1 + age3 + edu3 + race_w1)
        saveRDS(msm1a, paste0(models, "msm1a.RDS")) 


# MSM 1B: Age (7 categories) (all ages)
msm1b <- msm ( alc6 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q_ageCat,
              center=FALSE, control = list(trace=1, maxit=1000, fnscale = 3000000),
              covariates = ~ female_w1 + age7 + edu3 + race_w1)
        saveRDS(msm1b, paste0(models, "msm1b.RDS"))
        
    
# MSM 1: Ages less than 90years (le90) *********************************************************************************      
# MSM 1A: Age (3 categories) 
msm1a_le90 <- msm (alc6 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCat, 
              center=FALSE, control = list(trace=1, maxit=1000, fnscale = 3000000),
              covariates = ~ female_w1 + age3 + edu3 + race_w1)
      saveRDS(msm1a_le90, paste0(models, "msm1a_le90.RDS")) 


# MSM 1B: Age (7 categories)
msm1b_le90 <- msm ( alc6 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCat,
              center=FALSE, control = list(trace=1, maxit=1000, fnscale = 3000000),
              covariates = ~ female_w1 + age7 + edu3 + race_w1)
      saveRDS(msm1b_le90, paste0(models, "msm1b_le90.RDS"))
        
        
# MSM 1C: Age (continuous)
msm1c_le90 <- msm(alc6 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCont,
            center=FALSE, control = list(trace=1, maxit=5000, fnscale = 750000),
            covariates = ~ female_w1 + age.c + edu3 + race_w1)
        saveRDS(msm1c_le90, paste0(models, "msm1c_le90.RDS")) 


# MSM 1D: Age (squared)
msm1d_le90 <- msm(alc6 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCont,
            center=FALSE, control = list(trace=1, maxit=5000, fnscale = 750000),
            covariates = ~ female_w1 + age.c + age_sq.c + edu3 + race_w1)
        saveRDS(msm1d_le90, paste0(models, "msm1d_le90.RDS")) 





# MSM 2: AlcUse (5levels) ----------------------------------------------------------------------------------------

# Specify allowed transitions
# Only allow adjacent transitions, except for transitions back to lifetime abstainers, and abstainer->former drinker

Q <- rbind ( c(0,     0,    0.25,  0,    0),
             c(0,     0,    0.25,  0,    0),
             c(0,     0.25, 0,     0.25, 0),
             c(0,     0,    0.25,  0,    0.25),
             c(0,     0,    0,     0.25, 0))

# Specifies initial values
Q_ageCat  <- crudeinits.msm(alc5 ~ years, idnum, data=nesarc_expanded, qmatrix=Q)      # When age is categorical
Q_ageCont <- crudeinits.msm(alc5 ~ years, idnum, data=nesarc_expanded_ageC, qmatrix=Q) # When age is continuous 



# MSM 2: All ages **************************************************************************************
# MSM 2A: Age (3 categories)
msm2a <- msm (alc5 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q_ageCat, 
             center=FALSE, control = list(trace=1, maxit=500, fnscale = 3000000),
              covariates = ~ female_w1 + age3 + edu3 + race_w1)
        saveRDS(msm2a, paste0(models, "msm2a.RDS")) 


# MSM 2B: Age (7 categories)
msm2b <- msm (alc5 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q_ageCat, 
              center=FALSE,control = list(trace=1, maxit=600, fnscale = 3000000),
              covariates = ~ female_w1 + age7 + edu3 + race_w1)
        saveRDS(msm2b, paste0(models, "msm2b.RDS")) 

        
        
# MSM 2: Ages less than 90years (le90) *********************************************************************************      
# MSM 2A: Age (3 categories)
msm2a_le90 <- msm (alc5 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCat, 
              center=FALSE, control = list(trace=1, maxit=500, fnscale = 3000000),
              covariates = ~ female_w1 + age3 + edu3 + race_w1)
      saveRDS(msm2a_le90, paste0(models, "msm2a_le90.RDS")) 


# MSM 2B: Age (7 categories)
msm2b_le90 <- msm (alc5 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCat, 
            center=FALSE,control = list(trace=1, maxit=600, fnscale = 3000000),
            covariates = ~ female_w1 + age7 + edu3 + race_w1)
      saveRDS(msm2b_le90, paste0(models, "msm2b_le90.RDS")) 

# MSM 2C: Age (continuous)
msm2c_le90 <- msm(alc5 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCont,
            center=FALSE, control = list(trace=1, maxit=5000, fnscale = 750000),
            covariates = ~ female_w1 + age.c + edu3 + race_w1)
        saveRDS(msm2c_le90, paste0(models, "msm2c_le90.RDS")) 


# MSM 2D: Age (squared)
msm2d_le90 <- msm(alc5 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCont,
            center=FALSE, control = list(trace=1, maxit=5000, fnscale = 750000),
            covariates = ~ female_w1 + age.c + age_sq.c + edu3 + race_w1)
        saveRDS(msm2d_le90, paste0(models, "msm2d_le90.RDS")) 




# MSM 3: AlcUse (4levels) ---------------------------------------------------------------------------------------------

# Specify allowed transitions; only allow adjacent transitions
Q <- rbind ( c(0,    0.25,  0,    0),
             c(0.25, 0,     0.25, 0),
             c(0,    0.25,  0,    0.25),
             c(0,    0,     0.25, 0))

# Specifies initial values
Q_ageCat  <- crudeinits.msm(alc4 ~ years, idnum, data=nesarc_expanded, qmatrix=Q)      # When age is categorical
Q_ageCont <- crudeinits.msm(alc4 ~ years, idnum, data=nesarc_expanded_ageC, qmatrix=Q) # When age is continuous 



# MSM 3: All ages **************************************************************************************
# MSM 3A: Age (3 categories)
msm3a <- msm (alc4 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q_ageCat, 
              center=FALSE, control = list(trace=1, maxit=600, fnscale = 3000000),
                  covariates = ~ female_w1 + age3 + edu3 + race_w1)
        saveRDS(msm3a, paste0(models, "msm3a.RDS")) 

   
# MSM 3B: Age (7 categories)
msm3b <- msm (alc4 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q_ageCat,
              center=FALSE, control = list(trace=1, maxit=600, fnscale = 3000000),
              covariates = ~ female_w1 + age7 + edu3 + race_w1)
        saveRDS(msm3b, paste0(models, "msm3b.RDS")) # Save Results

        
        
        
# MSM 3: Ages less than 90years (le90) *********************************************************************************    
# MSM 3A: Age (3 categories)
msm3a_le90 <- msm (alc4 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCat, 
              center=FALSE, control = list(trace=1, maxit=600, fnscale = 3000000),
              covariates = ~ female_w1 + age3 + edu3 + race_w1)
      saveRDS(msm3a_le90, paste0(models, "msm3a_le90.RDS")) 


# MSM 3B: Age (7 categories)
msm3b_le90 <- msm (alc4 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCat,
              center=FALSE, control = list(trace=1, maxit=600, fnscale = 3000000),
              covariates = ~ female_w1 + age7 + edu3 + race_w1)
      saveRDS(msm3b_le90, paste0(models, "msm3b_le90.RDS")) # Save Results
        

# MSM 3C: Age (continuous)
msm3c_le90 <- msm(alc4 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCont,
              center=FALSE, control = list(trace=1, maxit=5000, fnscale = 750000),
              covariates = ~ female_w1 + age.c + edu3 + race_w1)
       saveRDS(msm3c_le90, paste0(models, "msm3c_le90.RDS")) 

      
# MSM 3D_le90: Age (squared)
msm3d_le90 <- msm(alc4 ~ years, subject=idnum, data = nesarc_expanded_le90, qmatrix = Q_ageCont,
            center=FALSE, control = list(trace=1, maxit=5000, fnscale = 750000),
            covariates = ~ female_w1 + age.c + age_sq.c + edu3 + race_w1)
        saveRDS(msm3d_le90, paste0(models, "msm3d_le90.RDS")) # Save Results

     


# Compare the MSM models (with adjacent transitions) -------------------------------------------------------

# Models with participants aged less than 90 years
msm1a_le90 <- readRDS(paste0(models, "msm1a_le90.RDS"))
msm1b_le90 <- readRDS(paste0(models, "msm1b_le90.RDS"))
msm1c_le90 <- readRDS(paste0(models, "msm1c_le90.RDS"))
msm1d_le90 <- readRDS(paste0(models, "msm1d_le90.RDS"))
msm2a_le90 <- readRDS(paste0(models, "msm2a_le90.RDS"))
msm2b_le90 <- readRDS(paste0(models, "msm2b_le90.RDS"))
msm2c_le90 <- readRDS(paste0(models, "msm2c_le90.RDS"))
msm2d_le90 <- readRDS(paste0(models, "msm2d_le90.RDS"))
msm3a_le90 <- readRDS(paste0(models, "msm3a_le90.RDS"))
msm3b_le90 <- readRDS(paste0(models, "msm3b_le90.RDS"))
msm3c_le90 <- readRDS(paste0(models, "msm3c_le90.RDS"))
msm3d_le90 <- readRDS(paste0(models, "msm3d_le90.RDS"))

AIC(msm1a_le90, msm1b_le90, msm1c_le90, msm1d_le90, 
    msm2a_le90, msm2b_le90, msm2c_le90, msm2d_le90, 
    msm3a_le90, msm3b_le90, msm3c_le90, msm3d_le90) %>% 
  arrange(AIC)

# Best model: AlcUse 7 categotries, age 7 categories



# Models with all participants
msm1a <- readRDS(paste0(models, "msm1a.RDS"))
msm1b <- readRDS(paste0(models, "msm1b.RDS"))
msm2a <- readRDS(paste0(models, "msm2a.RDS"))
msm2b <- readRDS(paste0(models, "msm2b.RDS"))
msm3a <- readRDS(paste0(models, "msm3a.RDS"))
msm3b <- readRDS(paste0(models, "msm3b.RDS"))

AIC(msm1a, msm1b, msm2a, msm2b, msm3a, msm3b) %>% arrange(AIC)

# Best model: msm3b -> AlcUse 7 categotries, age 7 categories



# Non-adjacent transitions (AlcUse 4 levels; Age 7 levels) -----------------------------------------------------------

# Specify allowed transitions
# Only allow adjacent transitions, and transitions back to non-drinker
Q <- rbind (c(0,    0.25,  0,    0),
            c(0.25, 0,     0.25, 0),
            c(0.25, 0.25,  0,    0.25),
            c(0.25, 0,     0.25, 0))


# Specifies initial values
Q <- crudeinits.msm(alc4 ~ years, idnum, data=nesarc_expanded, qmatrix=Q)

# MSM 3E: Age (7 categories); non-adjacent transitions 
msm3e <- msm (alc4 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q,
              center=FALSE, control = list(trace=1, maxit=600, fnscale = 3000000),
              covariates = ~ female_w1 + age7 + edu3 + race_w1)
      saveRDS(msm3e, paste0(models, "msm3e.RDS")) 

# Compare the models
msm3b <- readRDS(paste0(models, "msm3b.RDS"))
msm3c <- readRDS(paste0(models, "msm3c.RDS"))

AIC(msm3b, msm3c)

# Best model: msm3b -> AlcUse 7 categotries, age 7 categories, with adjacent transitions only  



# Unadjusted version of optimal MSM modle

# Specify allowed transitions; only allow adjacent transitions
Q <- rbind (c(0,    0.25,  0,    0),
            c(0.25, 0,     0.25, 0),
            c(0,    0.25,  0,    0.25),
            c(0,    0,     0.25, 0))

# Specifies initial values
Q <- crudeinits.msm(alc4 ~ years, idnum, data=nesarc_expanded, qmatrix=Q)

# Unadjusted Model; Age (7 categories) -----------------------------------------------------------
msm3b_crude <- msm (alc4 ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q,
                    center=FALSE, control = list(trace=1, maxit=600, fnscale = 3000000))
saveRDS(msm3b_crude, paste0(models, "msm3b_crude.RDS")) 



# HED Model -----------------------------------------------------------------------------------

# Specify allowed transitions
# only allow adjacent transitions
Q <- rbind (c(0,     0.25,    0,     0,    0),
  c(0.25,  0,       0.25,  0,    0),
  c(0,     0.25,    0,     0.25, 0),
  c(0,     0,       0.25,  0,    0.25),
  c(0,     0,       0,     0.25, 0))

# Specify initial values 
Q <- crudeinits.msm(hed ~ years, idnum, data=nesarc_expanded, qmatrix=Q)

# Run MSM model (unadjusted)
hed_crude.msm <- msm (hed ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q, 
                      center=FALSE, control = list(trace=1, maxit=500, fnscale = 3000000))
saveRDS(hed_crude.msm, paste0(models, "hed_crude.msm.RDS"))


# Run MSM model (adjusted for covariates)
hed.msm <- msm (hed ~ years, subject=idnum, data = nesarc_expanded, qmatrix = Q, 
                center=FALSE, control = list(trace=1, maxit=500, fnscale = 3000000),
                covariates = ~ female_w1 + age7 + edu3 + race_w1)  # For functions to work, the order of covariates should be: sex, age, edu, race
saveRDS(hed.msm, paste0(models, "hed.msm.RDS"))

