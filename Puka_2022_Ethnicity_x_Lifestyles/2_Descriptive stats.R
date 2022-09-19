# Race x Lifestyle Differential Vulnerability & Exposure Project
# Descriptive Statistics

# LOAD DATA AND SET FILE LOCATIONS

# load libraries
library(tidyverse)  # data management
library(skimr)      # descriptive statistics
library(gmodels)    # CrossTable command
library(tableone)   # create table one
library(survival)   # surivval analyses
library(survminer)  # surivval analyses
library(timereg)    # additive survival models
library(survey)     # for survey weighted cox model
library(biostat3)   # survRate command
library(broom)      # data/model management
library(fmsb)       # Radar plots


# Specify file locations
data    <- "C:/.../Processed data/"  # Location of data
output  <- "C:/.../Output/"          # Location of figures/tables

# Load data (participants aged 25-85 years)
nhis_all    <- readRDS (paste0(data, "nhis_all18_85.rds"))
nhis        <- readRDS (paste0(data, "nhis18_85.rds"))
nhis_male   <- filter(nhis, female==0)
nhis_female <- filter(nhis, female==1)


# DESCRIPTIVE STATISTICS 

# Table 1: Participant characteristics - STRATIFIED BY SEX

tab1 <-CreateTableOne(vars= c("age", "yrs_followup","allcause_death.factor", "alcohol5v2.factor","smoking4.factor",  
                               "bmi_cat.factor", "phy_act3.factor", "edu.factor",  "married.factor", "us_born"), 
                      factorVars = c("allcause_death.factor",  "alcohol5v2.factor", "smoking4.factor", "bmi_cat.factor", 
                                      "phy_act3.factor", "edu.factor",  "married.factor", "us_born"), 
                      strata= c("ethnicity.factor", "female.factor"), addOverall = TRUE, data=nhis)
  table1_v1 <- print(tab1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE)  # Shows sample size and %
  table1_v2 <- print(tab1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE, format="p") # shows % only
  write.csv(table1_v1, file = file.path(output, "Table1 Demographics_V1.csv"))
  write.csv(table1_v2, file = file.path(output, "Table1 Demographics_V2.csv"))
  kableone(table1_v1)

  
# Person years and death rate: tstop = person years; event = # events; rate = events per person year 
survRate(Surv(yrs_followup, allcause_death)~female.factor+ethnicity.factor, data=nhis) %>% view
  mutate(group = paste(female.factor, ethnicity.factor, sep=", "),
         rate_10000py = rate * 10000) %>% remove_rownames() %>%
  dplyr::select (group, tstop, rate_10000py) %>% t()%>%
  write.csv(file = file.path(output, "Table1 Demographics_V3.csv"))
  
survRate(Surv(yrs_followup, allcause_death)~ethnicity.factor, data=nhis) %>%view
  mutate(rate_10000py = round(rate * 10000, 0)) %>% 
  dplyr::select (ethnicity.factor, rate_10000py)



# Race/ethnicity proportions
nhis %>% group_by (ethnicity.factor) %>%
  summarise (n = n()) %>%
  mutate (perc = round(n / sum(n) *100 ,0))

nhis %>% filter (ethnicity.factor == "Other") %>% 
  group_by (ethnicity_detail) %>%
  summarise (n = n()) %>%
  mutate (perc = round(n / sum(n) *100 ,0))



# eTable 1: Attrition - Participant characteristics - STRATIFIED BY SEX
nhis_all_male <- filter(nhis_all, female==0)
nhis_all_female <- filter(nhis_all, female==1)

all_vars <- c("age", "alcohol5v2.factor","smoking4.factor", "bmi_cat.factor", "phy_act3.factor", "edu.factor",  "married.factor")
factor_vars <- c("alcohol5v2.factor", "smoking4.factor", "bmi_cat.factor", "phy_act3.factor", "edu.factor",  "married.factor")

tab_e1 <-CreateTableOne(all_vars, factor_vars, strata="lost", data=nhis_all_male)
table_e1 <- print(tab_e1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE, smd=TRUE, format="p") # shows % only
write.csv(table_e1, file = file.path(output, "Table e1 (men) Attrition.csv"))
kableone(table_e1)

tab_e1 <-CreateTableOne(all_vars, factor_vars, strata="lost", data=nhis_all_female)
table_e1 <- print(tab_e1, noSpaces = TRUE, catDigits = 0, contDigits = 1, printToggle = FALSE, test=FALSE, smd=TRUE, format="p") # shows % only
write.csv(table_e1, file = file.path(output, "Table e1 (women) Attrition.csv"))
kableone(table_e1)




# Figure 2: Survival plot 
ggsurvplot_facet(fit = survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity.factor, data = nhis), 
  data=nhis, facet.by="female.factor", censor = FALSE,xlim = c(18, 100), 
  conf.int = TRUE, 
  xlab = "Age (years)", 
  ylab = "Overall survival probability")

      # Age Medium Survival:
      survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity.factor, data = nhis)
      survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity.factor, data = nhis_female)
      survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity.factor, data = nhis_male)

      
# eFigure 1: Survival plot 
ggsurvplot_facet(fit = survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity_detail, data = nhis), 
  data=nhis, facet.by="female.factor", censor = FALSE,xlim = c(18, 100), 
  conf.int = TRUE, 
  xlab = "Age (years)", 
  ylab = "Overall survival probability") 
      
      # Age Medium Survival:
      survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity_detail, data = nhis)
      survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity_detail, data = nhis_female)
      survfit(Surv(bl_age, end_age, allcause_death) ~ ethnicity_detail, data = nhis_male)
    
      
    
# Radar plots
radar <- nhis %>%
  mutate(high_risk_drinker = ifelse(alcohol5v2.factor %in% c("Category II", "Category III"), 1, 0),
         everyday_smoker = ifelse(smoking4.factor %in% c("Current everyday smoker"), 1, 0),
         sedentary = ifelse(phy_act3.factor %in% c("Sedentary"), 1, 0),
         obese = ifelse(bmi_cat.factor %in% c("Obese"), 1, 0)) %>%
  group_by(ethnicity.factor, female.factor) %>%
  summarise(across(c(high_risk_drinker, everyday_smoker, sedentary, obese), mean)) %>%
  mutate (across(c(high_risk_drinker, everyday_smoker, sedentary, obese), ~ round(.x*100, 0)))  %>% 
  ungroup()


# Stratified by sex
my_radarchart <- function(data, ethnicity, title = NULL, legend=FALSE){
      
      min <- 0
      max <- 60
      
      data <- data %>%
        filter(ethnicity.factor==ethnicity) %>% 
        dplyr::select(-ethnicity.factor) %>%
        mutate(female.factor = recode(female.factor, "Male" = "Men", "Female" = "Women")) %>%
        column_to_rownames(var = "female.factor") %>%
        add_row(high_risk_drinker=min, everyday_smoker=min, sedentary=min, obese=min, .before=1) %>%
        add_row(high_risk_drinker=max, everyday_smoker=max, sedentary=max, obese=max, .before=1) 
      
      
      radarchart(data, axistype = 1,
                # Customize the polygon
                pcol = c("#00AFBB", "#FC4E07"), pfcol = scales::alpha(c("#00AFBB", "#FC4E07"), 0.15), 
                plwd = 2, plty = 1,
                # Customize the grid
                cglcol = "grey", cglty = 1, cglwd = 0.8,
                # Customize the axis
                axislabcol = "grey", 
                # Variable labels
                vlcex = 1, vlabels = c("Category II/III \n drinker", "Current \n everyday \n smoker", "Sedentary", "Obese"),
                caxislabels = seq(0,60,15), title = title)
      
      if(legend==TRUE){
        legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col= c("#00AFBB", "#FC4E07") , text.col = "black", cex=1.2, pt.cex=3)
      }
    }
    
    tiff(file.path(output, "Figure 2 - Prevalence of risk factors.tiff"), width=2200, height = 2000, res=300)
    op <- par(mar = c(1, 1, 1, 1))
    par(mfrow=c(2,2))
    my_radarchart(radar, "White", title="Non-Hispanic White")
    my_radarchart(radar, "Black", title="Non-Hispanic Black", legend=TRUE)
    my_radarchart(radar, "Hispanic", title="Hispanic/Latinx")
    my_radarchart(radar, "Other", title="Non-Hispanic Other")
    par(mfrow=c(1,1))
    dev.off()
    
    
    
# Stratified by ethnicity    
my_radarchart2 <- function(data, sex, title = NULL, legend=FALSE){
  
  min <- 0
  max <- 60
  
  data <- data %>%
    filter(female.factor==sex) %>% 
    filter(ethnicity.factor!="Other") %>% 
    dplyr::select(-female.factor) %>%
    #mutate(female.factor = recode(female.factor, "Male" = "Men", "Female" = "Women")) %>%
    column_to_rownames(var = "ethnicity.factor") %>%
    add_row(high_risk_drinker=min, everyday_smoker=min, sedentary=min, obese=min, .before=1) %>%
    add_row(high_risk_drinker=max, everyday_smoker=max, sedentary=max, obese=max, .before=1) 
  
  radarchart(data, axistype = 1,
    # Customize the polygon
    pcol = c("#00AFBB", "#FC4E07", "#E7B800", "#99999980"), 
            pfcol = scales::alpha(c("#00AFBB", "#FC4E07", "#E7B800"), 0.15), 
    plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = 1, vlabels = c("Category II/III \n drinker", "Current \n everyday \n smoker", "Sedentary", "Obese"),
    caxislabels = seq(0,60,15), title = title)
  
  if(legend==TRUE){
    legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , 
          col= c("#00AFBB", "#FC4E07", "#E7B800") , text.col = "black", cex=1.2, pt.cex=3)
  }
}


    tiff(file.path(output, "Figure 2_v2 - Prevalence of risk factors.tiff"), width=1400, height = 2500, res=300)
    op <- par(mar = c(1, 1, 1, 1))
    par(mfrow=c(2,1))
    my_radarchart2(radar, "Male", title="Men", legend=TRUE)
    my_radarchart2(radar, "Female", title="Women")
    par(mfrow=c(1,1))
    dev.off()
    


    
# Race by education counts
library(plotly)   
library(htmlwidgets)
race_edu <- nhis %>%
  dplyr::select(female.factor, ethnicity.factor, edu.factor) %>% 
  count(female.factor, ethnicity.factor, edu.factor) %>%
  rename(sex = female.factor, ethnicity = ethnicity.factor, edu = edu.factor) %>% 
  mutate(edu = factor(edu, levels=c("Highschool", "Some college", "Bachelors"))) %>% 
  ggplot(aes(ethnicity, edu, fill= n)) + geom_tile() + facet_wrap(~sex)

ggplotly(race_edu)  
ggplotly(race_edu) %>% saveWidget(file=paste0(output, "Race by education counts.html"))
    

  
  
  