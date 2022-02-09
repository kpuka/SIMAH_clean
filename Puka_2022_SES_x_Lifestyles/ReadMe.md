# About

This is the statistical code for:
Puka, Buckley, Mulia, Lasserre, Rehm, and Probst (2022). Educational Attainment and Lifestyle Risk Factors Associated with All-Cause Mortality in the United States: Decomposing differential exposure and vulnerability. JAMA Health Forum.

This repo is hereby licensed for use under the GNU GPL version 3.



The data used is publicly available from the National Center for Health Statistics.

- The 1997-2014 National Health Interview Survey (NHIS) survey are available here: https://www.cdc.gov/nchs/nhis/1997-2018.htm

- The 2015 Public-Use Linked Mortality Files are available here: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

The data files were extracted and initially manipulated using SAS. All analyses were completed in R. The SAS and R code are available in this repository. The key variables included in the analyses are listed below.

# Key Variables:

| Variable           	| Label                                 	| Values (Value Labels)                                                                                   	|
|-------------------	|---------------------------------------	|----------------------------------------------------------------------------------------------------------	|
| allcause_death    	| Status at last follow-up              	| 0 (Alive); 1 (Deceased)                                                                                   |
| alcohol5v2        	| Alcohol use                           	| 1 (Never Drinker); 2 (Former Drinker); 3 (Category I; reference); 4 (Category II); 5 (Category III)      	|
| smoking4           	| Smoking status                        	| 1 (Never smoker; reference); 2 (Former smoker); 3 (Current some day smoker); 4 (Current everyday smoker) 	|
| bmi_cat           	| Body mass index                       	| 1 (Underweight); 2 (Healthy weight; reference); 3 (Overweight); 4 (Obese)                                	|
| phy_act3          	| Physical activity                     	| 1 (Sedentary); 2 (Somewhat active); 3 (Active; reference)                                                	|
| edu            	    | Educational attainment                	| 1 (High school); 2 (Some college); 3 (Bachelors; reference)                                               	|
| female            	| Sex                                   	| 0 (Men); 1 (Women)                                                                                       	|
| bl_age         	    | Age (years) at baseline 	              | Continuous                                                                                               	|
| end_age        	    | Age (years) at last follow-up/death   	| Continuous                                                                                              	|
| ethnicity      	    | Ethnicity                             	| 1 (Non-Hispanic White; reference); 2 (Non-Hispanic Black); 3 (Hispanic); 4 (Other)                       	|
| married        	    | Marital status                        	| 0 (Not married/living together); 1 (Married/cohabitating)                                                	|
| srvy_yr        	    | Survey year                           	| Integer: 1997 to 2014                                                                                    	|
| hed            	    | Heavy episodic drinking               	| 1 (No HED); 2 (HED <1/month); 3 (HED >1/month, <1/week); 4 (HED >=1/week)                                	|

*Note: in the R scripts variables with the suffix ".factor" use the value label noted above, whereas variables without the suffix use the numeric value noted above*
