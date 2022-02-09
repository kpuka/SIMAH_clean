

## Functions to extract formatted results



# Function to extract results from Aalen model, multiple by 10,000py and round/format results
aalen_10000py <- function(model, x) {
  mu <- model$gamma[x]       
  var <-  model$var.gamma[x,x]
  confint.lower <- round((mu - (1.96 * sqrt(var)))*10000,1) # CI * 10,000 to get result per 10,000py
  confint.upper <- round((mu + (1.96 * sqrt(var)))*10000,1) # CI * 10,000 to get result per 10,000py
  mu <- round(mu*10000,1)                                   # mu * 10,000 to get result per 10,000py 
  output<-paste0(mu, " (",confint.lower,", ", confint.upper, ")")
  return(cat(output, "\n"))}   #cat() returns the text without quotes and without the leading numbers [1], [2]...



# Function to extract and format HR results
cox_HR <- function(model, x) {
  mu <- model$coefficients[x]       
  var <-  model$var[x,x]
  confint.lower <- round((exp(mu - (1.96 * sqrt(var)))),2) # Exp() lower CI and round it to 2 decimals
  confint.upper <- round((exp(mu + (1.96 * sqrt(var)))),2) # Exp() upper CI and round it to 2 decimals
  mu <- round(exp(mu),2)                                   # Exp() estimate and round it to 2 decimals
  output<-paste0(mu, " (",confint.lower,", ", confint.upper, ")")
  return(cat(output, "\n"))}   #cat() returns the text without quotes and without the leading numbers [1], [2]...

