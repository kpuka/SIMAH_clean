

## Functions for Causal Mediation 


# Direct, indirect and mediated interactive effects and standard errors are derived directly from the summary() command below
# Total effect is obtained by the sum of the three separate effects
# Confidence intervals for total effects and mediated proportions are computed using the code below:

# Function to get the (Not-Robust) total effect and proportion mediated
getTE_NotRobust <- function(CMed_model, v){
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



        # # Function to get the (Robust) total effect and proportion mediated
        getTE_Robust <- function(CMed_model, v){
          TE <- sum(CMed_model$gamma[v])
          mu <- CMed_model$gamma[v]
          Omega <- CMed_model$robvar.gamma[v,v] # To obtain robust estimates
          temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
          temp_TE <- apply(temp,1,sum)
          med_prop <- c(mu/TE,1)
          med_prop_CI <- rbind(t(apply(temp/temp_TE, 2, quantile, c(0.025, 0.975))), c(1,1))
          output <- cbind(c(mu,TE), c(apply(temp,2,sd),sd(temp_TE)), med_prop, med_prop_CI)
          colnames(output) <- c("Est.", "SE", "med_prop", "lowerCI", "UpperCI")
          rownames(output) <- c(rownames(CMed_model$gamma)[v],"TE")
          return(output)}





# Function to get the (Not-Robust) total combined indirect effect  
getIE_NotRobust <- function(CMed_model, v){
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



        # Function to get the (Robust) total combined indirect effect  
        getIE_Robust <- function(CMed_model, v){
          IE <- sum(CMed_model$gamma[v])
          mu <- CMed_model$gamma[v]
          Omega <- CMed_model$robvar.gamma[v,v] # To obtain robust estimates
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
getTE_IE_NotRobust <- function(CMed_model, v, z){
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


      

      # Function to get the (Not Robust)  proportion mediated of the combined indirect effect
      getTE_IE_Robust <- function(CMed_model, v, z){
        #total effect
        TE <- sum(CMed_model$gamma[v])
        mu <- CMed_model$gamma[v]
        Omega <- CMed_model$robvar.gamma[v,v]  # To obtain robust estimates
        require(MASS)
        temp <- mvrnorm(n=10^4, mu=mu, Sigma=Omega)
        temp_TE <- apply(temp,1,sum)
        IE <- sum(CMed_model$gamma[z])
        muIE <- CMed_model$gamma[z]
        OmegaIE <- CMed_model$robvar.gamma[z,z] # To obtain robust estimates
        require(MASS)
        tempIE <- mvrnorm(n=10^4, mu=muIE, Sigma=OmegaIE)
        temp_IE <- apply(tempIE,1,sum)
        med_prop <- c(IE/TE,1)
        med_prop_CI <- (temp_IE/temp_TE)
        output <- cbind(IE, med_prop, quantile)
        quantile <- quantile(med_prop_CI, c(0.025, 0.975))
        output <- cbind(IE, med_prop, quantile)
        return(output)}
