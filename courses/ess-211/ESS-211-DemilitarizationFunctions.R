# ESS-211 Final Project
# Modeling de-militarization in the US
# Christopher Anderson
# 12-2016

#############################
# below are the functions used in the ESS-211-DemilitarizationProject.R script
library(car)

#############################
# first, functions for modeling federal spending

# set up a function for exponential growth
#  based on fitting a polynomial using R's built in functions
#  returns a model
growth.exponential <- function(dframe, degree=2){
  model <- lm(FedSpending.BUSD ~ poly(Year, degree = degree), data = dframe)
  return(model)
}

# set up a function for logistic growth
#  based on the the logistic growth function:
#  y = phi1/(1 + exp(-(phi2 + phi3 * x)))
#  where y is the population, phi1 is the asymptote (or carrying capacity),
#  i don't remember phi2, phi3 is the growth parameter, and 
#  x is the response variable, in this case, time
growth.logistic <- function(dframe, phi1, phiIndex){
  # we're going to use the index provided to insert the phi1 value
  #  so the logistic growth fits to the asymptote
  dframe$FedSpending.BUSD[phiIndex:nrow(dframe)] <- phi1
  
  # we're going to have to guess the initial parameters for 
  #  phi2 and phi3, which we will do using a logit transform
  coefs <- coef(lm(logit(FedSpending.BUSD / phi1) ~ Year, data = dframe))
  phi2 <- coefs[1]
  phi3 <- coefs[2]
  
  # then we'll perform a non-linear least-squares fit to the data
  model <- nls(FedSpending.BUSD ~ (phi1 / (1 + exp(-(phi2 + phi3 * Year)))), 
               start = list(phi1 = phi1, phi2 = phi2, phi3 = phi3),
               data = dframe)
  
  return(model)
}

#############################
# we're going to create similar models, exponential and logistic,
#  to model the growth in military spending

# first, exponential spending
milGrowth.exponential <- function(dframe, degree=2){
  model <- lm(MilitarySpending.BUSD ~ poly(Year, degree = degree), data = dframe)
  return(model)
}

# then, logistic spending
milGrowth.logistic <- function(dframe, phi1, phiIndex){
  dframe$MilitarySpending.BUSD[phiIndex:nrow(dframe)] <- phi1
  
  # we're going to have to guess the initial parameters for 
  #  phi2 and phi3, which we will do using a logit transform
  coefs <- coef(lm(logit(MilitarySpending.BUSD / phi1) ~ Year, data = dframe))
  phi2 <- coefs[1]
  phi3 <- coefs[2]
  
  # then we'll perform a non-linear least-squares fit to the data
  model <- nls(MilitarySpending.BUSD ~ (phi1 / (1 + exp(-(phi2 + phi3 * Year)))), 
               start = list(phi1 = phi1, phi2 = phi2, phi3 = phi3),
               data = dframe)
  
  return(model)
}

#############################
# now, to model the growth in veteran spending
# first, exponential spending
vetGrowth.exponential <- function(dframe, degree=2){
  model <- lm(VeteranSpending.BUSD ~ poly(Year, degree = degree), data = dframe)
  return(model)
}

# then, logistic spending
vetGrowth.logistic <- function(dframe, phi1, phiIndex){
  dframe$VeteranSpending.BUSD[phiIndex:nrow(dframe)] <- phi1
  
  # we're going to have to guess the initial parameters for 
  #  phi2 and phi3, which we will do using a logit transform
  coefs <- coef(lm(logit(VeteranSpending.BUSD / phi1) ~ Year, data = dframe))
  phi2 <- coefs[1]
  phi3 <- coefs[2]
  
  # then we'll perform a non-linear least-squares fit to the data
  model <- nls(VeteranSpending.BUSD ~ (phi1 / (1 + exp(-(phi2 + phi3 * Year)))), 
               start = list(phi1 = phi1, phi2 = phi2, phi3 = phi3),
               data = dframe)
  
  return(model)
}

#############################
# unemployment modeling

# framework 1
#  finds the mean of the unemployment data, then adds normally distributed noise based on the 
#  standard deviation provided
unemployment.framework1 <- function(dframe, sd){
  
  # get basic unemployment info
  unemployment.years <- which(dframe$UnemploymentRate < 1.)
  unemployment.yearsToModel <- is.na(dframe$UnemploymentRate)
  unemployment.mean <- mean(dframe$UnemploymentRate[unemployment.years])
  
  # generate some noise
  noise <- rnorm(sum(unemployment.yearsToModel), sd = sd)
  
  # add the noise to the mean
  framework1 <- dframe$UnemploymentRate
  framework1[unemployment.yearsToModel] = unemployment.mean + noise
  return(framework1)
}

# framework 2
#  finds the trendline of the unemployment data, then adds normally distributed noise
#  based on the standard deviation provided
unemployment.framework2 <- function(dframe, sd){
  
  # get basic unemployment info
  unemployment.years <- which(dframe$UnemploymentRate < 1.)
  unemployment.yearsToModel <- is.na(dframe$UnemploymentRate)
  
  # find trendline in the unemployment data
  linearModel <- lm(dframe$UnemploymentRate[unemployment.years] ~ yearly$Year[unemployment.years])
  coefs <- coef(linearModel)
  
  # generate some noise
  noise <- rnorm(sum(unemployment.yearsToModel), sd = sd)
  
  # add the noise to the linear trendline
  framework2 <- dframe$UnemploymentRate
  framework2[unemployment.yearsToModel] <- (yearly$Year[unemployment.yearsToModel] * coefs[2] + coefs[1]) + noise
}
