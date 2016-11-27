# ESS-211 Final Project
# Modeling de-militarization in the US
# Christopher Anderson
# 12-2016

#############################
# below are the functions used in the ESS-211-DemilitarizationProject.R script
library(car)

# set up a function for exponential growth
#  based on fitting a polynomial using R's built in functions
#  returns a model
growth.exponential <- function(x, y, degree=2){
  model <- lm(y ~ poly(x, degree = degree))
  return(model)
}

# set up a function for logistic growth
#  based on the the logistic growth function:
#  y = phi1/(1 + exp(-(phi2 + phi3 * x)))
#  where y is the population, phi1 is the asymptote (or carrying capacity),
#  i don't remember phi2, phi3 is the growth parameter, and 
#  x is the response variable, in this case, time
growth.logistic <- function(x, y, phi1){
  # create a data frame for these functions
  dframe <- data.frame(y, x)
  
  # we're going to have to guess the initial parameters for 
  #  phi2 and phi3, which we will do using a logit transform
  coefs <- coef(lm(logit((y / phi1) ~ x, data = dframe)))
  phi2 <- coefs[1]
  phi3 <- coefs[2]
  
  # then we'll perform a non-linear least-squares fit to the data
  model <- nls(y ~ (phi1 / (1 + exp(-(phi2 + phi3 * x)))), 
               start = list(phi1 = phi1, phi2 = phi2, phi3 = phi3),
               data = dframe)
  
  return(model)
}