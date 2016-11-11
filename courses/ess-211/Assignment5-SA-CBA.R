# ESS-211 assignment 5, due 11-11-16
#  cba 11/2016

# set current working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# source the written functions
source("ESS-211-Functions.R")

#############################
# set up cba-written functions used in this assignment 

# to assess how mu.star and sigma change as a function of changing ranges and midpoints,
#  we'll create a function and shift the min/max values a few times
script.morris <- function(mins, maxs, names, title, k, r, p, delta, fun){
  
  # set midpoints for each parameter
  par.mids <- 0.5 * (mins + maxs)
  
  # set n parameters
  npars <- length(mins)
  
  # run morris method for the function
  effects <- morris(k, r, p, delta, mins, maxs, fun)
  
  # calculate mean effect
  mu <- apply(effects, 1, mean)
  
  # calculate sd effect
  sigma <- apply(effects, 1, sd)
  
  # calculate mean absolute effect
  mu.star <- apply(abs(effects), 1, mean)
  
  # set up some plotting params
  pch <- rep(19, npars)
  cols <- rainbow(npars)
  
  # plot the results
  plot(mu.star, sigma, pch = pch, main = title, col = cols)
  legend("topleft", legend = names, pch = pch, col = cols)
  
  # return the effects
  return(list(mu, sigma, mu.star))
}

# i have to create a modified version of the indian summer monsoon model i wrote for assignment 2
#  to include # years as a parameter
ism.modified <- function(inputParameters, nYears=200){
  
  # create an output vector to save results
  output <- rep(0, nYears)
  
  # run model over the total number of years
  for (i in 1:nYears){
    output[i] <- ism.rainfall(inputParameters)
  }
  
  # return the output
  return(output)
}

#############################
# task 1 - running the morris method with different parameter ranges

# define original bungee parameters to go +/- 10 on each of the parameters
par.mins <- c(40, 67, 20)
par.maxs <- c(60, 74, 40)

# define the parameter names
par.names <- c("Height", "Mass", "# Strands")

# set n parameters
k <- length(par.mins)

# n times to compute effect for each parameter
r <- 10

# n possible levels for each parameter (should be an even #)
p <- 4

# the increment to adjust parameter values by
delta <- p / (2 * (p - 1))

# run the morris script and change the min/max values
script.morris(par.mins, par.maxs, par.names, "Default parameters", k, r, p, delta, hmin)
script.morris(par.mins-5, par.maxs, par.names, "Min values - 5", k, r, p, delta, hmin)
script.morris(par.mins, par.maxs+5, par.names, "Max values + 5", k, r, p, delta, hmin)
script.morris(par.mins+5, par.maxs+5, par.names, "Min/max values + 5", k, r, p, delta, hmin)

#############################
# task 2 - running sensitivity analysis on the indian summer monsoon model

# set names and ranges
par.names <- c("Pstrong", "Pweak", "tau", "prmax", "P_init")
par.mins <- c(8, 0, 14, .7, .65)
par.maxs <- c(10, 2, 21, .9, .85)

# set the morris parameters per assignment instructions
k <- length(par.mins)
p <- 4
r <- 20
delta <- 2/3

# run the morris script
script.morris(par.mins, par.maxs, par.names, "Schewe-Levermann rainfall SA", k, r, p, delta, ism.modified)
