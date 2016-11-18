# ESS-211 assignment 5, due 11-11-16
#  cba 11/2016

# set current working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# source the written functions
source("ESS-211-Functions.R")

#############################
# test local SA using provided code

# define bungee parameters and names
par.names <- c("Height", "Mass", "# Strands")
par.mins <- c(40 , 67 , 20)
par.maxs <- c(60 , 74 , 40)

# define midpoints for parameters
par.mids <- 0.5 * (par.mins + par.maxs)

# set n parameters
npars <- length(par.mins)

# create array to hold results
difs <- numeric(length = npars)

# loop through each parameter set and
for (par in 1:npars){
  
  # create two input vectors to assign min/max values
  input1 <- par.mids 
  input2 <- par.mids
  
  # assign min value to first vector
  input1[par] <- par.mins[par]
  
  # assign max value to second vector
  input2[par] <- par.maxs[par]
  
  # compute differences between model outputs using different min/max vals
  difs[par] = hmin(input2) - hmin(input1)
}

# set up plotting parameters
cols <- c("salmon", "light green", "orange")
ylim <- c(min(difs), max(difs))

# plot the output, showing the difference in model outputs based on varying min/max
#  for each parameter in an OAT fashion.
barplot(difs, xlab = "Parameter", ylab = "Max - Min", names.arg = par.names, ylim = ylim, col = cols)

#############################
# test one-at-a-time global SA

# define n parameters
k <- npars

# n times to compute effect for each parameter
r <- 10

# n possible levels for each parameter (should be an even #)
p <- 4

# the increment to adjust parameter values by
delta <- p / (2 * (p - 1))

# define range for parameters
par.mins <- c(40, 67, 20)
par.maxs <- c(60, 74, 40)

# run morris method for hmin
effects <- morris(k, r, p, delta, par.mins, par.maxs, hmin)

# calculate mean effect
mu <- apply(effects, 1, mean)

# calculate sd effect
sigma <- apply(effects, 1, sd)

# calculate mean absolute effect
mu.star <- apply(abs(effects), 1, mean)

# plot the results
plot(mu.star, sigma, pch = c("H", "M", "s"))