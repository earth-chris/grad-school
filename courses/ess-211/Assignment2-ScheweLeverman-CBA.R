# ESS-211 assignment 2, due 10-19-16
#  cba 10/2016
#
# Question 1:
#  First year default precip : 3.8 mm/day"
#  100th year default precip : 6.06666666666667 mm/day"

#############################
# below is the function used for this assignment

# set up the ISM rainfall function as translated from matlab
ismRainfall <- function(inputVector, lSeason=135){
  # the input vector should be a 5-element vector, in this order:
  # [1] pStrong  : precipitation in wet state (mm/day)
  # [2] pWeak    : precipitation in dry state (mm/day)
  # [3[ tau      : memory length in time steps (days)
  # [4] prMax    : maximum probability of either state.
  # [5] pInit    : initial probability of strong state.
  #                must be <= prMax
  # lSeason  : length of the wet season (days). default = 135
  
  # split the input vector for legibility
  pStrong <- inputVector[1]
  pWeak <- inputVector[2]
  tau <- inputVector[3]
  prMax <- inputVector[4]
  prInit <- inputVector[5]
  
  ###
  # create a vector to store the season's precipitation values
  P <- rep(0, lSeason)
  
  ###
  # loop through days in a season to determine rainfall
  #  sadly, I could not find a way to vectorize this without breaking something
  for (n in 1:lSeason) {
    
    # set a random number
    pr <- runif(1)
    
    # determine probability based on memory effect
    if (n > tau){
      p <- (sum(P[(n-tau):(n-1)]) / tau - pWeak) / (pStrong - pWeak)
	  } else {
	  	p <- prInit
	  }
	
	  # limit high and low probabilities based on max probability
	  if (p > prMax) {
	  	p <- prMax
	  } else if (p < (1-prMax)) {
	  	p <- 1-prMax
	  }

  	# assign precip values based on probability deviation from random
  	if (pr < p) { 
  		P[n] <- pStrong
  	} else {
  		P[n] <- pWeak
  	}
  }
  
  # calculate the mean precipitation
  pMean <- mean(P)
  
  # return mean precip for the wet season
  return(pMean)
}
#############################
# PROBLEM 1
# below is the script to generate the model results

# set default parameters to test
pStrong <- 9.
pWeak <- 0.
tau <- 17
prMax <- 0.8
prInit <- 0.75
testVector <- c(pStrong, pWeak, tau, prMax, prInit)

# calculate the mean precipitation for a number of years to run
pastYears <- 6030
futureYears <- 250

# create a vector to store the outputs for the first model runs
defaultVector <- rep(NA, pastYears)

# set the seed
set.seed(5489)

# loop through each year and calculate Pmean
for (i in 1:pastYears){
  defaultVector[i] <- ismRainfall(testVector)
}

# report the first and 100th year's precipitation
print(paste("First year default precip :", defaultVector[1], "mm/day"))
print(paste("100th year default precip :", defaultVector[100], "mm/day"))

#############################
# PROBLEM 2
# set up custom vectors and matrices for four different scenarios

# 1. current climate
s1vec <- c(9.0, 0.0, 17, 0.8, 0.75)
s1 <- rep(0, pastYears)

# 2. 2150-2200 climate
s2vec <- c(10.9, 1.9, 17, 0.82, 0.2)
s2 <- rep(0, futureYears)

# 3. high atmospheric saturation
s3vec <- c(10.9, 1.9, 17, 0.8, 0.75)
s3 <- rep(0, pastYears)

# 4. changing sea level pressure
s4vec <- c(9.0, 0.0, 17, 0.8, 0.2)
s4 <- rep(0, pastYears)

# run each of the four scenarios 
for (i in 1:pastYears){
  s1[i] <- ismRainfall(s1vec)
  s3[i] <- ismRainfall(s3vec)
  s4[i] <- ismRainfall(s4vec)
}
for (i in 1:futureYears){
  s2[i] <- ismRainfall(s2vec)
}

#############################
# create the four-panel figures

# create a four-panel layout, leaving room for a title
par(mfcol = c(2,2), oma = c(0,0,2,0))

# we'll set the min/max bounds for the plots based on the range of all outputs
xmin <- min(c(min(s1), min(s2), min(s3), min(s4)))
xmax <- max(c(max(s1), max(s2), max(s3), max(s4)))

# create a vector for bins to use in the histograms
breaks <- c((xmin*5):((xmax+0.2)*5) * 0.2)

# set up variables for the plot labels
xlab <- "mm / day"
ylabPast <- paste("days out of", pastYears)
ylabFuture <- paste("days out of", futureYears)

# set the colors for each plot
s1color <- 'light blue'
s2color <- 'orange'
s3color <- 'salmon'
s4color <- 'green'

# set the titles for each plot
s1title <- "Current Climate"
s2title <- "2150-2200 Climate"
s3title <- "High Atmospheric Sat."
s4title <- "Delta Sea Level Pressure"

# plot each of the four histograms, adding quartile lines and the means
s1plot <- hist(s1, breaks = breaks, xlab = xlab, ylab = ylabPast, main = s1title, col = s1color)
abline(v = quantile(s1)[2:4], lty = 2)
abline(v = mean(s1), lwd = 3)

s3plot <- hist(s3, breaks = breaks, xlab = xlab, ylab = ylabPast, main = s3title, col = s3color)
abline(v = quantile(s3)[2:4], lty = 2)
abline(v = mean(s3), lwd = 3)

s2plot <- hist(s2, breaks = breaks, xlab = xlab, ylab = ylabFuture, main = s2title, col = s2color)
abline(v = quantile(s2)[2:4], lty = 2)
abline(v = mean(s2), lwd = 3)

s4plot <- hist(s4, breaks = breaks, xlab = xlab, ylab = ylabPast, main = s4title, col = s4color)
abline(v = quantile(s4)[2:4], lty = 2)
abline(v = mean(s4), lwd = 3)

# and add a title
title('Schewe-Levermann Indian Summer Monsoon Rainfall', outer = TRUE)

# that's all, folks!