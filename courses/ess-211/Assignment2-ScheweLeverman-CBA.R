# ESS-211 assignment 2, due 10-19-16
#  cba 10/2016
#

# clear the workspace
rm(ls = ls())

#############################
# below are the functions used for this assignment

# set up the ISM rainfall function as translated from matlab
ismRainfall <- function(inputVector, lSeason=135){
  # the input vector should be a 5-element vector, in this order:
  # [1] pStrong  : precipitation in wet state (mm/day)
  # [2] pWeak    : precipitation in dry state (mm/day)
  # [3[ tau      : memory length in time steps (days)
  # [4] prMax    : maximum probability of either state.
  # [5] pInit    : initial probability of strong state.
  #            must be <= prMax
  # lSeason  : length of the wet season (days). default = 135
  
  # split the input vector for legibility
  pStrong <- inputVector[1]
  pWeak <- inputVector[2]
  tau <- inputVector[3]
  prMax <- inputVector[4]
  prInit >- inputVector[5]
	
  ###
  # create the vectors to store outputs
  
  # create a vector to store the season's precipitatoin values
  P <- rep(0, lSeason)
  
  # create a vector to store random probabilities
  pr <- runif(lSeason)
  
  # create a vector to store modeled wet/dry probabilities
  p <- rep(0, lSeason)
  
  ###
  # determine probabilities based on memory effect
  #  sadly, I could not find a way to vectorize this D:!
  for (i in 1:lSeason) {
    
	if (i > tau){
		p[i] <- (sum(P[(i-tau):i-1]) / tau - pWeak) / (pStrong - pWeak)
	} else {
		p[i] <- prInit
	}
	
	# limit high and low probabilities based on max probability
	if (p[i] > prMax) {
		p[i] <- prMax
	} else if (p[i] < (1-prMax)) {
		p[i] <- 1-prMax
	}

	# assign precip values based on probability deviation from random
	if (pr[i] < p[i]) { 
		P[i] <- pStrong
	} else {
		P[i] <- pWeak
	}
  }
  
  # calculate the mean precipitation
  pMean <- mean(P)
  
  # return mean precip for the wet season
  return(pMean)
  
  # return the full season's precip vector
  #return(pSeason)
}
#############################
# below is the script to generate the model results

# set the seed so model outputs are deterministic
set.seed(5489)

# set default parameters to test
pStrong <- 9.
pWeak <- 0.
tau <- 17
prMax <- 0.8
prInit <- 0.75
testVector = c(pStrong, pWeak, tau, prMax, prInit)

# calculate the mean precipitation for a number of years to run
nYears <- 6030

# we'll apply this as a matrix, so create a 5 x nYears matrix to apply to
test <- apply(t(matrix(rep(testVector, nYears), nrow = length(testVector))), 1, ismRainfall)

# report the first and 100th year's precipitation
print(paste0("First year default precip : ", test[1], "mm/day"))
print(paste0("100th year default precip : ", test[100], "mm/day"))

# set up custom vectors for four different scenarios
# current climate
s1vec <- c(9.0, 0.0, 17, 0.8, 0.75)

# 2150-2200 climate
s2vec <- c(10.9, 1.9, 17, 0.82, 0.2)

# high atmospheric saturation
s3vec <- c(10.9, 1.9, 17, 0.8, 0.75)

# changing sea level pressure
s4vec <- c(9.0, 0.0, 0.8, 0.2)

# run each of the four scenarios using apply
s1 <- apply(rep(s1vec, nYears), 1, ismRainfall)
s2 <- apply(rep(s2vec, nYears), 1, ismRainfall)
s3 <- apply(rep(s3vec, nYears), 1, ismRainfall)
s4 <- apply(rep(s4vec, nYears), 1, ismRainfall)

#############################
# create the four-panel figures

# create a four-panel layout
par(mfcol = c(2,2))

# create a vector for bins to use in the histograms
breaks = c(0:50 * 0.2)

# set up variables for the plot labels
xlab <- "mm / day"
ylab <- paste0("days out of", nYears)

# plot each of the four histograms
s1plot <- hist(s1, breaks = breaks, xlab = xlab, ylab = ylab)
s2plot <- hist(s2, breaks = breaks, xlab = xlab, ylab = ylab)
s3plot <- hist(s3, breaks = breaks, xlab = xlab, ylab = ylab)
s4plot <- hist(s4, breaks = breaks, xlab = xlab, ylab = ylab)