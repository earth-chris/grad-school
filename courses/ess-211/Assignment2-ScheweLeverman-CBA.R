# ESS-211 assignment 2, due 10-19-16
#  cba 10/2016
#

# clear the workspace
rm(ls = ls()

#############################
# below are the functions used for this assignment

# set up the ISM rainfall function as translated from matlab
ismRainfall <- function(pStrong, pWeak, tau, prMax, prInit,
                        lSeason=135){
  # pStrong  : precipitation in wet state (mm/day)
  # pWeak    : precipitation in dry state (mm/day)
  # tau      : memory length in time steps (days)
  # prMax    : maximum probability of either state.
  # pInit    : initial probability of strong state.
  # lSeason  : length of the wet season (days). default = 135
  #            must be <= prMax
  
  ###
  # we're going to try vectorizing this whole code and use which()
  #  statements so no for/while/ifs are necessary
  
  # create a vector to store the season's precipitatoin values
  pSeason <- rep(NA, lSeason)
  
  # create a vector to store random probabilities
  pr <- runif(lSeason)
  
  # create a vector to store modeled wet/dry probabilities
  p <- rep(NA, lSeason)
  
  ###
  # determine probabilities based on memory effect
  
  # set [tau] number of days to initial probability of strong state
  p[1:tau] <- pInit
  
  # loop through all values after tau to calculate memory-based probability
  for (i in (tau+1):lSeason){
    p[i] <- ((sum(p[(i-tau):i-1]) / tau) - pWeak) / (pStrong - pWeak)
  }
  
  ###
  # limit the high and low probabilities based on maximum probability threshold
  p[p > prMax] <- prMax
  p[p < (1-prMax)] <- 1-prMax
  
  ###
  # assign precipitation values
  pSeason[pr < p] <- pStrong
  pSeason[pr >= p] <- pWeak
  
  # calculate the mean precipitation
  pMean <- sum(pSeason) / lSeason
  
  # return mean precip for the wet season
  return(pMean)
  
  # return the full season's precip vector
  #return(pSeason)
}
#############################
# below is the script to generate the model results

# set the seed so model outputs are deterministic
set.seed(5489)

# set default parameters to use
pStrong <- 9.
pWeak <- 0.
tau <- 17
prMax <- 0.8
prInit <- 0.75

# calculate the mean precipitation for a number of years to run
nYears <- 6030

test <- apply(c(pStrong, pWeak, tau, prMax, prInit), 1, ismRainfall)

# report the first and 100th year's precipitation
print(paste0("First year default precip : ", test[1]))
print(paste0("100th year default precip : ", test[100]))

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