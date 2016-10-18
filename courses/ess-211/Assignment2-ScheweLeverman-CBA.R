# ESS-211 assignment 2, due 10-19-16
#  cba 10/2016
#

#############################
# below are the functions used for this assignment

# set up the ISM rainfall function as translated from matlab
ismRainfall <- function(lSeason, pStrong, pWeak,
                        tau, prMax, pInit){
  # lSeason  : length of the season
  # pStrong  : precipitation in wet state (mm/day)
  # pWeak    : precipitation in dry state (mm/day)
  # tau      : memory length in time steps (days)
  # prMax    : maximum probability of either state
  # pInit    : initial probability of strong state.
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
}
#############################
# below is the script to generate the model results

# set the seed so model outputs are deterministic
set.seed(5489)
