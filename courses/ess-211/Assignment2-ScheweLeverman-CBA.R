# ESS-211 assignment 2, due 10-19-16
#  cba 10/2016
#

#############################
# below are the functions used for this assignment

# set up a function to calculate memory effect using time from tau
ismMemory <- function(day, tau, pSeason)

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
  # we're going to try vectorizing this whole code and use which()
  #  statements so no for/while/ifs are necessary
  
  # create a vector that will store the daily precip. for the year
  pSeason <- rep(NA, lSeason)
  
  # create a vector to store the random numbers
  pr <- runif(lSeason)
  
  # return cumulative precip
  return(pSum)
}
#############################
# below is the script to generate the model results

# set the seed so model outputs are deterministic
set.seed(5489)