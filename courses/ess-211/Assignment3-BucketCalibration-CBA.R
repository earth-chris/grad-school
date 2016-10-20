# ESS-211 assignment 3, due 10-26-16
#  cba 10/2016

#############################
# functions provided by teaching group for assignment 3

# 1-parameter model (bucket capacity)
bucket1 <- function(pars, P, PET) {
  # in this model:
  # pars : the input parameter values - one element - the bucket's maximum capacity
  # P    : precipitation vector (i.e. water in)
  # PET  : potential evapotranspiration vector (i.e. water out)
  
  # set max capacity for the bucket
  full.bucket <- pars[1]
  
  # create a vector to store runoff based on the # of precip values we have
  n <- length(P)
  runoff <- rep(0, n)
  
  # set the starting state of the bucket as the lower of a) the bucket's max capacity or b) the 
  #  difference between precip and PET (unless thats < 0)
  bucket <- min(full.bucket, max(P[1] - PET[1], 0))
  
  # loop through each time step in the precip/PET vectors
  for (m in 2:n) {
    
    # I don't know why these lines exist since they aren't used elsewhere
    3
    infilt <- P[m]
    AET <- PET[m]
    
    # fill/drain the bucket based on the difference in precip and PET
    bucket <- bucket + P[m] - PET[m]
    
    # send bucket overflow to runoff, or to 0 if it empties
    if (bucket > full.bucket) {
      runoff[m] <- bucket - full.bucket
      bucket <- full.bucket
    } else if (bucket<0) {
      bucket <- 0
    }
  }
  
  # return our output runoff values
  return(runoff)
}

# 3-parameter model (bucket capacity, varying P scaling factor, varying PET scaling factor)
bucket3 = function(pars,P,PET) {
  # in this model:
  # pars : the input parameter values - 3 elements 
  #        [1] the bucket's maximum capacity
  #        [2] the minimum infiltration (i.e. the fraction of precip that doesn't infiltrate the soil)
  #        [3] the wilting point (i.e. when soil dries out)
  # P    : precipitation vector (i.e. water in)
  # PET  : potential evapotranspiration vector (i.e. water out)
  
  # set variables from the input parameters
  full.bucket <- pars[1]
  min.infilt <- pars[2]
  wilt.point <- pars[3]
  
  # initialize the bucket values
  bucket <- min(full.bucket, max(P[1] - PET[1], 0)) 
  
  # create an actual evapotranspiration vector. but why? it gets redefined later...
  n <- length(P)
  AET <- runoff <- rep(0,n)
  
  # loop through each time step in the precip/PET vectors
  for (m in 2:n) {
    
    # set the saturation level
    beta = bucket/full.bucket 
    
    # update infiltration fraction based on how full the bucket is
    infilt.frac <- (1 - min.infilt) + beta * (2 * min.infilt - 1)
    
    # update actual ET fraction based on how full the bucket is
    aet.frac <- max(0, min(1, (beta - wilt.point) / (1 - 2 * wilt.point)))
    
    # set actual infiltration based on precip and the fraction of infiltration
    infilt <- P[m] * infilt.frac
    
    # set runoff as difference between precip and infiltration
    runoff[m] <- P[m] - infilt
    
    # set actual ET as product of potential ET and the fraction of actual ET
    AET <- PET[m] * aet.frac
    
    # fill/drain the bucket
    bucket <- bucket + infilt - AET
    
    # send bucket overflow to runoff, or to 0 if it empties
    if (bucket > full.bucket) {
      runoff[m] <- runoff[m] + bucket - full.bucket
      bucket <- full.bucket
    } else if (bucket < 0) {
      bucket <- 0
    }
  }
  
  # return our output runoff values
  return(runoff)
}