# contains frequently used functions for ESS-211

#############################
# functions provided by teaching group for assignment 1 (PET)

# saturated vapor pressure for a given temperature
pet.e0 <- function(temp){
  0.6108*exp(17.27*temp/(temp+237.3))
}

# this is how FAO suggests one do it for daily esat b/c of nonlinearity of e0 function 
pet.esat <- function(tmax,tmin){ 
  (pet.e0(tmax)+pet.e0(tmin))/2  
}

# slope of saturation curve
pet.s <- function(temp){
  4098*(0.6108*exp(17.27*temp/(temp+237.3)))/(temp+237.3)^2    
}

# extraterrestrial net radiation (MJ / m^2 / day)
pet.Ra <- function(doy,lat){   #lat in decimal degrees       #p79 FAO
  psi=(pi/180)*lat
  delta=0.409*sin((2*pi*doy/365)-1.39)
  ws=acos(-tan(psi)*tan(delta))
  dr=1+0.033*cos(2*pi*doy/365)
  Ra=(24*60/pi)*0.082*dr*(ws*sin(psi)*sin(delta)+cos(psi)*cos(delta)*sin(ws))
  Ra
}

# psychrometric constant [kPa C-1]
pet.psychro <- function(elev){             #elev in m
  press=101.3*((293-0.0065*elev)/293)^5.26
  0.665e-3*press
}

# daylength (returns # hours in the day with sunlight. only for areas between 65N-6)
pet.daylength <-function(lat, doy){
  psi=(pi/180)*lat
  delta=0.409*sin((2*pi*doy/365)-1.39)
  ws=acos(-tan(psi)*tan(delta))
  d=24*ws/pi
  d
}

# net surface radiation
pet.Rn <- function(lat, elev, doy, dew, tmax, tmin, Rs){
  ea=pet.e0(dew)
  Ra=pet.Ra(doy,lat)
  Rso=(0.75+(2e-5)*elev)*Ra          #p85
  Rnl=4.903e-9*(((tmax+273)^4+(tmin+273)^4)/2)*(0.34-0.14*sqrt(ea))*(1.35*(Rs/Rso)-0.35)  #p86
  Rn=(1-0.23)*Rs-Rnl             #(Rns=(1-0.23)Rs) ; Rn=Rns-Rnl p87
  Rn
}

#############################
# functions written by cba for assignment 1 (PET)

# method 1: priestly taylor
pet.priestlyTaylor <- function(doy, tMax, tMin, tMean, RH, tDew, Rs, elev, lat){
  # define constant 'a'
  a <- 1.26
  
  # calculate slope of vapor pressure curve using temp above
  vpSlope <- pet.s(tMean)
  
  # calculate the psychometric constant
  gamma <- pet.psychro(elev)
  
  # set the ground heat flux to zero
  groundHeatFlux <- 0
  
  # calculate net radiation
  netRadiation <- pet.Rn(lat, elev, doy, tDew, tMax, tMin, Rs)
  
  # set lambda
  lambda <- 2.501 - (0.002361 * tMean)
  
  # run the PET calculation
  PET <- a * ((vpSlope * (netRadiation - groundHeatFlux)) / 
                (lambda * (vpSlope + gamma)))
  return(PET)
}

# method 2: modified priestly taylor
pet.modifiedPriestlyTaylor <- function(tMax, tMin, Rs){
  # define constant 'albedo'
  albedo <- 0.23
  
  # calculate constant TD
  TD <- (0.6 * tMax) + (0.4 * tMin)
  
  # calculate EEQ
  EEQ <- Rs * (4.88e-3 - (4.37e-3 * albedo)) * (TD + 29)
  
  # PET is calculated differently based on tMax
  if (tMax < 5){
    PET <- EEQ * 0.01 * exp(0.18 * (tMax + 20))
  } else if (tMax > 24){
    PET <- EEQ * ((tMax -24) * 0.05 + 1.1)
  } else {
    PET <- EEQ * 1.1
  }
  return(PET)
}

# method 3: hammon
pet.hammon <- function(doy, tMax, tMin, tMean, lat){
  # calculate the day length
  dl <- pet.daylength(lat, doy)
  
  # transform day length to fraction of a day (i.e. 0-1 with 0 as no sunlight, 1 as all-day light
  dl <- dl / 24
  
  # calculate the saturated vapor pressure
  vpSat <- pet.esat(tMax, tMin)
  
  # run the PET calculation
  PET <- 715.5 * dl * (vpSat / (tMean + 273.2))
  return(PET)
}

# method 4: hargreaves
pet.hargreaves <- function(doy, tMax, tMin, tMean, lat){
  # calculate lambda
  lambda <- 2.501 - (0.002361 * tMean)
  
  # calculate atmospheric net radiation
  atmRad <- pet.Ra(doy, lat)
  
  # run the PET calculation
  PET <- (0.0023 * (tMean + 17.8) * ((tMax - tMin)^0.5) * atmRad) / lambda
  return(PET)
}

# method 5: linacre
pet.linacre <- function(tMean, elev, lat, tDew){
  # calculate Tm
  Tm <- tMean + (0.006 * elev)
  
  # run the PET calculation
  PET <- (500 * (Tm / (100 - lat)) + (15 * (tMean - tDew))) / (80 - tMean)
  return(PET)
}

# method 6: turc
pet.turc <- function(tMean, RH, Rs){
  # convert Rs from units of MJ/m^2 to cal/cm^2
  RsCal <- Rs * 23.9
  
  # convert RH to 1-100 scale
  #RH = RH * 100
  
  # PET is calculated differently based on relative humidity
  if (RH < 50){
    PET <- (1 + ((50 - RH) / 70)) * ((0.013 * tMean * (RsCal + 50)) / (tMean + 15))
  } else {
    PET <- (0.013 * tMean * (RsCal + 50)) / (tMean + 15) 
  }
  return(PET)
}

#############################
# functions written by cba for assignment 2

# set up the ISM rainfall function as translated from matlab
ism.rainfall <- function(inputVector, lSeason=135){
  # the input vector should be a 5-element vector, in this order:
  # [1] pStrong  : precipitation in wet state (mm/day)
  # [2] pWeak    : precipitation in dry state (mm/day)
  # [3] tau      : memory length in time steps (days)
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

#############################
# functions provided by teaching group for assignment 4

# calculate cross validation rmse
rmse <- function(train, test) {
  
  # perform a linear regression
  fit <- lm(y~., data = train) 
  
  # calculate rmse of training fit
  train.err <- sqrt(mean((fit$resid)^2))
  
  # calculate the rmse of training vs test data
  test.err = sqrt(mean((test$y - predict(fit, test))^2)) #RMS CV RESIDUALS
  
  # return the training and test error
  return(c(train.err, test.err))
}

# and one written by cba for assignment 4
# a function to calculate rmse for k splits
cv.err <- function(data.frame, K){
  # split the data into ~equal chunks
  #  first, find the length of the data frame to split
  nr <- nrow(data.frame)
  
  # split the data using the cut funtion
  groups <- cut(1:nr, K, label = FALSE)
  
  # create a matrix that will contain the 2 x K training and test errors
  errorMatrix <- matrix(ncol = 2, nrow = K)
  
  # loop through each split and calculate training and test errors
  for (i in 1:K){
    errorMatrix[i,] <- rmse(data.frame[groups != i,], data.frame[groups == i,])
  }
  
  # return the matrix
  return(errorMatrix)
}

#############################
# functions provided by the teaching group for assignment 5

# define the linear oscillator equation for determining bungee height
hmin <- function (x){
  return(x[1] - 2 * x[2] * 9.8 / (1.5 * x[3]))
}

# define function to sample from set {0, 1/(p -1), ... 1}
base.samp <- function(p , n){
  # limit to range of (0, 1-delta)
  x <- (p - 1 - p / 2) * runif(n)
  return(round(x) / (p -1))
}

# define the morris function
morris <- function(k, r, p, delta, par.mins, par.maxs, fun){
  
  # set the range for parameters
  par.range <- par.maxs - par.mins
  
  # create an array to save the effects
  effects <- array(dim = c(k, r))
  
  # loop through each computed effect
  for (r.ind in 1:r){
    
    # create an array of ones
    J <- array(1, dim = c(k +1 , k)) 
    
    # set a lower triangular array of ones
    B <- lower.tri(J) * 1 
    
    # set base vector
    xstar <- base.samp(p, k)
    
    # set array with 1 or -1 in diagonal
    D <- array(0, dim = c(k, k)) 
    diag(D) <- 1.0 - 2 * (runif(k) < .5)
    
    # set random permutation array with 1 in diagonal
    P <- array(0, dim = c(k, k))
    diag(P) <- 1
    P <- P[, sample(c(1:k), k, replace = FALSE)]
    Bstar <- (t(xstar * t(J)) + (delta / 2) * ((2*B - J)%*%D + J)) %*% P
    
    # set vector for model outputs
    y <- numeric(length = (k +1))
    
	# run model on each parameter permutation
	for (i in 1:(k+1)){
      y[i] <- fun(par.mins + par.range * Bstar [i,])
    } 
	
	# calculate effect based on which parameter changed
    for (i in 1:k){
      par.change <- Bstar[i+1,] - Bstar[i,]
      i2 <- which(par.change != 0)
      effects[i2, r.ind] <- (y[i +1] - y[i]) * (1 - 2 * (par.change[ i2 ] < 0))
    }
  }
  
  # return k by r array of computed effects
  return(effects) 
}

# define variance-based global sensitivity analysis (vsa)
vsa <- function(fun, par.mins, par.maxs, nrun){
  
  # set ranges
  par.range <- par.maxs - par.mins
  k <- length(par.range)
  
  # create U(0,1) n by k matrix
  M <- array(runif(nrun * k), dim = c(nrun, k))
  
  # transform to fit range of each parameter
  for (i in 1:k) {
    M[, i] <- M[, i] * par.range[i] + par.mins[i] 
  }
  
  # make a different input matrix, M2, the same way
  M2 <- array(runif(nrun * k), dim = c(nrun, k))
  for (i in 1:k) {
    M2[, i] = M2[, i] * par.range[i] + par.mins[i] 
  }
  
  # make K different matrices, where in the jth matrix, all parameters are taken from M2 except parameter j is taken from M
  NJ.list <- list()
  for (j in 1:k){
    temp <- M2
    temp[, j] <- M[, j]
    NJ.list[[j]] <- temp
  }
  
  # compute the first-order effects (S) of factor J, we need to compute the values UJ, as discussed in class
  S <- numeric(length = k)
  ST <- numeric(length = k)
  y1 <- y2 <- y3 <- numeric(length = nrun)
  for (i in 1:nrun){
    y1[i] <- fun(M[i,])
    y2[i] <- fun(M2[i,])
  }
  
  # compute the expected value (EY) and total variance (VY) of Y, for which we will use M2
  EY <- mean(y1)
  
  # will estimate EY^2 using both M and M2
  VY <- sum(y1 * y1) / (nrun - 1) - EY^2
  
  for (j in 1:k){
    NJ <- NJ.list[[j]]
    for (i in 1:nrun) {
      y3[i] = fun(NJ[i,])
    }
    
    # everything but factor j resampled
    UJ <- sum(y1 * y3) / (nrun - 1)
    
    # only factor j resampled
    UJ2 = sum(y2 * y3) / (nrun - 1)
    
    S[j] = (UJ - EY^2) / VY
    ST[j] = 1.0 - (UJ2 - EY^2) / VY
  }
  return(list(S,ST))
}

# a way to script the morris calculations
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
  cex <- rep(3, npars)
  cols <- rainbow(npars)
  xlab <- "Mean absolute effect"
  ylab <- "SD effect"
  
  # plot the results
  plot(mu.star, sigma, pch = pch, main = title, col = cols, xlab = xlab, ylab = ylab, cex=cex)
  legend("topleft", legend = names, pch = pch, col = cols)
  
  # return the effects
  return(list(mu, sigma, mu.star))
}

#############################
# functions for assignment 6

# create a function to index a matrix or a vector using another vector
whichVec <- function(toIndex, indexVec){
  nIndices <- length(indexVec)
  outputIndex <- c()
  for (i in seq(1, nIndices)){
    tmpIndices <- which(toIndex == indexVec[i], arr.ind = TRUE)
    if (length(tmpIndices) > 0) outputIndex <- append(outputIndex, tmpIndices)
  }
  return(outputIndex)
}

#############################
# methods to calculate loss functions

# root mean-squared error
loss.rmse <- function(y, yhat){
	sqrt(mean((y-yhat)^2, na.rm=TRUE))
} 

# mean absolute error
loss.abs <- function(y, yhat){
	mean(abs(y-yhat), na.rm=TRUE)
} 
# create a function to calculate loss function
minfunc <- function(pars, loss, model, y, P, PET){
	loss(y, model(pars, P, PET))
}

# create a method to calculate the modal value of a vector
mode <- function(x) {
  xunique <- unique(x)
  xunique[which.max(tabulate(match(x, xunique)))]
}

# create a method to calculate the standard error of a vector
se <- function(x) {
  sd(x) / sqrt(length(x))
}

#############################
# plotting utilities

# add an alpha value to a colour (from https://gist.githubusercontent.com/mages/5339689/raw/2aaa482dfbbecbfcb726525a3d81661f9d802a8e/add.alpha.R)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}