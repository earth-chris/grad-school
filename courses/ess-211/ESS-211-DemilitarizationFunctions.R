# ESS-211 Final Project
# Modeling de-militarization in the US
# Christopher Anderson
# 12-2016

#############################
# below are the functions used in the ESS-211-DemilitarizationProject.R script
library(car)

#############################
# first, functions for modeling federal spending

# set up a function for polynomial growth
#  based on fitting a polynomial using R's built in functions
#  returns a model
growth.polynomial <- function(x, y, degree=2){
  model <- lm(y ~ poly(x, degree = degree))
  return(model)
}

# set up a function for exponential growth using
#  the p = p0 * e ^ r*t model
growth.exponential <- function(y0, r, t){
  y <- y0 * exp(r * t)
  return(y)
}

# and a function to apply it
growth.exponential.apply <- function(y0, r, nYears){
  
  # create an output storage vector
  y <- rep(0, nYears)
  for (i in 1:nYears){
    y[i] <- growth.exponential(y0, r, i)
  }
  return(y)
}

# set up function to derive the growth rate from data
growth.calc.r <- function(y0, y, t){
  r <- log(y / y0) / t
  return(r)
}

# set up a loss function to find best fit of r
growth.minfunc <- function(lossfunc, model, y, y0, r, nYears){
  lossfunc(y, model(y0, r, nYears))
}

# set up a function to calibrate exponential growth by looking at all
#  growth rates based on each starting point and each year
#  to get the range of growth rates to fit with
calibrate.growth <- function(years, population){
  
  # determine number of starting years we'll use
  nYears <- length(years) - 1
  
  # create an output vector to store results
  r.guess <- rep(0, sum(seq(1:nYears)))
  
  # set up a counter to use for each iteration
  counter = 1
  
  # loop through each start year
  for (i in 1:nYears){
    
	  # find the starting year for this iteration
	  y0 <- years[i]
	
	  # find the index to use
	  y0ind = which(years == y0)
	
	  # loop through each combination of years and get r
	  for (j in (y0ind+1):(nYears+1)){
	    r.guess[counter] <- growth.calc.r(population[y0ind], population[j], j-y0ind)
	    counter = counter + 1
	  }
  }
  
  # run the optimization using min and max r, using the mean as the starting guess
  opt <- optim(mean(r.guess), growth.minfunc, loss = loss.rmse, model = growth.exponential.apply,
               y = population[2:(nYears+1)], y0 = population[1], nYears = nYears, method = 'L-BFGS-B',
               lower = min(r.guess), upper=max(r.guess))
  r.best <- opt$par
  return(r.best)
}


# set up a function for logistic growth
#  based on the the logistic growth function:
#  y = phi1/(1 + exp(-(phi2 + phi3 * x)))
#  where y is the population, phi1 is the asymptote (or carrying capacity),
#  i don't remember phi2, phi3 is the growth parameter, and 
#  x is the response variable, in this case, time
growth.logistic <- function(x, y, phi1, phiIndex){
  # we're going to use the index provided to insert the phi1 value
  #  so the logistic growth fits to the asymptote
  y[phiIndex:length(y)] <- phi1
  
  # we're going to have to guess the initial parameters for 
  #  phi2 and phi3, which we will do using a logit transform
  coefs <- coef(lm(logit(y / phi1) ~ x))
  phi2 <- coefs[1]
  phi3 <- coefs[2]
  
  # then we'll perform a non-linear least-squares fit to the data
  model <- nls(y ~ (phi1 / (1 + exp(-(phi2 + phi3 * x)))), 
               start = list(phi1 = phi1, phi2 = phi2, phi3 = phi3))
  
  return(model)
}

#############################
# unemployment modeling

# framework 1
#  finds the mean of the unemployment data, then adds normally distributed noise based on the 
#  standard deviation provided
unemployment.framework1 <- function(dframe, sd){
  
  # get basic unemployment info
  unemployment.years <- which(dframe$UnemploymentRate < 1.)
  unemployment.mean <- mean(dframe$UnemploymentRate[unemployment.years])
  
  # generate some noise
  noise <- rnorm(nrow(dframe), sd = sd)
  
  # add the noise to the mean
  framework1 = unemployment.mean + noise
  return(framework1)
}

# framework 2
#  finds the trendline of the unemployment data, then adds normally distributed noise
#  based on the standard deviation provided
unemployment.framework2 <- function(dframe, sd){
  
  # get basic unemployment info
  unemployment.years <- which(dframe$UnemploymentRate < 1.)
  
  # set x and y vars
  x <- dframe$Year[unemployment.years]
  y <- dframe$UnemploymentRate[unemployment.years]
  
  # find trendline in the unemployment data
  linearModel <- lm(y ~ x)
  
  # generate some noise
  noise <- rnorm(nrow(dframe), sd = sd)
  
  # apply the fit to the data and add noise
  predicted <- predict(linearModel, data.frame(x = dframe$Year, y = dframe$UnemploymentRate))
  framework2 <- predicted + noise
  return(framework2)
}

# framework 3 
#  fits a sin curve to unemployment data
#  from this stack exchange: http://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
unemployment.framework3 <- function(dframe, sd, per = NA){
  
  # get basic unemployment info
  unemployment.years <- which(dframe$UnemploymentRate < 1.)
  
  # set x and y vars
  x <- dframe$Year[unemployment.years]
  y <- dframe$UnemploymentRate[unemployment.years]
  
  # calculate the spectrum and period of the data, if not set
  ssp <- spectrum(y, plot=FALSE)
  if (is.na(per)){
    per <- 1 / ssp$freq[ssp$spec == max(ssp$spec)]
  }
  
  # least squares fit the data
  linearModel <- lm(y ~ sin(2 * pi / per * x) + cos(2 * pi / per * x))
  
  # generate some noise
  noise <- rnorm(nrow(dframe), sd = sd)
  
  # apply the fit to the data and add noise
  predicted <- predict(linearModel, data.frame(x = dframe$Year, y = dframe$UnemploymentRate))
  framework3 <- predicted + noise
  return(framework3)
}

#############################
# modeling veteran populations based on death / de-enlistment rates. 
#  will calibrate this model to derive the first parameter set
veterans.deathDelistmentRates <- function(deathDelistmentRates, years=c(2013), yearlyData=yearly){
  
  # the deathDelistment rates contain the vectors for 1) the exponential rate of death and 
  #  2) the linear delistment rate, and 3) the veteran population at t0
  deathRate <- deathDelistmentRates[1]
  delistmentRate <- deathDelistmentRates[2]
  vetPopulationT0 <- deathDelistmentRates[3]
  
  # set up the output vector based on the number of years we are analyzing
  nVeterans <- rep(0, length(years))
  
  # loop through each year we're trying to get veteran population for
  for (i in seq(1,length(years))){
    
    # get the right indices
    yearInds <- which(yearlyData$Year == years[i])
    nYears <- years[i] - yearlyData$Year[1]
    
    # calculate the delistment per year
    delistment <- delistmentRate * yearlyData$TotalEnlistment[1:yearInds]
    
    # set the veteran population to add to/subtract from each year
    vetPopulation <- vetPopulationT0
    
    # loop through each year and add/subtract vets based on delistment/death
    for (j in seq(1, nYears)){
      vetPopulation <- growth.exponential(vetPopulation, deathRate, 1)
      vetPopulation <- vetPopulation + delistment[j]
    }
    
    # calculate the number of veteran as (delistment rate * yearly enlisted population) + t0population * e^rt
    nVeterans[i] <- vetPopulation
  }
  return(nVeterans)
}

# set up a loss function for the veterans model
veterans.minfunc <- function(pars, loss, model, y, years, yearlyData){
  loss(y, model(pars, years, yearlyData))
}

# set up a function to calibrate the veterans model
veterans.calibrateRates <- function(mins, maxs, nGuesses, nFolds, yearlyData){
  
  # create a list that stores the output data
  optList = list()
   
  # create the cross-folds
  yIndices <- which(is.na(yearlyData$TotalVeterans) == FALSE)
  yGroups <- cut(1:length(yIndices), nFolds, label = FALSE)
  
  # loop through each fold and calibrate for each year
  for (i in 1:nFolds){
    
    # create an array to store the outputs
    optMatrix <- array(dim = c(nGuesses, 3))
    
    # create a series of random pulls from the min/max of each parameter
    random.par1 <- runif(nGuesses, min = mins[1], max = maxs[1])
    random.par2 <- runif(nGuesses, min = mins[2], max = maxs[2])
    random.par3 <- runif(nGuesses, min = mins[3], max = maxs[3])
    
    # loop through each year and calibrate the model
    for (j in 1:nGuesses){
      
      # set the parameter guesses
      par0 <- c(random.par1[j], random.par2[j], random.par3[j])
      
      # run the optimization
      opt <- optim(par0, veterans.minfunc, loss = loss.rmse, model = veterans.deathDelistmentRates,
                   years = yearlyData$Year[yIndices[yGroups != i]], yearlyData = yearlyData, 
                   y = yearlyData$TotalVeterans[yIndices[yGroups != i]], method='L-BFGS-B',
                   lower = mins, upper = maxs)
      
      # save the results to our output matrix
      optMatrix[j,] <- opt$par
    }
    
    # then save these to the output list
    optList[[i]] <- optMatrix
  }
  return(optList)
}


#############################
# basic model functions

# a function to get the indices for years to model
#seq(year.start, year.end-year.demil, by=year.demil)