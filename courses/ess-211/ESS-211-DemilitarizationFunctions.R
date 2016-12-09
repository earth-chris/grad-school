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

# set up a sinusoidal growth model
growth.sin <- function(x, y, sd, per = NA){
  
  # get indices for good data
  growth.inds <- which(is.na(y) == FALSE)
  
  # generate some noise
  noise <- rnorm(length(x), sd = sd)
  
  # set x and y vars
  xOrig <- x
  yOrig <- y
  x <- x[growth.inds]
  y <- y[growth.inds]
  
  # calculate the spectrum and period of the data, if not set
  ssp <- spectrum(y, plot=FALSE)
  if (is.na(per)){
    per <- 1 / ssp$freq[ssp$spec == max(ssp$spec)]
  }
  
  # least squares fit the data
  linearModel <- lm(y ~ sin(2 * pi / per * x) + cos(2 * pi / per * x))
  
  # apply the fit to the data and add noise
  predicted <- predict(linearModel, data.frame(x = xOrig, y = yOrig))
  sinModel <- predicted + noise
  return(sinModel)
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
# inflation rate functions

# takes historic spending data and adjusts it to 2016 USD
inflation.adjust <- function(values, rates){
  
  # find the data to use
  inds <- which(is.na(rates) == FALSE)
  
  # loop through each year and adjust
  for (i in 1:(length(inds)-1)){
    for(j in (i+1):length(inds)){
      values[i] <- growth.exponential(values[i], rates[j], 1)
    }
  }
  return(values)
}

#############################
# functions to calculate demilitarization rates

# linear demilitarization rate
demil.linear <- function(enlistmentT0, minEnlistment, nYears){
  
  # get the linear rate
  demilRate <- (enlistmentT0 - minEnlistment) / nYears
  
  # calculate the output enlistment numbers
  enlistment <- seq((nYears-1), 0) * demilRate
  return(enlistment)
}

# exponential demilitarization rate
demil.exp <- function(enlistmentT0, minEnlistment, nYears){
  
  # get the exponential rate
  demilRate <- growth.calc.r(enlistmentT0, minEnlistment, nYears + 1)
  
  # calculate the output enlistment numbers
  enlistment <- growth.exponential.apply(enlistmentT0, demilRate, nYears)
}

#############################
# ok, so now that we have veteran population data to analyze, we'll start running our actual model, where we vary
#  a) when demilitarization occurs
#  b) over how man years it occurs
#  c) whether we use historic or modeled data
#
# the idea is that we'll have two potential model outputs to measure the effect of demilitarization
#  a) the highest ratio of people added to the labor force times the amount of money freed from demilitarization, and 
#  b) the total money saved on veteran's costs (for n years beyond the start of demilitarization)

# create a function to calculate labor index
labor.index <- function(newLabor, population, laborParticipation, unemployment, newSpending, fedSpending){
  
  # calculate number of unemployed
  nUnemployed <- population * laborParticipation * unemployment
  
  # calculate the labor index
  laborIndex <- (newSpending / fedSpending) / (newLabor / nUnemployed) 
  return(laborIndex)
}

# create a version for sensitivity analysis
labor.index.sa <- function(params){
  
  # extrapolate the values based on a single vector input
  newLabor <- params[1] 
  population <- params[2] 
  laborParticipation <- params[3] 
  unemployment <- params[4] 
  newSpending <- params[5] 
  fedSpending <- params[6]
  
  # calculate number of unemployed
  nUnemployed <- population * laborParticipation * unemployment
  
  # calculate the labor index
  laborIndex <- (newSpending / fedSpending) / (newLabor / nUnemployed) 
  return(laborIndex)
}

# this is the big model
demil.model <- function(nYears, nVetYears, demilFunction, minEnlistment, deathDelistmentRates, data, returnParameter, plot = FALSE){
  # where:
  #  nYears = the number of years over which demilitarization occurrs
  #  demilFuntion = the method used to calculate demilitarization (e.g. linear, exponential)
  #  minEnlistment = the minimum number of service members to keep (can be set to ~0)
  #  data = the data frame with expected variable types (e.g. $GDP.BUSD, etc.)
  #  returnParameter specifies whether to return veterans spending data
  #   or laborIndex data
  
  # get years to process
  nRows <- length(data$Year)
  startYear <- data$Year[1]
  endYear <- data$Year[nRows-nVetYears]
  startInd <- 1
  endInd <- which(data$Year == endYear)
  
  # calculate the baseline veteran population
  data$TotalVeterans[1] <- deathDelistmentRates[3]
  data$TotalVeterans[2:nRows] <- veterans.deathDelistmentRates(deathDelistmentRates, years = (startYear+1):(startYear+nRows-1), yearlyData = data)
  
  # convert populations from millions of people to just # people
  data$TotalPopulation <- data$TotalPopulation * 10e6
  
  # calculate the business as usual per-capita veteran spending
  data['PerCapitaVeteranSpending'] = data$VeteranSpending.BUSD / data$TotalVeterans
  
  # create a new column for business as usual veteran populations
  data['OriginalVeterans'] <- data$TotalVeterans
  
  # create a matrix to return of mean/stdev labor indices
  DLI <- matrix(nrow = length(startInd:endInd), ncol = 2)
  
  # and a vector to return veteran's costs
  veteranCosts <- vector(length = length(startInd:endInd))
  
  # loop through each year and calculate metrics
  for (yearInd in startInd:endInd){
    
    # create a new data frame to store outputs
    dataModel <- data
    
    # calculate the number of soldiers delisting
    delistSoldiers <- demilFunction(data$TotalEnlistment[yearInd], minEnlistment, nYears)
    
    # don't consider fractions of humans, because that is bad.
    delistSoldiers <- ceiling(delistSoldiers)
    
    # set the delistment data in the data frame, and set future enlistment to the minimum enlistment
    dataModel$TotalEnlistment[yearInd:(yearInd + nYears - 1)] <- delistSoldiers
    dataModel$TotalEnlistment[(yearInd + nYears):nRows] <- minEnlistment
    
    # do the same for military spending, but first, calculate per-capita mil spending at the time
    #  to determine what minimum spending will be at minimum enlistment
    minSpending <- data$MilitarySpending.BUSD[yearInd] / data$TotalEnlistment[yearInd]
    delistSpending <- demilFunction(data$MilitarySpending.BUSD[yearInd], minSpending, nYears)
    dataModel$MilitarySpending.BUSD[yearInd:(yearInd + nYears - 1)] <- delistSpending
    dataModel$MilitarySpending.BUSD[(yearInd + nYears):nRows] <- minSpending
    
    # calculate the number of soldiers entering the economy per year
    enteringSoldiers <- vector(length = nYears)
    enteringSoldiers[1] <- data$TotalEnlistment[yearInd] - delistSoldiers[1]
    for (i in 2:nYears){
      enteringSoldiers[i] <- delistSoldiers[i-1] - delistSoldiers[i]
    }
    
    # and calculate the amount of money entering the economy
    enteringSpending <- vector(length = nYears)
    enteringSpending[1] <- data$MilitarySpending.BUSD[yearInd] - delistSpending[1]
    for (i in 2:nYears){
      enteringSpending[i] <- delistSpending[i-1] - delistSpending[i]
    }
    
    # then calculate the labor index for each year
    laborIndex <- vector(length = nYears)
      for (i in 1:nYears){
      laborIndex[i] <- labor.index(enteringSoldiers[i], data$TotalPopulation[i], data$LaborParticipationRate[i],
                                   data$UnemploymentRate[i], enteringSpending[i], data$FedSpending.BUSD[i])
    }
    
    # return mean/sd DLI
    DLI[yearInd, 1] <- mean(laborIndex)
    DLI[yearInd, 2] <- sd(laborIndex)
    
    # the number of service members entering the economy is also the number of people that enter the 
    #  veteran population. add them here
    dataModel$TotalVeterans[yearInd:(yearInd + nYears - 1)] <- dataModel$TotalVeterans[yearInd:(yearInd + nYears - 1)] + enteringSoldiers
    
    # then apply the death rate to forecast future numbers
    dataModel$TotalVeterans[(yearInd + nYears):nRows] <- growth.exponential(dataModel$TotalVeterans[(yearInd + nYears)], deathDelistmentRates[1], length((yearInd + nYears):nRows))
    
    # calculate the difference in veteran spending
    vetDifferences <- data$PerCapitaVeteranSpending * (data$TotalVeterans - dataModel$TotalVeterans)
    
    # and report the total costs for n years ahead
    veteranCosts[yearInd] <- sum(vetDifferences[yearInd:(yearInd + nVetYears)])
  }
  
  # return the parameter set
  if (returnParameter == 'DLI'){
    
    # plot the output
    if (plot == TRUE){
      title <- paste0("Demilitarized Labor Index - ", nYears, " Year Demilitarization")
      xlab <- "Year"
      ylab <- "DLI"
      pch <- 19
      col <- add.alpha("Black", 0.8)
      cex <- 2
      plot(data$Year[startInd:endInd], DLI[startInd:endInd], xlab = xlab, ylab = ylab, main = title, col = col, pch = pch, cex = cex)
    }
    # return the parameters
    return(DLI)
  } else {
    
    # plot the output
    if (plot == TRUE){
      title <- paste0("Total Veterans Savings After ", nVetYears, " Years")
      xlab <- "Year"
      ylab <- "Billions of Dollars"
      pch <- 19
      col <- add.alpha("Dark Green", 0.8)
      cex <- 2
      plot(data$Year[startInd:endInd], veteranCosts[startInd:endInd], xlab = xlab, ylab = ylab, main = title, col = col, pch = pch, cex = cex)
    }
    return(veteranCosts)
  }
}