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

# set up function to derive the growth rate from data
growth.calc.r <- function(y0, y, t){
  r <- ln(y0 / y) / t
  return(r)
}

# set up a loss function to find best fit of r
growth.minfunc <- function(lossfunc, model, y, y0, r, t){
  lossfunc(y, model(y0, r, t))
}

# set up a function to calibrate exponential growth by looking at all
#  growth rates based on each starting point and each year
#  to get the range of growth rates to fit with
calibrate.growth <- function(years, population,){
  
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
	for (j in (y0ind+1):(nYears)){
	  r.guess[counter] <- growth.calc.r(population[y0ind], population[j], j-y0ind)
	  counter = counter + 1
	}
	
	# run the optimization using min and max r, using the mean as the starting guess
	opt <- optim(mean(r.guess), growth.minfunc, loss = loss.rmse, model = growth.exponential,
			y = population[2:(nYears+1)], y0 = population[1:nYears], t = 1, method = 'L-BFGS-B',
			lower = min(r.guess), upper=max(r.guess))
	r.best <- opt$par
  }
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

