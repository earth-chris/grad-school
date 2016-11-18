# ESS-211 assignment 3, due 10-26-16
#  cba 10/2016

# change to the working directory and read the rdata
#setwd("~/cba/aei-grad-school/courses/ess-211/")
setwd("~/src/aei-grad-school/courses/ess-211/")

# this functions file contains bucket1, bucket3, and the loss functions
source("ESS-211-Functions.R")

# load the provided P, PET, and runoff data
load('course_material/runoff.Rdata')

#############################
# task 1 - finding parameters for bucket3 using runoff data

# first, set the upper and lower bounds of the parameters for infilt and wilt
#  bucket size range 0-200 based on in-class examples
#  infiltration ranges 0-0.5 
#  wilting point ranges 0-0.5
bounds.lower <- c(0.001, 0, 0)
bounds.upper <- c(200, 0.5, 0.5)

# in order to find our parameters, we'll randomly set par0 by sampling
#  the upper and lower bounds using a uniform distribution

# set the number of iterations to try
n_searches <- 100

# create a vector for bucket size, min infilt and wilt point
random.bucket <- runif(n_searches, min = bounds.lower[1], max = bounds.upper[1])
random.infilt <- runif(n_searches, min = bounds.lower[2], max = bounds.upper[2])
random.wilt <- runif(n_searches, min = bounds.lower[3], max = bounds.upper[3])

# create a matrix to store the outputs from each parameter run
opt_matrix <- matrix(nrow = n_searches, ncol = 3)

# report starting
print("Beginning parameter search for bucket3 function")

# run the optimization on series of initial parameter guesses
for (i in seq(1,n_searches)){
  
  # set par0 equal to the randomly assigned starting parameters
  par0 <- c(random.bucket[i], random.infilt[i], random.wilt[i])
  
  # run the optimization for these parameters
  opt <- optim(par0, minfunc, loss=loss.rmse, model=bucket3, y=runoff, P=P, 
              PET=PET, method='L-BFGS-B', lower=bounds.lower ,upper=bounds.upper)
  
  # assign the output params to the storage matrix
  opt_matrix[i,] <- opt$par
  
  # and report the values for each iteration
  print(opt$par)
}

# plot the outputs to examine
par(mfcol = c(1,3))
hist(opt_matrix[,1], breaks = seq(0,200, by=20), xlab = "Retrieved value of full.bucket", main = "full.bucket", col='orange', las=1)
hist(opt_matrix[,2], breaks = seq(0,0.5, by=0.05), xlab = "Retrieved value of min.infilt", main = "min.infilt", col='light green', las=1)
hist(opt_matrix[,3], breaks = seq(0,0.5, by=0.05), xlab = "Retrieved value of wilt.point", main = "wilt.point", col='salmon', las=1)

# the values retrieved from this optimization vary substantially, but each
#  optimized parameter has an approximately normal distribution. will round
#  the outputs to a few digits and look for the modal value returned
modal.bucket.rmse <- mode(round(opt_matrix[,1], digits = 1))
modal.infilt.rmse <- mode(round(opt_matrix[,2], digits = 2))
modal.wilt.rmse <- mode(round(opt_matrix[,3], digits = 2))

# report the retrieved optimized parameters
#  should be ~ 112, 0.23, 0.15
print(paste("Optimized RMSE parameters for bucket 3:", modal.bucket.rmse, modal.infilt.rmse, modal.wilt.rmse))
opt_rmse <- c(modal.bucket.rmse, modal.infilt.rmse, modal.wilt.rmse)

#############################
# task 2, comparing root mean squared and mean of absolute errors

# perform the same calibration process, but using mean absolute error for the loss function
for (i in seq(1,n_searches)){
  
  # set par0 equal to the randomly assigned starting parameters
  par0 <- c(random.bucket[i], random.infilt[i], random.wilt[i])
  
  # run the optimization for these parameters
  opt <- optim(par0, minfunc, loss=loss.abs, model=bucket3, y=runoff, P=P, 
               PET=PET, method='L-BFGS-B', lower=bounds.lower ,upper=bounds.upper)
  
  # assign the output params to the storage matrix
  opt_matrix[i,] <- opt$par
  
  # and report the values for each iteration
  print(opt$par)
}

# find the modal values from using mean absolute error loss function
modal.bucket.abs <- mode(round(opt_matrix[,1], digits = 1))
modal.infilt.abs <- mode(round(opt_matrix[,2], digits = 2))
modal.wilt.abs <- mode(round(opt_matrix[,3], digits = 2))

# report parameters
print(paste("Optimized MAE parameters for bucket 3:", modal.bucket.abs, modal.infilt.abs, modal.wilt.abs))
opt_abs <- c(modal.bucket.abs, modal.infilt.abs, modal.wilt.abs)

# calculate the runoff values for the different optimizations
runoff_rmse <- bucket3(opt_rmse, P, PET)
runoff_abs <- bucket3(opt_abs, P, PET)

# if I got the parameters correct above, the rmse and mean of absolute errors should be 0
print(paste("Root mean squared error:", loss.rmse(runoff, runoff_rmse)))
print(paste("Mean of absolute error :", loss.abs(runoff, runoff_abs)))

# set the y axis min/max
ymin <- min(runoff, runoff_rmse, runoff_abs)
ymax <- max(runoff, runoff_rmse, runoff_abs)
ylim <- c(ymin, ymax)

# set the labels
title <- "Comparing runoff with optimized parameters"
xlab <- "Time step"
ylab <- "Runoff"

# set line and color params
pch_runoff <- 19
lty_rmse <- 2
lty_abs <- 1
color_runoff <- 'green'
color_rmse <- 'black'
color_abs <- 'orange'

# plot the runoff values for each simulation
plot(runoff, ylim = ylim, main = title, xlab = xlab, ylab = ylab, pch = pch_runoff, col = color_runoff)
lines(runoff_abs, lty = lty_abs, col = color_abs)
lines(runoff_rmse, lty = lty_rmse, col = color_rmse)

# add the legend
legend('topright', legend=c('"True" runoff', 'RMSE-loss runoff', 'MAE-loss runoff'),
       col=c(color_runoff, color_rmse, color_abs), lty=c(0, lty_rmse, lty_abs), 
       pch=c(pch_runoff, NA, NA))

#############################
# task 3 - adding noise to estimates of P and PET

# want to add noise from a random distribution
stdev.P <- sd(P)
stdev.PET <- sd(PET)

# we want to iteratively add noise to P and PET from 10-200%, run opt, save the values, then repeat 10 times
noise_fraction <- seq(0.1, 2.0, by = 0.1)
n_iter <- 10

# create the output matrix to store outputs
opt_array_noise <- array(0, dim = c(3, length(noise_fraction), n_iter))

for (i in seq(1,length(noise_fraction))){
	
	# get your normal distribution of noise for each iteration
	noise.P <- rnorm(n_iter, sd = (stdev.P * noise_fraction[i]))
	noise.PET <- rnorm(n_iter, sd = (stdev.PET * noise_fraction[i]))
	
	# loop through the n iterations and try retrieving the parameters for each noise addition
	for (j in seq(1, n_iter)){
	  
	  # create new P and PET vectors with noise added
	  P_temp <- P + noise.P[j]
	  PET_temp <- PET + noise.PET[j]
	  
	  # set negative values to zero
	  P_temp[P_temp < 0] = 0
	  PET_temp[PET_temp < 0] = 0
	  
	  # run optim using the optimized parameters retrieved from the RMSE loss function.
	  #  that should retrieve the correct parameters under noiseless conditions, but here 
	  #  demonstrates how noise affects retrieved parameters
	  opt <- optim(opt_rmse, minfunc, loss=loss.rmse, model=bucket3, y=runoff, P=P_temp, 
	               PET=PET_temp, method='L-BFGS-B', lower=bounds.lower ,upper=bounds.upper)
	  
	  # add the returned params to the matrices
	  opt_array_noise[, i, j] <- opt$par
	}
}

# now that we have generated info for n iterations at varying levels of noise, plot the
#  output parameters (so, for bucket3, each of the 3 parameters)

# first set up plot titles
title.param1 <- "Bucket size param retrieved with added noise"
title.param2 <- "Infiltration param retrieved with added noise"
title.param3 <- "Wilt point param retrieved with added noise"
ylab.param1 <- "Bucket size"
ylab.param2 <- "Infiltration"
ylab.param3 <- "Wilt point"
xlab <- "% added noise"

# lets set up a data frame so we can create a boxplot to represent the retrieved parameters
df.param1 <- data.frame(opt_array_noise[1,,], row.names = paste0(noise_fraction * 100, "%"))
df.param2 <- data.frame(opt_array_noise[2,,], row.names = paste0(noise_fraction * 100, "%"))
df.param3 <- data.frame(opt_array_noise[3,,], row.names = paste0(noise_fraction * 100, "%"))

# create box plots for each parameter
par(mfcol=c(1,3))
boxplot(t(df.param1), xlab = xlab, ylab = ylab.param1, main = title.param1, las = 2, col = 'orange')
boxplot(t(df.param2), xlab = xlab, ylab = ylab.param2, main = title.param2, las = 2, col = 'light green')
boxplot(t(df.param3), xlab = xlab, ylab = ylab.param3, main = title.param3, las = 2, col = 'salmon')

#############################
# task 4 - adding noise to runoff

# add noise based on standard deviation
stdev.runoff <- sd(runoff)

for (i in seq(1,length(noise_fraction))){
  
  # get your normal distribution of noise for each iteration
  noise.runoff <- rnorm(n_iter, sd = (stdev.runoff * noise_fraction[i]))
  
  # loop through the n iterations and try retrieving the parameters for each noise addition
  for (j in seq(1, n_iter)){
    
    # create new P and PET vectors with noise added
    runoff_temp <- runoff + noise.runoff[j]
    
    # set negative values to zero
    runoff_temp[runoff_temp < 0] = 0
    
    # run optim using the optimized parameters retrieved from the RMSE loss function.
    #  that should retrieve the correct parameters under noiseless conditions, but here 
    #  demonstrates how noise affects retrieved parameters. This time, noisy 'truth' data
    opt <- optim(opt_rmse, minfunc, loss=loss.rmse, model=bucket3, y=runoff_temp, P=P, 
                 PET=PET, method='L-BFGS-B', lower=bounds.lower ,upper=bounds.upper)
    
    # add the returned params to the matrices
    opt_array_noise[, i, j] <- opt$par
  }
}

# lets set up a data frame so we can create a boxplot to represent the retrieved parameters
df.param1 <- data.frame(opt_array_noise[1,,], row.names = paste0(noise_fraction * 100, "%"))
df.param2 <- data.frame(opt_array_noise[2,,], row.names = paste0(noise_fraction * 100, "%"))
df.param3 <- data.frame(opt_array_noise[3,,], row.names = paste0(noise_fraction * 100, "%"))

# create box plots for each parameter
par(mfcol=c(1,3))
boxplot(t(df.param1), xlab = xlab, ylab = ylab.param1, main = title.param1, las = 2, col = 'orange')
boxplot(t(df.param2), xlab = xlab, ylab = ylab.param2, main = title.param2, las = 2, col = 'light green')
boxplot(t(df.param3), xlab = xlab, ylab = ylab.param3, main = title.param3, las = 2, col = 'salmon')

# that's all, folks!