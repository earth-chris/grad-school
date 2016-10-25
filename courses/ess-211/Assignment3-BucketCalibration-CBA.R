# ESS-211 assignment 3, due 10-26-16
#  cba 10/2016

# change to the working directory and read the rdata
setwd("~/cba/aei-grad-school/courses/ess-211/")

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
random.bucket <- runif(n_searches, min=bounds.lower[1], max=bounds.upper[1])
random.infilt <- runif(n_searches, min=bounds.lower[2], max=bounds.upper[2])
random.wilt <- runif(n_searches, min=bounds.lower[3], max=bounds.upper[3])

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

# the values retrieved from this optimization vary substantially, but each
#  optimized parameter has an approximately normal distribution. will round
#  the outputs to a few digits and look for the modal value returned
modal.bucket <- mode(round(opt_matrix[,1], digits = 1))
modal.infilt <- mode(round(opt_matrix[,2], digits = 2))
modal.wilt <- mode(round(opt_matrix[,3], digits = 2))

# report the retrieved optimized parameters
#  should be ~ 112, 0.23, 0.15
print(paste("Optimized parameters for bucket 3:", modal.bucket, modal.infilt, modal.wilt))
opt_rmse <- c(modal.bucket, modal.infilt, modal.wilt)

#############################
# task 2, comparing root mean squared and mean of absolute errors

# perform the same calibration process, bus using mean absolute error for loss function
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

# if I got the parameters correct above, the rmse and mean of absolute errors should be 0
print(paste("Root mean squared error:", loss.rmse(runoff, c(modal.bucket, modal.infilt, modal.wilt))))
print(paste("Mean of absolute error :", loss.abs(runoff, c(modal.bucket, modal.infilt, modal.wilt))))