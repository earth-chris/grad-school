# ESS-211 assignment 3, due 10-26-16
#  cba 10/2016

# change to the working directory and read the rdata
setwd("~/Downloads/source/aei-grad-school/courses/ess-211/")

# this functions file contains bucket1, bucket3, and the loss functions
source("ESS-211-Functions.R")

# load the provided P, PET, and runoff data
load('course_material/runoff.Rdata')

#############################
# task 1 - finding parameters for bucket3 using runoff data

# first, set the upper and lower bounds of the parameters for infilt and wilt
bounds.upper <- c()
bounds.lower <- c()

# run the optimization
opt = optim(par0, minfunc, loss=loss.rmse, model=bucket3, y=runoff, P=P, 
		PET=PET, method='L-BFGS-B', lower=bounds.lower ,upper=bounds.upper)

#############################

# run it with optim function, which calculates the loss function over a series of parameters
bucket.optim <- optim(par0, minfunc, loss=loss.rmse, model=bucket3, y=classbucket3, P=P, PET=PET, method = 'L-BFGS-B')