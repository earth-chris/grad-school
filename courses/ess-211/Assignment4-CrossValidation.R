# ESS-211 assignment 4, due 11-4-16
#  cba 11/2016

# set current working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# source the written functions
source("ESS-211-Functions.R")
library(gplots)

# read the input data for assignment 4
data <- read.csv(file = "course_material/exercise3.data.csv", header = TRUE)

#############################
# task 1 - cross validation error on k splits

# write a function to calculate rmse for k splits
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
# task 2 - varying k

# set up a vector for the different values of k we'll test, from 2 (50% hold-out) to
#  nrows/2 (two samples per hold-out)
K <- seq(2, nrow(data)/2)

# then append a final value to the vector, the number of rows, to use for n-fold x-val
K <- append(K, nrow(data))

# create a matrix to store outputs of mean/stdev for test/training results 
testMatrix <- matrix(ncol = 2, nrow = length(K))
trainMatrix <- matrix(ncol = 2, nrow = length(K))

# loop through each of the k-fold methods and calculate mean/stdev for test/training error
for (i in 1:length(K)){
  # run the k-fold xval procedure on the i'th value of K
  ifoldMatrix <- cv.err(data, K[i])
  
  # assign mean and stdev of each xval procedure to storage matrices
  trainMatrix[i,1] <- mean(ifoldMatrix[,1])
  trainMatrix[i,2] <- sd(ifoldMatrix[,1])
  testMatrix[i,1] <- mean(ifoldMatrix[,2])
  testMatrix[i,2] <- sd(ifoldMatrix[,2])
}

# replace the stdev values with standard error values
testMatrix[,2] <- testMatrix[,2] / sqrt(K)
trainMatrix[,2] <- trainMatrix[,2] / sqrt(K)

# find the min and max ranges for all data to plot on the same axes
ymin <- min(c(trainMatrix[,1]-trainMatrix[,2], testMatrix[,1]-testMatrix[,2]))
ymax <- max(c(trainMatrix[,1]+trainMatrix[,2], testMatrix[,1]+testMatrix[,2]))
buffer <- (ymax - ymin) * 0.05
ylim <- c(ymin-buffer, ymax+buffer)

# plot the means and standard errors of the training and test error with varying K
#  first, set up titles, etc.
mainTitle <- "Training and test error as a function of the number of folds"
trainTitle <- "Training error"
testTitle <- "Testing error"
ylab <- "RMSE"
xlab <- "N folds"

# set up x axis as the number of training samples per fold
xvals <- nrow(data)/K

# set up plot structure
par(mfcol=c(1,2), oma = c(0,0,2,0))

# set up the training data plot
plotCI(K, trainMatrix[,1], type = 'p', uiw = trainMatrix[,2], liw = trainMatrix[,2], ylim = ylim,
       pch = 19, ylab = ylab, xlab = xlab, main = trainTitle, col = 'black', barcol = 'orange')

# add a legend
legend('topright', legend=c('Mean', 'Standard error'),
       col=c('black', 'orange'), lty=c(NA, 1), pch=c(19,NA))

# set up the test data plot
plotCI(K, testMatrix[,1], type = 'p', uiw = testMatrix[,2], liw = testMatrix[,2], ylim = ylim,
       pch = 19, ylab = ylab, xlab = xlab, main = testTitle, col = 'black', barcol = 'aquamarine')

# add a legend
legend('topright', legend=c('Mean', 'Standard error'),
       col=c('black', 'aquamarine'), lty=c(NA, 1), pch=c(19,NA))

# and add a title
title(mainTitle, outer = TRUE)

#############################
# task 3 - adding noise

# add three columnds of noise to the data
data[,5:7] <- rnorm(3 * nrow(data), sd = 4)

# create matrices to store outputs
noiseTrainMatrix <- matrix(ncol = 2, nrow = ncol(data)-1)
noiseTestMatrix <- matrix(ncol = 2, nrow = ncol(data)-1)

# iteratively add columns to the linear regression to see how cross validation errors change
#  with additional variables and noise
for (i in 2:ncol(data)){
  
  # we're going to use the leave-one-out method to assess the contribution of additional variables/noise
  columnMatrix <- cv.err(data[,1:i], nrow(data))
  
  # assign mean/stdev to matrix
  noiseTrainMatrix[i-1,1] <- mean(columnMatrix[,1])
  noiseTrainMatrix[i-1,2] <- sd(columnMatrix[,1])
  noiseTestMatrix[i-1,1] <- mean(columnMatrix[,2])
  noiseTestMatrix[i-1,2] <- sd(columnMatrix[,2])
}

# replace the stdev values with standard error values
noiseTestMatrix[,2] <- noiseTestMatrix[,2] / sqrt(nrow(data))
noiseTrainMatrix[,2] <- noiseTrainMatrix[,2] / sqrt(nrow(data))

# plot the means and standard errors of the training and test error with varying number of terms to fit
#  first, set up titles, etc.
mainTitle <- "Training and test error as a function of the number of terms (SD = 4)"
trainTitle <- "Training error"
testTitle <- "Testing error"
ylab <- "RMSE"
xlab <- "N terms fit"

# find the min and max ranges for all data to plot on the same axes
ymin <- min(c(noiseTrainMatrix[,1]-noiseTrainMatrix[,2], noiseTestMatrix[,1]-noiseTestMatrix[,2]))
ymax <- max(c(noiseTrainMatrix[,1]+noiseTrainMatrix[,2], noiseTestMatrix[,1]+noiseTestMatrix[,2]))
buffer <- (ymax - ymin) * 0.05
ylim <- c(ymin-buffer, ymax+buffer)

# set up plot structure
par(mfcol=c(1,2), oma = c(0,0,2,0))

# set up the training data plot
plotCI(1:(ncol(data)-1), noiseTrainMatrix[,1], type = 'p', uiw = noiseTrainMatrix[,2], liw = noiseTrainMatrix[,2], ylim = ylim,
       pch = 19, ylab = ylab, xlab = xlab, main = trainTitle, col = 'black', barcol = 'orange')

# add a legend
legend('topright', legend=c('Mean', 'Standard error'),
       col=c('black', 'orange'), lty=c(NA, 1), pch=c(19,NA))

# set up the test data plot
plotCI(1:(ncol(data)-1), noiseTestMatrix[,1], type = 'p', uiw = noiseTestMatrix[,2], liw = noiseTestMatrix[,2], ylim = ylim,
       pch = 19, ylab = ylab, xlab = xlab, main = testTitle, col = 'black', barcol = 'aquamarine')

# add a legend
legend('topright', legend=c('Mean', 'Standard error'),
       col=c('black', 'aquamarine'), lty=c(NA, 1), pch=c(19,NA))

# and add a title
title(mainTitle, outer = TRUE)