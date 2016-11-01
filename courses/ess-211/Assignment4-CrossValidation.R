# ESS-211 assignment 4, due 11-4-16
#  cba 11/2016

# set current working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# source the written functions
source("ESS-211-Functions.R")

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
# task 3 - varying k

# set up a vector for the different values of k we'll test, from 2 (50% hold-out) to
#  nrows/2 (two samples per hold-out)
K <- seq(2, nrow(data)/2)

# then append a final value to the vector, the number of rows, to use for n-fold x-val
K <- append(K, nrow(data))

# create a matrix to store outputs of mean/stdev for test/training results 
testMatrix <- matrix(ncol = 2, nrow = length(K))
trainMatrix <- matrix(ncol = 2, nrow = length(K))

# loop through each of the k-fold methods and calculate mean/stdev test/training error
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

# plot the means and confidence intervals of the training and test error with varying K
#  first, set up titles, etc.
trainTitle <- "Training error"
testTitle <- "Testing error"
ylab <- "RMSE"
xlab <- "N folds"

# set up x axis as the number of training samples per fold
xvals <- nrow(data)/K

par(mfcol=c(1,2))

# set up the training data plot
plotCI(K, trainMatrix[,1], type = 'p', uiw = trainMatrix[,2], liw = trainMatrix[,2], 
       pch = 19, ylab = ylab, xlab = xlab, main = trainTitle, col = 'black', barcol = 'orange')

# set up the test data plot
plotCI(K, testMatrix[,1], type = 'p', uiw = testMatrix[,2], liw = testMatrix[,2], 
       pch = 19, ylab = ylab, xlab = xlab, main = trainTitle, col = 'black', barcol = 'aquamarine')
