#ESS-211 Final Project
# Modeling de-militarization in the US
# Christopher Anderson
# 12-2016

#############################
# set up the environment

# set working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# load ESS-211 functions
source("ESS-211-Functions.R")
source("ESS-211-DemilitarizationFunctions.R")

# load the raw data
enlistment <- read.csv('demilitarization_project/Enlistment.csv', header = TRUE)
demographics <- read.csv('demilitarization_project/Demographics.csv', row.names = 1, header = TRUE)
yearly <- read.csv('demilitarization_project/YearlyData.csv', header = TRUE)

#############################
# set parameters to run with

# set the years to model
year.start <- 1954
year.end <- 2016

# set the number of years for demilitarization
year.demil <- 10

model.demil <- function(population, enlistment, fedSpending, milSpending, vetSpending, unemployment, laborParticipation, demographics){
  
}