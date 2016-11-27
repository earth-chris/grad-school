# ESS-211 Final Project
# Modeling de-militarization in the US
# Christopher Anderson
# 12-2016

#############################
# set up the environment

# set working directory
setwd("~/src/aei-grad-school/courses/ess-211/")

# load ESS-211 functions
source("ESS-211-Functions.R")

# load the raw data to a data frame
enlistment <- read.csv('demilitarization_project/Enlistment.csv', header = TRUE)
demographics <- read.csv('demilitarization_project/Demographics.csv', row.names = 1, header = TRUE)
yearly <- read.csv('demilitarization_project/YearlyData.csv', header = TRUE)

# append no-data to the data frames for years we will model
finalYear = 2100
startYear = max(yearly$Year) + 1

# loop through each year
for (year in seq(startYear, finalYear)){
  
  # create the NA vectors to append
  yearVec <- c(year, rep(NA, ncol(yearly) - 1))
  enlistVec <- c(year, rep(NA, ncol(enlistment) - 1))
  
  # append the new rows
  yearly <- rbind(yearly, yearVec)
  enlistment <- rbind(enlistment, enlistVec)
}

#############################
# create some exploratory plots

# number of unique conflicts per year
cols <- rep("", nrow(yearly))
cols[yearly$PresidentialParty == "Republican"] = "Red"
cols[yearly$PresidentialParty == "Democrat"] = "Blue"
ylab <- "Number of Unique Conflicts"
xlab <- "Year"
title <- "Unique Military Conflicts Per Year"
barplot(yearly$nConflicts, col = add.alpha(cols, 0.6), names.arg = yearly$Year, ylab = ylab, xlab = xlab, main = title)

# average conflicts by party
cols <- c("Red", "Blue")
xlab <- c("Republican", "Democrat")
title <- "Unique Military Conflicts by Presidential Party"
conflictList <- list(yearly$nConflicts[yearly$PresidentialParty == "Republican"], yearly$nConflicts[yearly$PresidentialParty == "Democrat"])
boxplot(conflictList, col = add.alpha(cols, 0.6), names = xlab, main = title, ylab = ylab)
