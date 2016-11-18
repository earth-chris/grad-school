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
enlistmint <- read.csv('demilitarization_project/Enlistment.csv', header = TRUE)
demographics <- read.csv('demilitarization_project/Demographics.csv', row.names = 1, header = TRUE)
yearly <- read.csv('demilitarization_project/YearlyData.csv', header = TRUE)

# exploratory plots
cols <- rep("", nrow(yearly))
cols[yearly$PresidentialParty == "Republican"] = "Red"
cols[yearly$PresidentialParty == "Democrat"] = "Blue"
ylab <- "Number of Unique Conflicts"
xlab <- "Year"
title <- "Unique Military Conflicts Per Year"

barplot(yearly$nConflicts, col = cols, names.arg = yearly$Year, ylab = ylab, xlab = xlab, main = title)
