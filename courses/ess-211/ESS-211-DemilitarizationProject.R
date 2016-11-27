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
source("ESS-211-DemilitarizationFunctions.R")

# load the raw data to a data frame
enlistment <- read.csv('demilitarization_project/Enlistment.csv', header = TRUE)
demographics <- read.csv('demilitarization_project/Demographics.csv', row.names = 1, header = TRUE)
yearly <- read.csv('demilitarization_project/YearlyData.csv', header = TRUE)

# set parameter and variable values we'll use later in modeling
finalYear = 2050
startYear = max(yearly$Year) + 1

# set plotting parameters
colDem <- add.alpha("Blue", 0.6)
colRep <- add.alpha("Red", 0.6)
colFed <- add.alpha("Purple", 0.8)
colMil <- add.alpha("Orange", 0.8)
colPred <- add.alpha("Dark Green", 0.8)
colReal <- add.alpha("Black", 0.8)

#############################
# create some exploratory plots

# number of unique conflicts per year
cols <- rep("", nrow(yearly))
cols[yearly$PresidentialParty == "Republican"] = colRep
cols[yearly$PresidentialParty == "Democrat"] = colDem
ylab <- "Number of Unique Conflicts"
xlab <- "Year"
title <- "Unique Military Conflicts per Year"
barplot(yearly$nConflicts, col = add.alpha(cols, 0.6), names.arg = yearly$Year, ylab = ylab, xlab = xlab, main = title)

# average conflicts by party
cols <- c(colRep, colDem)
xlab <- c("Republican", "Democrat")
ylab <- "Number of Unique Conflicts per Year"
title <- "Unique Military Conflicts per Year by Presidential Party"
conflictList <- list(yearly$nConflicts[yearly$PresidentialParty == "Republican"], yearly$nConflicts[yearly$PresidentialParty == "Democrat"])
boxplot(conflictList, col = add.alpha(cols, 0.6), names = xlab, main = title, ylab = ylab)

# comparing enlistment and total US population
cols <- c(colFed, colMil)
xlab <- "Year"
ylabPop <- "US Population (millions of people)"
ylabEnlist <- "Military Enlistment (millions of people)"
title <- "US and Military Populations"
lwd <- c(3, 3)
par(mar = c(5,5,2,5))
with(yearly, plot(Year, TotalPopulation, type = 'l', col = cols[1], xlab = xlab, ylab = ylabPop, lwd = lwd[1], main = title))
par(new = TRUE)
with(yearly, plot(Year, (TotalEnlistment / 1e6), type = 'l', col = cols[2], lwd = lwd[2], xlab = NA, ylab = NA, axes = FALSE))
axis(side = 4)
mtext(ylabEnlist, side = 4, line = 3)

# comparing total federal spending and military spending
cols <- c(colFed, colMil)
xlab <- "Year"
ylabFed <- "Total Federal Spending (B $USD)"
ylabMil <- "Total Military Spending (B $USD)"
title <- "Federal and Military Spending"
lwd <- c(3, 3)
par(mar = c(5,5,2,5))
with(yearly, plot(Year, FedSpending.BUSD, type = 'l', col = cols[1], xlab = xlab, ylab = ylabFed, lwd = lwd[1], main = title))
par(new = TRUE)
with(yearly, plot(Year, MilitarySpending.BUSD, type = 'l', col = cols[2], xlab = NA, ylab = NA, lwd = lwd[2], main = title, axes = FALSE))
axis(side = 4)
mtext(ylabMil, side = 4, line = 3)

#############################
# format the data frames so they contain all historic data
#  and NAs for data we need to model

# loop through each year and append no-data to the data frames for years we will model
for (year in seq(startYear, finalYear)){
  
  # create the NA vectors to append
  yearVec <- c(year, rep(NA, ncol(yearly) - 1))
  enlistVec <- c(year, rep(NA, ncol(enlistment) - 1))
  
  # append the new rows
  yearly <- rbind(yearly, yearVec)
  enlistment <- rbind(enlistment, enlistVec)
}

#############################
# fill the missing data for 'business-as-usual' scenarios, wherein we assume what 
#  various predictor values will be based on historic data

# first, we'll look at what exponential federal spending looks like when fit to all data.
#  the exponential growth model is our null model because exponential growth underlies
#  the policies set by the federal reserve regarding growth of the money supply.
#  growth in money supply does not necessarily mean growth in fed spending, but it historically has.
fitYears <- which(yearly$FedSpending.BUSD > 0)
fedSpending.exponential <- growth.exponential(x = yearly$Year[fitYears], y = yearly$FedSpending.BUSD[fitYears], degree = 2)
fedSpending.predicted.exp <- predict(fedSpending.exponential, x = yearly$Year)

# plot the output
ylab <- "Total Federal Spending (B $USD)"
xlab <- "Year"
title <- "Actual and Modeled (exponential) Federal Spending"
legend <- c("Actual Federal Spending", "Modeled Federal Spending")
lwd <- c(3, 3)
plot(yearly$Year, yearly$FedSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1])
par(new = TRUE)
plot(yearly$Year, fedSpending.predicted.exp, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred), lwd = lwd)
