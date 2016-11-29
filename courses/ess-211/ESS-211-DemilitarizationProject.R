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

#############################
# set parameter and variable values we'll use later in modeling

# set the years we'll model until
endModelYear <- 2050
startModelYear <- max(yearly$Year) + 1
startYear <- min(yearly$Year)
endYear <- max(yearly$Year)

# set the limits on spending
maxSpending.fed <- 8000 # in billions of USD
maxSpending.mil <- 800 # in billions of USD
maxSpending.vet <- maxSpending.mil # assumes max budget for vets will be at least how much we spend on active military

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
legend("topleft", col = cols, lwd = lwd, legend = c(ylabFed, ylabMil))

#############################
# format the data frames so they contain all historic data
#  and NAs for data we need to model

# loop through each year and append no-data to the data frames for years we will model
for (year in seq(startModelYear, endModelYear)){
  
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

# set years we have historic data from
fitYears <- which(yearly$FedSpending.BUSD > 0)

# first, we'll look at what exponential federal spending looks like when fit to all data.
#  the exponential growth model is our null model because exponential growth underlies
#  the policies set by the federal reserve regarding growth of the money supply.
#  growth in money supply does not necessarily mean growth in fed spending, but it historically has.
fedSpending.exponential <- growth.exponential(yearly[fitYears,], degree = 2)
fedSpending.predicted.exp <- predict(fedSpending.exponential, yearly)

# calculate rmse for exponential growth
rmse.fedSpending.exp <- sqrt(mean(residuals(fedSpending.exponential)^2))

# plot the output
ylab <- "Total Federal Spending (B $USD)"
xlab <- "Year"
title <- "Actual and Modeled (exponential) Federal Spending"
legend <- c("Actual Federal Spending", "Modeled Federal Spending", paste("RMSE:", format(rmse.fedSpending.exp, nsmall = 2)))
ylim <- c(min(yearly$FedSpending.BUSD, fedSpending.predicted.exp, na.rm = TRUE), 
          max(yearly$FedSpending.BUSD, fedSpending.predicted.exp, na.rm = TRUE))
lwd <- c(3, 3, NA)
plot(yearly$Year, yearly$FedSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, fedSpending.predicted.exp, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

# next, we'll look at what logistic growth in federal spending looks like.
#  we're going to set the ceiling of federal spending to that
#  specified in the parameters above, suggesting that at the final year we model,
#  we will be spending the max amount. we'll change this as we do sensitivity analysis, etc.
fedSpending.logistic <- growth.logistic(yearly, maxSpending.fed, nrow(yearly))
fedSpending.predicted.log <- predict(fedSpending.logistic, yearly)

# calculate rmse for logistic growth
rmse.fedSpending.log <- sqrt(mean(residuals(fedSpending.logistic)^2))

# plot the output
ylab <- "Total Federal Spending (B $USD)"
xlab <- "Year"
title <- "Actual and Modeled (logistic) Federal Spending"
legend <- c("Actual Federal Spending", "Modeled Federal Spending", paste("RMSE:", format(rmse.fedSpending.log, nsmall = 2)))
ylim <- c(min(yearly$FedSpending.BUSD, fedSpending.predicted.log, na.rm = TRUE), 
          max(yearly$FedSpending.BUSD, fedSpending.predicted.log, na.rm = TRUE))
lwd <- c(3, 3, NA)
plot(yearly$Year, yearly$FedSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, fedSpending.predicted.log, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

#############################
# compare exponential vs logistic growth for military spending

# first, exponential growth
milSpending.exponential <- milGrowth.exponential(yearly[fitYears,], degree = 2)
milSpending.predicted.exp <- predict(milSpending.exponential, yearly)

# next, logistic growth
milSpending.logistic <- milGrowth.logistic(yearly, maxSpending.mil, nrow(yearly))
milSpending.predicted.log <- predict(milSpending.logistic, yearly)

# get rmse for each
rmse.milSpending.exp <- sqrt(mean(residuals(milSpending.exponential)^2))
rmse.milSpending.log <- sqrt(mean(residuals(milSpending.logistic)^2))

# then plot the outputs
ylab <- "Total Military Spending (B $USD)"
xlab <- "Year"
title <- "Actual and Modeled (exponential) Military Spending"
legend <- c("Actual Military Spending", "Modeled Military Spending", paste("RMSE:", format(rmse.milSpending.exp, nsmall = 2)))
ylim <- c(min(yearly$MilitarySpending.BUSD, milSpending.predicted.exp, na.rm = TRUE), 
          max(yearly$MilitarySpending.BUSD, milSpending.predicted.exp, na.rm = TRUE))
lwd <- c(3, 3, NA)
plot(yearly$Year, yearly$MilitarySpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, milSpending.predicted.exp, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

legend <- c("Actual Military Spending", "Modeled Military Spending", paste("RMSE:", format(rmse.milSpending.log, nsmall = 2)))
ylim <- c(min(yearly$MilitarySpending.BUSD, milSpending.predicted.log, na.rm = TRUE), 
          max(yearly$MilitarySpending.BUSD, milSpending.predicted.log, na.rm = TRUE))
lwd <- c(3, 3, NA)
plot(yearly$Year, yearly$MilitarySpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, milSpending.predicted.log, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

#############################
# compare exponential vs logistic growth for veteran spending

# first, exponential growth
vetSpending.exponential <- vetGrowth.exponential(yearly[fitYears,], degree = 2)
vetSpending.predicted.exp <- predict(vetSpending.exponential, yearly)

# next, logistic growth
vetSpending.logistic <- vetGrowth.logistic(yearly, maxSpending.vet, nrow(yearly))
vetSpending.predicted.log <- predict(vetSpending.logistic, yearly)

# get rmse for each
rmse.vetSpending.exp <- sqrt(mean(residuals(vetSpending.exponential)^2))
rmse.vetSpending.log <- sqrt(mean(residuals(vetSpending.logistic)^2))

# then plot the outputs
ylab <- "Total Veteran Spending (B $USD)"
xlab <- "Year"
title <- "Actual and Modeled (exponential) Veteran Spending"
legend <- c("Actual Veteran Spending", "Modeled Veteran Spending", paste("RMSE:", format(rmse.vetSpending.exp, nsmall = 2)))
ylim <- c(min(yearly$VeteranSpending.BUSD, vetSpending.predicted.exp, na.rm = TRUE), 
          max(yearly$VeteranSpending.BUSD, vetSpending.predicted.exp, na.rm = TRUE))
lwd <- c(3, 3, NA)
plot(yearly$Year, yearly$VeteranSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, vetSpending.predicted.exp, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

legend <- c("Actual Military Spending", "Modeled Military Spending", paste("RMSE:", format(rmse.vetSpending.log, nsmall = 2)))
ylim <- c(min(yearly$VeteranSpending.BUSD, vetSpending.predicted.log, na.rm = TRUE), 
          max(yearly$VeteranSpending.BUSD, vetSpending.predicted.log, na.rm = TRUE))
lwd <- c(3, 3, NA)
plot(yearly$Year, yearly$VeteranSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, vetSpending.predicted.log, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

#############################
# modeling future unemployment rates

# unemployment is difficult to model, and there are many high level economics papers
#  devoted to doing so. That is beyond the scope of this project, for now.
#  I'm planning to model future unemployment using three basic frameworks:
#  1) using the mean of historical unemployment and adding noise based on a normal distribution,
#  2) using a linear trendline of historical unemployment and adding noise based on a normal distribution, and
#  3) fitting a sin curve to historical unemployment and adding noise based on a normal distribution.

# get basic unemployment data
unemployment.mean <- mean(yearly$UnemploymentRate[unemployment.years])
unemployment.sd <- sd(yearly$UnemploymentRate[unemployment.years])

# run framework 1
framework1 <- unemployment.framework1(yearly, unemployment.sd)

# run framework 2
framework2 <- unemployment.framework2(yearly, unemployment.sd)