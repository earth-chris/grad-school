# ESS-211 Final Project
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
maxSpending.vet <- maxSpending.mil / 2 # assumes max budget for vets will be half we spend on active military

# set plotting parameters
colDem <- add.alpha("Blue", 0.6)
colRep <- add.alpha("Red", 0.6)
colFed <- add.alpha("Purple", 0.8)
colMil <- add.alpha("Orange", 0.8)
colPred <- add.alpha("Dark Green", 0.8)
colReal <- add.alpha("Black", 0.8)

#############################
# adjust all USD measurements for inflation
yearly$GDP.BUSD <- inflation.adjust(yearly$GDP.BUSD, yearly$InflationRate)
yearly$FedSpending.BUSD <- inflation.adjust(yearly$FedSpending.BUSD, yearly$InflationRate)
yearly$MilitarySpending.BUSD <- inflation.adjust(yearly$MilitarySpending.BUSD, yearly$InflationRate)
yearly$VeteranSpending.BUSD <- inflation.adjust(yearly$VeteranSpending.BUSD, yearly$InflationRate)

#############################
# create some exploratory plots

# number of unique conflicts per year
cols <- rep(NA, nrow(yearly))
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

# calibrate the growth rate
fedSpending.growth.rate <- calibrate.growth(yearly$Year[fitYears], yearly$FedSpending.BUSD[fitYears])

# apply it to the data
fedSpending.predicted.exp <- growth.exponential.apply(yearly$FedSpending.BUSD[1], fedSpending.growth.rate, nYears = nrow(yearly))

# calculate rmse for exponential growth
rmse.fedSpending.exp <- sqrt(mean((yearly$FedSpending.BUSD[fitYears] - fedSpending.predicted.exp[fitYears])^2))

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

# next, we'll look at a 2-term, polynomial growth function
#  instead of using the PERT-based exponential growth model
fedSpending.polynomial <- growth.polynomial(yearly$Year[fitYears], yearly$FedSpending.BUSD[fitYears], degree = 2)
fedSpending.predicted.poly <- predict(fedSpending.polynomial, data.frame(x = yearly$Year))

# calculate rmse for polynomial growth
rmse.fedSpending.poly <- sqrt(mean(residuals(fedSpending.polynomial)^2))

# plot the output
title <- "Actual and Modeled (polynomial) Federal Spending"
legend <- c("Actual Federal Spending", "Modeled Federal Spending", paste("RMSE:", format(rmse.fedSpending.poly, nsmall = 2)))
ylim <- c(min(yearly$FedSpending.BUSD, fedSpending.predicted.poly, na.rm = TRUE), 
          max(yearly$FedSpending.BUSD, fedSpending.predicted.poly, na.rm = TRUE))
plot(yearly$Year, yearly$FedSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, fedSpending.predicted.poly, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

# next, we'll look at what logistic growth in federal spending looks like.
#  we're going to set the ceiling of federal spending to that
#  specified in the parameters above, suggesting that at the final year we model,
#  we will be spending the max amount. we'll change this as we do sensitivity analysis, etc.
fedSpending.logistic <- growth.logistic(yearly$Year, yearly$FedSpending.BUSD, maxSpending.fed, nrow(yearly))
fedSpending.predicted.log <- predict(fedSpending.logistic, data.frame(x = yearly$Year))

# calculate rmse for logistic growth
rmse.fedSpending.log <- sqrt(mean(residuals(fedSpending.logistic)^2))

# plot the output
title <- "Actual and Modeled (logistic) Federal Spending"
legend <- c("Actual Federal Spending", "Modeled Federal Spending", paste("RMSE:", format(rmse.fedSpending.log, nsmall = 2)))
ylim <- c(min(yearly$FedSpending.BUSD, fedSpending.predicted.log, na.rm = TRUE), 
          max(yearly$FedSpending.BUSD, fedSpending.predicted.log, na.rm = TRUE))
plot(yearly$Year, yearly$FedSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, fedSpending.predicted.log, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

# and try a third model that uses GDP for predicting federal spending
fedSpending.gdp <- growth.polynomial(yearly$GDP.BUSD[fitYears], yearly$FedSpending.BUSD[fitYears], degree = 1)
fedSpending.predicted.gdp <- predict(fedSpending.gdp, data.frame(x = yearly$GDP.BUSD))

# calculate rmse for gdp-based growth
rmse.fedSpending.gdp <- sqrt(mean(residuals(fedSpending.gdp)^2))

# plot the output
title <- "Actual and Modeled (GDP-Based polynomial) Federal Spending"
legend <- c("Actual Federal Spending", "Modeled Federal Spending", paste("RMSE:", format(rmse.fedSpending.log, nsmall = 2)))
ylim <- c(min(yearly$FedSpending.BUSD, fedSpending.predicted.gdp, na.rm = TRUE), 
          max(yearly$FedSpending.BUSD, fedSpending.predicted.gdp, na.rm = TRUE))
plot(yearly$Year, yearly$FedSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, fedSpending.predicted.gdp, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

#############################
# compare exponential vs logistic growth for military spending

# first, exponential growth
milSpending.growth.rate <- calibrate.growth(yearly$Year[fitYears], yearly$MilitarySpending.BUSD[fitYears])
milSpending.predicted.exp <- growth.exponential.apply(yearly$MilitarySpending.BUSD[1], milSpending.growth.rate, nYears = nrow(yearly))

# next, polynomial growth
milSpending.polynomial <- growth.polynomial(yearly$Year[fitYears], yearly$MilitarySpending.BUSD[fitYears], degree = 2)
milSpending.predicted.poly <- predict(milSpending.polynomial, data.frame(x = yearly$Year))

# next, logistic growth
#  just kidding, when adjusted for inflation this fails
#milSpending.logistic <- growth.logistic(yearly$Year, yearly$MilitarySpending.BUSD, maxSpending.mil, nrow(yearly))
#milSpending.predicted.log <- predict(milSpending.logistic, data.frame(x = yearly$Year))

# next, gdp growth
milSpending.gdp <- growth.polynomial(yearly$GDP.BUSD[fitYears], yearly$MilitarySpending.BUSD[fitYears], degree = 2)
milSpending.predicted.gdp <- predict(milSpending.gdp, data.frame(x = yearly$GDP.BUSD))

# get rmse for each
rmse.milSpending.exp <- sqrt(mean((yearly$MilitarySpending.BUSD[fitYears] - milSpending.predicted.exp[fitYears])^2))
rmse.milSpending.poly <- sqrt(mean(residuals(milSpending.polynomial)^2))
#rmse.milSpending.log <- sqrt(mean(residuals(milSpending.logistic)^2))
rmse.milSpending.gdp <- sqrt(mean(residuals(milSpending.gdp)^2))

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

title <- "Actual and Modeled (polynomial) Military Spending"
legend <- c("Actual Military Spending", "Modeled Military Spending", paste("RMSE:", format(rmse.milSpending.poly, nsmall = 2)))
ylim <- c(min(yearly$MilitarySpending.BUSD, milSpending.predicted.poly, na.rm = TRUE), 
          max(yearly$MilitarySpending.BUSD, milSpending.predicted.poly, na.rm = TRUE))
plot(yearly$Year, yearly$MilitarySpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, milSpending.predicted.poly, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

#title <- "Actual and Modeled (logistic) Military Spending"
#legend <- c("Actual Military Spending", "Modeled Military Spending", paste("RMSE:", format(rmse.milSpending.log, nsmall = 2)))
#ylim <- c(min(yearly$MilitarySpending.BUSD, milSpending.predicted.log, na.rm = TRUE), 
#          max(yearly$MilitarySpending.BUSD, milSpending.predicted.log, na.rm = TRUE))
#plot(yearly$Year, yearly$MilitarySpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
#par(new = TRUE)
#plot(yearly$Year, milSpending.predicted.log, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
#legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

title <- "Actual and Modeled (GDP-based polynomial) Military Spending"
legend <- c("Actual Military Spending", "Modeled Military Spending", paste("RMSE:", format(rmse.milSpending.gdp, nsmall = 2)))
ylim <- c(min(yearly$MilitarySpending.BUSD, milSpending.predicted.gdp, na.rm = TRUE), 
          max(yearly$MilitarySpending.BUSD, milSpending.predicted.gdp, na.rm = TRUE))
plot(yearly$Year, yearly$MilitarySpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, milSpending.predicted.gdp, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

#############################
# compare exponential vs logistic growth for veteran spending

# first, exponential growth
vetSpending.growth.rate <- calibrate.growth(yearly$Year[fitYears], yearly$VeteranSpending.BUSD[fitYears])
vetSpending.predicted.exp <- growth.exponential.apply(yearly$VeteranSpending.BUSD[1], vetSpending.growth.rate, nYears = nrow(yearly))

# next, polynomial growth
vetSpending.polynomial <- growth.polynomial(yearly$Year[fitYears], yearly$VeteranSpending.BUSD[fitYears], degree = 2)
vetSpending.predicted.poly <- predict(vetSpending.polynomial, data.frame(x = yearly$Year))

# next, logistic growth
vetSpending.logistic <- growth.logistic(yearly$Year, yearly$VeteranSpending.BUSD, maxSpending.vet, nrow(yearly))
vetSpending.predicted.log <- predict(vetSpending.logistic, data.frame(x = yearly$Year))

# get rmse for each
rmse.vetSpending.exp <- sqrt(mean((yearly$VeteranSpending.BUSD[fitYears] - vetSpending.predicted.exp[fitYears])^2))
rmse.vetSpending.poly <- sqrt(mean(residuals(vetSpending.polynomial)^2))
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

ylab <- "Total Veteran Spending (B $USD)"
xlab <- "Year"
title <- "Actual and Modeled (polynomial) Veteran Spending"
legend <- c("Actual Veteran Spending", "Modeled Veteran Spending", paste("RMSE:", format(rmse.vetSpending.poly, nsmall = 2)))
ylim <- c(min(yearly$VeteranSpending.BUSD, vetSpending.predicted.poly, na.rm = TRUE), 
          max(yearly$VeteranSpending.BUSD, vetSpending.predicted.poly, na.rm = TRUE))
lwd <- c(3, 3, NA)
plot(yearly$Year, yearly$VeteranSpending.BUSD, type = 'l', col = colReal, xlab = xlab, ylab = ylab, lwd = lwd[1], ylim = ylim, main = title)
par(new = TRUE)
plot(yearly$Year, vetSpending.predicted.poly, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim, axes = FALSE)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

title <- "Actual and Modeled (logistic) Veteran Spending"
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
unemployment.years <- which(yearly$UnemploymentRate < 1.)
unemployment.yearsToModel <- is.na(yearly$UnemploymentRate)
unemployment.mean <- mean(yearly$UnemploymentRate[unemployment.years])
unemployment.sd <- sd(yearly$UnemploymentRate[unemployment.years])

# run framework 1
framework1 <- unemployment.framework1(yearly, unemployment.sd)

# run framework 2
framework2 <- unemployment.framework2(yearly, unemployment.sd)

# run framework 3
framework3 <- unemployment.framework3(yearly, unemployment.sd)

# calculate rmses
rmse.framework1 <- sqrt(mean((yearly$UnemploymentRate[unemployment.years] - framework1[unemployment.years])^2))
rmse.framework2 <- sqrt(mean((yearly$UnemploymentRate[unemployment.years] - framework2[unemployment.years])^2))
rmse.framework3 <- sqrt(mean((yearly$UnemploymentRate[unemployment.years] - framework3[unemployment.years])^2))

# plot em
ylim <- c(min(yearly$UnemploymentRate, framework1, na.rm = TRUE), 
          max(yearly$UnemploymentRate, framework1, na.rm = TRUE))
xlab = "Year"
ylab = "Unemployment Rate"
title = "Yearly US Unemployment"
legend = c("Real Unemployment Rate", "Mean + Noise Unemployment", paste("RMSE:", format(rmse.framework1, nsmall = 2)))
plot(yearly$Year, yearly$UnemploymentRate, type = 'l', xlab = xlab, ylab = xlab, main = title, ylim = ylim, lwd=lwd[2])
par(new=TRUE)
plot(yearly$Year, framework1, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

legend = c("Real Unemployment Rate", "Trendline + Noise Unemployment", paste("RMSE:", format(rmse.framework2, nsmall = 2)))
plot(yearly$Year, yearly$UnemploymentRate, type = 'l', xlab = xlab, ylab = xlab, main = title, ylim = ylim, lwd=lwd[2])
par(new=TRUE)
plot(yearly$Year, framework2, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

legend = c("Real Unemployment Rate", "Sin Curve + Noise Unemployment", paste("RMSE:", format(rmse.framework3, nsmall = 2)))
plot(yearly$Year, yearly$UnemploymentRate, type = 'l', xlab = xlab, ylab = xlab, main = title, ylim = ylim, lwd=lwd[2])
par(new=TRUE)
plot(yearly$Year, framework3, type = 'l', col = colPred, xlab = NA, ylab = NA, lwd = lwd[2], ylim = ylim)
legend("topleft", legend = legend, col = c(colReal, colPred, NA), lwd = lwd)

#############################
# modeling veteran populations

# we don't have veteran population data for our historical records, but we want it for the future and
#  historical estimates to find the number of enlistments prevented, so we can project the difference 
#  in what would have been spent on veteran's costs.
#  we have other data starting 1954, but we want to know the veteran population then.
#  world wars I & II were the last wars where the majority of enlistees would have gone from the enlistmnet
#  to the veterans communities, but we don't have any death rate esimates for WWI. We do, however, 
#  for WWII. so we'll model the number of veterans based on their veteran populations

# total service members serving world war two minus battle and non-theater deaths
ww2ServiceMembers <- 16112566 - 291557 - 113842

# and the living veterans
ww2LivingVeterans <- 1711000

# set start and end years for calculation
ww2t0 <- 1945
ww2t1 <- 2016

# find the growth rate
deathRate.ww2 <- growth.calc.r(ww2ServiceMembers, ww2LivingVeterans, (ww2t1 - ww2t0))

# then estimate the number of veterans at 1954
vetPopulationT0 <- growth.exponential(ww2ServiceMembers, deathRate.ww2, (yearly$Year[1] - ww2t0))

# the above was done to get a general estimate for the total number of WWII veterans that would have been 
#  alive at our t0 (1954). I then took a range around that number, and growth rates, and decided to perform
#  a model calibration to get historic veteran population data

# the model takes 3 args we'll vary: 1) the veteran death rate, 2) the average delistment rate, and 
#  3) the veteran population at t0. 
vetMins <- c(-0.2, 0.05, vetPopulationT0 * 0.8)
vetMaxs <- c(-0.005, 0.15, vetPopulationT0 * 1.2)

# get indices for years beyond start date
enlistmentIndices = which(yearly$Year >= 2013)

# we have to set unknown enlistment numbers. based on the data, I will forecast future enlistment to be the mean
#  of the bush/obama era.
yearly$TotalEnlistment[enlistmentIndices] = mean(yearly$TotalEnlistment[(which(yearly$Year==2000)):(which(yearly$Year==2014))])

# or, set enlistment to 0 to try and get just the rates based on the data's assumption of no new vets
yearly$TotalEnlistment[enlistmentIndices] = 0

# set parameters for cross-calibration
nGuesses = 20
nFolds = 5

# now run the k-fold calibration - warning, takes a bit
vetParams <- veterans.calibrateRates(vetMins, vetMaxs, nGuesses, nFolds, yearly)

# with these parameters, we can assess k-fold cross-validation error
yIndices <- which(is.na(yearlyData$TotalVeterans) == FALSE)
yGroups <- cut(1:length(yIndices), nFolds, label = FALSE)

# calculate the calibration and cross-validation r^2 and rmse
vetCorr <- matrix(nrow = nFolds, ncol = 2)
vetRMSE <- matrix(nrow = nFolds, ncol = 2)

for (i in 1:nFolds){
  
  # apply the model
  vetApplied = veterans.deathDelistmentRates(colMeans(vetParams[[i]]), yearly$Year[2:nrow(yearly)], yearly)
  
  # add the t0 data (since we model from year 2+ using t0 as a parameter)
  vetApplied = append(colMeans(vetParams[[i]])[3], vetApplied)
  
  # set calibration and validation r^2
  vetCorr[i, 1] = cor(yearly$TotalVeterans[yIndices[yGroups == i]], vetApplied[yIndices[yGroups == i]])
  vetCorr[i, 2] = cor(yearly$TotalVeterans[yIndices[yGroups != i]], vetApplied[yIndices[yGroups != i]])
  
  # set calibration and validation RMSE
  vetRMSE[i, 1] <- sqrt(mean((yearly$TotalVeterans[yIndices[yGroups == i]] - vetApplied[yIndices[yGroups == i]])^2))
  vetRMSE[i, 2] <- sqrt(mean((yearly$TotalVeterans[yIndices[yGroups != i]] - vetApplied[yIndices[yGroups != i]])^2))
}

# the correlation coefficients are all 0.98 and above, so we won't bother plotting, but we'll plot the cal/val rmses
ylim = c(min(vetRMSE), max(vetRMSE))
xlab = 'Fold'
ylab = 'RMSE'
title = "Veteran population RMSE by fold"
legend = c("Calibration", "Validation")
barplot(vetRMSE[,1], col = colReal, ylim = ylim, main = title, xlab = xlab, ylab = ylab)
par(new=TRUE)
barplot(vetRMSE[,2], col = colPred, ylim = ylim, main = title, xlab = xlab, ylab = ylab)
legend("top", legend = legend, col = c(colReal, colPred), pch = c(19,19))

# let's perform some sensitivity analysis of this model

# define the parameter names
vetNames <- c("DeathRate", "DelistmentRate", "t0 Veteran Population")

# set n parameters
k <- length(vetMins)

# n times to compute effect for each parameter
r <- 10

# n possible levels for each parameter (should be an even #)
p <- 4

# the increment to adjust parameter values by
delta <- p / (2 * (p - 1))

# run the morris script, which produces a plot
script.morris(vetMins, vetMaxs, vetNames, "Morris Sensitivity Analysis\nVeteran Populations", k, r, p, delta, veterans.deathDelistmentRates)

# and run the vsa script
nRandom = 500
vetVsa <- vsa(veterans.deathDelistmentRates, vetMins, vetMaxs, nRandom)

# then produce our own plot
xlab <- "Main effects"
ylab <- "Total effects"
title <- "Variance-based Sensitivity Analysis\nVeteran Populations"
pch <- rep(19, k)
cols <- rainbow(k)
cex <- rep(3, k)
plot(vetVsa[[1]], vetVsa[[2]], pch = pch, main = title, col = cols, xlab = xlab, ylab = ylab)
legend("topleft", legend = vetNames, pch = pch, col = cols)

# then we'll plot some of the outputs
foldToUse = 5
yearly$TotalEnlistment[enlistmentIndices] = 0
vetApplied = veterans.deathDelistmentRates(colMeans(vetParams[[foldToUse]]), yearly$Year[2:nrow(yearly)], yearly)
vetApplied = append(colMeans(vetParams[[foldToUse]])[3], vetApplied)
vetAppliedOrig = vetApplied
yearly$TotalEnlistment[enlistmentIndices] = mean(yearly$TotalEnlistment[(which(yearly$Year==2000)):(which(yearly$Year==2014))])
vetApplied = veterans.deathDelistmentRates(colMeans(vetParams[[foldToUse]]), yearly$Year[2:nrow(yearly)], yearly)
vetApplied = append(colMeans(vetParams[[foldToUse]])[3], vetApplied)
ylim = c(min(append(vetApplied, append(vetAppliedOrig, yearly$TotalVeterans)), na.rm=TRUE), 
         max(append(vetApplied, append(vetAppliedOrig, yearly$TotalVeterans)), na.rm=TRUE))
xlab = "Year"
ylab = "Veteran Population"
title = "Modeled vs Real US Veteran Population"
legend = c("No new veterans", "With new veterans", "VA Projections")
cols = c(colMil, colFed, colReal)
lwd = c(2, 2, 2)
lty = c(1, 1, 2)
plot(yearly$Year, vetApplied, type='l', xlab = xlab, ylab = ylab, main = title, lwd = lwd[2], col = cols[2], ylim=ylim)
par(new = TRUE)
plot(yearly$Year, vetAppliedOrig, type='l', xlab = xlab, ylab = ylab, main = title, lwd = lwd[1], col = cols[1], ylim=ylim)
par(new = TRUE)
plot(yearly$Year, yearly$TotalVeterans, type='l', xlab = xlab, ylab = ylab, main = title, lwd = lwd[3], col = cols[3], ylim=ylim, lty=lty[3])
legend("bottom", legend=legend, lty=lty, col=cols, lwd=lwd)

#############################
# FINAL MODELING EXERCISE

# ok, so now that we have veteran population data to analyze, we'll start running our actual model, where we vary
#  a) when demilitarization occurs
#  b) over how man years it occurs
#  c) whether we use historic or modeled data
#
# the idea is that we'll have two potential model outputs to measure the effect of demilitarization
#  a) the highest ratio of people added to the labor force times the amount of money freed from demilitarization, and 
#  b) the total money saved on veteran's costs (for n years beyond the start of demilitarization)

# perform SA on the DLI (labor index) model

# define the parameter names
pNames <- c("newLabor", "Population", "laborParticipation", "Unemployment", "newSpending", "fedSpending")

# set mins/maxs
minDelist <- 0.1 * min(yearly$TotalEnlistment, na.rm = TRUE)
maxDelist <- 0.1 * max(yearly$TotalEnlistment, na.rm = TRUE)
minPop <- min(yearly$TotalPopulation, na.rm = TRUE)
maxPop <- max(yearly$TotalPopulation, na.rm = TRUE)
minLabor <- min(yearly$LaborParticipationRate, na.rm = TRUE)
maxLabor <- max(yearly$LaborParticipationRate, na.rm = TRUE)
minUnemp <- min(yearly$UnemploymentRate, na.rm = TRUE)
maxUnemp <- max(yearly$UnemploymentRate, na.rm = TRUE)
minMil <- 0.1 * min(yearly$MilitarySpending.BUSD, na.rm = TRUE)
maxMil <- 0.1 * max(yearly$MilitarySpending.BUSD, na.rm = TRUE)
minFed <- min(yearly$FedSpending.BUSD, na.rm = TRUE)
maxFed <- max(yearly$FedSpending.BUSD, na.rm = TRUE)
pMins <- c(minDelist, minPop, minLabor, minUnemp, minMil, minFed)
pMaxs <- c(maxDelist, maxPop, maxLabor, maxUnemp, maxMil, maxFed)

# set n parameters
k <- length(pMins)

# n times to compute effect for each parameter
r <- 100

# n possible levels for each parameter (should be an even #)
p <- 10

# the increment to adjust parameter values by
delta <- p / (2 * (p - 1))

# run the morris script, which produces a plot
script.morris(pMins, pMaxs, pNames, "Morris Sensitivity Analysis\nInflation Adjusted Demilitarized Labor Index", k, r, p, delta, labor.index.sa)

# then we'll run VSA
nRandom = 500
dliVsa <- vsa(labor.index.sa, pMins, pMaxs, nRandom)

# then produce our own plot
xlab <- "Main effects"
ylab <- "Total effects"
title <- "Variance-based Sensitivity Analysis\nInflation Adjusted Demilitarized Labor Index"
pch <- rep(19, k)
cols <- rainbow(k)
cex <- rep(3, k)
plot(dliVsa[[1]], dliVsa[[2]], pch = pch, main = title, col = cols, xlab = xlab, ylab = ylab, cex=cex)
legend("topleft", legend = pNames, pch = pch, col = cols)

