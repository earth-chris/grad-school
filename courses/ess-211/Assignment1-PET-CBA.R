# ESS-211 assignment 1, due 10-12-16
#  cba 10/2016
# 
# below are the retrieved values from the 6 PET models. These can be retrieved from
#  sourcing this script
#
#  "Priestly Taylor: 4.51250168465471"
#  "Modified Priestly Taylor: 5.706963375975"
#  "Hammon: 3.5098169976678"
#  "Hargreaves: 4.92657369213604"
#  "Linacre: 4.44736768602802"
#  "Turc: 7.52164961720322"

# clear the workspace to remove any environmental conflicts
rm(list = ls())

# load the associated helper functions for this assignment
source("Assignment1-PET-Functions.R")

# write a function to calculate the mean temperature as it is not given
meanTemp <- function(tMin, tMax){
    tMean <- (tMin + tMax)/2
}

###
# below are functions for the 6 methods for calculating PET
###

# method 1: priestly taylor
priestlyTaylor <- function(doy, tMax, tMin, tMean, RH, tDew, Rs, elev, lat){
    # define constant 'a'
    a <- 1.26
 
    # calculate slope of vapor pressure curve using temp above
    vpSlope <- s(tMean)
 
    # calculate the psychometric constant
    gamma <- psychro(elev)
    
    # set the ground heat flux to zero
    groundHeatFlux <- 0
    
    # calculate net radiation
    netRadiation <- Rn(lat, elev, doy, tDew, tMax, tMin, Rs)
 
    # set lambda
    lambda <- 2.501 - (0.002361 * tMean)
    
    # run the PET calculation
    PET <- a * ((vpSlope * (netRadiation - groundHeatFlux)) / 
      (lambda * (vpSlope + gamma)))
    return(PET)
}

# method 2: modified priestly taylor
modifiedPriestlyTaylor <- function(tMax, tMin, Rs){
    # define constant 'albedo'
    albedo <- 0.23
    
    # calculate constant TD
    TD <- (0.6 * tMax) + (0.4 * tMin)
    
    # calculate EEQ
    EEQ <- Rs * (4.88e-3 - (4.37e-3 * albedo)) * (TD + 29)
    
    # PET is calculated differently based on tMax
    if (tMax < 5){
       PET <- EEQ * 0.01 * exp(0.18 * (tMax + 20))
    } else if (tMax > 24){
        PET <- EEQ * ((tMax -24) * 0.05 + 1.1)
    } else {
        PET <- EEQ * 1.1
    }
    return(PET)
}

# method 3: hammon
hammon <- function(doy, tMax, tMin, tMean, lat){
    # calculate the day length
    dl <- daylength(lat, doy)
	  
    # transform day length to fraction of a day (i.e. 0-1 with 0 as no sunlight, 1 as all-day light
	  dl <- dl / 24
    
    # calculate the saturated vapor pressure
    vpSat <- esat(tMax, tMin)
    
    # run the PET calculation
    PET <- 715.5 * dl * (vpSat / (tMean + 273.2))
    return(PET)
}

# method 4: hargreaves
hargreaves <- function(doy, tMax, tMin, tMean){
    # calculate lambda
    lambda <- 2.501 - (0.002361 * tMean)
    
    # calculate atmospheric net radiation
    atmRad <- Ra(doy, lat)
    
    # run the PET calculation
    PET <- (0.0023 * (tMean + 17.8) * ((tMax - tMin)^0.5) * atmRad) / lambda
    return(PET)
}

# method 5: linacre
linacre <- function(tMean, elev, lat){
    # calculate Tm
    Tm <- tMean + (0.006 * elev)
    
    # run the PET calculation
    PET <- (500 * (Tm / (100 - lat)) + (15 * (tMean - tDew))) / (80 - tMean)
    return(PET)
}

# method 6: turc
turc <- function(tMean, RH, Rs){
    # convert Rs from units of MJ/m^2 to cal/cm^2
	  RsCal <- Rs * 23.9
    
    # PET is calculated differently based on relative humidity
    if (RH < 50){
        PET <- (1 + ((50 - RH) / 70)) * ((0.013 * tMean * (RsCal + 50)) / (tMean + 15))
    } else {
       PET <- (0.013 * tMean * (RsCal + 50)) / (tMean + 15) 
    }
    return(PET)
}

###
# below are variables and commands to run each function using the following info
#  Report PET according to each method for:
#  date : Aug 21, 2013
#  Tmax : 27.7 C
#  Tmin : 13.3 C
#  RH   : 67%
#  Tdew : 13.9C
#  Rs   : 22.5 MJ/m^2
###

# set day of year using julian calendar
doyString <- "August 21, 2013"
doy <- 233

# set max temp
tMax <- 27.7

# set min temp
tMin <- 13.3

# set dew temp
tDew <- 13.9

# set relative humidity
RH <- 0.67

# set the input solar radiation
Rs <- 22.5

# we need elevation, lat lon, and mean temp, which are not provided in the assignment
# set elevation and lat/lon from this website
#  http://www.latlong.net/place/gilroy-ca-usa-6030.html
elev <- 61
lat <- 37.005783
lon <- -121.568275

# calculate the mean temperature since it is used in several functions
tMean <- meanTemp(tMin, tMax)

# report the parameters to stdOut
print("ESS-211 Assignment 1 - Christopher Anderson")
print("Running with the following parameters:")
print(paste0("Day of the year : ", doyString))
print(paste0("Julian date     : ", doy))
print(paste0("Min temp (deg C): ", tMin))
print(paste0("Max temp (deg C): ", tMax))
print(paste0("Dew temp (dec C): ", tDew))
print(paste0("Rel humidity (%): ", RH))
print(paste0("Solar radiation (MJ/m^2) : ", Rs))
print("----------")

# run each PET model and report results
print(paste0("Priestly Taylor: ", priestlyTaylor(doy, tMax, tMin, tMean, RH, tDew, Rs, elev, lat)))
print(paste0("Modified Priestly Taylor: ", modifiedPriestlyTaylor(tMax, tMin, Rs)))
print(paste0("Hammon: ", hammon(doy, tMax, tMin, tMean, lat)))
print(paste0("Hargreaves: ", hargreaves(doy, tMax, tMin, tMean)))
print(paste0("Linacre: ", linacre(tMean, elev, lat)))
print(paste0("Turc: ", turc(tMean, RH, Rs)))