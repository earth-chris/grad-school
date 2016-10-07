# ESS-211 assignment 1, due 10-12-16
#  cba 10/2016

###
# below are the functions provided in PET.helper.functions.r
###

# saturated vapor pressure for a given temperature
e0 <- function(temp){
    0.6108*exp(17.27*temp/(temp+237.3))
}

# this is how FAO suggests one do it for daily esat b/c of nonlinearity of e0 function 
esat <- function(tmax,tmin){ 
    (e0(tmax)+e0(tmin))/2  
}

# slope of saturation curve
s <- function(temp){
    4098*(0.6108*exp(17.27*temp/(temp+237.3)))/(temp+237.3)^2    
}

# extraterrestrial net radiation (MJ / m^2 / day)
Ra <- function(doy,lat){   #lat in decimal degrees       #p79 FAO
  psi=(pi/180)*lat
  delta=0.409*sin((2*pi*doy/365)-1.39)
  ws=acos(-tan(psi)*tan(delta))
  dr=1+0.033*cos(2*pi*doy/365)
  Ra=(24*60/pi)*0.082*dr*(ws*sin(psi)*sin(delta)+cos(psi)*cos(delta)*sin(ws))
  Ra
}

# psychrometric constant [kPa C-1]
psychro <- function(elev){             #elev in m
  press=101.3*((293-0.0065*elev)/293)^5.26
  0.665e-3*press
}

# daylength (returns # hours in the day with sunlight. only for areas between 65N-6)
daylength <-function(lat, doy){
  psi=(pi/180)*lat
  delta=0.409*sin((2*pi*doy/365)-1.39)
  ws=acos(-tan(psi)*tan(delta))
  d=24*ws/pi
  d
}
 
# net surface radiation
Rn <- function(lat, elev, doy, dew, tmax, tmin, Rs){
  ea=e0(dew)
  Ra=Ra(doy,lat)
  Rso=(0.75+(2e-5)*elev)*Ra          #p85
  Rnl=4.903e-9*(((tmax+273)^4+(tmin+273)^4)/2)*(0.34-0.14*sqrt(ea))*(1.35*(Rs/Rso)-0.35)  #p86
  Rn=(1-0.23)*Rs-Rnl             #(Rns=(1-0.23)Rs) ; Rn=Rns-Rnl p87
  Rn
}

###
# below are various functions by cba
###

tMean <- function(tMin, tMax){
    tMean <- (tMin + tMax)/2
}

###
# below are the 6 methods for calculating PET
###

# method 1: priestly taylor
priestlyTaylor <- function(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon){
    # define constant 'a'
    a <- 1.26
 
    # calculate the mean temperature for the day
    meanTemp <- tMean(tMin, tMax)
 
    # calculate slope of vapor pressure curve using temp above
    vpSlope <- s(meanTemp)
 
    # calculate the psychometric constant
    gamma <- psychro(elev)
    
    # set the ground heat flux to zero
    groundHeatFlux <- 0
    
    # calculate net radiation
    netRadiation <- Rn(lat, elev, doy, tDew, tMax, tMin, Rs)
 
    # set lambda
    lambda <- 2.501 - (0.002361 * meanTemp)
    
    # run the PET calculation
    PET <- a * ((vpSlope * (netRadiation - groundHeatFlux)) / 
      (lambda * (vpSlope + gamma)))
    return(PET)
}

# method 2: modified priestly taylor
modifiedPriestlyTaylor <- function(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon){
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
        EEQ * ((tMax -24) * 0.05 + 1.1)
    } else {
        PET <- EEQ * 1.1
    }
    return(PET)
}

# method 3: hammon
hammon <- function(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon){
    # calculate the day length
    dl <- daylength(lat, doy)
	
	# transform day length to fraction of a day (i.e. 0-1 with 0 as no sunlight, 1 as all-day light
	dl <- dl / 24
    
    # calculate the saturated vapor pressure
    vpSat <- esat(tMax, tMin)
    
    # calculate the mean temperature
    meanTemp <- tMean(tMin, tMax)
    
    # run the PET calculation
    PET <- 715.5 * dl * (vpSat / (meanTemp + 273.2))
    return(PET)
}

# method 4: hargreaves
hargreaves <- function(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon){
    # calculate mean temperature
    meanTemp <- tMean(tMin, tMax)
    
    # calculate lambda
    lambda <- 2.501 - (0.002361 * meanTemp)
    
    # calculate atmospheric net radiation
    atmRad <- Ra(doy, lat)
    
    # run the PET calculation
    PET <- (0.0023 * (meanTemp + 17.8) * ((tMax - tMin)^0.5) * atmRad) / lambda
    return(PET)
}

# method 5: linacre
linacre <- function(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon){
    # calculate mean temp
    meanTemp <- tMean(tMin, tMax)
    
    # calculate Tm
    Tm <- meanTemp + (0.006 * elev)
    
    # run the PET calculation
    PET <- (500 * (Tm / (100 - lat)) + (15 * (meanTemp - tDew))) / (80 - meanTemp)
    return(PET)
}

# method 6: turc
turc <- function(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon){
    # calculate mean temp
    meanTemp <- tMean(tMin, tMax)
	
	# convert Rs from units of MJ/m^2 to cal/cm^2
	RsCal <- Rs * 23.9
    
    # PET is calculated differently based on relative humidity
    if (RH < 50){
        PET <- (1 + ((50 - RH) / 70)) * ((0.013 * meanTemp * (RsCal + 50)) / (meanTemp + 15))
    } else {
       PET <- (0.013 * meanTemp * (RsCal + 50)) / (meanTemp + 15) 
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

# set elevation and lat/lon from this website
#  http://www.latlong.net/place/gilroy-ca-usa-6030.html
elev <- 61
lat <- 37.005783
lon <- -121.568275

# calculate the mean temperature since it is used in several functions
meanTemp <- tMean(tMin, tMax)

# calculate lambda since it is used in several functions
lambda <- 2.501 - (0.002361 * meanTemp)

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
print(paste0("Priestly Taylor: ", priestlyTaylor(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon)))
print(paste0("Modified Priestly Taylor: ", modifiedPriestlyTaylor(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon)))
print(paste0("Hammon: ", hammon(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon)))
print(paste0("Hargreaves: ", hargreaves(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon)))
print(paste0("Linacre: ", linacre(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon)))
print(paste0("Turc: ", turc(doy, tMax, tMin, RH, tDew, Rs, elev, lat, lon)))