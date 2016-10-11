# contains frequently used functions for ESS-211

#############################
# functions provided by teaching group for assignment 1 (PET)

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

#############################
# functions written by cba for assignment 1 (PET)

# method 1: priestly taylor
petPriestlyTaylor <- function(doy, tMax, tMin, tMean, RH, tDew, Rs, elev, lat){
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
petModifiedPriestlyTaylor <- function(tMax, tMin, Rs){
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
petHammon <- function(doy, tMax, tMin, tMean, lat){
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
petHargreaves <- function(doy, tMax, tMin, tMean){
  # calculate lambda
  lambda <- 2.501 - (0.002361 * tMean)
  
  # calculate atmospheric net radiation
  atmRad <- Ra(doy, lat)
  
  # run the PET calculation
  PET <- (0.0023 * (tMean + 17.8) * ((tMax - tMin)^0.5) * atmRad) / lambda
  return(PET)
}

# method 5: linacre
petLinacre <- function(tMean, elev, lat){
  # calculate Tm
  Tm <- tMean + (0.006 * elev)
  
  # run the PET calculation
  PET <- (500 * (Tm / (100 - lat)) + (15 * (tMean - tDew))) / (80 - tMean)
  return(PET)
}

# method 6: turc
petTurc <- function(tMean, RH, Rs){
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

#############################