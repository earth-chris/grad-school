#############################
# functions provided by the teaching group for assignment 1 (PET)

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