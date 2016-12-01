# ESS-211 assignment 5, due 12-2-16
#  cba 12/2016

# set current working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# source the course functions
source("ESS-211-Functions.R")

#############################
# task 1 - downloading data from nasa-larc

# set the names of the sites
siteNames <- c("MtTamalpais", "LasCruces", "Tambopata", "Bakossi", "Hainan")
siteList <- list()
siteElev <- rep(0, length(siteNames))

# set output file names to download to our working directory
weatherFiles <- paste0(siteNames, "WeatherData.txt")

# set the lat/lon coordinates for each site
lat.lon <- matrix(nrow = length(siteNames), ncol = 2)
lat.lon[1,] <- c(37.923992, -122.596555)
lat.lon[2,] <- c(8.790813, -82.954302)
lat.lon[3,] <- c(-13.075157, -69.471705)
lat.lon[4,] <- c(5.052105, 9.522459)
lat.lon[5,] <- c(19.170315, 109.765304)

# set the start and end dates
monthStart <- 12
monthEnd <- monthStart
dayStart <- 1
dayEnd <- dayStart
yearStart <- 2006
yearEnd <- 2016

# set names to use in the output data frames
dfNames <- c("WEYR", "WEDAY",  "SRAD", "TMAX", "TMIN", "RAIN", "WIND", "TDEW", "T2M", "RH2M")

# create the strings necessary to download from a url
baseUrl <- "https://power.larc.nasa.gov/cgi-bin/agro.cgi?email=&step=1"
urlDataRequest <- "&submit=Yes&p=toa_dwn&p=swv_dwn&p=lwv_dwn&p=T2M&p=T2MN&p=T2MX&p=RH2M&p=DFP2M&p=RAIN&p=WS10M"

# loop through each of our sites and download the text file from a URL
for (i in seq(1, length(siteNames))){
  urlLat <- paste0("&lat=", lat.lon[i,1])
  urlLon <- paste0("&lon=", lat.lon[i,2])
  urlMs <- paste0("&ms=", monthStart)
  urlDs <- paste0("&ds=", dayStart)
  urlYs <- paste0("&ys=", yearStart)
  urlMe <- paste0("&me=", monthEnd)
  urlDe <- paste0("&de=", dayEnd)
  urlYe <- paste0("&ye=", yearEnd)
  url <- paste0(baseUrl, urlLat, urlLon, urlMs, urlDs, urlYs, urlMe, urlDe, urlYe, urlDataRequest)
  
  # downlaod the file from the url
  download.file(url, weatherFiles[i])
  
  # read the file into a data frame, then add it to our list
  dfSite <- read.table(weatherFiles[i], skip=14)
  siteList[[i]] <- dfSite
  
  # read the elevation values for the extra credit
  fileRef <- file(weatherFiles[i])
  fileLines <- readLines(fileRef)
  for (j in seq(1, length(fileLines))){
    lineStr <- strsplit(fileLines[j], "\\s+")[[1]]
    if (length(lineStr) > 5){
      if (lineStr[5] == "WELEV"){
        mdStr <- strsplit(fileLines[j+1], "\\s+")[[1]]
        siteElev[i] <- mdStr[5]
      }
    }
  }
  close(fileRef)
}

#############################
# problem 2 - world map plotting