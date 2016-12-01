# ESS-211 assignment 5, due 12-2-16
#  cba 12/2016

# set current working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# source the course functions
source("ESS-211-Functions.R")
library(maps)
library(dplyr)
library(zoo)

#############################
# task 1 - downloading data from nasa-larc

# set the names of the sites
siteNames <- c("MtTamalpais", "LasCruces", "Tambopata", "Bakossi", "Hainan")
siteList <- list()
siteElev <- rep(0, length(siteNames))
siteColors <- rainbow(5)

# set output file names to download to our working directory
weatherFiles <- paste0(siteNames, "WeatherData.txt")

# set the lat/lon coordinates for each site
lat.lon <- matrix(nrow = length(siteNames), ncol = 2)
lat.lon[1,] <- c(37.923992, -122.596555)
lat.lon[2,] <- c(8.790813, -82.954302)
lat.lon[3,] <- c(-13.075157, -69.471705)
lat.lon[4,] <- c(5.052105, 9.522459)
lat.lon[5,] <- c(19.170315, 109.765304)

# set summer dates
summer.N <- seq(171, 265) # june 20-sept 22
summer.S <- c(seq(1, 79), seq(355, 366)) # dec 21-mar 20

# enter the no-data value
noData <- -99

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
  names(dfSite) <- dfNames
  
  # set nodata values to NA, then fix them
  dfSite[which(dfSite == noData, arr.ind = TRUE)] <- NA
  dfSite <- na.approx(dfSite)
  
  # add the data frame to the list
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

# loop through each site and get the mean summer temp and precip
meanPrecip <- c()
meanTemp <- c()
for (i in seq(1,length(siteNames))){
  
  # find the indices of summer values based on latitude
  if (lat.lon[i,1] >= 0){
    summerVec <- whichVec(siteList[[i]][,"WEDAY"], summer.N)
  } else {
    summerVec <- whichVec(siteList[[i]][,"WEDAY"], summer.S)
  }
  
  # find mean temp and precip values
  meanPrecip[i] <- mean(siteList[[i]][summerVec,"RAIN"], na.rm = TRUE)
  meanTemp[i] <- mean((siteList[[i]][summerVec, "TMAX"] + siteList[[i]][summerVec, "TMIN"])/2, na.rm = TRUE)
  
  # add a new column to the data frame here to specify that we will be comparing precip and temp during summer
  nSiteRows <- nrow(siteList[[i]])
  inSeasonVec <- rep(FALSE, nSiteRows)
  inSeasonVec[summerVec] <- TRUE
  siteList[[i]]["inSeason"] <- inSeasonVec
}

# set min/max temp values to scale the plot symbols
cex.min <- 1
cex.max <- 3
tmp.min <- min(meanTemp)
tmp.max <- max(meanTemp)
cex.size <- ((meanTemp - tmp.min) / (tmp.max - tmp.min)) * (cex.max - cex.min) + cex.min

# set up the base world map
map(fill = TRUE, col = add.alpha("Black", 0.075))

# add a legend
legend("bottomleft", col = siteColors, legend = siteNames, cex=0.9, pch = rep(19, length(siteNames)))

# loop through each of the sites and plot the names, coordinates, mean precip and temp for the sites
for (i in seq(1, length(siteNames))){
  points(lat.lon[i,2], lat.lon[i,1], cex = cex.size[i], pch = 19, col = siteColors[i])
  text(lat.lon[i,2], lat.lon[i,1], labels = paste("MeanP:", format(round(meanPrecip[i], 1), nsmall = 1), 
      "\nMeanT:", format(round(meanTemp[i], 1), nsmall=1)), pos = 2, cex=0.8, bg="White")
  #legend(lat.lon[i,2], lat.lon[i,1], legend=c(paste("MeanP:", format(round(meanPrecip[i], 1), nsmall = 1)), 
  #       paste("nMeanT:", format(round(meanTemp[i], 1), nsmall=1))), cex=0.8)
}

#############################
# problem 3 computing seasonal PET for each year






