# ESS-211 assignment 5, due 12-2-16
#  cba 12/2016

# set current working directory
setwd("~/cba/aei-grad-school/courses/ess-211/")

# source the course functions
source("ESS-211-Functions.R")
library(maps)
library(dplyr)
library(zoo)
library(gplots)

#############################
# task 1 - downloading data from nasa-larc

# set the names of the sites
siteNames <- c("MtTamalpais", "LasCruces", "Tambopata", "Bakossi", "Hainan")
siteList <- list()
nSites <- length(siteNames)
siteElev <- rep(0, nSites)
siteColors <- rainbow(5)

# set output file names to download to our working directory
weatherFiles <- paste0(siteNames, "WeatherData.txt")

# set the lat/lon coordinates for each site
lat.lon <- matrix(nrow = nSites, ncol = 2)
lat.lon[1,] <- c(37.923992, -122.596555)
lat.lon[2,] <- c(8.790813, -82.954302)
lat.lon[3,] <- c(-13.075157, -69.471705)
lat.lon[4,] <- c(5.052105, 9.522459)
lat.lon[5,] <- c(19.170315, 109.765304)

# set summer dates
summer.N <- seq(171, 265) # june 20-sept 22
summer.S <- c(seq(1, 79), seq(355, 366)) # dec 21-mar 20

# set growing season dates
growingSeason <- seq(60, 120) # mar 1 - apr 30

# enter the no-data value
noData <- -99

# set the start and end dates
monthStart <- 12
monthEnd <- monthStart
dayStart <- 1
dayEnd <- dayStart
yearStart <- 2006
yearEnd <- 2016
years <- yearStart:yearEnd
nYears <- yearEnd - yearStart

# set names to use in the output data frames
dfNames <- c("WEYR", "WEDAY",  "SRAD", "TMAX", "TMIN", "RAIN", "WIND", "TDEW", "T2M", "RH2M")

# create the strings necessary to download from a url
baseUrl <- "https://power.larc.nasa.gov/cgi-bin/agro.cgi?email=&step=1"
urlDataRequest <- "&submit=Yes&p=toa_dwn&p=swv_dwn&p=lwv_dwn&p=T2M&p=T2MN&p=T2MX&p=RH2M&p=DFP2M&p=RAIN&p=WS10M"

# loop through each of our sites and download the text file from a URL
for (i in seq(1, nSites)){
  
  # set up urls
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
  dfSite <- data.frame(na.approx(dfSite, na.rm = TRUE))
  
  # add a column for mean temp
  dfSite["TMEAN"] = (dfSite[, "TMAX"] + dfSite[, "TMIN"]) / 2
  
  # add new columns for lat/lon
  dfSite["LAT"] <- rep(lat.lon[i,1], nrow(dfSite))
  dfSite["LON"] <- rep(lat.lon[i,2], nrow(dfSite))
  
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
    
    # add elevation column
    dfSite["ELEV"] <- rep(as.numeric(siteElev[i]), nrow(dfSite))
    
    # add the data frame to the list
    siteList[[i]] <- dfSite
  }
  close(fileRef)
}

#############################
# problem 2 - world map plotting

# loop through each site and get the mean summer temp and precip
meanPrecip <- c()
meanTemp <- c()
for (i in seq(1, nSites)){
  
  # find the indices of summer values based on latitude
  if (lat.lon[i,1] >= 0){
    summerVec <- whichVec(siteList[[i]][,"WEDAY"], summer.N)
  } else {
    summerVec <- whichVec(siteList[[i]][,"WEDAY"], summer.S)
  }
  
  # find mean temp and precip values
  meanPrecip[i] <- mean(siteList[[i]][summerVec,"RAIN"], na.rm = TRUE)
  meanTemp[i] <- mean(siteList[[i]][summerVec, "TMEAN"], na.rm = TRUE)
  
  # add a new column for growing season
  nSiteRows <- nrow(siteList[[i]])
  inSeasonVec <- rep(FALSE, nSiteRows)
  growingVec <- whichVec(siteList[[i]][,"WEDAY"], growingSeason)
  inSeasonVec[growingVec] <- TRUE
  siteList[[i]]["inSeason"] <- inSeasonVec
}

# set min/max temp values to scale the plot symbols
cex.min <- 1
cex.max <- 3
tmp.min <- min(meanTemp)
tmp.max <- max(meanTemp)
cex.size <- ((meanTemp - tmp.min) / (tmp.max - tmp.min)) * (cex.max - cex.min) + cex.min

# set up the base world map
map(fill = TRUE, col = add.alpha("Black", 0.075), bg = "Grey")

# add a legend
legend("bottomleft", col = siteColors, legend = siteNames, cex=0.9, pch = rep(19, length(siteNames)))

# loop through each of the sites and plot the names, coordinates, mean precip and temp for the sites
for (i in seq(1, nSites)){
  points(lat.lon[i,2], lat.lon[i,1], cex = cex.size[i], pch = 19, col = siteColors[i])
  text(lat.lon[i,2], lat.lon[i,1], labels = paste("MeanP:", format(round(meanPrecip[i], 1), nsmall = 1), 
      "\nMeanT:", format(round(meanTemp[i], 1), nsmall=1)), pos = 2, cex=0.8, bg="White", col="White")
  #legend(lat.lon[i,2], lat.lon[i,1], legend=c(paste("MeanP:", format(round(meanPrecip[i], 1), nsmall = 1)), 
  #       paste("nMeanT:", format(round(meanTemp[i], 1), nsmall=1))), cex=0.8)
}

#############################
# problem 3 computing seasonal PET for each year

# this comment is where i register my displeasure with having to use dplyr here. 

# create vectors/arrays to store outputs
petMethods <- c("PriestlyTaylor", "ModPriesTaylor", "Hamon", "Hargreaves", "Linacre", "Turc")
nMethods <- length(petMethods)
petMeans <- array(dim=c(nSites, nMethods, nYears))
petStdvs <- array(dim=c(nSites, nMethods, nYears))
meanPrecip <- array(dim=c(nSites, nMethods))

for (i in seq(1, nSites)){
  petGroup <- group_by(filter(siteList[[i]], inSeason==TRUE), WEYR)
  petSummary <- summarize(petGroup,
    mean.PriestlyTaylor = mean(pet.priestlyTaylor(WEDAY, TMAX, TMIN, TMEAN, RH2M, TDEW, SRAD, ELEV, LAT), na.rm=TRUE),
    mean.ModifiedPriestlyTaylor = mean(pet.modifiedPriestlyTaylor(TMAX, TMIN, SRAD), na.rm=TRUE),
    mean.Hamon = mean(pet.hammon(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    mean.Hargreaves = mean(pet.hargreaves(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    mean.Linacre = mean(pet.linacre(TMEAN, ELEV, LAT, TDEW), na.rm=TRUE),
    mean.Turc = mean(pet.turc(TMEAN, RH2M, SRAD), na.rm=TRUE),
    mean.Precip = mean(RAIN),
    sd.PriestlyTaylor = sd(pet.priestlyTaylor(WEDAY, TMAX, TMIN, TMEAN, RH2M, TDEW, SRAD, ELEV, LAT), na.rm=TRUE),
    sd.ModifiedPriestlyTaylor = sd(pet.modifiedPriestlyTaylor(TMAX, TMIN, SRAD), na.rm=TRUE),
    sd.Hamon = sd(pet.hammon(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    sd.Hargreaves = sd(pet.hargreaves(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    sd.Linacre = sd(pet.linacre(TMEAN, ELEV, LAT, TDEW), na.rm=TRUE),
    sd.Turc = sd(pet.turc(TMEAN, RH2M, SRAD), na.rm=TRUE))
  
  # assign this absurd stuff to an output matrix
  petMeans[i, 1,] <- petSummary$mean.PriestlyTaylor
  petMeans[i, 2,] <- petSummary$mean.ModifiedPriestlyTaylor
  petMeans[i, 3,] <- petSummary$mean.Hamon
  petMeans[i, 4,] <- petSummary$mean.Hargreaves
  petMeans[i, 5,] <- petSummary$mean.Linacre
  petMeans[i, 6,] <- petSummary$mean.Turc
  petStdvs[i, 1,] <- petSummary$sd.PriestlyTaylor
  petStdvs[i, 2,] <- petSummary$sd.ModifiedPriestlyTaylor
  petStdvs[i, 3,] <- petSummary$sd.Hamon
  petStdvs[i, 4,] <- petSummary$sd.Hargreaves
  petStdvs[i, 5,] <- petSummary$sd.Linacre
  petStdvs[i, 6,] <- petSummary$sd.Turc
  meanPrecip[i,] <- petSummary$mean.Precip
  
  # create new columns to analyze in problem 5
  #siteList[[i]]["meanPriestlyTaylor"] = rep(petSummary$mean.PriestlyTaylor, nrow(siteList[[i]]))
  #siteList[[i]]["meanModifiedPriestlyTaylor"] = rep(petSummary$mean.ModifiedPriestlyTaylor, nrow(siteList[[i]]))
  #siteList[[i]]["meanHamon"] = rep(petSummary$mean.Hamon, nrow(siteList[[i]]))
  #siteList[[i]]["meanHargreaves"] = rep(petSummary$mean.Hargreaves, nrow(siteList[[i]]))
  #siteList[[i]]["meanLinacre"] = rep(petSummary$mean.Linacre, nrow(siteList[[i]]))
  #siteList[[i]]["meanTurc"] = rep(petSummary$mean.Turc, nrow(siteList[[i]]))
  #siteList[[i]]["meanPrecip"] = rep(petSummary$mean.Precip, nrow(siteList[[i]]))
}

# prepare the bar plot for these values
colVec <- rep(siteColors, each = nMethods)
xlabVec <- rep(petMethods, nSites)
nVals <- length(meanVec)
title <- "Potential Evapotranspiration by Site and by Method \n Growing Season: March to April"
ylab <- "PET"

# create a 10-panel layout to plot with room for a title
par(mfrow = c(2,5), oma = c(0,0,2,0))

for (i in seq(1,nYears)){
  meanVec <- c(t(petMeans[,,i]))
  stdVec <- c(t(petStdvs[,,i]))
  
  # set up the bar plot
  barplot2(meanVec, names.arg = xlabVec, las = 2, ylab = ylab, col = colVec, plot.ci = TRUE, 
           ci.l = (meanVec - stdVec), ci.u = (meanVec + stdVec), main = years[i])
}

# add a title to the final plot
title(title, outer = TRUE)

#############################
# problem 4 - same as above, but for P - PET

pmPetMeans <- array(dim=c(nSites, nMethods, nYears))
pmPetStdvs <- array(dim=c(nSites, nMethods, nYears))

for (i in seq(1, nSites)){
  petGroup <- group_by(filter(siteList[[i]], inSeason==TRUE), WEYR)
  petSummary <- summarize(petGroup,
    mean.PriestlyTaylor = mean(RAIN - pet.priestlyTaylor(WEDAY, TMAX, TMIN, TMEAN, RH2M, TDEW, SRAD, ELEV, LAT), na.rm=TRUE),
    mean.ModifiedPriestlyTaylor = mean(RAIN - pet.modifiedPriestlyTaylor(TMAX, TMIN, SRAD), na.rm=TRUE),
    mean.Hamon = mean(RAIN - pet.hammon(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    mean.Hargreaves = mean(RAIN - pet.hargreaves(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    mean.Linacre = mean(RAIN - pet.linacre(TMEAN, ELEV, LAT, TDEW), na.rm=TRUE),
    mean.Turc = mean(RAIN - pet.turc(TMEAN, RH2M, SRAD), na.rm=TRUE),
    sd.PriestlyTaylor = sd(RAIN - pet.priestlyTaylor(WEDAY, TMAX, TMIN, TMEAN, RH2M, TDEW, SRAD, ELEV, LAT), na.rm=TRUE),
    sd.ModifiedPriestlyTaylor = sd(RAIN - pet.modifiedPriestlyTaylor(TMAX, TMIN, SRAD), na.rm=TRUE),
    sd.Hamon = sd(RAIN - pet.hammon(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    sd.Hargreaves = sd(RAIN - pet.hargreaves(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    sd.Linacre = sd(RAIN - pet.linacre(TMEAN, ELEV, LAT, TDEW), na.rm=TRUE),
    sd.Turc = sd(RAIN - pet.turc(TMEAN, RH2M, SRAD), na.rm=TRUE))
  
  # assign this absurd stuff to an output matrix
  pmPetMeans[i, 1,] <- petSummary$mean.PriestlyTaylor
  pmPetMeans[i, 2,] <- petSummary$mean.ModifiedPriestlyTaylor
  pmPetMeans[i, 3,] <- petSummary$mean.Hamon
  pmPetMeans[i, 4,] <- petSummary$mean.Hargreaves
  pmPetMeans[i, 5,] <- petSummary$mean.Linacre
  pmPetMeans[i, 6,] <- petSummary$mean.Turc
  pmPetStdvs[i, 1,] <- petSummary$sd.PriestlyTaylor
  pmPetStdvs[i, 2,] <- petSummary$sd.ModifiedPriestlyTaylor
  pmPetStdvs[i, 3,] <- petSummary$sd.Hamon
  pmPetStdvs[i, 4,] <- petSummary$sd.Hargreaves
  pmPetStdvs[i, 5,] <- petSummary$sd.Linacre
  pmPetStdvs[i, 6,] <- petSummary$sd.Turc
}

# prepare the bar plot for these values
colVec <- rep(siteColors, each = nMethods)
xlabVec <- rep(petMethods, nSites)
nVals <- length(meanVec)
title <- "P - PET by Site and by Method \n Growing Season: March to April"
ylab <- "P - PET"

# create a 10-panel layout to plot with room for a title
par(mfrow = c(2,5), oma = c(0,0,2,0))

for (i in seq(1,nYears)){
  meanVec <- c(t(pmPetMeans[,,i]))
  stdVec <- c(t(pmPetStdvs[,,i]))
  
  # set up the bar plot
  barplot2(meanVec, names.arg = xlabVec, las = 2, ylab = ylab, col = colVec, plot.ci = TRUE, 
           ci.l = (meanVec - stdVec), ci.u = (meanVec + stdVec), main = years[i])
}

# add title to the final plot
title(title, outer = TRUE)

#############################
# problem 5 - deriving a metric to asses importance of P vs PET in driving yearly P-PET

# I am choosing to assess the relative importance of P vs PET in driving
#  variability in P - PET by calculating P(yearlyMean) - PET(everyValue)
#  and P(everyValue) - PET(yearlyMean), then subtract each of these from the 
#  original P(everyValue)-PET(everyValue) calculations. smaller residuals
#  should mean the variance in whatever parameter used every value
#  was more important in driving the P-PET signal

# so calculate P(yearlyMean) - PET(everyValue) and P(everyValue) - PET(yearlyMean)
pmeanmPETMeans <- array(dim=c(nSites, nMethods, nYears))
pmPETmeanMeans <- array(dim=c(nSites, nMethods, nYears))

for (i in seq(1, nSites)){
  petGroup <- group_by(filter(siteList[[i]], inSeason==TRUE), WEYR)
  petSummary <- summarize(petGroup,
    mean.PriestlyTaylor = mean(mean(RAIN, na.rm=TRUE) - pet.priestlyTaylor(WEDAY, TMAX, TMIN, TMEAN, RH2M, TDEW, SRAD, ELEV, LAT), na.rm=TRUE),
    mean.ModifiedPriestlyTaylor = mean(mean(RAIN, na.rm=TRUE) - pet.modifiedPriestlyTaylor(TMAX, TMIN, SRAD), na.rm=TRUE),
    mean.Hamon = mean(mean(RAIN, na.rm=TRUE) - pet.hammon(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    mean.Hargreaves = mean(mean(RAIN, na.rm=TRUE) - pet.hargreaves(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE),
    mean.Linacre = mean(mean(RAIN, na.rm=TRUE) - pet.linacre(TMEAN, ELEV, LAT, TDEW), na.rm=TRUE),
    mean.Turc = mean(mean(RAIN, na.rm=TRUE) - pet.turc(TMEAN, RH2M, SRAD), na.rm=TRUE),
    mean2.PriestlyTaylor = mean(RAIN - mean(pet.priestlyTaylor(WEDAY, TMAX, TMIN, TMEAN, RH2M, TDEW, SRAD, ELEV, LAT), na.rm=TRUE), na.rm=TRUE),
    mean2.ModifiedPriestlyTaylor = mean(RAIN - mean(pet.modifiedPriestlyTaylor(TMAX, TMIN, SRAD), na.rm=TRUE), na.rm=TRUE),
    mean2.Hamon = mean(RAIN - mean(pet.hammon(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE), na.rm=TRUE),
    mean2.Hargreaves = mean(RAIN - mean(pet.hargreaves(WEDAY, TMAX, TMIN, TMEAN, LAT), na.rm=TRUE), na.rm=TRUE),
    mean2.Linacre = mean(RAIN - mean(pet.linacre(TMEAN, ELEV, LAT, TDEW), na.rm=TRUE), na.rm=TRUE),
    mean2.Turc = mean(RAIN - mean(pet.turc(TMEAN, RH2M, SRAD), na.rm=TRUE), na.rm=TRUE))
  
  # assign this absurd stuff to an output matrix
  pmeanmPETMeans[i, 1,] <- petSummary$mean.PriestlyTaylor
  pmeanmPETMeans[i, 2,] <- petSummary$mean.ModifiedPriestlyTaylor
  pmeanmPETMeans[i, 3,] <- petSummary$mean.Hamon
  pmeanmPETMeans[i, 4,] <- petSummary$mean.Hargreaves
  pmeanmPETMeans[i, 5,] <- petSummary$mean.Linacre
  pmeanmPETMeans[i, 6,] <- petSummary$mean.Turc
  pmPETmeanMeans[i, 1,] <- petSummary$mean2.PriestlyTaylor
  pmPETmeanMeans[i, 2,] <- petSummary$mean2.ModifiedPriestlyTaylor
  pmPETmeanMeans[i, 3,] <- petSummary$mean2.Hamon
  pmPETmeanMeans[i, 4,] <- petSummary$mean2.Hargreaves
  pmPETmeanMeans[i, 5,] <- petSummary$mean2.Linacre
  pmPETmeanMeans[i, 6,] <- petSummary$mean2.Turc
}

# subtract the real P-PET
pmeanmPETMeans <- pmeanmPETMeans - pmPetMeans
pmPETmeanMeans <- pmPETmeanMeans - pmPetMeans

# loop through and average by year
averageEveryP <- matrix(nrow = nSites, ncol = nMethods)
averageEveryPET <- matrix(nrow = nSites, ncol = nMethods)
for (i in 1:nSites){
  averageEveryP[i,] <- rowMeans(pmeanmPETMeans[i,,], na.rm = TRUE)
  averageEveryPET[i,] <- rowMeans(pmeanmPETMeans[i,,], na.rm = TRUE)
}