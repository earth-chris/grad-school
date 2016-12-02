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
monthStart <- 1
monthEnd <- monthStart
dayStart <- 1
dayEnd <- dayStart
yearStart <- 2006
yearEnd <- 2016
years <- yearStart:yearEnd
nYears <- yearEnd - yearStart

# set the temp and precip change values
tempLow <- 1
tempMean <- 2
tempHigh <- 3
precipLow <- -0.1
precipMean <- 0.0
precipHigh <- 0.1

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
pmPetMeans <- array(dim=c(nSites, nMethods, nYears))
pmPetStdvs <- array(dim=c(nSites, nMethods, nYears))
pmeanmPETMeans <- array(dim=c(nSites, nMethods, nYears))
pmPETmeanMeans <- array(dim=c(nSites, nMethods, nYears))
lowTChangePET <- array(dim=c(nSites, nYears))
lowTChangePmPET <- array(dim=c(nSites, nYears))
highTChangePET <- array(dim=c(nSites, nYears))
highTChangePmPET <- array(dim=c(nSites, nYears))
lowPChangePET <- array(dim=c(nSites, nYears))
lowPChangePmPET <- array(dim=c(nSites, nYears))
highPChangePET <- array(dim=c(nSites, nYears))
highPChangePmPET <- array(dim=c(nSites, nYears))

#meanPrecip <- array(dim=c(nSites, nYears))

for (i in seq(1, nSites)){
  # group by year an add new columns for the PET calculations
  petGroup <- group_by(filter(siteList[[i]], inSeason==TRUE), WEYR)
  petGroup <- mutate(petGroup, PriestlyTaylor = pet.priestlyTaylor(WEDAY, TMAX, TMIN, TMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     ModifiedPriestlyTaylor = pet.modifiedPriestlyTaylor(TMAX, TMIN, SRAD),
                     Hamon = pet.hammon(WEDAY, TMAX, TMIN, TMEAN, LAT),
                     Hargreaves = pet.hargreaves(WEDAY, TMAX, TMIN, TMEAN, LAT),
                     Linacre = pet.linacre(TMEAN, ELEV, LAT, TDEW),
                     Turc = pet.turc(TMEAN, RH2M, SRAD),
                     meanPrecip = mean(RAIN),
                     lowTMAX = TMAX + tempLow,
                     lowTMIN = TMIN + tempLow,
                     lowTMEAN = TMEAN + tempLow,
                     meanTMAX = TMAX + tempMean,
                     meanTMIN = TMIN + tempMean,
                     meanTMEAN = TMEAN + tempMean,
                     maxTMAX = TMAX + tempHigh,
                     maxTMIN = TMIN + tempHigh,
                     maxTMEAN = TMEAN + tempMean,
                     lowRAIN = RAIN + (RAIN * precipLow),
                     meanRAIN = RAIN + (RAIN * precipMean),
                     highRAIN = RAIN + (RAIN * precipHigh))
  
  # report the summary statistics for PET variables
  petSummary <- summarize(petGroup,
    mean.PriestlyTaylor = mean(PriestlyTaylor, na.rm=TRUE),
    mean.ModifiedPriestlyTaylor = mean(ModifiedPriestlyTaylor, na.rm=TRUE),
    mean.Hamon = mean(Hamon, na.rm=TRUE),
    mean.Hargreaves = mean(Hargreaves, na.rm=TRUE),
    mean.Linacre = mean(Linacre, na.rm=TRUE),
    mean.Turc = mean(Turc, na.rm=TRUE),
    mean.Precip = mean(RAIN),
    sd.PriestlyTaylor = sd(PriestlyTaylor, na.rm=TRUE),
    sd.ModifiedPriestlyTaylor = sd(ModifiedPriestlyTaylor, na.rm=TRUE),
    sd.Hamon = sd(Hamon, na.rm=TRUE),
    sd.Hargreaves = sd(Hargreaves, na.rm=TRUE),
    sd.Linacre = sd(Linacre, na.rm=TRUE),
    sd.Turc = sd(Turc, na.rm=TRUE))
  
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
  #meanPrecip[i,] <- petSummary$mean.Precip
  
  # create new columns to analyze in problems 4 and 5
  petGroup <- mutate(petGroup, MeanPriestlyTaylor = mean(PriestlyTaylor, na.rm = TRUE),
                     MeanModifiedPriestlyTaylor = mean(ModifiedPriestlyTaylor, na.rm = TRUE),
                     MeanHamon = mean(Hamon, na.rm = TRUE),
                     MeanHargreaves = mean(Hargreaves, na.rm = TRUE),
                     MeanLinacre = mean(Linacre, na.rm = TRUE),
                     MeanTurc = mean(Turc, na.rm = TRUE),
                     pmPriestlyTaylor = RAIN - PriestlyTaylor,
                     pmModifiedPriestlyTaylor = RAIN - ModifiedPriestlyTaylor,
                     pmHamon = RAIN - Hamon,
                     pmHargreaves = RAIN - Hargreaves,
                     pmLinacre = RAIN - Linacre,
                     pmTurc = RAIN - Turc)
  
  # report the summary statistics for P-PET variables
  petSummary <- summarize(petGroup,
                          mean.pmPriestlyTaylor = mean(pmPriestlyTaylor, na.rm=TRUE),
                          mean.pmModifiedPriestlyTaylor = mean(pmModifiedPriestlyTaylor, na.rm=TRUE),
                          mean.pmHamon = mean(pmHamon, na.rm=TRUE),
                          mean.pmHargreaves = mean(pmHargreaves, na.rm=TRUE),
                          mean.pmLinacre = mean(pmLinacre, na.rm=TRUE),
                          mean.pmTurc = mean(pmTurc, na.rm=TRUE),
                          sd.pmPriestlyTaylor = sd(pmPriestlyTaylor, na.rm=TRUE),
                          sd.pmModifiedPriestlyTaylor = sd(pmModifiedPriestlyTaylor, na.rm=TRUE),
                          sd.pmHamon = sd(pmHamon, na.rm=TRUE),
                          sd.pmHargreaves = sd(pmHargreaves, na.rm=TRUE),
                          sd.pmLinacre = sd(pmLinacre, na.rm=TRUE),
                          sd.pmTurc = sd(pmTurc, na.rm=TRUE))
  
  # assign this absurd stuff to an output matrix
  pmPetMeans[i, 1,] <- petSummary$mean.pmPriestlyTaylor
  pmPetMeans[i, 2,] <- petSummary$mean.pmModifiedPriestlyTaylor
  pmPetMeans[i, 3,] <- petSummary$mean.pmHamon
  pmPetMeans[i, 4,] <- petSummary$mean.pmHargreaves
  pmPetMeans[i, 5,] <- petSummary$mean.pmLinacre
  pmPetMeans[i, 6,] <- petSummary$mean.pmTurc
  pmPetStdvs[i, 1,] <- petSummary$sd.pmPriestlyTaylor
  pmPetStdvs[i, 2,] <- petSummary$sd.pmModifiedPriestlyTaylor
  pmPetStdvs[i, 3,] <- petSummary$sd.pmHamon
  pmPetStdvs[i, 4,] <- petSummary$sd.pmHargreaves
  pmPetStdvs[i, 5,] <- petSummary$sd.pmLinacre
  pmPetStdvs[i, 6,] <- petSummary$sd.pmTurc
  
  # then go through another round to add more columns to asses for problem 5
  petGroup <- mutate(petGroup, pmMeanPriestlyTaylor = RAIN - MeanPriestlyTaylor,
                     pmMeanModifiedPriestlyTaylor = RAIN - MeanModifiedPriestlyTaylor,
                     pmMeanHamon = RAIN - MeanHamon,
                     pmMeanHargreaves = RAIN - MeanHargreaves,
                     pmMeanLinacre = RAIN - MeanLinacre,
                     pmMeanTurc = RAIN - MeanTurc,
                     meanPmPriestlyTaylor = meanPrecip - PriestlyTaylor,
                     meanPmModifiedPriestlyTaylor = meanPrecip - ModifiedPriestlyTaylor,
                     meanPmHamon = meanPrecip - Hamon,
                     meanPmHargreaves = meanPrecip - Hargreaves,
                     meanPmLinacre = meanPrecip - Linacre,
                     meanPmTurc = meanPrecip - Turc)
  
  # report the summary statistics for the P(yearlyMean) - PET(everyValue) and P(everyValue) - PET(yearlyMean) variables
  petSummary <- summarize(petGroup,
                          mean.pmMeanPriestlyTaylor = mean(pmMeanPriestlyTaylor, na.rm=TRUE),
                          mean.pmMeanModifiedPriestlyTaylor = mean(pmMeanModifiedPriestlyTaylor, na.rm=TRUE),
                          mean.pmMeanHamon = mean(pmMeanHamon, na.rm=TRUE),
                          mean.pmMeanHargreaves = mean(pmMeanHargreaves, na.rm=TRUE),
                          mean.pmMeanLinacre = mean(pmMeanLinacre, na.rm=TRUE),
                          mean.pmMeanTurc = mean(pmMeanTurc, na.rm=TRUE),
                          mean.meanPmPriestlyTaylor = mean(meanPmPriestlyTaylor, na.rm=TRUE),
                          mean.meanPmModifiedPriestlyTaylor = mean(meanPmModifiedPriestlyTaylor, na.rm=TRUE),
                          mean.meanPmHamon = mean(meanPmHamon, na.rm=TRUE),
                          mean.meanPmHargreaves = mean(meanPmHargreaves, na.rm=TRUE),
                          mean.meanPmLinacre = mean(meanPmLinacre, na.rm=TRUE),
                          mean.meanPmTurc = mean(meanPmTurc, na.rm=TRUE))
  
  # assign this absurd stuff to an output matrix
  pmeanmPETMeans[i, 1,] <- petSummary$mean.meanPmPriestlyTaylor
  pmeanmPETMeans[i, 2,] <- petSummary$mean.meanPmModifiedPriestlyTaylor
  pmeanmPETMeans[i, 3,] <- petSummary$mean.meanPmHamon
  pmeanmPETMeans[i, 4,] <- petSummary$mean.meanPmHargreaves
  pmeanmPETMeans[i, 5,] <- petSummary$mean.meanPmLinacre
  pmeanmPETMeans[i, 6,] <- petSummary$mean.meanPmTurc
  pmPETmeanMeans[i, 1,] <- petSummary$mean.pmMeanPriestlyTaylor
  pmPETmeanMeans[i, 2,] <- petSummary$mean.pmMeanModifiedPriestlyTaylor
  pmPETmeanMeans[i, 3,] <- petSummary$mean.pmMeanHamon
  pmPETmeanMeans[i, 4,] <- petSummary$mean.pmMeanHargreaves
  pmPETmeanMeans[i, 5,] <- petSummary$mean.pmMeanLinacre
  pmPETmeanMeans[i, 6,] <- petSummary$mean.pmMeanTurc
  
  # now mutate for problem 6 climate modeling
  petGroup <- mutate(petGroup, lowTChangePET = pet.priestlyTaylor(WEDAY, lowTMAX, lowTMIN, lowTMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     lowTChangePmPET = meanRAIN - pet.priestlyTaylor(WEDAY, lowTMAX, lowTMIN, lowTMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     highTChangePET = pet.priestlyTaylor(WEDAY, maxTMAX, maxTMIN, maxTMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     highTChangePmPET = meanRAIN - pet.priestlyTaylor(WEDAY, maxTMAX, maxTMIN, maxTMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     lowPChangePET = pet.priestlyTaylor(WEDAY, meanTMAX, meanTMIN, meanTMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     lowPChangePmPET = lowRAIN - pet.priestlyTaylor(WEDAY, meanTMAX, meanTMIN, meanTMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     highPChangePET = pet.priestlyTaylor(WEDAY, meanTMAX, meanTMIN, meanTMEAN, RH2M, TDEW, SRAD, ELEV, LAT),
                     highPChangePmPET = highRAIN - pet.priestlyTaylor(WEDAY, meanTMAX, meanTMIN, meanTMEAN, RH2M, TDEW, SRAD, ELEV, LAT))
  
  # summarize the changes in climate
  petSummary <- summarize(petGroup,
                          mean.lowTChangePET = mean(lowTChangePET, na.rm = TRUE),
                          mean.lowTChangePmPET = mean(lowTChangePmPET, na.rm = TRUE),
                          mean.highTChangePET = mean(highTChangePET, na.rm = TRUE),
                          mean.highTChangePmPET = mean(highTChangePmPET, na.rm = TRUE),
                          mean.lowPChangePET = mean(lowPChangePET, na.rm = TRUE),
                          mean.lowPChangePmPET = mean(lowPChangePmPET, na.rm = TRUE),
                          mean.highPChangePET = mean(highPChangePET, na.rm = TRUE),
                          mean.highPChangePmPET = mean(highPChangePmPET, na.rm = TRUE))
  
  # and assign to output variables
  lowTChangePET[i,] <- petSummary$mean.lowTChangePET
  lowTChangePmPET[i,] <- petSummary$mean.lowTChangePmPET
  highTChangePET[i,] <- petSummary$mean.highTChangePET
  highTChangePmPET[i,] <- petSummary$mean.highTChangePmPET
  lowPChangePET[i,] <- petSummary$mean.lowPChangePET
  lowPChangePmPET[i,] <- petSummary$mean.lowPChangePmPET
  highPChangePET[i,] <- petSummary$mean.highPChangePET
  highPChangePmPET[i,] <- petSummary$mean.highPChangePmPET
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

# just kidding, dplyr doesn't work how I think it does. cutting bait now to work on problem 6

#############################
# problem 6 - climate modeling

lowTChangePET[i,] <- petSummary$mean.lowTChangePET
lowTChangePmPET[i,] <- petSummary$mean.lowTChangePmPET
highTChangePET[i,] <- petSummary$mean.highTChangePET
highTChangePmPET[i,] <- petSummary$mean.highTChangePmPET
lowPChangePET[i,] <- petSummary$mean.lowPChangePET
lowPChangePmPET[i,] <- petSummary$mean.lowPChangePmPET
highPChangePET[i,] <- petSummary$mean.highPChangePET
highPChangePmPET[i,] <- petSummary$mean.highPChangePmPET

# set up the plot utility
par(mfrow = c(2,4), oma = c(0,0,2,0))

# set up labels
title <- "Forecasting Effects of Climate Change on PET, P - PET"
ylabPmPET <- "P - PET"
ylabPET <- "PET"

# plot the changes
barplot2(rowMeans(lowTChangePET), names.arg = siteNames, las = 2, ylab = ylabPET, col = siteColors, main = "Low T Change")
barplot2(rowMeans(highTChangePET), names.arg = siteNames, las = 2, ylab = ylabPET, col = siteColors, main = "High T Change")
barplot2(rowMeans(lowPChangePET), names.arg = siteNames, las = 2, ylab = ylabPET, col = siteColors, main = "Low P Change")
barplot2(rowMeans(highPChangePET), names.arg = siteNames, las = 2, ylab = ylabPET, col = siteColors, main = "High P Change")
barplot2(rowMeans(lowTChangePmPET), names.arg = siteNames, las = 2, ylab = ylabPmPET, col = siteColors, main = "Low T Change")
barplot2(rowMeans(highTChangePmPET), names.arg = siteNames, las = 2, ylab = ylabPmPET, col = siteColors, main = "High T Change")
barplot2(rowMeans(lowPChangePmPET), names.arg = siteNames, las = 2, ylab = ylabPmPET, col = siteColors, main = "Low P Change")
barplot2(rowMeans(highPChangePmPET), names.arg = siteNames, las = 2, ylab = ylabPmPET, col = siteColors, main = "High P Change")

# add title to the final plot
title(title, outer = TRUE)
