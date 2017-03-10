###3#
# function to create a ternary plot of remote sensing data sources for ecological modeling
#
# started 3/7/17 cba
#####

library(ggtern)

# load the csv file 
infile <- "~/cba/aei-grad-school/scratch/satellite_scales.csv"
satData <- read.csv(file = infile, header = TRUE)

# transform the data to get a good distribution of values
satData["Frequency"] = 1/satData$Frequency
satData["logExtent"] = log10(satData$Extent)
satData["lnExtent"] = log(satData$Extent)
satData["logScale"] = log10(satData$UseScale)
satData['lnScale'] = log(satData$UseScale)
satData["logFreq"] = log10(satData$Frequency)
satData['lnFreq'] = log(satData$Frequency)
satData["logSwath"] = log10(satData$Swath)
satData['lnSwath'] = log(satData$Swath)

# set some custom jiggering for overlapping plots
hyp = which(satData$Satellite == "Hyperion")
lst = which(satData$Satellite == "Landsat")
ast = which(satData$Satellite == "ASTER")
ali = which(satData$Satellite == "ALI")
pmx = which(satData$Satellite == "CBERS-PANMUX")
spt = which(satData$Satellite == "SPOT")

timeJitter <- 0.004
scaleJitter <- 150
satData$Frequency[hyp] = satData$Frequency[hyp] + timeJitter
satData$UseScale[hyp] = satData$UseScale[hyp] - scaleJitter
satData$Frequency[lst] = satData$Frequency[lst] - timeJitter
satData$UseScale[lst] = satData$UseScale[lst] - scaleJitter
satData$Frequency[ast] = satData$Frequency[ast] - timeJitter
satData$UseScale[ast] = satData$UseScale[ast] + scaleJitter
satData$Frequency[ali] = satData$Frequency[ali] + timeJitter
satData$UseScale[ali] = satData$UseScale[ali] + scaleJitter
satData$Frequency[pmx] = satData$Frequency[pmx] - timeJitter/2
satData$Frequency[spt] = satData$Frequency[spt] + timeJitter/2

# filter for data
gd <- which(satData$UseScale != 900 & satData$UseScale < 1e6 & satData$Satellite != "ResourceSat-1")
jd <- which(satData$UseScale == 900)
ad <- which(satData$UseScale < 1e5 & satData$Satellite != "ResourceSat-1")
ad <- which(satData$Satellite != "ResourceSat-1" & satData$Satellite != "ASCAT")

# set up plot formats
ylab <- "Spatial Resolution (m^2)"
ydata <- satData$logScale

xlab <- "Temporal Resolution (freq.)"
xdata <- satData$Frequency

# set the data used for plot shapes
shapeValues = seq(length(unique(satData$EBV.Class))) + 20
shapeValues = rep(21, length(unique(satData$EBV.Class)))
EBV.Class <- satData$EBV.Class

# set the data used for plot colors
cbBase <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- colorRampPalette(cbBase)(length(unique(satData$EBV)))
cbBase <- c("#CC79A7", "#0072B2", "#D55E00", "#009E73", "#56B4E9", "#F0E442", "#E69F00", "#000000")
cbPalette <- colorRampPalette(cbBase)(length(unique(satData$EBV)))

# set up breaks for log scale plotting
xbr = seq(0, max(satData$Frequency) * 1.5, by = 0.05)
ybr = seq(0, max(satData$UseScale) * 1.5, by = 100)

# set up the EBV Class plot function
title <- "Remote Sensing of Biodiversity Classes"
base = ggplot(data = satData[ad,], aes(Frequency, UseScale))
base = base + geom_point(aes(color = EBV.Class), size = 10) +
  geom_label_repel(aes(Frequency, UseScale, label = Satellite), segment.color = "black", 
                  fontface = 'bold', color = 'black', box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.7, "lines"), alpha = 0.8) +
  labs(x = xlab, y = ylab, title = title) +
  scale_y_log10(minor_breaks = ybr) +
  scale_x_log10(minor_breaks = xbr) +
  theme_showticks() + theme_linedraw() + theme_showsecondary() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cbBase)
base

# set up the EBV type function
title <- "Remote Sensing of Essential Biodiversity Variables"
base = ggplot(data = satData[ad,], aes(Frequency, UseScale))
base = base + geom_point(aes(color = EBV), size = 10) +
  geom_label_repel(aes(Frequency, UseScale, label = Satellite), segment.color = "black", 
                  fontface = 'bold', color = 'black', box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.7, "lines"), alpha = 0.8) +
  labs(x = xlab, y = ylab, title = title) +
  scale_y_log10(minor_breaks = ybr) +
  scale_x_log10(minor_breaks = xbr) +
  theme_showticks() + theme_linedraw() + theme_showsecondary() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cbPalette)
base

# set up the RS Class function
title <- "Sensor Class"
base = ggplot(data = satData[ad,], aes(Frequency, UseScale))
base = base + geom_point(aes(color = Sensor.Class), size = 10) +
  geom_label_repel(aes(Frequency, UseScale, label = Satellite), segment.color = "black", 
                  fontface = 'bold', color = 'black', box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.7, "lines"), alpha = 0.8) +
  labs(x = xlab, y = ylab, title = title) +
  scale_y_log10(minor_breaks = ybr) +
  scale_x_log10(minor_breaks = xbr) +
  theme_showticks() + theme_linedraw() + theme_showsecondary() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cbBase)
base

# set up the RS Class function
title <- "Sensor Type"
base = ggplot(data = satData[ad,], aes(Frequency, UseScale))
base = base + geom_point(aes(color = Sensor.Type), size = 10) +
  geom_label_repel(aes(Frequency, UseScale, label = Satellite), segment.color = "black", 
                  fontface = 'bold', color = 'black', box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.7, "lines"), alpha = 0.8) +
  labs(x = xlab, y = ylab, title = title) +
  scale_y_log10(minor_breaks = ybr) +
  scale_x_log10(minor_breaks = xbr) +
  theme_showticks() + theme_linedraw() + theme_showsecondary() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cbBase)
base