###3#
# function to create a ternary plot of remote sensing data sources for ecological modeling
#
# started 3/7/17 cba
#####

library(ggtern)

# load the csv file 
infile <- "~/src/aei-grad-school/scratch/satellite_scales.csv"
satData <- read.csv(file = infile, header = TRUE)

# transform the data to get a good distribution of values
satData["Frequency"] = 1/satData$Frequency
satData["sqrtFreq"] = sqrt(satData$Frequency)
satData["logExtent"] = log10(satData$Extent)
satData["lnExtent"] = log(satData$Extent)
satData["logScale"] = log10(satData$UseScale)
satData['lnScale'] = log(satData$UseScale)
satData["logFreq"] = log10(satData$Frequency)
satData['lnFreq'] = log(satData$Frequency)
satData["logSwath"] = log10(satData$Swath)
satData['lnSwath'] = log(satData$Swath)

# filter for data
gd <- which(satData$UseScale < 1e6)

# set up plot formats
xlab <- "Sampling Unit (m^2)"
xdata <- satData$logScale

ylab <- "Extent (km^2)"
ydata <- satData$logExtent

zlab <- "Frequency (1/days)"
zdata <- satData$Frequency

# set the data used for plot shapes
shapeValues = seq(length(unique(satData$EBV.Class))) + 20
EBV.Class <- satData$EBV.Class

# set the data used for plot colors
cbBase <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- colorRampPalette(cbbBase)(length(unique(satData$EBV.Class)))

title <- "Satellite Scales"

# set up the basic plot function
base = ggtern(data = satData[gd,], aes(logScale, EBVScale, sqrtFreq))
base = base + geom_point(aes(fill = EBV.Class, shape = EBV.Class, size = 1)) +
  scale_shape_manual(values = shapeValues) + 
  scale_fill_manual(values = cbBase) +
  theme_nomask() +
  labs(x = xlab, y = ylab, z = zlab) + 
  theme_showarrows() +
  scale_x_log10()# + scale_y_log10()
base = base + theme_notitles() + theme_nolabels()
base

base = ggtern(data = satData[ad,], aes(logScale, EBVScale, sqrtFreq))
base = base + geom_point(aes(color = Sensor.Class), size = 10) +
  #geom_label_repel(aes(Frequency, UseScale, label = Satellite), segment.color = "black", 
  #                 fontface = 'bold', color = 'black', box.padding = unit(0.4, "lines"),
  #                 point.padding = unit(0.7, "lines"), alpha = 0.8) +
  labs(x = xlab, y = ylab, z = zlab) + 
  theme_showarrows() +
  scale_color_manual(values = cbBase) +
  theme_notitles() + theme_nolabels()
base
