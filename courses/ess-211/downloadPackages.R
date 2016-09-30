#####
# this script is used to download and install packages to get around stupid stuff with my chromebook
#  cba 2016
#####

# set up the mirror site
mirror <- "http://cran.us.r-project.org"

# set up the library path
library <- "/home/cba/Downloads/source/R/"

# create a list for the packages to download
packages <- c("rgdal", "randomForest", "raster")

for (package in packages){
	install.packages(package, lib = library, repos = mirror)
}
