# be sure to install rgdal, raster, sp, and maptools packages

# missed early part of discussion since I was playing around with installing packages.

# raster package uses column major orientation for indexing, which is the opposite of how they deal with matrices.
#  ...wtf

# they have an xyFromCell and cellFromXY function

# extract() is the principal function used to get cell and layer values for any set of xy coordinates
library("rgdal")
library("raster")
library("sp") 
library("maptools")

# still having trouble loading libraries, and they aren't keeping functions up long enough to follow along

# will deal with namespace conflicts, e.g. the 'select' function from the raster package and from dplyr
#  can use things like raster::extract(yatta, blatta) to specify which namespace is necessary to use

# can use projection() function to learn the projection info for a raster or vector object