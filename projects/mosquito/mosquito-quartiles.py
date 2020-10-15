#!/usr/bin/python
#####
# splits the predictor data sets into quartiles and maps the unique combinations
#####

import aei
import gdal
import numpy as np

# set base directory for files
base = '/home/cba/Downloads/mosquito/'

# file for CR-wide data
#cr_file = base + 'CCB-LC-Predictors-masked.tif'
#cr_ref = gdal.Open(cr_file)

# file for Coto Brus data
cb_file = base + 'CCB-LC-Predictors-southern.tif'
cb_ref = gdal.Open(cb_file)
nx = cb_ref.RasterXSize
ny = cb_ref.RasterYSize
geo = cb_ref.GetGeoTransform()
prj = cb_ref.GetProjection()

# file for current sampling data
fs_cb_file = base + 'CR-plot-locations-southern.tif'
fs_cb_ref = gdal.Open(fs_cb_file)

# file for Aedes aegyptii points
ae_file = base + 'CR-mosquito-locations-50m.tif'
ae_ref = gdal.Open(ae_file)

# set the colors and band names to plot
plt_bands = ["Tree Cover (%)", "Soil Cover (%)", "Veg. Cover (%)", "Impervious Cover (%)", 
    "Min. Temp (C)", "Median Temp (C)", "Max. Temp (C)"]
bands = ["TreeCover", "Soil-Cover", "Veg-Cover", "Impervious-Cover", "Temp-Min", 
    "Temp-Median", "Temp-Max"]
    
# create a key for the number of factors
factor = []
for i in range(len(bands)):
    factor.append(1 * 10 ** i)
    
# create an output array to store the final unique values
uniq = np.zeros((ny, nx))

# set the percentiles to use
percentiles = [0, 25, 50, 75, 100]

# loop through each predictor, split into quartiles, and output into a new image
for i in range(len(bands)):
    print("Starting {} quartile mapping".format(bands[i]))
    
    # get the band reference
    bref = cb_ref.GetRasterBand(i+1)
    
    # create an output array for this band
    oarr = np.zeros((ny, nx))
    
    # read the data
    band = bref.ReadAsArray()
    
    # get only good data
    ndval = bref.GetNoDataValue()
    gd = np.where(band != ndval)
    
    # subset the data to good indices
    data = band[gd[0], gd[1]]
    
    # get the percentiles
    pvals = []
    for pct in percentiles:
        pvals.append(np.percentile(data, pct))
        
    # loop through and get the indices for where these distributions fall
    for j in range(len(percentiles)-1):
        ll = np.where(data >= pvals[j])
        ul = np.where(data < pvals[j + 1])
        ol = np.intersect1d(ll, ul)
        oarr[gd[0][ol], gd[1][ol]] = j + 1
        uniq[gd[0][ol], gd[1][ol]] += (j + 1) * factor[i]
        
    # write the output to a new file
    ofile = "{}CR-quartiles-{}.tif".format(base, bands[i])
    oref = gdal.GetDriverByName("GTiff").Create(ofile, nx, ny, 1, gdal.GDT_Byte)
    oref.SetGeoTransform(geo)
    oref.SetProjection(prj)
    oband = oref.GetRasterBand(1)
    oband.WriteArray(oarr)
    oband = None
    oarr = None
    oref = None
    
# write the full unique quartile band
ofile = "{}CR-quartiles-unique.tif".format(base)
oref = gdal.GetDriverByName("GTiff").Create(ofile, nx, ny, 1, gdal.GDT_Int32)
oref.SetGeoTransform(geo)
oref.SetProjection(prj)
oband = oref.GetRasterBand(1)
oband.WriteArray(uniq)
oband = None
uniq = None
oref = None