# 
import aei
import gdal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

%matplotlib tk

# set the working directories
base = '/home/cba/Downloads/scale-conceptual/'
plots = base + 'plots/'
tch_file = base + 'marin_tch_byte.tif'
gd_file = base + 'marin_tch_mask.tif'

# set the number of resolutions to assess
res = [20, 30, 50, 100, 250, 500, 1000]

# open the suitability raster and read in good indices
gdref = gdal.Open(gd_file)
gd = gdref.ReadAsArray()
gd_ext1 = gd == 1
bd = gd == 0
#gd_ext2 = np.where(gd > 1)
gd = None
gdref = None

# read the tch data into memory 
tchref = gdal.Open(tch_file)
#tchb1 = tchref.GetRasterBand(1)
#ndval = tchb1.GetNoDataValue()
tch = tchref.ReadAsArray()
tchb1 = None
tchref = None

# set the nodata values to nan
tch[bd] = np.nan

# set the number of pixels to sample 
fullres = 2 * 2
testres = min(res) * min(res)
n_var = testres / fullres
n_loc = n_var
n_rep = 1000

# create a set of arrays to store the data
wg_var = np.zeros((len(res), n_rep, n_loc))
bg_var = np.zeros((len(res), n_rep, n_loc))

# loop through each resolution, and calculate within/between grain variance
for i in range(len(res)):
    rdif = res[i] / 2
    for j in range(n_rep)):
        
        # find the random pixel locations to select
        loc_ext1 = np.random.choice(len(gd_ext1[0]), n_loc)
        #loc_ext2 = np.random.choice(len(gd_ext2[0]), n_loc)
        
        # loop through each random location
        for k in range(n_loc):
            # sample the x/y space around this pixel
            xmin = loc_ext1[k]- rdif
            xmax = loc_ext1[k] + rdif
            ymin = loc_ext1[k] - rdif
            ymax = loc_ext1[k] + rdif
            subset = tch[xmin:xmax, ymin:ymax]
            avg_ext1 = np.nanmean(subset)
            var_ext1 = np.nanvar(subset)
            wg_var[i,j,k] = var_ext1
            bg_var[i,j,k] = avg_ext1
            
            # then do it for the second extent
            #xmin = gd[0][loc_ext2[k]] - rdif
            #xmax = gd[0][loc_ext2[k]] + rdif
            #ymin = gd[1][loc_ext2[k]] - rdif
            #ymax = gd[1][loc_ext2[k]] + rdif
            #subset = tch[xmin:xmax, ymin:ymax]
            #avg_ext2 = np.nanmean(subset)
            #var_ext2 = np.nanvar(subset)