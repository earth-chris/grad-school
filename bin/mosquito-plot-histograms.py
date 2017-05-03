#!/usr/bin/python
#####
# plots distributions of predictor variables 
#####

import aei
import gdal
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# set base directory for files
base = '/home/cba/Downloads/mosquito-proposal/'

# file for CR-wide data
cr_file = base + 'CR-aligned.tif'
cr_ref = gdal.Open(cr_file)

# file for Coto Brus data
cb_file = base + 'CR-aligned-coto-brus.tif'
cb_ref = gdal.Open(cb_file)

# file for current sampling data
fs_cr_file = base + 'CR-plot-locations.tif'
fs_cb_file = base + 'CR-plot-locations-coto-brus.tif'

# set the colors and band names to plot
plt_bands = ["Temperature (C)", "Soil Cover (%)", "Veg. Cover (%)", "Impervious Cover (%)"]
bands = ["Temperature", "Soil-Cover", "Veg-Cover", "Impervious-Cover"]
output_files = []
for band in bands:
    output_files.append("{base}CR-{band}-plot.png".format(base = base, band = band))
cols = aei.color.color_blind(3)

# read plot sampling locations into memory to get raster indices
#fs_cr_ref = gdal.Open(fs_cr_file)
#fs_cr_band = fs_cr_ref.ReadAsArray()
#fs_cr = np.where(fs_cr_band == 1)
#fs_cr_band = None

fs_cb_ref = gdal.Open(fs_cb_file)
fs_cb_band = fs_cb_ref.ReadAsArray()
fs_cb = np.where(fs_cb_band == 1)
fs_cb_band = None

# count how many pixels in coto brus field sites to use as a limit on the
#  number of random samples to pull from the national- and coto brus-wide data
ns = fs_cb[0].shape[0]
nscr = 2e5

# loop through each band and plot histograms of nation-wide data,
#  coto brus-wide data, and the data from all plots
for i in range(cr_ref.RasterCount):
    cr_band = cr_ref.GetRasterBand(i+1)
    cb_band = cb_ref.GetRasterBand(i+1)
    cr_nodata = cr_band.GetNoDataValue()
    cb_nodata = cb_band.GetNoDataValue()
    
    # read the national data into memory and extract the values to plot
    cr_arr = cr_band.ReadAsArray()
    
    # if this is the first band, select a random subset of pixels to sample and use that throughout
    if i == 0:
        cr_ind = np.where(cr_arr != cr_nodata)
        cr_rnd = np.int16(np.floor(np.random.uniform(0, cr_ind[0].shape[0], ns)))
    
    # subset the data and clear the rest out of memory    
    cr_data = cr_arr[cr_ind[0][cr_rnd], cr_ind[1][cr_rnd]]
    cr_arr = None
    
    # read the coto brus data into memory and extract the values to plot
    cb_arr = cb_band.ReadAsArray()
    
    # pull the field plot data
    fs_data = cb_arr[fs_cb[0], fs_cb[1]]
    
    if i == 0:
        cb_ind = np.where(cb_arr != cb_nodata)
        cb_rnd = np.int16(np.floor(np.random.uniform(0, cb_ind[0].shape[0], ns)))
        
    cb_data = cb_arr[cb_ind[0][cb_rnd], cb_ind[1][cb_rnd]]
    cb_arr = None
    
    # find the min and max for each data set to set bounds for plot
    xmin = min(cr_data.min(), cb_data.min(), fs_data.min())
    xmax = max(cr_data.max(), cb_data.max(), fs_data.max())
    
    # get a density distribution for each 
    cr_dns = gaussian_kde(cr_data)
    cb_dns = gaussian_kde(cb_data)
    fs_dns = gaussian_kde(fs_data)
    
    # set custom covariance to smooth multiple peaks
    covar = 0.25
    cr_dns.covariance_factor = lambda : covar
    cb_dns.covariance_factor = lambda : covar
    fs_dns.covariance_factor = lambda : covar
    cr_dns._compute_covariance()
    cb_dns._compute_covariance()
    fs_dns._compute_covariance()
    
    # set an x scale to plot smooth density distribution
    xs = np.linspace(xmin, xmax, 200)
    
    # plot each function in a single plot
    plt.figure()
    plt.plot(xs, cr_dns(xs), label = "Costa Rica", color = cols[0])
    plt.plot(xs, cb_dns(xs), label = "Coto Brus", color = cols[1])
    plt.plot(xs, fs_dns(xs), label = "Field Plots", color = cols[2])
    plt.xlabel(plt_bands[i])
    plt.ylabel("Density")
    plt.title("Costa Rica {var} Distributions".format(bands[i]))
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_files[i], dpi = 200)
    plt.close()

# that's all, folks!
