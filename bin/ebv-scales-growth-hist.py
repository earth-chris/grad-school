#!/usr/bin/python
#####
# plots distributions of growth EBV variables 
#####

import aei
import gdal
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# set base directory for files
base = '/home/cba/Downloads/scale-figures/'

# file to read
file = base + 'CR-GROWTH-stack4.tif'

# set the output file
ofile = base + 'CR-histograms.tif'

# band names
names = ['TreeCover', 'NIRv', 'Temperature', 'Insolation', 'CloudCover']
titles = ['Tree Cover (%)', 'NIRv (unitless)', 'Temperature (C)',
    'Solar Insolation (kWh / m2)', 'Annual Cloud Cover (%)']
    
# read the data reference
ref = gdal.Open(file)
ndval = 0
data = ref.ReadAsArray()

# get the no data locations
gd = np.where(data[0,:,:] != ndval)

# loop through each band and plot a density distribution
plt.figure(2)
counter = 0
for i in [3,0,1,2]:
    
    # subset to a 1d array
    band_data = data[i,gd[0], gd[1]]
    
    # find min/max for plot bounds
    xmin = np.percentile(band_data, 2)
    xmax = np.percentile(band_data, 98)
    
    # set custom covariance to smooth peaks
    covar = 0.15
    dns = gaussian_kde(band_data)
    dns.covariance_factor = lambda : covar
    
    # set xscale to plot
    xs = np.linspace(xmin, xmax, 100)
    
    # set a vertical structure
    plt.subplot(411 + counter)
    
    # plot the function
    plt.plot(xs, dns(xs), color = 'black')
    plt.xlabel(titles[i])
    plt.yticks([], [])
    plt.ylabel('%')
    
    counter += 1
    
# save the final figure
plt.tight_layout()
plt.savefig(ofile, dpi = 300)