#!/usr/bin/python
#####
# plots distributions of growth EBV variables 
#####

import aei
import gdal
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc
from matplotlib.colors import ColorConverter as cc
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
cor_titles = ['Tree Cover\n(%)', 'NIRv\n(unitless)', 'Temperature\n(C)',
    'Solar Insolation\n(kWh / m2)', 'Annual Cloud Cover (%)']

# create color maps for each plot
tree_cmap = lsc.from_list('tree_cmap', [(0, cc.to_rgba('#cc78a7')), (1, cc.to_rgba('#06de37'))])
nirv_cmap = lsc.from_list('nirv_cmap', [(0, cc.to_rgba('#e74726')), (1, cc.to_rgba('#5db3e5'))])
temp_cmap = lsc.from_list('temp_cmap', [(0, cc.to_rgba('#0773b3')), (1, cc.to_rgba('#f1e545'))])
insl_cmap = lsc.from_list('insl_cmap', [(0, cc.to_rgba('#0d0887')), 
    (0.5, cc.to_rgba('#da5b69')), (1, cc.to_rgba('#f0f921'))])
    
cmap = [tree_cmap, nirv_cmap, temp_cmap, insl_cmap]
    
# read the data reference
ref = gdal.Open(file)
ndval = 0
data = ref.ReadAsArray()

# get the no data locations
gd = np.where(data[0,:,:] != ndval)

# loop through each band and plot a density distribution
plt.figure(1)
counter = 0
for i in [3,0,1,2]:
    
    # subset to a 1d array
    band_data = data[i,gd[0], gd[1]]
    
    # find min/max for plot bounds
    xmin = np.percentile(band_data, 2)
    xmax = np.percentile(band_data, 98)
    
    # create ticks to plot
    xticks = np.around(np.arange(0,1.25, 0.25) * (xmax - xmin) + xmin, 2)
    
    # set custom covariance to smooth peaks
    covar = 0.5
    dns = gaussian_kde(band_data)
    if i == 2:
        dns.covariance_factor = lambda : covar
        dns._compute_covariance()
    
    # set xscale to plot
    npts = 100
    xs = np.linspace(xmin, xmax, npts)
    
    # create a normalizing function for color filling
    norm = matplotlib.colors.Normalize(vmin=xmin, vmax=xmax)
    
    # set a vertical structure for this plot format
    plt.subplot(411 + counter)
    
    # calculate the lines to plot
    ydata = dns(xs)
    
    # plot the function
    plt.plot(xs, dns(xs), color = 'black', lw = 3)
    plt.xlabel(titles[i])
    plt.yticks([], [])
    plt.xticks(xticks)
    plt.ylabel('%')
    
    # add the fill point by point
    for j in range(npts-1):
        plt.fill_between([xs[j], xs[j+1]], [ydata[j], ydata[j+1]], color = cmap[i](norm(xs[j])))
    
    counter += 1
    
# save the final figure
plt.tight_layout()
plt.savefig(ofile, dpi = 300)

# calculate correlation coefficients between variables
cor = np.corrcoef(data[0:4, gd[0], gd[1]])
arr = data[0:4, gd[0], gd[1]]

# create heat maps using hexbin to plot correlations between each variable
ctable = 'plasma'
#ctable = 'copper'
stanford_cmap = temp_cmap = lsc.from_list('stanford_cmap', [(0, cc.to_rgba("#000000")), (1, cc.to_rgba("#8C1515"))])

plt.figure(2)

# set up the plots to use
#plots = [[0,0], [0,1], [0,2], [0,3], [1,1], [1,2], [1,3], [2,2], [2,3], [3,3]]
plots = [[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]]
for plot in plots:
    #plt.subplot(441)
    #plt.subplot(4, 4, 1 + plot[0] + (plot[1] * 4))
    plt.subplot(3, 3, 1 + plot[0] + ((plot[1]-1) * 3))
    plt.hexbin(arr[plot[0]], arr[plot[1]], gridsize = 20, cmap = ctable,
        alpha = 0.95, linewidths = 0, bins = 'log')
    
    # label the correlation
    plt.title("pearson's r: {:0.2f}".format(cor[plot[0], plot[1]]))
        
    # label the axes
    if plot[0] == 0:
        plt.ylabel(cor_titles[plot[1]])
    if plot[1] == 3:
        plt.xlabel(cor_titles[plot[0]])
    