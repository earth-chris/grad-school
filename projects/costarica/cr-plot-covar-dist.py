# script to plot the density distributions of covariates at the field sites and for costa rica
import aei
import gdal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

%matplotlib tk

# set the paths to use
base = '/home/salo/Downloads/costa-rica/'
plots = base + 'plots/'
sp_file = base + 'sp-data/sites-covariates.csv'

# set the raster paths
biomass = base + 'raster/biomass.tif'
fcover = base + 'raster/fcover.tif'
radar = base + 'raster/radar.tif'
tassledcap = base + 'raster/tassled-cap.tif'
temperature = base + 'raster/temperature.tif'
treecover = base + 'raster/tree-cover.tif'

# read the species data into pandas
sp = pd.read_csv(sp_file)

# create a list of raster files, bands, and names to iterate through
rasters = [biomass, fcover, fcover, fcover, radar, radar, 
           tassledcap, tassledcap, tassledcap, temperature, 
           temperature, temperature, treecover]
bands = [0,0,1,2,0,1,0,1,2,0,1,2,0]
names = ['Biomass', 'SoilPct', 'VegPct', 'Impervious', 'HH', 'HV',
        'tcBrt', 'tcGrn', 'tcWet', 'TMin', 'TMed', 'TMax', 'TreeCover']
pnames = ['Biomass', 'Soil cover', 'Vegetation cover', 'Impervious cover', 'Radar-HH', 
          'Radar-HV','Tassled cap brightness', 'Tassled cap greenness', 
          'Tassled cap wetness', 'Min. annual temp', 'Median annual temp', 
          'Max. annual temp', 'Tree cover']
xunit = ['Mg C/ha', '%', '%', '%', 'dB', 'dB', 'unitless', 
         'unitless', 'unitless', 'C', 'C', 'C', '%']

# loop through each covariate and plot background vs field plots
for i in range(len(names)):

    print('[ STATUS ]: Plotting {} data'.format(names[i]))
    
    # read the full raster data into memory
    tref = gdal.Open(rasters[i])
    bref = tref.GetRasterBand(bands[i]+1)
    ndval = bref.GetNoDataValue()
    band = bref.ReadAsArray()
    gd = band != ndval
    bkgr = band[gd]
    
    # clear memory
    band = None
    bref = None
    tref = None
    
    # to speed things up a bit, sample just 1 mil. background points
    n_rnd = int(1e6)
    rnd = np.random.choice(len(bkgr), n_rnd, replace=False)
    bkgr = bkgr[rnd]
    
    # subset the data from the field plots
    plts = np.array(sp[names[i]])
    
    # create the plot device
    plt.figure(figsize=(5,5), dpi=150)
    
    # set the colors to use
    colors = aei.color.color_blind()
    dns = aei.plot.density_dist([bkgr, plts], plot=plt, 
        label=['Background', 'Plots'],
        color=[colors[0], colors[5]],
        title='{}'.format(pnames[i]),
        xlabel = '({})'.format(xunit[i]),
        ylabel = '', fill_alpha = 0.4)
    
    # prep to save the figure    
    dns.tight_layout()
    dns.savefig(plots + '{}-density-dist.png'.format(names[i]))
    dns.close()
    