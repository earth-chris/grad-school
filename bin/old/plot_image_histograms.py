#####
# plots histograms of raster data
#####
import os as os
import aei as aei
import numpy as np
import gdal as gdal
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

infiles = ['/home/cba/cba/costa_rica/Costa_Rica_Biomass_Baccini.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_insol_total_winter.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Inter_Annual_Cloud_SD.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Intra_Annual_Cloud_SD.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Mean_Annual_Cloud.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_SRTM_aspect.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_SRTM_slope.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_SRTM.tif',
           '/home/cba/cba/costa_rica/scratch/CR_forest_cover_proximity_resized.tif',
           '/home/cba/cba/costa_rica/scratch/Costa_Rica_inverted_forest_cover_proximity_resized.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Indigenous_Areas_distance.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Protected_Areas_distance.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Protected_and_Indigenous_Areas_distance.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_TOA_refl_median_frac.tif']
           
band_title = ['Biomass',
         'Total solar insolation',
         'Inter Annual Cloudiness',
         'Intra Annual Cloudiness',
         'Mean Annual Cloudiness',
         'Aspect',
         'Slope',
         'Elevation',
         'Distance to Forest',
         'Distance from Forest',
         'Distance to Indigenous Area',
         'Distance to Protected Area',
         'Distance to Protected or Indigenous Area',
         'Bare soil',
         'Photosynthetic Vegetation',
         'Non-photosynthetic Vegetation',
         'STDEV Bare',
         'STDEV PV',
         'STDEV NPV',
         'Bare-PV-NPV RMSE']

units = ['Mg C',
         'w / m^2',
         '% * 10000',
         '% * 10000',
         '% * 10000',
         'Degrees',
         'Degrees',
         'Meters ASL',
         'Meters',
         'Meters',
         'Meters',
         'Meters',
         'Meters',
         '%',
         '%',
         '%',
         '%',
         '%',
         '%',
         '%']

# set output directory for files
outdir = '/home/cba/cba/scratch/'

# set the title to append
append = 'All Plots 150m Buffer'

# if you only want to review data from certain region
#  (e.g. coto-brus, around field plots, etc.)
#  set an input mask file to use and set flag to True
subset = True
#subset_file = '/home/cba/cba/costa_rica/scratch/CotoBrus_merged_average_30m_resized.tif'
subset_file = '/home/cba/cba/costa_rica/scratch/Costa_Rica_All_Plots_150m_buffer.tif'
#subset_file = '/home/cba/cba/costa_rica/CotoBrus_mask.tif'

if subset:
    subset_file_ref = gdal.Open(subset_file)
    subset_band = subset_file_ref.GetRasterBand(1)
    subset_arr = subset_band.ReadAsArray()
    subset_bd = np.where(subset_arr <= 0.0)
    
    # clear memory
    subset_file_ref = None
    subset_band = None
    subset_arr = None
           
for i in range(len(infiles)):
    
    # open file
    file = gdal.Open(infiles[i])
    
    # loop through each band
    for j in range(file.RasterCount):
        
        # set output title
        title = os.path.basename(infiles[i])[:-4] + ' band %s ' % (j+1)
        
        # get band object
        band = file.GetRasterBand(j + 1)
        
        # get no-data
        nd = band.GetNoDataValue()
        
        # if it doesn't exist, use the min value in the band for nd
        if nd is None:
            nd = band.ComputeRasterMinMax[0]
        
        # read band into memory
        arr = band.ReadAsArray()
        
        # set areas outside of subset to no-data
        if subset:
            arr[subset_bd[0],subset_bd[1]] = nd
        
        # find indices for good data
        gd = np.where(arr != nd)
        
        # subset to 1-D array
        arr = arr[gd[0],gd[1]]
        
        # get min/max for plotting
        xmin = arr.min()
        xmax = arr.max()
        
        # plot the data
        plt.hist(arr, 50, facecolor='green', alpha=0.75)
        plt.title(append + ' ' + band_title[i+j])
        plt.ylabel('Frequency')
        plt.xlabel(units[i+j])
        plt.xlim([xmin,xmax])
        plt.grid(True)
        plt.show()
        
        # export the data
        plt.savefig(outdir + title + append + '.png', format='png')
        
        # close the plot
        plt.clf()
        plt.close()
        
        # clear memory
        arr = None
        band = None
        
    file = None
        
