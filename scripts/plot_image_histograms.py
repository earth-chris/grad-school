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
           '/home/cba/cba/costa_rica/scratch/CR_forest_cover_proximity.tif',
           '/home/cba/cba/costa_rica/scratch/Costa_Rica_inverted_forest_cover_proximity.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Indigenous_Areas_distance.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Protected_Areas_distance.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_Protected_and_Indigenous_Areas_distance.tif',
           '/home/cba/cba/costa_rica/Costa_Rica_TOA_refl_median_frac.tif']
           
band_title = ['Biomass (Mg C)',
         'Total solar insolation (w / m^2)',
         'Inter Annual Cloudiness',
         'Intra Annual Cloudiness',
         'Mean Annual Cloudiness',
         'Aspect (degrees)',
         'Slope (degrees)',
         'Elevation (m ASL)',
         'Distance to Forest (m)',
         'Distance from Forest (m)',
         'Distance to Indigenous Area (m)',
         'Distance to Protected Area (m)',
         'Distance to Protected or Indigenous Area (m)',
         'Bare soil (%)',
         'Photosynthetic Vegetation (%)',
         'Non-photosynthetic Vegetation (%)',
         'STDEV Bare (%)',
         'STDEV PV (%)',
         'STDEV NPV (%)',
         'Bare-PV-NPV RMSE (%)']

# set output directory for files
outdir = '/home/cba/cba/scratch/'
           
for i in range(len(infiles)):
    
    # open file
    file = gdal.Open(infiles[i])
    
    # loop through each band
    for j in range(file.RasterCount):
        
        # set output title
        title = os.path.basename(infiles[i])[:-4] + ' band %s' % (j+1)
        
        # get band object
        band = file.GetRasterBand(j + 1)
        
        # get no-data
        nd = band.GetNoDataValue()
        
        # if it doesn't exist, use the min value in the band for nd
        if nd is None:
            nd = band.ComputeRasterMinMax[0]
        
        # read band into memory
        arr = band.ReadAsArray()
        
        # find indices for good data
        gd = np.where(arr != nd)
        
        # subset to 1-D array
        arr = arr[gd[0],gd[1]]
        
        # get min/max for plotting
        xmin = arr.min()
        xmax = arr.max()
        
        # plot the data
        plt.hist(arr, 50, facecolor='green', alpha=0.75)
        plt.title(band_title[i+j])
        plt.ylabel('Frequency')
        plt.xlim([xmin,xmax])
        plt.grid(True)
        plt.show()
        
        # export the data
        plt.savefig(outdir + title + '.png', format='png')
        
        # close the plot
        plt.clf()
        plt.close()
        
        # clear memory
        arr = None
        band = None
        
    file = None
        
