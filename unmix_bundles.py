###
# comments
###
import pysptools.abundance_maps as abm
import spectral as spectral
import numpy as np
import gdal as gdal
import aei as aei

# set the input reflectance file
refl_file = '/home/cba/cba/costa_rica/scratch/LC80140542015350LGN00/raw/Test_TOA_envi'

# set the input spectral library
spec_file = '/home/cba/cba/aei-grad-school/data/tropical_atmosphere_landsat8_bundles_1.sli'

# and set up the base name for the output file (to have *_mean.tif,
#  *_stdev.tif, *_median.tif appended
out_file = '/home/cba/cba/costa_rica/scratch/LC80140542015350LGN00/raw/cba_unmix_testing'
out_ends = ['_mean.tif', '_median.tif', '_stdev.tif']

# stat number of endmembers of each type in spectral library
n_endmembers = 20
n_bundles = 3

# open the refl/spectral library files for reading and get params
refl_data = spectral.open_image(refl_file+'.hdr')
spec_data = spectral.envi.open(spec_file[:-4]+'.hdr', spec_file)

# create an output array for the results of the unmixing
#  we are going to output the mean, median and stdev of solutions
#  the shape of nl * ns * n_bundles
mean_arr = np.zeros((refl_data.nrows, refl_data.ncols, n_bundles))
medi_arr = np.zeros((refl_data.nrows, refl_data.ncols, n_bundles))
stdv_arr = np.zeros((refl_data.nrows, refl_data.ncols, n_bundles))

# set up unmixing class
unmixer = abm.FCLS()

# we'll read through each line in the file, and unmix each
# of the endmembers present
for j in range(refl_data.nrows):
    
    # report on the status every 200 lines
    if j in range(0,refl_data.nrows,200): 
        print("[ STATUS ]: Unmixing line %s") % (j+1)
    
    # create a temp array to 
    temp_arr = np.empty((n_endmembers, refl_data.ncols, n_bundles))
    
    # read the line
    line = refl_data.read_subimage([j], range(refl_data.ncols))
    
    # find the bad data to mask
    gd = np.where(line[0,:,0] != -0.1)[0]
    
    if len(gd) == 0:
        continue
    
    # loop through each endmember
    for i in range(n_endmembers):
        
        # set up n_endmember x n_band unmixing array
        # in the order of bare, pv, npv for viewing
        bundles = np.array([spec_data.spectra[(i+1)+(n_endmembers*2)], 
            spec_data.spectra[i+1], spec_data.spectra[(i+1)+n_endmembers]])
        
        # perform the fully constrained least squares unmixing    
        temp_arr[i, gd] = unmixer.map(line[:,gd], bundles).squeeze()
        
    # collapse the unmixed data to mean, median and mode per pixel
    mean_arr[j] = temp_arr.mean(0)
    medi_arr[j] = np.median(temp_arr, axis=0)
    stdv_arr[j] = temp_arr.std(0)
    
# now that we've finished calculating our arrays, let's output them to
#  nice, fun geo-tiffs
orig = gdal.Open(refl_file)
proj = orig.GetProjection()
geot = orig.GetGeoTransform()

for i in range(3):
    # create the output file
    outRaster = gdal.GetDriverByName("GTiff").Create(
        out_file+out_ends[i], refl_data.ncols, refl_data.nrows, n_bundles,
        gdal.GDT_Float32)
    
    outRaster.SetGeoTransform(geot)
    outRaster.SetProjection(proj)
    
    # determine output array by loop
    if i == 0:
        out_arr = mean_arr
    elif i == 1:
        out_arr = medi_arr
    elif i == 2:
        out_arr = stdv_arr
        
    # loop through each band and write out the array
    for j in range(3):
        outBand = outRaster.GetRasterBand(j+1)
        outBand.WriteArray(out_arr[:,:,j])
        outBand.FlushCache()