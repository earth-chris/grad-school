#!/usr/bin/python
#####
# assesses drivers of variance of nir-v using PROSAIL canopy modeling
#####

import numpy as np
import pyprosail
import random
import scipy
import math
import spectral as spectral
import matplotlib.pyplot as plt
from sklearn import tree
from sklearn import metrics
from sklearn.model_selection import train_test_split

#####
# set up output files and processing parameters
#####

# set number of random veg, bundles to simulate
n_bundles = 2000

# set the number of output bands (default prosail is 2101)
nb = 2101

# set a color palette
palette = ["#CC79A7", "#0072B2", "#D55E00", "#009E73", "#56B4E9", "#F0E442", "#E69F00", "#000000"]

# set output files to store the results
output_csv = 'prosail_spectra.csv'
output_sli = 'prosail_spectra.sli'
output_hdr = 'prosail_spectra.hdr'
output_spec = []

#####
# set up the leaf and canopy modeling parameters
#####

N = []
chloro = []
caroten = []
brown = []
EWT = []
LMA = []
soil_reflectance = []
LAI = []
hot_spot = []
LAD_inclination = []
LAD_bimodality = []
s_az = []
s_za = []
v_az = []
v_za = []
nir_v = []
brightness = []

# set the wavelengths to use in calculating nir_v
red_wl = 0.680
nir_wl = 0.800

# set a dummy run of prosail to find the index for each wavelength
spec = pyprosail.run(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
red_band = np.where(spec[:,0] == red_wl)[0][0]
nir_band = np.where(spec[:,0] == nir_wl)[0][0]

# create a function to calculate nir_v
def calc_nir_v(red, nir):
    return nir * (nir - red) / (nir + red)
    
# and some functions to optimize linear and logistic fits for plotting later
def func_sigmoid(x, x0, k, a, c):
    return a / (1 + np.exp(-k * (x - x0))) + c

def func_linear(x, m, b):
    return m * x + b
    
def func_fit(x, y, function):
    opt, cov = scipy.optimize.curve_fit(function, x, y)
    y_fit = function(np.array(x), *opt)
    rsq = metrics.r2_score(y, y_fit)
    rms = np.sqrt(metrics.mean_squared_error(y, y_fit))
    return [y_fit, rsq, rms]
    
def bn(spectrum):
    return np.sqrt((spectrum ** 2).sum())
    
# find the wavelengths where water absorption does its thing
#  and remove them when calculating brightness
water_bands = [[1.350, 1.460], [1.790, 1.960]]

# start with nir-swir1 transition
gt = np.where(spec[:,0] > water_bands[0][0])
lt = np.where(spec[:,0] < water_bands[0][1])
nd1 = np.intersect1d(gt[0], lt[0])

# then swir1-swir2 transition
gt = np.where(spec[:,0] > water_bands[1][0])
lt = np.where(spec[:,0] < water_bands[1][1])
nd2 = np.intersect1d(gt[0], lt[0])

# concatenate into a single list of indices to remove
water_inds = np.concatenate((nd1, nd2))

# loop through the bundles and generate random canopy parameters
for i in range(n_bundles):
    
    # structural coefficient (arbitrary units)
    #  range 1.3 - 2.5 from Rivera et al. 2013 http://dx.doi.org/10.3390/rs5073280
    N.append(random.uniform(1.3,2.5))

    # total chlorophyll content (ug/cm^2)
    #  range ~ 5 - 75 from Rivera et al. 2013, but i'll set a little more conservative
    chloro.append(random.gauss(35, 30))
    while chloro[-1] < 10 or chloro[-1] > 60:
        chloro[-1] = random.gauss(35, 30)

    # total carotenoid content (ug/cm^2)
    #  kinda fudged this to be like 1/4 of total chl
    caroten.append(random.gauss(8.75, 7.5))
    while caroten[-1] < 2 or caroten[-1] > 15:
        caroten[-1] = random.gauss(8.75, 7.5)

    # brown pigment content (arbitrary units) - not gonna mess with this
    brown.append(0)

    # equivalent water thickness (cm)
    #  range 0.002 - 0.05 from Rivera et al. 2013
    EWT.append(random.uniform(0.002, 0.05))

    # leaf mass per area (g/cm^2)
    #  global range 0.0022 - 0.0365 (median 0.01)
    #  from Asner et al. 2011 http://dx.doi.org/10.1016/j.rse.2011.08.020
    # gonna go a little more conservative
    LMA.append(random.gauss(0.012, 0.005))
    while LMA[-1] < 0.005 or LMA[-1] > 0.0250:
        LMA[-1] = random.gauss(0.012, 0.005)

    # soil reflectance metric (wet soil = 0, dry soil = 1)
    soil_reflectance.append(random.uniform(0,1))

    # leaf area index (unitless, cm^2 leaf area/cm^2 ground area)
    #  range 0.01 - 18.0 (5.5 mean) globally
    #  range 0.2 - 8.7 (3.6 mean) for crops
    #  range 0.6 - 2.8 (1.3 mean) for desert plants
    #  range 0.5 - 6.2 (2.6 mean) for boreal broadleaf forest
    #  range 0.5 - 8.5 (4.6 mean) for boreal needle forest
    #  range 0.8 - 11.6 (5.1 mean) for temperate broadleaf forest
    #  range 0.01 - 15.0 (5.5 mean) for temperate needle forest
    #  range 0.6 - 8.9 (4.8 mean) for tropical broadleaf forest
    #  range 0.3 - 5.0 (1.7 mean) for grasslands
    #  range 1.6 - 18.0 (8.7 mean) for plantations
    #  range 0.4 - 4.5 (2.1 mean) for shrublands
    #  range 0.2 - 5.3 (1.9 mean) for tundra
    #  range 2.5 - 8.4 (6.3 mean) for wetlands
    #  from Asner, Scurlock and Hicke 2003 http://dx.doi.org/10.1046/j.1466-822X.2003.00026.x
    LAI.append(random.gauss(3,2))
    while LAI[-1] < 0.5 or LAI[-1] > 15:
        LAI[-1] = random.gauss(3,2)

    # hot spot parameter (derived from brdf model)
    #  range 0.05-0.5 from Rivera et al. 2013
    hot_spot.append(random.uniform(0.05, 0.5))

    # leaf distribution function parameter.
    #  range LAD_inc -0.4 -  0.4, LAD_bim -0.1 - 0.2 for trees
    #  range LAD_inc -0.1 -  0.3, LAD_bim  0.3 - 0.5 for lianas
    #  range LAD_inc -0.8 - -0.2, LAD_bim -0.1 - 0.3 for Palms
    #  from Asner et al. 2011
    LAD_inclination.append(random.uniform(-0.4, 0.4))
    LAD_bimodality.append(random.uniform(-0.1, 0.2))

    # viewing and solar angle parameters
    #  solar zenith ranges cludged from http://gis.stackexchange.com/questions/191692/maximum-solar-zenith-angle-for-landsat-8-images
    #  I couldn't find good data on the range of possible solar or viewing azimuth.
    #  I decided to set view parameters to 0 to assume nice, clean nadir viewing, and let the sun vary.
    s_za.append(random.uniform(20, 70))
    s_az.append(random.uniform(0,360))
    v_az.append(0)
    v_za.append(0)

#####
# set up the loop for each atmosphere/canopy model
#####

# first create the output array that will contain all the resulting spectra
output_array = np.zeros([nb, (n_bundles) + 1])

# loop through each veg / wood / soil bundle
for j in range(n_bundles):
    
    # load prosail and run the canopy model
    LIDF = (LAD_inclination[j], LAD_bimodality[j])
    spectrum = pyprosail.run(N[j], chloro[j], caroten[j],  
                brown[j], EWT[j], LMA[j], soil_reflectance[j], 
                LAI[j], hot_spot[j], s_za[i], s_az[i],
                v_za[i], v_az[i], LIDF)

    # add the modeled spectrum to the output array
    output_array[:, (j+1)] = spectrum[:,1]
    
    # add a new name to label in the output spectral library
    output_spec.append('veg bundle ' + str(j+1))
    
    # calculate nirv for this spectrum
    nir_v.append(calc_nir_v(spectrum[red_band, 1], spectrum[nir_band, 1]))
    
    # calculate the brightness scalar for this spectrum
    brightness.append(bn(np.delete(spectrum[:,1], water_inds)))

# now that the loop has finished we can export our results to a csv file
output_array[:, 0] = spectrum[:,0]
np.savetxt(output_csv, output_array.transpose(), delimiter=",", fmt = '%.3f')
    
# output a spectral library
with open(output_sli, 'w') as f: 
    output_array[:,1:].transpose().tofile(f)

# write the ENVI header file for the spectral library    
metadata = {
    'samples' : nb,
    'lines' : n_bundles,
    'bands' : 1,
    'data type' : 5,
    'header offset' : 0,
    'interleave' : 'bsq',
    'byte order' : 0,
    'sensor type' : 'prosail',
    'spectra names' : output_spec,
    'wavelength units' : 'micrometers',
    'wavelength' : output_array[:,0]
    }
spectral.envi.write_envi_header(output_hdr, metadata, is_library=True)

# start running some analysis
y = nir_v
x = []
for j in range(n_bundles):
    x.append([N[j], chloro[j], caroten[j], LMA[j], soil_reflectance[j], LAI[j], hot_spot[j], 
        LAD_inclination[j], LAD_bimodality[j], s_az[j], s_za[j], brightness[j]])
        
# split in to train/test splits to eval regression
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = 0.3)
    
# train a few models
mod_simple = tree.DecisionTreeRegressor(max_depth = 2)
mod_medium = tree.DecisionTreeRegressor(max_depth = 5)
mod_large  = tree.DecisionTreeRegressor()

mod_simple.fit(x_train, y_train)
mod_medium.fit(x_train, y_train)
mod_large.fit(x_train, y_train)

# run some predictions
y_test_simple = mod_simple.predict(x_test)
y_test_medium = mod_medium.predict(x_test)
y_test_large = mod_large.predict(x_test)
y_test_list = [y_test_simple, y_test_medium, y_test_large]

# print some outputs
mod = [mod_simple, mod_medium, mod_large]
names = ['simple', 'medium', 'large']

for i in range(len(mod)):
    print("{mod} model".format(mod = names[i]))
    print("r-squared    : {:0.3f}".format(metrics.r2_score(y_test, y_test_list[i])))
    print("explained var: {:0.3f}".format(metrics.explained_variance_score(y_test, y_test_list[i])))
    print("rmse         : {:0.3f}".format(np.sqrt(metrics.mean_squared_error(y_test, y_test_list[i]))))
    print("mean abs err : {:0.3f}".format(metrics.mean_absolute_error(y_test, y_test_list[i])))
    print("feature importance")
    print("N            : {:0.3f}".format(mod[i].feature_importances_[0]))
    print("chlorophyll  : {:0.3f}".format(mod[i].feature_importances_[1]))
    print("carotenoids  : {:0.3f}".format(mod[i].feature_importances_[2]))
    print("LMA          : {:0.3f}".format(mod[i].feature_importances_[3]))
    print("soil_refl    : {:0.3f}".format(mod[i].feature_importances_[4]))
    print("LAI          : {:0.3f}".format(mod[i].feature_importances_[5]))
    print("hot spot     : {:0.3f}".format(mod[i].feature_importances_[6]))
    print("LAD_incl     : {:0.3f}".format(mod[i].feature_importances_[7]))
    print("LAD_biomod   : {:0.3f}".format(mod[i].feature_importances_[8]))
    print("solar azi    : {:0.3f}".format(mod[i].feature_importances_[9]))
    print("solar zen    : {:0.3f}".format(mod[i].feature_importances_[10]))
    print("brightness   : {:0.3f}".format(mod[i].feature_importances_[11]))
    print("-----")
    print("")
    
# start plotting some outputs
plt.figure(1)
plt.suptitle("Canopy trait drivers of NIRv\nLAI {:0.1f} - {:0.1f}".format(min(LAI), max(LAI)))

# first, LAI and nir_v
plt.subplot(321)
lai_fit = func_fit(np.array(LAI), nir_v, func_linear)
plt.scatter(LAI, nir_v, c = palette[0], label = "rmse: {:0.3f}".format(lai_fit[2]))
plt.plot(LAI, lai_fit[0], c = 'black', label = "r-squared: {:0.3f}".format(lai_fit[1]))
plt.title("LAI")
plt.ylabel("NIRv")
plt.legend()

# then, LMA and nir_v
plt.subplot(322)
lma_fit = func_fit(LMA, nir_v, func_linear)
plt.scatter(LMA, nir_v, c = palette[1], label = "rmse: {:0.3f}".format(lma_fit[2]))
plt.plot(LMA, lma_fit[0], c = 'black', label = "r-squared: {:0.3f}".format(lma_fit[1]))
plt.title("LMA")
plt.legend()

# then leaf angle distribution
plt.subplot(323)
lad_fit = func_fit(LAD_inclination, nir_v, func_linear)
plt.scatter(LAD_inclination, nir_v, c = palette[2], label = "rmse: {:0.3f}".format(lad_fit[2]))
plt.plot(LAD_inclination, lad_fit[0], c = "black", label = "r-squared: {:0.3f}".format(lad_fit[1]))
plt.title("LAD")
plt.ylabel("NIRv")
plt.legend()

# then chlorophyll
plt.subplot(324)
chl_fit = func_fit(chloro, nir_v, func_linear)
plt.scatter(chloro, nir_v, c = palette[3], label = "rmse: {:0.3f}".format(chl_fit[2]))
plt.plot(chloro, chl_fit[0], c = "black", label = "r-squared: {:0.3f}".format(chl_fit[1]))
plt.title("CHL")
plt.legend()

# then brightness
plt.subplot(325)
brt_fit = func_fit(brightness, nir_v, func_linear)
plt.scatter(brightness, nir_v, c = palette[4], label = "rmse: {:0.3f}".format(brt_fit[2]))
plt.plot(brightness, brt_fit[0], c = "black", label = "r-squared: {:0.3f}".format(brt_fit[1]))
plt.title("Brightness")
plt.ylabel("NIRv")
plt.legend()

# save the plot, fool
plt.savefig("nirv-params.png")
plt.show()
plt.close()

# plot some canopy spectra
spec_min = output_array[:, np.where(nir_v == min(nir_v))[0] + 1]
spec_max = output_array[:, np.where(nir_v == max(nir_v))[0] + 1]
spec_mean = np.mean(np.array(output_array[:,1:]), axis = 1)
plt.plot(output_array[:,0], spec_min, c = palette[6], label = 'Min NIRv: {:0.3f}'.format(min(nir_v)))
plt.plot(output_array[:,0], spec_mean, c = 'black', label = 'Mean reflectance')
plt.plot(output_array[:,0], spec_max, c = palette[3], label = 'Max NIRv: {:0.3f}'.format(max(nir_v)))
plt.xlabel("Wavelength (um)")
plt.ylabel("Reflectance (%)")
plt.title("Simulated canopy reflectance")
plt.legend()

# then, prediction results from these parameters
plt.figure()
plt.scatter(y_test_large, y_test, c = palette[4], label = "r-squared: {:0.3f}".format(metrics.r2_score(y_test, y_test_large)))
plt.xlabel("predicted")
plt.ylabel("test data")
plt.title("Predicted NIR-v from PROSAIL parameters")
plt.legend()
plt.savefig("nirv-predicted.png")
plt.show()
plt.close()