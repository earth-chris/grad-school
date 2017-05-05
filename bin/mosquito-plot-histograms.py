#!/usr/bin/python
#####
# plots distributions of predictor variables 
#####

import aei
import gdal
import numpy as np
from sklearn import tree
from sklearn import ensemble
from sklearn import svm
from sklearn import metrics
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# set base directory for files
base = '/home/cba/Downloads/mosquito-proposal/'

# file for CR-wide data
cr_file = base + 'CR-aligned-byte.tif'
cr_ref = gdal.Open(cr_file)

# file for Coto Brus data
cb_file = base + 'CR-aligned-southern.tif'
cb_ref = gdal.Open(cb_file)

# file for current sampling data
fs_cb_file = base + 'CR-plot-locations-southern.tif'
fs_cb_ref = gdal.Open(fs_cb_file)

# file for Aedes aegyptii points
ae_file = base + 'CR-mosquito-locations-90m.tif'
ae_ref = gdal.Open(ae_file)

# set the colors and band names to plot
plt_bands = ["Temperature (C)", "Soil Cover (%)", "Veg. Cover (%)", "Impervious Cover (%)", "Cloud Cover (%)"]
bands = ["Temperature", "Soil-Cover", "Veg-Cover", "Impervious-Cover", "Cloud-Cover"]
output_files = []
for band in bands:
    output_files.append("{base}CR-{band}-plot.png".format(base = base, band = band))
cols = aei.color.color_blind(len(bands))

# set the models to apply
models = ["DecisionTree", "SVM", "RandomForest", "GradientBoosting"]
mods = [tree.DecisionTreeClassifier(), svm.SVC(), ensemble.RandomForestClassifier(),
    ensemble.AdaBoostClassifier()]
output_models = []
for model in models:
    output_models.append("{base}CR-southern-prediction-{model}.tif".format(base = base, model = model))

# grab the indices for the aedes aegyptii and field plot data
fs_cb_band = fs_cb_ref.ReadAsArray()
fs_cb = np.where(fs_cb_band == 1)
fs_cb_band = None
ae_band = ae_ref.ReadAsArray()
ae_ind = np.where(ae_band == 1)
ae_band = None

# count how many pixels in coto brus field sites to use as a limit on the
#  number of random samples to pull from the national- and coto brus-wide data
nscr = int(2e5)
nfs = fs_cb[0].shape[0]
nae = ae_ind[0].shape[0]
print("Number of field plot samples: {}".format(nfs))
print("Number of aedes aegyptii samples: {}".format(nae))

# set up a variable to hold data to predict aedes occurrences
y = np.zeros(nae * 4)
y[0:nae] = 1
x = np.zeros((nae * 4, len(bands)))

# loop through each band and plot histograms of nation-wide data,
#  coto brus-wide data, and the data from all plots
for i in range(cr_ref.RasterCount):
    print("Processing band: {band}".format(band = bands[i]))
    cr_band = cr_ref.GetRasterBand(i+1)
    cb_band = cb_ref.GetRasterBand(i+1)
    cr_nodata = cr_band.GetNoDataValue()
    cb_nodata = cb_band.GetNoDataValue()
    
    # read the national data into memory and extract the values to plot
    cr_arr = cr_band.ReadAsArray()
    
    # if this is the first band, select a random subset of pixels to sample and use that throughout
    if i == 0:
        cr_ind = np.where(cr_arr != cr_nodata)
        cr_rnd = np.int32(np.floor(np.random.uniform(0, cr_ind[0].shape[0], nscr)))
        cr_ae = np.int32(np.floor(np.random.uniform(0, cr_ind[0].shape[0], nae*3)))
    
    # subset the national and aedes data and clear the rest out of memory    
    cr_data = cr_arr[cr_ind[0][cr_rnd], cr_ind[1][cr_rnd]]
    ae_data = cr_arr[ae_ind[0], ae_ind[1]]
    cr_ae_data = cr_arr[cr_ind[0][cr_ae], cr_ind[1][cr_ae]]
    cr_arr = None
    
    # read the coto brus data into memory and extract the values to plot
    cb_arr = cb_band.ReadAsArray()
    if i == 0:
        cb_ind = np.where(cb_arr != cb_nodata)
    
    # pull the field plot and coto brus data
    fs_data = cb_arr[fs_cb[0], fs_cb[1]]
    cb_rnd = np.int32(np.floor(np.random.uniform(0, cb_ind[0].shape[0], nscr)))
    cb_data = cb_arr[cb_ind[0][cb_rnd], cb_ind[1][cb_rnd]]
    cb_arr = None
    
    # find the min and max for each data set to set bounds for plot
    xmin = min(np.percentile(cr_data, 2), np.percentile(cb_data, 2), 
        np.percentile(fs_data, 2), np.percentile(ae_data, 2))
    xmax = max(np.percentile(cr_data, 98), np.percentile(cb_data, 98), 
        np.percentile(fs_data, 98), np.percentile(ae_data, 98))
    
    # get a density distribution for each 
    cr_dns = gaussian_kde(cr_data)
    cb_dns = gaussian_kde(cb_data)
    fs_dns = gaussian_kde(fs_data)
    ae_dns = gaussian_kde(ae_data)
    
    # set custom covariance to smooth multiple peaks
    covar = 0.25
    cr_dns.covariance_factor = lambda : covar
    cb_dns.covariance_factor = lambda : covar
    fs_dns.covariance_factor = lambda : covar
    ae_dns.covariance_factor = lambda : covar
    cr_dns._compute_covariance()
    cb_dns._compute_covariance()
    fs_dns._compute_covariance()
    ae_dns._compute_covariance()
    
    # set an x scale to plot smooth density distribution
    xs = np.linspace(xmin, xmax, 200)
    
    # plot each function in a single plot
    plt.figure()
    plt.plot(xs, cr_dns(xs), label = "Costa Rica", color = cols[0])
    plt.plot(xs, cb_dns(xs), label = "Southern CR", color = cols[1])
    plt.plot(xs, fs_dns(xs), label = "Field Plots", color = cols[2])
    plt.plot(xs, ae_dns(xs), label = "Aedes aegyptii", color = cols[3])
    plt.xlabel(plt_bands[i])
    plt.ylabel("Density")
    plt.title("Costa Rica {var} Distributions".format(var = bands[i]))
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_files[i], dpi = 200)
    plt.close()
    
    # add data to the predictor array
    x[0:nae, i] = ae_data
    x[nae:, i] = cr_ae_data

    # that's all, folks!
    
# now, try and predict aedes occurrence!
# first, split into test/train data
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = 0.3)
    
# then loop through the classifiers
for i in range(len(mods)):
    print("----------")
    print("Classifying using: {mod}".format(mod = models[i]))
    
    # just runnin' on defaults here
    model = mods[i]
    model.fit(x_train[:,:-1], y_train)
    
    # predict on the training data and report results
    y_pred = model.predict(x_test[:,:-1])
    report = metrics.classification_report(y_test, y_pred, 
        target_names = ['random sampling', 'aedes aegyptii'])
    print(report)
    if i in [0, 2, 3]:
        print("Feature importance:")
        for j in range(len(bands)-1):
            print("{}: {:0.3f}%".format(bands[j], model.feature_importances_[j]))
        
    # apply the prediction to the southern costa rica region
    scarr = cb_ref.ReadAsArray()
    scarr = scarr[:-1, cb_ind[0], cb_ind[1]]
    newpred = model.predict(scarr.transpose())
    newarr = np.zeros((cb_ref.RasterYSize, cb_ref.RasterXSize), np.byte)
    newarr[cb_ind[0], cb_ind[1]] = newpred
    
    # write to an output file
    newref = gdal.GetDriverByName("GTiff").Create(output_models[i], cb_ref.RasterXSize, 
        cb_ref.RasterYSize, 1, gdal.GDT_Byte)
    newref.SetGeoTransform(cb_ref.GetGeoTransform())
    newref.SetProjection(cb_ref.GetProjection())
    newband = newref.GetRasterBand(1)
    newband.WriteArray(newarr)
    newband = None
    newarr = None
    newref = None