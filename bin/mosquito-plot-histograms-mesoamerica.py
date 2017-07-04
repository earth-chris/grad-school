#!/usr/bin/python
#####
# plots distributions of predictor variables 
#####

import aei
import gdal
import pandas as pd
import numpy as np
from sklearn import tree
from sklearn import ensemble
from sklearn import svm
from sklearn import metrics
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# set base directory for files
base = '/home/cba/cba/aei-grad-school/'
figs = base + 'figures/'

# file for the GBIF points/random samples
ma_file = base + 'data/mesoamerica_predictor_data.csv'
cr_file = base + 'data/costa_rica_predictor_data.csv'
ma = pd.read_csv(ma_file)
cr = pd.read_csv(cr_file)

# filter some bad data
ma = ma[ma['TempMax'] != 255]
cr = cr[cr['TempMax'] != 255]

# file for field sample data
cb_file = base + 'data/aedes-sites-extract.csv'
sites = pd.read_csv(cb_file)

# set the colors and band names to plot
plt_bands = ["Soil Cover (%)", "Veg. Cover (%)", "Impervious Cover (%)", 
    "Min. Temp (C)", "Median Temp (C)", "Max. Temp (C)"]
bands = ["Soil-Cover", "Veg-Cover", "Impervious-Cover", "Temp-Min", 
    "Temp-Median", "Temp-Max"]
columns = ['Soil', 'Veg', 'Imp', 'TempMin', 'TempMed', 'TempMax']
labels = ['Mesoamerica', 'Aedes-aegypti', 'Aedes-albopictus', 'Field-sites', 'Costa Rica']
output_files = []
for band in bands:
    output_files.append("{base}Mesoamerica-{band}-plot.png".format(base = figs, band = band))
cols = aei.color.color_blind(len(bands))

# set the models to apply
models = ["DecisionTree", "SVM", "RandomForest", "GradientBoosting"]
mods = [tree.DecisionTreeClassifier(), svm.SVC(), ensemble.RandomForestClassifier(),
    ensemble.AdaBoostClassifier()]
output_models = []
for model in models:
    output_models.append("{base}Mesoamerica-prediction-{model}.tif".format(
        base = base + 'scratch/', model = model))

# count how many pixels in coto brus field sites to use as a limit on the
#  number of random samples to pull from the national- and coto brus-wide data
nrnd = ma[ma['class'] == 0].count()[0]
naeg = ma[ma['class'] == 1].count()[0]
nalb = ma[ma['class'] == 2].count()[0]
nplt = sites.count()[0]

# report
print("Number of field plot sites: {}".format(nplt))
print("Number of aedes aegyptii samples  : {}".format(naeg))
print("Number of aedes albopictus samples: {}".format(nalb))

# loop through each band and plot histograms of nation-wide data,
#  coto brus-wide data, and the data from all plots
for i in range(len(columns)):
    
    # subset data to just what we're working with
    m0_data = ma[ma['class'] == 0][columns[i]]
    m1_data = ma[ma['class'] == 1][columns[i]]
    m2_data = ma[ma['class'] == 2][columns[i]]
    st_data = sites[columns[i]]
    cr_data = cr[cr['class'] == 1][columns[i]]
    
    # find the min and max for each data set to set bounds for plot
    xmin = min(np.percentile(m0_data, 2), np.percentile(m1_data, 2), 
        np.percentile(m2_data, 2), np.percentile(st_data, 2))
    xmax = max(np.percentile(m0_data, 98), np.percentile(m1_data, 98), 
        np.percentile(m2_data, 98), np.percentile(st_data, 98))
    
    # get a density distribution for each 
    m0_dns = gaussian_kde(m0_data)
    m1_dns = gaussian_kde(m1_data)
    m2_dns = gaussian_kde(m2_data)
    st_dns = gaussian_kde(st_data)
    cr_dns = gaussian_kde(cr_data)
    
    # set custom covariance to smooth multiple peaks
    covar = 0.25
    m0_dns.covariance_factor = lambda : covar
    m1_dns.covariance_factor = lambda : covar
    m2_dns.covariance_factor = lambda : covar
    st_dns.covariance_factor = lambda : covar
    cr_dns.covariance_factor = lambda : covar
    m0_dns._compute_covariance()
    m1_dns._compute_covariance()
    m2_dns._compute_covariance()
    st_dns._compute_covariance()
    cr_dns._compute_covariance()
    
    # set an x scale to plot smooth density distribution
    xs = np.linspace(xmin, xmax, 200)
    
    # plot each function in a single plot
    plt.figure()
    #plt.plot(xs, m0_dns(xs), label = labels[0], color = cols[0])
    plt.plot(xs, m1_dns(xs), label = labels[1], color = cols[1])
    plt.plot(xs, m2_dns(xs), label = labels[2], color = cols[2])
    plt.plot(xs, st_dns(xs), label = labels[3], color = cols[3])
    #plt.plot(xs, cr_dns(xs), label = labels[4], color = cols[4])
    plt.xlabel(plt_bands[i])
    plt.ylabel("Density")
    plt.title("{var} Distributions".format(var = bands[i]))
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_files[i], dpi = 200)
    plt.close()

    # that's all, folks!
    
# now, try and predict aedes occurrence!
# first, split into test/train data
y = ma.pop('class')
x = ma
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = 0.3)
    
# then loop through the classifiers
for i in range(len(mods)):
    print("----------")
    print("Classifying using: {mod}".format(mod = models[i]))
    
    # just runnin' on defaults here
    model = mods[i]
    model.fit(x_train, y_train)
    
    # predict on the training data and report results
    y_pred = model.predict(x_test)
    report = metrics.classification_report(y_test, y_pred, 
        target_names = ['random sampling', 'aedes aegyptii', 'aedes albopictus'])
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