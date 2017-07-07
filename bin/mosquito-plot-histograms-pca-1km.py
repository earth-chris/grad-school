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
from sklearn import linear_model
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# set base directory for files
base = aei.params.environ['AEI_GS'] + '/'
figs = base + 'figures/'
fig_base = "Mesoamerica_pca_1km_"

# set file for predicting distributions
pred_file = "/home/cba/Downloads/mesoamerica/CCB-LC-Mesoamerica-stacked-0000080640-0000107520.tif"

# open the file, get the metadata, read it into memory, and subset it
pref = gdal.Open(pred_file)
prj = pref.GetProjection()
geo = pref.GetGeoTransform()
nx = pref.RasterXSize
ny = pref.RasterYSize
b1 = pref.GetRasterBand(1).ReadAsArray()
gd = np.where(b1 != 255)
arr = pref.ReadAsArray()
pred = arr[:, gd[0], gd[1]].transpose()
arr = None
pref = None

# set this to false to forego plotting
build_plots = False

# file for the GBIF points/random samples
pc_file = base + 'data/mesoamerica_1km_pca.csv'
ma_file = base + 'data/mesoamerica_predictor_data.csv'
pc = pd.read_csv(pc_file)
ma = pd.read_csv(ma_file)

# filter some bad data
pc = pc[pc['value_0'] != 255]
ma = ma[ma['TempMin'] != 255]

# determine which features to plot/predict
f = ma

# split into the background, aegypti, and albopictus portions
labels = ['background', 'aegypti', 'albopictus']
bk = f[f['class'] == 0]
ae = f[f['class'] == 1]
aa = f[f['class'] == 2]
band_names = list(bk.columns)
band_names.remove('class')

# set the models to apply
models = ["DecisionTree", "SVM", "RandomForest", "GradientBoosting", "MaxEnt"]
mods = [tree.DecisionTreeClassifier(), svm.SVC(), ensemble.RandomForestClassifier(),
    ensemble.AdaBoostClassifier(), linear_model.LogisticRegression()]
output_models = []
for model in models:
    output_models.append("{base}Mesoamerica-prediction-{model}.tif".format(
        base = base + 'scratch/', model = model))

# plot the data, if you'd like
if build_plots:
    # for each pc, plot a density distribution
    for i in range(len(band_names)):
        
        # set output figure name
        fig_name = "{}{}{}.png".format(figs, fig_base, band_names[i])
        
        # create a list for plotting
        vec = [np.array(bk[band_names[i]]), np.array(ae[band_names[i]]),
            np.array(aa[band_names[i]])]
        
        # plot the distributions
        plot = aei.plot.density_dist(vec, label = labels, 
            xlabel = band_names[i], aei_color = aei.color.stanford_bright)
        
        # save the figure
        plot.savefig(fig_name, dpi = 200)
    
#####
# start classifying distributions

# set albopictus class to 1 to run as binary classification
aa['class'] = 1

# build models for bk/ae and bk/aa based on balanced classes
nbk = len(bk)
nae = len(ae)
naa = len(aa)

# indices for random background sampling
rae = np.random.choice(nbk, nae)
raa = np.random.choice(nbk, naa)

# subset the random samples and append each species to a single data frame
bkae = bk.iloc[rae].append(ae)
bkaa = bk.iloc[raa].append(aa)

# remove the last few pc's
#for key in band_names[2:]:
#    bkae.pop(key)
#    bkaa.pop(key)

# convert to numpy arrays, and get y and xdata for each species
yae = np.array(bkae.pop('class'))
xae = np.array(bkae)
yaa = np.array(bkaa.pop('class'))
xaa = np.array(bkaa)

# get global range of background data
bk.pop('class')
bk_min = np.array(bk.min(axis = 0))
bk_max = np.array(bk.max(axis = 0))

# scale data from 0-1
xae_scaled = (xae - bk_min) / (bk_max - bk_min)
xaa_scaled = (xaa - bk_min) / (bk_max - bk_min)

# split into train/test groups
xae_train, xae_test, yae_train, yae_test = train_test_split(
    xae, yae, test_size = 0.25)
xaa_train, xaa_test, yaa_train, yaa_test = train_test_split(
    xaa, yaa, test_size = 0.25)

# iterate over each model
for i in range(len(mods)):
    
    # report starting
    print("----------")
    print("Classifying using: {mod}".format(mod = models[i]))
    
    # just runnin' on defaults here
    model = mods[i]
    
    # fit for aegypti first
    model.fit(xae_train, yae_train)
    
    # get the calibration accuracies
    yae_calib = model.predict(xae_train)
    report = metrics.classification_report(yae_train, yae_calib, 
        target_names = ['background', 'aedes aegyptii'])
    print("Aedes aegypti")
    print("Training results")
    print(report)
    
    # predict on the training data and report results
    yae_pred = model.predict(xae_test)
    report = metrics.classification_report(yae_test, yae_pred, 
        target_names = ['background', 'aedes aegyptii'])
    print("Test results")
    print(report)
    #if i in [0, 2, 3]:
    #    print("Feature importance:")
    #    for j in range(len(band_names)-1):
    #        print("{}: {:0.3f}%".format(band_names[j], model.feature_importances_[j]))
    
    # fit for albopictus next
    model.fit(xaa_train, yaa_train)
    
    # get the calibration accuracies
    yaa_calib = model.predict(xaa_train)
    report = metrics.classification_report(yaa_train, yaa_calib,
        target_names = ['background', 'aedes albopictus'])
    print("-----")
    print("Aedes albopictus")
    print("Training results")
    print(report)
    
    # predict on the training data and report results
    yaa_pred = model.predict(xaa_test)
    report = metrics.classification_report(yaa_test, yaa_pred, 
        target_names = ['background', 'aedes albopictus'])
    print("Test results")
    print(report)
    #if i in [0, 2, 3]:
    #    print("Feature importance:")
    #    for j in range(len(band_names)-1):
    #        print("{}: {:0.3f}%".format(band_names[j], model.feature_importances_[j]))
    
    #####
    # predict on the full image extent