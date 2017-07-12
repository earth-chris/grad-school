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
figs = base + 'scratch/'
fig_base = "Mesoamerica_pca_1km_"

# set file for predicting distributions
pred_file = "/home/cba/Downloads/mesoamerica/CCB-LC-Mesoamerica-stacked-0000080640-0000107520_southern.tif"

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
arr = np.zeros((ny, nx), np.byte) + 255
pref = None

# set this to false to forego plotting/predicting on images
build_plots = False
build_scatter = False
write_data = True
tuning = True

# set this to true to subset the predictor data
subset_pcs = False
last_pc = 3

# set the color palette
palette = aei.color.color_blind

# file for the GBIF points/random samples
pc_file = base + 'data/mesoamerica_1km_pca.csv'
ma_file = base + 'data/mesoamerica_predictor_data.csv'
pc = pd.read_csv(pc_file)
ma = pd.read_csv(ma_file)

# filter some bad data
pc = pc[pc['PCA 1'] != 255]
ma = ma[ma['TempMin (C)'] != 255]

# determine which features to plot/predict
ftype = "Biophysical"
f = ma

# split into the background, aegypti, and albopictus portions
labels = ['background', 'aegypti', 'albopictus']
bk = f[f['class'] == 0]
ae = f[f['class'] == 1]
aa = f[f['class'] == 2]
band_names = list(bk.columns)
band_names.remove('class')

# set the models to apply
models = ["DecisionTree", "SVM", "RandomForest", "AdaBoosting",
    "GradientBoosting", "MaxEnt"]
mods = [tree.DecisionTreeClassifier(), svm.SVC(), ensemble.RandomForestClassifier(n_jobs = -1),
    ensemble.GradientBoostingClassifier(), ensemble.AdaBoostClassifier(), linear_model.LogisticRegression()]
output_models = []
for model in models:
    output_models.append("{base}Mesoamerica-prediction-{model}".format(
        base = base + 'scratch/', model = model))

#####
# plot the density distributions, if you'd like
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
            xlabel = band_names[i], aei_color = palette)
        
        # save the figure
        plot.savefig(fig_name, dpi = 200)
        plot.close()
        
#####
# plot 2-d scatter plots for each pair-wise combination of predictors
if build_scatter:
    
    # set up list of all unique pairs
    pairs = aei.fn.pairwise(band_names)
    pair_loc = aei.fn.pairwise(range(1, len(band_names) + 1))
    pair_new = []
    for item in pair_loc:
        pair_new.append([item[0], item[1] - 2])
    nplots = max(pair_new)[0]
    
    # set the color palette
    colors = palette(2)
    
    # set up random sampling to balance the background/aegypti with
    #  the smallest class, albopictus
    nbk = len(bk)
    nae = len(ae)
    naa = len(aa)
    rbk = np.random.choice(nbk, naa*2)
    rae = np.random.choice(nae, naa)
    bkaa = bk.iloc[rbk]
    aeaa = ae.iloc[rae]
    
    # set the figure number
    plt.figure(1)
    
    for i in range(len(pairs)):
        
        # set the subplot
        xloc = pair_new[i][0]
        yloc = pair_new[i][1]
        plt.subplot(nplots, nplots, (yloc * nplots) + xloc)
        
        # plot each pairwise set by grouping
        plt.scatter(bkaa[pairs[i][0]], bkaa[pairs[i][1]], color = colors[0],
            label = 'background', alpha = 0.3)
        plt.scatter(aeaa[pairs[i][0]], aeaa[pairs[i][1]], color = colors[1],
            label = 'aegypti', alpha = 0.3)
        plt.scatter(aa[pairs[i][0]], aa[pairs[i][1]], color = colors[1],
            label = 'albopictus', alpha = 0.3)
        
        # plot metadata
        if yloc == (nplots - 1):
            plt.xlabel(pairs[i][0])
        if xloc == 1:
            plt.ylabel(pairs[i][1])
            if yloc == 0:
                plt.legend(bbox_to_anchor = (1, 0.5), loc = 'center left')
    
    # set up meta plot parameters
    plt.suptitle("{} feature covariance".format(ftype))
    plt.tight_layout()
    
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
if subset_pcs:
    for key in band_names[last_pc:]:
        bkae.pop(key)
        bkaa.pop(key)
        band_names.remove(key)
        pred = pred[:,0:last_pc]

# convert to numpy arrays, and get y and xdata for each species
yae = np.array(bkae.pop('class'))
xae = np.array(bkae)
yaa = np.array(bkaa.pop('class'))
xaa = np.array(bkaa)

# get global range of background data
bk.pop('class')
bk_min = np.array(bkae.min(axis = 0))
bk_max = np.array(bkae.max(axis = 0))

# scale data from 0-1
xae_scaled = (xae - bk_min) / (bk_max - bk_min)
xaa_scaled = (xaa - bk_min) / (bk_max - bk_min)
pred_scaled = (pred - bk_min) / (bk_max - bk_min)

# split into train/test groups
xae_train, xae_test, yae_train, yae_test = train_test_split(
    xae_scaled, yae, test_size = 0.25)
xaa_train, xaa_test, yaa_train, yaa_test = train_test_split(
    xaa_scaled, yaa, test_size = 0.25)

# set up arrays to store performance metrics from each model
nm = len(mods)
train_precision = np.zeros((nm,2))
train_recall = np.zeros((nm,2))
train_fscore = np.zeros((nm,2))
train_acc = np.zeros((nm, 2))
test_precision = np.zeros((nm,2))
test_recall = np.zeros((nm,2))
test_fscore = np.zeros((nm,2))
test_acc = np.zeros((nm, 2))
test_auc = np.zeros((nm,2))
test_cm_ae = []
test_cm_aa = []
model_prediction_ae = np.zeros((nm, len(yae_test)))
model_prediction_aa = np.zeros((nm, len(yaa_test)))

# create the model tuning object, and the list of models to tune
ae_tuner = aei.model.tune(xae_train, yae_train)
aa_tuner = aei.model.tune(xaa_train, yaa_train)

models = ["DecisionTree", "SVM", "RandomForest", "AdaBoosting",
"GradientBoosting", "MaxEnt"]

ae_model = [ae_tuner.DecisionTreeClassifier, ae_tuner.SVC, 
    ae_tuner.RandomForestClassifier, ae_tuner.AdaBoostClassifier, 
    ae_tuner.GradientBoostClassifier, ae_tuner.LogisticRegression]
aa_model = [aa_tuner.DecisionTreeClassifier, aa_tuner.SVC, 
    aa_tuner.RandomForestClassifier, aa_tuner.AdaBoostClassifier, 
    aa_tuner.GradientBoostClassifier, aa_tuner.LogisticRegression]

# iterate over each model
for i in range(len(mods)):
    
    # report starting
    print("----------")
    print("Classifying using: {mod}".format(mod = models[i]))
    
    # set up logic to perform model tuning, or if we just run on defaults
    if tuning:
        print("Performing model tuning")
        ae_tuner.param_grid = None
        tuner = ae_model[i]
        tuner(param_grid = None)
        model = ae_tuner.best_estimator
    else:
        # just runnin' on defaults here
        model = mods[i]
        model.fit(xae_train, yae_train)
    
    # get the calibration accuracies
    yae_calib = model.predict(xae_train)
    report = metrics.classification_report(yae_train, yae_calib, 
        target_names = ['background', 'aedes aegyptii'])
    print("Aedes aegypti")
    print("Training results")
    print(report)
    
    # assign training results to output arrays
    prfs = np.array(metrics.precision_recall_fscore_support(yae_train, yae_calib))
    train_precision[i,0] = prfs[0][1]
    train_recall[i,0] = prfs[1][1]
    train_fscore[i,0] = prfs[2][1]
    train_acc[i,0] = metrics.accuracy_score(yae_train, yae_calib)
    
    # predict on the training data and report results
    yae_pred = model.predict(xae_test)
    report = metrics.classification_report(yae_test, yae_pred, 
        target_names = ['background', 'aedes aegyptii'])
    print("Test results")
    print(report)
    if models[i] in ["DecisionTree", "RandomForest", "GradientBoosting"]:
        print("Feature importance:")
        for j in range(len(band_names)):
            print("{}: {:0.3f}%".format(band_names[j], model.feature_importances_[j]))
    
    # assign test results to output arrays
    prfs = np.array(metrics.precision_recall_fscore_support(yae_test, yae_pred))
    test_precision[i,0] = prfs[0][1]
    test_recall[i,0] = prfs[1][1]
    test_fscore[i,0] = prfs[2][1]
    test_acc[i,0] = metrics.accuracy_score(yae_test, yae_pred)
    test_auc[i,0] = metrics.roc_auc_score(yae_test, yae_pred)
    test_cm_ae.append(metrics.confusion_matrix(yae_test, yae_pred))
    model_prediction_ae[i] = yae_pred
    
    # predict on the full image extent
    if write_data:
        
        # run the prediction
        print("Predicting on image extent")
        im_pred = model.predict(pred_scaled)
        arr[gd[0], gd[1]] = im_pred
        outfile = '{}-aegypti-{}.tif'.format(output_models[i], ftype)
        
        # write the output file
        pref = gdal.GetDriverByName("GTiff").Create(outfile, nx, ny, 1, gdal.GDT_Byte)
        pref.SetProjection(prj)
        pref.SetGeoTransform(geo)
        bref = pref.GetRasterBand(1)
        bref.SetNoDataValue(255)
        bref.WriteArray(arr)
        bref = None
        pref = None
    
    #####
    # fit for albopictus next
    
    # set up logic to perform model tuning, or if we just run on defaults
    if tuning:
        print("Performing model tuning")
        aa_tuner.param_grid = None
        tuner = aa_model[i]
        tuner(param_grid = None)
        model = aa_tuner.best_estimator
    else:
        # just runnin' on defaults here
        model = mods[i]
        model.fit(xaa_train, yaa_train)
    
    # get the calibration accuracies
    yaa_calib = model.predict(xaa_train)
    report = metrics.classification_report(yaa_train, yaa_calib,
        target_names = ['background', 'aedes albopictus'])
    print("-----")
    print("Aedes albopictus")
    print("Training results")
    print(report)
    
    # assign training results to output arrays
    prfs = np.array(metrics.precision_recall_fscore_support(yaa_train, yaa_calib))
    train_precision[i,1] = prfs[0][1]
    train_recall[i,1] = prfs[1][1]
    train_fscore[i,1] = prfs[2][1]
    train_acc[i,1] = metrics.accuracy_score(yaa_train, yaa_calib)
    
    # predict on the training data and report results
    yaa_pred = model.predict(xaa_test)
    report = metrics.classification_report(yaa_test, yaa_pred, 
        target_names = ['background', 'aedes albopictus'])
    print("Test results")
    print(report)
    if models[i] in ["DecisionTree", "RandomForest", "GradientBoosting"]:
        print("Feature importance:")
        for j in range(len(band_names)):
            print("{}: {:0.3f}%".format(band_names[j], model.feature_importances_[j]))
    
    # assign test results to output arrays
    prfs = np.array(metrics.precision_recall_fscore_support(yaa_test, yaa_pred))
    test_precision[i,1] = prfs[0][1]
    test_recall[i,1] = prfs[1][1]
    test_fscore[i,1] = prfs[2][1]
    test_acc[i,1] = metrics.accuracy_score(yaa_test, yaa_pred)
    test_auc[i,1] = metrics.roc_auc_score(yaa_test, yaa_pred)
    test_cm_aa.append(metrics.confusion_matrix(yaa_test, yaa_pred))
    model_prediction_aa[i] = yaa_pred
    
    # predict on the full image extent
    if write_data:
        
        # run the prediction
        print("Predicting on image extent")
        im_pred = model.predict(pred_scaled)
        arr[gd[0], gd[1]] = im_pred
        outfile = '{}-albopictus-{}.tif'.format(output_models[i], ftype)
        
        # write the output file
        pref = gdal.GetDriverByName("GTiff").Create(outfile, nx, ny, 1, gdal.GDT_Byte)
        pref.SetProjection(prj)
        pref.SetGeoTransform(geo)
        bref = pref.GetRasterBand(1)
        bref.SetNoDataValue(255)
        bref.WriteArray(arr)
        bref = None
        pref = None
    
#####
# start plotting these meta model values

# plot the train vs test scores for each model method
colors = palette(4)
bar_width = 0.2
index = np.arange(nm)
alpha = 0.8
edgecolor = 'black'

metrics = ["Precision", "Recall", "F-Score"]
train_vars = [train_precision, train_recall, train_fscore]
test_vars = [test_precision, test_recall, test_fscore]

for i in range(len(metrics)):
    fig, ax = plt.subplots()
    plt.bar(index, train_vars[i][:,0], bar_width, color = colors[0], 
        label = 'aegypti training data', alpha = alpha, edgecolor = edgecolor)
    plt.bar(index + bar_width, test_vars[i][:,0], bar_width,
        color = colors[1], label = "aegypti test data", alpha = alpha, edgecolor = edgecolor)
    plt.bar(index + bar_width * 2, train_vars[i][:,1], bar_width, 
        color = colors[2], label = 'albopictus training data', alpha = alpha, edgecolor = edgecolor)
    plt.bar(index + bar_width * 3, test_vars[i][:,1], bar_width,
        color = colors[3], label = "albopictus test data", alpha = alpha, edgecolor = edgecolor)
    plt.xlabel("Model Type")
    plt.ylabel(metrics[i])
    plt.title("{} model {}".format(ftype, metrics[i]))
    plt.xticks(index + (1.5 * bar_width), models)
    plt.legend()
    plt.tight_layout()
    
# and plot the accuracy/AUC scores for each model
colors = aei.color.stanford_bright(4)
fig, ax = plt.subplots()
plt.bar(index, test_acc[:,0], bar_width, color = colors[0],
    label = 'aegypti accuracy', alpha = alpha, edgecolor = edgecolor)
plt.bar(index + bar_width, test_auc[:,0], bar_width, color = colors[1],
    label = 'aegypti AUC', alpha = alpha, edgecolor = edgecolor)
plt.bar(index + bar_width * 2, test_acc[:,1], bar_width, color = colors[2],
    label = 'albopictus accuracy', alpha = alpha, edgecolor = edgecolor)
plt.bar(index + bar_width * 3, test_auc[:,1], bar_width, color = colors[3],
    label = 'albopictus AUC', alpha = alpha, edgecolor = edgecolor)
plt.xlabel("Model Type")
plt.ylabel("Score")
plt.title("{} model accuracy scores".format(ftype))
plt.xticks(index + (1.5 * bar_width), models)
plt.legend()
plt.ylim(0.25, 0.9)
plt.tight_layout()