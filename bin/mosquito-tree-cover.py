#!/usr/bin/python
#####
# tries predicting tree cover in coto brus
#####

import gdal
import scipy
import numpy as np
from sklearn import tree
from sklearn import ensemble
from sklearn import linear_model
from sklearn import svm
from sklearn import metrics
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pickle

# some functions for plotting later
def func_linear(x, m, b):
    return m * x + b
    
def func_fit(x, y, function):
    opt, cov = scipy.optimize.curve_fit(function, x, y)
    y_fit = function(np.array(x), *opt)
    rsq = metrics.r2_score(y, y_fit)
    rms = np.sqrt(metrics.mean_squared_error(y, y_fit))
    return [y_fit, rsq, rms]

# set base directory for files
base = '/home/salo/Downloads/tree-cover/'

# predictor data
#pred_file = base + 'coto_brus_predictors_masked2.tif'
pred_file = base + 'coto_brus_despeckle.tif'
pred_ref = gdal.Open(pred_file)
pred_prj = pred_ref.GetProjection()
pred_geo = pred_ref.GetGeoTransform()
nx = pred_ref.RasterXSize
ny = pred_ref.RasterYSize
nb = pred_ref.RasterCount
#pred_names = ["Palsar HH", "Palsar HV", "% Soil", "% Veg", "% Urban", "NIRv"]
pred_names = ["Palsar HH", "Palsar HV", "NIRv", "% Soil", "% Veg", "% Urban"]

# training data
train_file = base + 'coto_brus_tree_cover.tif'
train_ref = gdal.Open(train_file)
train_nd = train_ref.GetRasterBand(1).GetNoDataValue()

# name of the output mode
model_file = base + 'coto_brus_tree_cover_model'

# read the data into memory and reduce dimensions
train_arr = train_ref.ReadAsArray()
gd = np.where(train_arr != train_nd)
train = train_arr[gd[0], gd[1]]
train_arr = None

pred_arr = pred_ref.ReadAsArray()
pred = pred_arr[:, gd[0], gd[1]]
pred_arr = None
pred = pred.transpose()

# transform the predictor data to 0-1 range
scaler = preprocessing.MinMaxScaler()
pred = scaler.fit_transform(pred)

# save the scaler to use when applying later
pickle.dump(scaler, open(base + 'tree_cover_scaler.sav', 'wb'))

# create the train/test split
x_train, x_test, y_train, y_test = train_test_split(
    pred, train, test_size = 0.8)

# create a more advanced train/test split based on even
#  sampling across tree cover bins (here, 10%)
binsize = 0.1
bins = np.arange(0, 1.00001, binsize)
#bins = np.array([0, 0.1, 1.0001])

# use 90% of the least common bin size as the number of per-bin samples
train_hist = np.histogram(train, bins = bins)
nrnd = int(np.floor(min(train_hist[0]) * 0.9))

# to get good representation of low tree cover, add a bunch more 
#  samples from that class
lt10 = 10 
nrnd_lt10 = lt10 * nrnd

# create arrays to store the train and test data
n_train = int(np.floor(nrnd * len(train_hist[0]) + (nrnd_lt10 - nrnd)))
n_test = train.shape[0] - n_train

x_train = np.zeros((n_train, nb))
y_train = np.zeros(n_train)
x_test = np.zeros((n_test, nb))
y_test = np.zeros(n_test)

# create a counter to store the indices for where the test data go in the test array
train_start = 0
test_start = 0

for i in range(len(bins)-1):
    # find the indices for this bin, then randomize them
    ll = np.where(train >= bins[i])
    if i != len(bins)-2:
        ul = np.where(train < bins[i+1])
    else:
        ul = np.where(train <= bins[i+1])
    intersect = np.intersect1d(ul, ll)
    np.random.shuffle(intersect)
    
    # separate the indices for each bin into train and test sets
    if i == 0:
        int_train = intersect[:nrnd_lt10]
        int_test = intersect[nrnd_lt10:]
    else:    
        int_train = intersect[:nrnd]
        int_test = intersect[nrnd:]
    
    # extract these data from the response and predictor data
    n_test = len(int_test)
    if i == 0:
        n_train = nrnd_lt10
    else:
        n_train = nrnd
    
    x_train[train_start:train_start+n_train, :] = pred[int_train, :]
    y_train[train_start:train_start+n_train] = train[int_train]
    train_start += n_train
    
    x_test[test_start:test_start+n_test] = pred[int_test, :]
    y_test[test_start:test_start+n_test] = train[int_test]
    test_start += n_test
    
# create the regression models
rf = ensemble.RandomForestRegressor(n_jobs = 7, verbose = 2, n_estimators=15)
gb = ensemble.GradientBoostingRegressor(learning_rate=0.1, min_samples_split=2,
    n_estimators=800, verbose = 2)
et = ensemble.ExtraTreesRegressor(n_jobs = 7, verbose = 1)
br = ensemble.BaggingRegressor(n_jobs = 7, verbose = 1)
ab = ensemble.AdaBoostRegressor()
ol = linear_model.LinearRegression(n_jobs = 7)

models = [rf, gb, et, br, ab, ol]
model_names = ['RandomForest', 'GradientBoosting', 'ExtraTrees',
    'Bagging', 'AdaBoost', 'LinearModel']
    
output_pred = {}
output_test = {}

for i in range(len(models)):
    
    # report starting
    print("-----")
    print("Processing {} model".format(model_names[i]))
    
    # train the model
    models[i].fit(x_train, y_train)
    
    model_file_name = "{}_{}.sav".format(model_file, model_names[i])
    
    # save the model
    pickle.dump(models[i], open(model_file_name, 'wb'))
    
    # test the model
    y_pred = models[i].predict(x_test)
    output_test[model_names[i]] = y_pred
    
    # output the full model to a dictionary
    opred = models[i].predict(pred)
    output_pred[model_names[i]] = opred
    
    # get metrics
    r_squared = metrics.r2_score(y_test, y_pred)
    explained_var = metrics.explained_variance_score(y_test, y_pred)
    rmse = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
    mae = metrics.mean_absolute_error(y_test, y_pred)
    
    # print model results
    print("r-squared    : {:0.3f}".format(r_squared))
    print("explained var: {:0.3f}".format(explained_var))
    print("rmse         : {:0.3f}".format(rmse))
    print("mean abs err : {:0.3f}".format(mae))
    #print("var. importance")
    #for j in range(len(pred_names)):
    #    print("{} : {:0.3f}".format(pred_names[j], rf_def.feature_importances_[j]))
    
    # write to an output raster
    output_file_name = "{}_{}.tif".format(model_file, model_names[i])
    oref = gdal.GetDriverByName("GTiff").Create(output_file_name,
        nx, ny, 1, gdal.GDT_Float32)
    oref.SetProjection(pred_prj)
    oref.SetGeoTransform(pred_geo)
    ob = oref.GetRasterBand(1)
    outarr = np.zeros((ny, nx))
    outarr[gd[0], gd[1]] = opred
    ob.WriteArray(outarr)
    outarr = None
    ob = None
    oref = None
    
    # plot the output
    #plt.figure(1)
    #lin_fit = func_fit(y_test, y_pred, func_linear)
    #plt.hexbin(y_pred, y_test, cmap = 'plasma', alpha = 0.95, bins = 'log',
    #    label = "rmse: {:0.3f}".format(rmse), gridsize = 20)
    #plt.plot(lin_fit[0], y_test, c = "black", label = "r-sqared: {:0.3f}".format(lin_fit[1]))
    #plt.xlabel("Predicted")
    #plt.ylabel("Actual")
    #plt.title("{} Tree Cover".format(model_names[i]))
    #plt.legend()
    #plt.show()
    #plt.tight_layout()
    #outfile = "{}_{}.png".format(model_file, model_names[i])
    #plt.savefig(outfile)
    #plt.close()
