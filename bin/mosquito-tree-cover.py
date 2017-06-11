#!/usr/bin/python
#####
# plots distributions of predictor variables 
#####

import gdal
import scipy
import numpy as np
from sklearn import tree
from sklearn import ensemble
from sklearn import svm
from sklearn import metrics
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
pred_file = base + 'coto_brus_predictors_masked2.tif'
pred_ref = gdal.Open(pred_file)
pred_names = ["Palsar HH", "Palsar HV", "% Soil", "% Veg", "% Urban", "NIRv"]

# training data
train_file = base + 'coto_brus_tree_cover.tif'
train_ref = gdal.Open(train_file)
train_nd = train_ref.GetRasterBand(1).GetNoDataValue()

# name of the output mode
model_file = base + 'coto_brus_tree_cover_model.sav'

# read the data into memory and reduce dimensions
train_arr = train_ref.ReadAsArray()
gd = np.where(train_arr != train_nd)
train = train_arr[gd[0], gd[1]]
train_arr = None

pred_arr = pred_ref.ReadAsArray()
pred = pred_arr[:, gd[0], gd[1]]
pred_arr = None

# create the train/test split
x_train, x_test, y_train, y_test = train_test_split(
    pred.transpose(), train, test_size = 0.5)
    
# create the regression model
rf_def = ensemble.RandomForestRegressor(n_jobs = 7, verbose = 1)

# train the model
rf_def.fit(x_train, y_train)

# save the model
pickle.dump(rf_def, open(model_file, 'wb'))

# test the model
y_pred = rf_def.predict(x_test)

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
print("var. importance")
for i in range(len(pred_names)):
    print("{} : {:0.3f}".format(pred_names[i], rf_def.feature_importances_[i]))

# plot the output
plt.figure(1)
lin_fit = func_fit(y_test, y_pred, func_linear)
plt.hexbin(y_pred, y_test, cmap = 'plasma', alpha = 0.95, bins = 'log',
    label = "rmse: {:0.3f}".format(rmse))
plt.plot(lin_fit[0], y_test, c = "black", label = "r-sqared: {:0.3f}".format(lin_fit[1]))
