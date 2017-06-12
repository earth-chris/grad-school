#!/usr/bin/python
#####
# applies model predictions to tiled CR predictor data
#####

import gdal
import scipy
import numpy as np
from sklearn import tree
from sklearn import ensemble
from sklearn import linear_model
from sklearn import svm
from sklearn import metrics
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pickle

# set base directory
base = '/home/cba/Downloads/tree-cover'

# list the tiles to compute over
tiles = [base + '/tiles/pred_001.tif', base + '/tiles/pred_002.tif', base + '/tiles/pred_003.tif', 
    base + '/tiles/pred_004.tif', base + '/tiles/pred_005.tif', base + '/tiles/pred_006.tif', 
    base + '/tiles/pred_008.tif', base + '/tiles/pred_009.tif']
    
# list the models to apply
model_files = ["coto_brus_tree_cover_model_AdaBoost.sav", 
          "coto_brus_tree_cover_model_Bagging.sav", 
          "coto_brus_tree_cover_model_ExtraTrees.sav", 
          "coto_brus_tree_cover_model_GradientBoosting.sav", 
          "coto_brus_tree_cover_model_LinearModel.sav", 
          "coto_brus_tree_cover_model_RandomForest.sav"]
          
# load the models to memory
models = []
for file in model_files:
    models.append(pickle.load(open(file, 'r')))
    
# loop through each tile, apply the model, and save the output
for i in range(len(tiles)):
    
    # report
    print("-----")
    print("Predicting file: {}".format(tiles[i]))
    
    # read input reference
    tref = gdal.Open(tiles[i])
    tgeo = tref.GetGeoTransform()
    tprj = tref.GetProjection()
    nx = tref.RasterXSize
    ny = tref.RasterYSize
    
    # get nodata value and indices
    bref = tref.GetRasterBand(1)
    ndval = bref.GetNoDataValue()
    band = bref.ReadAsArray()
    gd = np.where(band != ndval)
    band = None
    bref = None
    
    # read the predictor data into memory
    print("Reading data into memory")
    pred_arr = tref.ReadAsArray()
    pred = pred_arr[:, gd[0], gd[1]].transpose()
    pred_arr = None
    
    # predict each model
    print("Predicting each model")
    opred = np.zeros(len(gd[0]))
    for model in models:
        y_pred = model.predict(pred)
        opred += y_pred
        
    # take the average
    opred /= float(len(models))
    
    # create the square array to write as output
    print("Creating output array")
    outarr = np.zeros((ny, nx)) + 2.55
    outarr[gd[0], gd[1]] = opred
    
    # write to an output file
    outfile = base + "/tiles/predicted_{:03d}.tif".format(i+1)
    oref = gdal.GetDriverByName("GTiff").Create(outfile, nx, ny, 1, gdal.GDT_Float32)
    oref.SetGeoTransform(tgeo)
    oref.SetProjection(tprj)
    oband = oref.GetRasterBand(1)
    oband.WriteArray(outarr)
    oband = None
    outarr = None
    oref = None