#
# script to build species ID models from NEON reflectance data
#####

# package imports
import aei
import copy
#import plotly
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import svm
from sklearn import metrics
from sklearn import ensemble
from sklearn import multiclass
from sklearn import calibration
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
%matplotlib tk

# set the seed for reproducibility
np.random.seed(1985)

# set the path to the ECODSEdataset files
path_sep = '/'
path_base = '/home/cba/Downloads/ECODSEdataset'
path_train_hyp = path_sep.join([path_base, 'Task3', 'GroundData', 'hyper_bands_train.csv'])
path_train_bands = path_sep.join([path_base, 'Task3', 'hyper_bands.csv'])
path_train_spp = path_sep.join([path_base, 'Task3', 'GroundData', 'species_id_train.csv'])
path_test_hyp = path_sep.join([path_base, 'Task3', 'GroundData', 'hyper_bands_test.csv'])

# set paths to new output files
path_plots = path_sep.join([path_base, 'plots'])
path_scratch = path_sep.join([path_base, 'scratch'])
path_outputs = path_sep.join([path_base, 'outputs'])
path_ova = path_sep.join([path_outputs, 'one_vs_all'])
path_ovo = path_sep.join([path_outputs, 'one_vs_one'])

# set various options for saving / printing outputs
verbose = True
remove_outliers = True
plot_spectra = True
plot_all_sp = False
plot_ova_cv = True
plot_ovo_cv = True
reduce_dims = True
tune_params = False
tune_ovo = False
use_transformed = True
auto_multiclass = True
auto_species = True
auto_genus = True

##########
# 
# Base functions
#
##########

##########
# 
# Scripting
#
# Order of operations will be:
#  
#  1. Reading the data and deriving useful indices
#  2. Preprocessing the data (to remove outlier points)
#  3. Plotting the data to explore each species spectral response
#  4. Performing data normalization and dimensionality reduction
#
##########

#####
# 1. Reading data
train_hyp = pd.read_csv(path_train_hyp)
train_bands = pd.read_csv(path_train_bands)
train_spp = pd.read_csv(path_train_spp)
test_hyp = pd.read_csv(path_test_hyp)
refl = np.array(train_hyp.loc[:,'band_1':])
crid = np.array(train_hyp.loc[:,'crown_id'])
hght = np.array(train_hyp.loc[:,'chm'])

# sort all the crown IDs to align with the species IDs
spid = np.zeros(len(crid), dtype=np.int)
sp_unique = list(train_spp['species'].unique())
sp_unique.sort()
n_species = len(sp_unique)
for i in range(n_species):
    cr_unique = train_spp[train_spp['species'] == sp_unique[i]]['crown_id'].unique()
    for j in range(len(cr_unique)):
        cr_ind = np.where(crid == cr_unique[j])
        spid[cr_ind[0]] = i # now spid is labeled 0:n_species-1 to correspond with the sp_unique list

# do the same for genus ID
geid = np.zeros(len(crid), dtype=np.int)
ge_unique = list(train_spp['genus_id'].unique())
ge_unique.sort()
n_genera = len(ge_unique)
for i in range(n_genera):
    cr_unique = train_spp[train_spp['genus_id'] == ge_unique[i]]['crown_id'].unique()
    for j in range(len(cr_unique)):
        cr_ind = np.where(crid == cr_unique[j])
        geid[cr_ind[0]] = i # now geid is labeled 0:n_genera-1 to correspond with the ge_unique list

# set the unique colors
sp_colors = aei.color.color_blind(n_species)

# retrieve the wavelength data and flag additional bad bands
#  first, remove data in the blue/UV since it is usually noisy and not particularly helpful for species id
thresh_green = 0.5
wv = np.array(train_bands['Band_nanometers'])
bands = np.array(train_bands['Noise_flag'])
bands[wv < thresh_green] = 1

# set the indices for bands to use for analyses
gb = bands == 0
bb = bands == 1

# save the best bands to use

#####
# 2. Data preprocessing

# first, perform a dimensionality reduction if set
if remove_outliers:
    
    if verbose:
        print("[ STATUS ]: Removing outliers from dataset")
    
    # we'll perform a PCA, just look through the first 20 or so components, and remove
    #  data points more than 3 stdv outside the mean
    mask = np.repeat(True, len(crid))
    n_components = 20
    reducer = PCA(n_components = n_components, whiten = True)
    reducer.fit(refl[:,gb])
    
    # print some stats
    if verbose:
        print('Explained variance:')
        print('\n'.join('PC {:0.3f}: {}'.format(*k) for k in enumerate(reducer.explained_variance_)))
        print('Explained variance ratio:')
        print('\n'.join('PC {:0.3f}: {}'.format(*k) for k in enumerate(reducer.explained_variance_ratio_)))
    
    # transform the input data
    transformed = reducer.transform(refl[:,gb])
    
    # since the data was whitened to unit variance, the stdv will be just 3
    std = 1
    
    # loop through each PC and maks out values > 3x std
    for i in range(n_components):
        outliers = np.where(abs(transformed[:, i]) > 3 * std)
        mask[outliers[0]] = False
        
    # finally, cull the data
    refl = refl[mask, :]
    refl[:, bb] = 0.
    crid = crid[mask]
    spid = spid[mask]
    geid = geid[mask]
    hght = hght[mask]

#####
# 3. Plotting the data

# we'll loop through each species and plot the mean spectra +/- stdv and the number of pixels
#  but first we'll count the number of usable pixels per class
n_spec = list()
for i in range(n_species):
    sp_ind = np.where(spid == i)
    spec_mn = np.mean(refl[sp_ind[0], :], axis=0)
    spec_sd = np.std(refl[sp_ind[0], :], axis=0)
    
    spec_mn[bb] = np.nan
    
    # save the number of spectra for each class for later
    n_spec.append(len(sp_ind[0]))
    
    if plot_spectra:
        # plot the mean spectra with the stdev highlighted
        plt.figure(figsize = (3,3), dpi=200)
        plt.plot(wv, spec_mn, linewidth = 2, c = 'black')
        #plt.plot(wv, spec_mn, linewidth = 2, alpha=0.0, c = 'white', label="n = {:d}".format(len(sp_ind[0])))
        plt.fill_between(wv, spec_mn-spec_sd, spec_mn+spec_sd, 
            facecolor = sp_colors[i], alpha = 0.8, label='{}'.format(sp_unique[i]))
        plt.ylim(0, 0.6)
        
        # hide the labels
        axs = plt.gca()
        axs.axes.yaxis.set_ticklabels([])
        axs.axes.xaxis.set_ticklabels([])
                
                # label the data
                #plt.xlabel('Wavelength ({}m)'.format(r"$\mu$"))
                #plt.ylabel('Reflectance (%)')
                #plt.title('{}'.format(sp_unique[i]), style='italic')
                #plt.title('{}\nn = {:d}'.format(sp_unique[i], len(sp_ind[0])))
        plt.legend()
        plt.tight_layout()
        
        # save the figure for later
        path_fig = path_sep.join([path_plots, sp_unique[i]])
        plt.savefig(path_fig + '-no-ticks.png', dpi = 300)
        plt.savefig(path_fig + '-no-ticks.svg', dpi = 300)
        plt.close()
        
    if plot_all_sp:
        # create the figure params
        if i == 0:
            plt.figure(figsize = (3, 3), dpi=200)
        
        # plot the mean of each species in different colors
        plt.plot(wv, spec_mn, linewidth=1.5, c=sp_colors[i], alpha=0.9, label=sp_unique[i])
        
if plot_all_sp:
    # label the data
    plt.ylim(0,0.6)
    #plt.xlabel('Wavelength ({}m)'.format(r"$\mu$"))
    #plt.ylabel('Reflectance (%)')
    #plt.title('All species')
    
    # hide the labels
    axs = plt.gca()
    axs.axes.yaxis.set_ticklabels([])
    axs.axes.xaxis.set_ticklabels([])
    
    plt.text(0.95, 0.9, 'All species', horizontalalignment='right', transform=axs.transAxes)
    
    #plt.legend()
    plt.tight_layout()
        
    # save teh figure
    path_fig = path_sep.join([path_plots, 'All-species-reflectance'])
    plt.savefig(path_fig+'.png', dpi=300)
    plt.savefig(path_fig+'.svg')
    plt.close()
        
#####
# 4. Data dimensionality reduction and exploration

# Perform a full dimensionality reduction on the data with just the first 100 PCs
#  now that the outliers have been removed
if reduce_dims:
    n_components = 100
    #reducer = PCA(n_components = n_components, whiten = True)
    reducer = PCA(whiten = True)
    transformed = reducer.fit_transform(refl[:,gb])
    
    # save the reducer to apply to data later
reducer_output = path_sep.join([path_outputs, 'Data_transformer.pickle'])
with open(reducer_output, 'wb') as f:
    pickle.dump(reducer, f)
    
######
# 5a. One-vs-all model tuning

# we're going to tune a one-vs-all classification approach based on balanced 
#  training and testing data sets
svc_refl_outputs = path_sep.join([path_outputs, 'SVC_onevsall_params_refl.pickle'])
gbc_refl_outputs = path_sep.join([path_outputs, 'GBC_onevsall_params_refl.pickle'])
rfc_refl_outputs = path_sep.join([path_outputs, 'RFC_onevsall_params_refl.pickle'])
svc_pca_outputs = path_sep.join([path_outputs, 'SVC_onevsall_params_pca.pickle'])
gbc_pca_outputs = path_sep.join([path_outputs, 'GBC_onevsall_params_pca.pickle'])
rfc_pca_outputs = path_sep.join([path_outputs, 'RFC_onevsall_params_pca.pickle'])

# create variables to save the parameters and scores
best_params_svc, best_score_svc = [], []
best_params_gbc, best_score_gbc = [], []
best_params_rfc, best_score_rfc = [], []
mean_svc, stdv_svc = [], []
mean_gbc, stdv_gbc = [], []
mean_rfc, stdv_rfc = [], []
min_class = min(n_spec)

# set which data we'll be tuning / calibrating on
if use_transformed:
    target_data = transformed
else:
    target_data = refl
    
if verbose:
    print("[ STATUS ]: Beginning one-vs-all classification")

#n_pcs = [5, 10, 20, 50, 100, 200, 300, 345]
n_pcs = np.arange(10, transformed.shape[1], 10)
#n_pcs = [100]
n_itr = 50
n_samples = np.arange(50, 550, 50)

xval = len(n_pcs)
#xval = len(n_samples)

acc_gbc_list = np.zeros((xval, n_itr))
acc_gvr_list = np.zeros((xval, n_itr))
scr_gbc_list = np.zeros((xval, n_itr))
scr_gvr_list = np.zeros((xval, n_itr))

acc_rfc_list = np.zeros((xval, n_itr))
acc_rvr_list = np.zeros((xval, n_itr))
scr_rfc_list = np.zeros((xval, n_itr))
scr_rvr_list = np.zeros((xval, n_itr))

for k in range(xval):
    for j in range(n_itr):
        target_data = transformed[:,0:n_pcs[k]]
        
        if auto_multiclass:
            n_per_class = 400
            #n_per_class = n_samples[k]
            class_weight = n_species
            tuning_y = np.zeros(n_species * n_per_class, dtype=np.uint8)
            tuning_y_ge = np.zeros(n_species * n_per_class, dtype=np.uint8)
            tuning_x = np.zeros((n_species * n_per_class, target_data.shape[1]))
            sample_weight = np.ones(n_species * n_per_class)
            
            for i in range(n_species):
                ind_class = np.where(spid == i)
                ind_rnd = np.random.randint(0, high=ind_class[0].shape[0], size=n_per_class)
                tuning_x[i * n_per_class:(i+1) * n_per_class] = target_data[ind_class[0][ind_rnd]]
                tuning_y[i * n_per_class:(i+1) * n_per_class] = i
                tuning_y_ge[i * n_per_class:(i+1) * n_per_class] = geid[ind_class[0]].min()
            
            scoring = 'accuracy'
            
            # create the train/test splits
        xtrain, xtest, ytrain, ytest = model_selection.train_test_split(
            tuning_x, tuning_y, test_size=0.5, stratify=tuning_y)
        
        # and the calibration/test results
        xcalib, xctest, ycalib, yctest = model_selection.train_test_split(
            xtest, ytest, test_size=0.5, stratify=ytest)
        
        # classify at the species level
        #if auto_species:
        # gbc first
        
        # check model tuning
        if tune_params:
            model_tuning = aei.model.tune(tuning_x, tuning_y)
            model_tuning.GradientBoostClassifier(scoring='f1_weighted')
            params = model_tuning.best_params
        else:
            params = {'criterion': 'friedman_mse',
                      'learning_rate': 0.1,
                      'max_depth': 10,
                      'max_features': 'sqrt',
                      'min_impurity_split': 1e-06,
                      'min_samples_leaf': 1,
                      'min_samples_split': 0.1,
                      'min_weight_fraction_leaf': 0.0,
                      'n_estimators': 200}
        gbc = ensemble.GradientBoostingClassifier(**params)
        
        #gbc = ensemble.GradientBoostingClassifier(max_depth=None, max_features='sqrt', n_estimators=300,
        #    learning_rate=0.01)
        
        gbc.fit(xtrain, ytrain)
        ovr_gbc = calibration.CalibratedClassifierCV(gbc, method='sigmoid', cv='prefit')
        ovr_gbc.fit(xcalib, ycalib)
        prob_ovr = ovr_gbc.predict_proba(xctest)
        prob_gbc = gbc.predict_proba(xctest)
        pred_ovr = ovr_gbc.predict(xctest)
        pred_gbc = gbc.predict(xctest)
        
        score_ovr = metrics.log_loss(yctest, prob_ovr)
        score_gbc = metrics.log_loss(yctest, prob_gbc)
        acc_ovr = metrics.accuracy_score(yctest, pred_ovr)
        acc_gbc = metrics.accuracy_score(yctest, pred_gbc)
        
        acc_gbc_list[k,j] = acc_gbc
        acc_gvr_list[k,j] = acc_ovr
        scr_gbc_list[k,j] = score_gbc
        scr_gvr_list[k,j] = score_ovr
                
                #ovr_gbc = multiclass.OneVsRestClassifier(gbc, n_jobs=-2)
                #ovo_gbc = multiclass.OneVsOneClassifier(gbc, n_jobs=-2)
                #cv_score_gbc_ovr = model_selection.cross_val_score(ovr_gbc, tuning_x, tuning_y, scoring=scoring)
                #cv_score_gbc_ovo = model_selection.cross_val_score(ovo_gbc, tuning_x, tuning_y, scoring=scoring)
                #mean_gbc.append(cv_score_gbc_ovr.mean())
                #mean_gbc.append(cv_score_gbc_ovo.mean())
                #stdv_gbc.append(cv_score_gbc_ovr.std() * 2)
                #stdv_gbc.append(cv_score_gbc_ovo.std() * 2)
                
        # now, random forest
        #rfc = ensemble.RandomForestClassifier(max_depth=None, max_features='sqrt', n_estimators=300)
        if tune_params:
            model_tuning.param_grid = None
            model_tuning.RandomForestClassifier(scoring='f1_weighted')
            model_tuning = aei.model.tune(tuning_x, tuning_y)
            model_tuning.GradientBoostClassifier(scoring='f1_weighted')
            params = model_tuning.best_params
        else:
            params = {'criterion': 'gini',
                      'max_depth': None,
                      'max_features': 'sqrt',
                      'min_impurity_split': 1e-06,
                      'min_samples_leaf': 1,
                      'min_samples_split': 2,
                      'min_weight_fraction_leaf': 0.0,
                      'n_estimators': 200}
                      
        rfc = ensemble.RandomForestClassifier(**params)
        rfc.fit(xtrain, ytrain)
        ovr_rfc = calibration.CalibratedClassifierCV(rfc, method='sigmoid', cv='prefit')
        ovr_rfc.fit(xcalib, ycalib)
        prob_ovr = ovr_rfc.predict_proba(xctest)
        prob_rfc = rfc.predict_proba(xctest)
        pred_ovr = ovr_rfc.predict(xctest)
        pred_rfc = rfc.predict(xctest)
        
        score_ovr = metrics.log_loss(yctest, prob_ovr)
        score_rfc = metrics.log_loss(yctest, prob_rfc)
        acc_ovr = metrics.accuracy_score(yctest, pred_ovr)
        acc_rfc = metrics.accuracy_score(yctest, pred_rfc)
        
        acc_rfc_list[k,j] = acc_rfc
        acc_rvr_list[k,j] = acc_ovr
        scr_rfc_list[k,j] = score_rfc
        scr_rvr_list[k,j] = score_ovr

# calculate mean and stdev for metrics
acc_rvr_mean = np.mean(acc_rvr_list, axis=1)
acc_rvr_stdv = np.std(acc_rvr_list, axis=1)
acc_gvr_mean = np.mean(acc_gvr_list, axis=1)
acc_gvr_stdv = np.std(acc_gvr_list, axis=1)

scr_rvr_mean = np.mean(scr_rvr_list, axis=1)
scr_rvr_stdv = np.std(scr_rvr_list, axis=1)
scr_gvr_mean = np.mean(scr_gvr_list, axis=1)
scr_gvr_stdv = np.std(scr_gvr_list, axis=1)

# temp variables to store results from no outlier removal    
acc_rfc_nr = copy.copy(acc_rfc_list)
acc_rvr_nr = copy.copy(acc_rvr_list)
acc_gbc_nr = copy.copy(acc_gbc_list)
acc_gvr_nr = copy.copy(acc_gvr_list)

scr_rfc_nr = copy.copy(scr_rfc_list)
scr_rvr_nr = copy.copy(scr_rvr_list)
scr_gbc_nr = copy.copy(scr_gbc_list)
scr_gvr_nr = copy.copy(scr_gvr_list)
    
# lets plot some figs
plt.figure(figsize=(4,4), dpi=200)
#plt.fill_between(n_samples, acc_rvr_mean-acc_rvr_stdv, acc_rvr_mean+acc_rvr_stdv,
plt.fill_between(n_pcs, acc_rvr_mean-acc_rvr_stdv, acc_rvr_mean+acc_rvr_stdv,
    facecolor = 'Orange', alpha = 0.4)
#plt.fill_between(n_samples, acc_gvr_mean-acc_gvr_stdv, acc_gvr_mean+acc_gvr_stdv,
plt.fill_between(n_pcs, acc_gvr_mean-acc_gvr_stdv, acc_gvr_mean+acc_gvr_stdv,
    facecolor = 'Blue', alpha = 0.4)
plt.plot(n_pcs, acc_rvr_mean, label='RFC', c='Orange', linewidth=2, alpha=0.9)
plt.plot(n_pcs, acc_gvr_mean, label='GBC', c='Blue', linewidth=2, alpha=0.9)
#plt.plot(n_samples, acc_rvr_mean, label='RFC', c='Orange', linewidth=2, alpha=0.9)
#plt.plot(n_samples, acc_gvr_mean, label='GBC', c='Blue', linewidth=2, alpha=0.9)
plt.xlabel("Number of principal components used")
#plt.xlabel("Number of bootstrap samples")
plt.ylabel('Model accuracy')
plt.legend()
plt.tight_layout()
plt.savefig(path_sep.join([path_plots, 'NPCs-accuracy.svg']))
plt.savefig(path_sep.join([path_plots, 'NPCs-accuracy.png']))
#plt.savefig(path_sep.join([path_plots, 'N-Samples-accuracy.svg']))
#plt.savefig(path_sep.join([path_plots, 'N-Samples-accuracy.png']))
plt.close()

plt.figure(figsize=(4,4), dpi=200)
#plt.fill_between(n_samples, scr_rvr_mean-scr_rvr_stdv, scr_rvr_mean+scr_rvr_stdv, 
plt.fill_between(n_pcs, scr_rvr_mean-scr_rvr_stdv, scr_rvr_mean+scr_rvr_stdv, 
    facecolor = 'Orange', alpha = 0.4)
#plt.fill_between(n_samples, scr_gvr_mean-scr_gvr_stdv, scr_gvr_mean+scr_gvr_stdv, 
plt.fill_between(n_pcs, scr_gvr_mean-scr_gvr_stdv, scr_gvr_mean+scr_gvr_stdv, 
    facecolor = 'Blue', alpha = 0.4)
plt.plot(n_pcs, scr_rvr_mean, label='RFC', c='Orange', linewidth=2, alpha=0.9)
plt.plot(n_pcs, scr_gvr_mean, label='GBC', c='Blue', linewidth=2, alpha=0.9)
#plt.plot(n_samples, scr_rvr_mean, label='RFC', c='Orange', linewidth=2, alpha=0.9)
#plt.plot(n_samples, scr_gvr_mean, label='GBC', c='Blue', linewidth=2, alpha=0.9)
plt.xlabel("Number of principal components used")
#plt.xlabel("Number of bootstrap samples")
plt.ylabel('Log loss')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(path_sep.join([path_plots, 'NPCs-loss.svg']))
plt.savefig(path_sep.join([path_plots, 'NPCs-loss.png']))
#plt.savefig(path_sep.join([path_plots, 'N-Samples-loss.svg']))
#plt.savefig(path_sep.join([path_plots, 'N-Samples-loss.png']))
plt.close()

# plot feature importance
plt.figure(figsize=(4,4), dpi=200)
plt.plot(rfc.feature_importances_, c='Orange', label='RFC', linewidth=2, alpha=0.9)
plt.ylabel('RFC feature importance score')
plt.xlabel('Principal component')
        #ovr_rfc = multiclass.OneVsRestClassifier(rfc, n_jobs=-2)
        #ovo_rfc = multiclass.OneVsOneClassifier(rfc, n_jobs=-2)
        #cv_score_rfc_ovr = model_selection.cross_val_score(ovr_rfc, tuning_x, tuning_y, scoring=scoring)
        #cv_score_rfc_ovo = model_selection.cross_val_score(ovo_rfc, tuning_x, tuning_y, scoring=scoring)
        #mean_rfc.append(cv_score_rfc_ovr.mean())
        #mean_rfc.append(cv_score_rfc_ovo.mean())
        #stdv_rfc.append(cv_score_rfc_ovr.std() * 2)
        #stdv_rfc.append(cv_score_rfc_ovo.std() * 2)
        
        # fit the final models
        ovr_gbc.fit(xtest, ytest)
        #ovr_gbc.fit(tuning_x, tuning_y)
        #ovo_gbc.fit(tuning_x, tuning_y)
        ovr_rfc.fit(xtest, ytest)
        #ovr_rfc.fit(tuning_x, tuning_y)
        #ovo_rfc.fit(tuning_x, tuning_y)
        
        # then, save 'em
gbc_file = path_sep.join([path_ova, 'Multiclass-GBC.pickle'])
with open(gbc_file, 'wb') as f:
    pickle.dump(ovr_gbc, f)
    
with open(gbc_file, 'r') as f:
    ovr_gbc = pickle.load(f)
    
#gbc_file = path_sep.join([path_ovo, 'Multiclass-GBC.pickle'])
#with open(gbc_file, 'wb') as f:
#    pickle.dump(ovo_gbc, f)
   
rfc_file = path_sep.join([path_ova, 'Multiclass-RFC.pickle'])
with open(rfc_file, 'wb') as f:
    pickle.dump(ovr_rfc, f)
    
rfc_file = path_sep.join([path_ova, 'Multiclass-RFC.pickle'])
with open(rfc_file, 'r') as f:
    ovr_rfc = pickle.load(f)
    
        #rfc_file = path_sep.join([path_ovo, 'Multiclass-RFC.pickle'])
        #with open(rfc_file, 'wb') as f:
        #    pickle.dump(ovo_rfc, f)
        
    # ok, now genus level
    if auto_genus:
        
        # create the train/test splits
        xtrain, xtest, ytrain, ytest = model_selection.train_test_split(
            tuning_x, tuning_y_ge, test_size=0.5, stratify=tuning_y_ge)
        
        # and the calibration/test results
        xcalib, xtest, ycalib, ytest = model_selection.train_test_split(
            xtest, ytest, test_size=0.5, stratify=ytest)
        
        # gbc first
        model_tuning = aei.model.tune(tuning_x, tuning_y_ge)
        model_tuning.GradientBoostClassifier(scoring='f1_weighted')
        gbc = ensemble.GradientBoostingClassifier(**model_tuning.best_params)
        #gbc = ensemble.GradientBoostingClassifier(max_depth=None, max_features='sqrt', n_estimators=200,
        #    learning_rate=0.01)
        
        gbc.fit(xtrain, ytrain)
        ovr_gbc = calibration.CalibratedClassifierCV(gbc, method='sigmoid', cv='prefit')
        ovr_gbc.fit(xcalib, ycalib)
        prob_ovr = ovr_gbc.predict_proba(xtest)
        prob_gbc = gbc.predict_proba(xtest)
        pred_ovr = ovr_gbc.predict(xtest)
        pred_gbc = gbc.predict(xtest)
        
        score_ovr = metrics.log_loss(ytest, prob_ovr)
        score_gbc = metrics.log_loss(ytest, prob_gbc)
        acc_ovr = metrics.accuracy_score(ytest, pred_ovr)
        acc_gbc = metrics.accuracy_score(ytest, pred_gbc)
        
        #ovr_gbc = multiclass.OneVsRestClassifier(gbc, n_jobs=-2)
        #ovo_gbc = multiclass.OneVsOneClassifier(gbc, n_jobs=-2)
        #cv_score_gbc_ovr = model_selection.cross_val_score(ovr_gbc, tuning_x, tuning_y_ge, scoring=scoring)
        #cv_score_gbc_ovo = model_selection.cross_val_score(ovo_gbc, tuning_x, tuning_y_ge, scoring=scoring)
        #mean_gbc.append(cv_score_gbc_ovr.mean())
        #mean_gbc.append(cv_score_gbc_ovo.mean())
        #stdv_gbc.append(cv_score_gbc_ovr.std() * 2)
        #stdv_gbc.append(cv_score_gbc_ovo.std() * 2)
        
        model_tuning.param_grid = None
        model_tuning.RandomForestClassifier(scoring='f1_weighted')
        rfc = ensemble.RandomForestClassifier(**model_tuning.best_params)
        #rfc = ensemble.RandomForestClassifier(**best_params_rfc[i])
        #rfc = ensemble.RandomForestClassifier(max_depth=None, max_features='sqrt', n_estimators=200)
        
        rfc.fit(xtrain, ytrain)
        ovr_rfc = calibration.CalibratedClassifierCV(rfc, method='sigmoid', cv='prefit')
        ovr_rfc.fit(xcalib, ycalib)
        prob_ovr = ovr_rfc.predict_proba(xtest)
        prob_rfc = rfc.predict_proba(xtest)
        pred_ovr = ovr_rfc.predict(xtest)
        pred_rfc = rfc.predict(xtest)
        
        score_ovr = metrics.log_loss(ytest, prob_ovr)
        score_gbc = metrics.log_loss(ytest, prob_rfc)
        acc_ovr = metrics.accuracy_score(ytest, pred_ovr)
        acc_gbc = metrics.accuracy_score(ytest, pred_rfc)
        
        #ovr_rfc = multiclass.OneVsRestClassifier(rfc, n_jobs=-2)
        #ovo_rfc = multiclass.OneVsOneClassifier(rfc, n_jobs=-2)
        #cv_score_rfc_ovr = model_selection.cross_val_score(ovr_rfc, tuning_x, tuning_y_ge, scoring=scoring)
        #cv_score_rfc_ovo = model_selection.cross_val_score(ovo_rfc, tuning_x, tuning_y_ge, scoring=scoring)
        #mean_rfc.append(cv_score_rfc_ovr.mean())
        #mean_rfc.append(cv_score_rfc_ovo.mean())
        #stdv_rfc.append(cv_score_rfc_ovr.std() * 2)
        #stdv_rfc.append(cv_score_rfc_ovo.std() * 2)
        
        # fit the final models
        ovr_gbc.fit(tuning_x, tuning_y_ge)
        #ovo_gbc.fit(tuning_x, tuning_y_ge)
        ovr_rfc.fit(tuning_x, tuning_y_ge)
        #ovo_rfc.fit(tuning_x, tuning_y_ge)
        
        # then, save 'em
        gbc_file = path_sep.join([path_ova, 'Multiclass-GBC-genus.pickle'])
        with open(gbc_file, 'wb') as f:
            pickle.dump(ovr_gbc, f)
            
        #gbc_file = path_sep.join([path_ovo, 'Multiclass-GBC-genus.pickle'])
        #with open(gbc_file, 'wb') as f:
        #    pickle.dump(ovo_gbc, f)
           
        rfc_file = path_sep.join([path_ova, 'Multiclass-RFC-genus.pickle'])
        with open(rfc_file, 'wb') as f:
            pickle.dump(ovr_rfc, f)
            
        #rfc_file = path_sep.join([path_ovo, 'Multiclass-RFC-genus.pickle'])
        #with open(rfc_file, 'wb') as f:
        #    pickle.dump(ovo_rfc, f)
        
else:
    # loop through and classify as each species vs all    
    for i in range(n_species):
        if verbose:
            print("[ STATUS ]: Classifying species: {}".format(sp_unique[i]))
        ind_target = np.where(spid == i)
        ind_backgr = np.where(spid != i)
        
        #n_per_class = np.int(np.ceil(len(ind_target[0]) / (n_species - 1)))
        
        # standardize the minimum number of samples per class to pull
        #n_backgr = min(min_class, n_per_class) * (n_species - 1)
        #n_target = min(len(ind_target[0]), (n_species - 1) * n_backgr)
        
        # create the arrays to store the randomly sampled classes
        #tuning_x = np.zeros((n_backgr + n_target, n_components))
        #tuning_y = np.zeros(n_backgr + n_target, dtype = np.uint8)
        #tuning_y[n_backgr:] = 1
        
        ###
        # NEW SAMPLING - WEIGHTED CLASSES
        ###
        #n_target = len(ind_target[0])
        #n_backgr = len(ind_backgr[0])
        #class_weight = float(n_backgr) / n_target
        #sample_weight = np.ones(n_target+n_backgr)
        #sample_weight[ind_target[0]] = class_weight
        #fit_params = {'sample_weight': sample_weight}
        
        #tuning_y = np.zeros(n_backgr+n_target, dtype=np.uint8)
        #tuning_y[ind_target[0]] = 1
        #tuning_x = target_data
        
        ###
        # OK, NEW SAMPLING AGAIN - set a fixed number per class 
        ###
        n_per_class = 400
        class_weight = n_species
        tuning_y = np.zeros(n_species * n_per_class, dtype=np.uint8)
        tuning_x = np.zeros((n_species * n_per_class, n_components))
        tuning_y[-n_per_class:] = 1
        
        sample_weight = np.ones(n_species * n_per_class)
        sample_weight[-n_per_class:] = class_weight
        fit_params = {'sample_weight': sample_weight}
        
        # loop through each non-target species, randomly sample, and add to the tuning set
        non_target_species = range(n_species)
        non_target_species.pop(i)
        for j in range(len(non_target_species)):
            
            # get indices for the class, and pull a random subset of those
            ind_class = np.where(spid == non_target_species[j])
        #    ind_rnd = ind_class[0][np.random.choice(len(ind_class[0]), n_per_class)]
        #    tuning_x[j * n_per_class:(j + 1) * n_per_class] = target_data[ind_rnd]
            ind_rnd = np.random.randint(0, high=ind_class[0].shape[0], size=n_per_class)
            tuning_x[j * n_per_class:(j + 1) * n_per_class] = target_data[ind_class[0][ind_rnd]]
        
        # then, randomly sample the target class
        #tuning_x[n_backgr:] = target_data[ind_target[0][np.random.choice(len(ind_target[0]), n_target)]]
        ind_rnd = np.random.randint(0, high=ind_target[0].shape[0], size=n_per_class)
        tuning_x[-n_per_class:] = target_data[ind_target[0][ind_rnd]]
        
        # now we should have populated tuning x/y variables. now tune the models for each species
        #  using multiple classifiers
        
        if tune_params:
            print("[ STATUS ]: Tuning model parameters")
            # first, SVC
            model_tuning = aei.model.tune(tuning_x, tuning_y)
           # model_tuning.SVC(scoring='f1_weighted', class_weight={0: 1, 1: class_weight})
            model_tuning.SVC(scoring='f1_weighted')
            best_params_svc.append(model_tuning.best_params)
            best_score_svc.append(model_tuning.best_score)
            
            # second, gradient boosting
            model_tuning.param_grid = None
            #model_tuning.GradientBoostClassifier(scoring='f1_weighted', fit_params=fit_params)
            model_tuning.GradientBoostClassifier(scoring='f1_weighted')
            best_params_gbc.append(model_tuning.best_params)
            best_score_gbc.append(model_tuning.best_score)
            
            # finally, random forest
            model_tuning.param_grid = None
            #model_tuning.RandomForestClassifier(scoring='f1_weighted', class_weight={0: 1, 1: class_weight})
            model_tuning.RandomForestClassifier(scoring='f1_weighted')
            best_params_rfc.append(model_tuning.best_params)
            best_score_rfc.append(model_tuning.best_score)
            
            # save the outputs
            with open(svc_pca_outputs, 'wb') as f:
                pickle.dump(best_params_svc, f)
            with open(gbc_pca_outputs, 'wb') as f:
                pickle.dump(best_params_gbc, f)
            with open(rfc_pca_outputs, 'wb') as f:
                pickle.dump(best_params_rfc, f)
                    
        else:
            with open(svc_refl_outputs, 'rb') as f:
                best_params_svc.append(pickle.load(f))
            with open(gbc_refl_outputs, 'rb') as f:
                best_params_gbc.append(pickle.load(f))
            with open(rfc_refl_outputs, 'rb') as f:
                best_params_rfc.append(pickle.load(f))
                    
        # ok, now that we have the best parameters for each model type, let's perform some
        #  cross validation to assess model performance
        scoring = 'f1_weighted'
        
        if verbose:
            print("[ STATUS ]: Performing cross validation")
        
        # first, svm
        svc = svm.SVC(class_weight={0: 1, 1: class_weight}, **best_params_svc[i])
        cv_score_svc = model_selection.cross_val_score(svc, tuning_x, tuning_y, scoring = scoring)
        mean_svc.append(cv_score_svc.mean())
        stdv_svc.append(cv_score_svc.std() * 2)
        
        # next, gradient boosting
        gbc = ensemble.GradientBoostingClassifier(**best_params_gbc[i])
        cv_score_gbc = model_selection.cross_val_score(gbc, tuning_x, tuning_y, scoring = scoring)
        mean_gbc.append(cv_score_gbc.mean())
        stdv_gbc.append(cv_score_gbc.std() * 2)
        
        # finally, random forest
        rfc = ensemble.RandomForestClassifier(class_weight={0: 1, 1: class_weight}, **best_params_rfc[i])
        cv_score_rfc = model_selection.cross_val_score(rfc, tuning_x, tuning_y, scoring = scoring)
        mean_rfc.append(cv_score_rfc.mean())
        stdv_rfc.append(cv_score_rfc.std() * 2)
        
        if verbose:
            print("[ STATUS ]: Fitting and saving final models")
        
        # train and save the final models using all data
        svc.fit(tuning_x, tuning_y)
        svc_file = path_sep.join([path_ova, sp_unique[i] + ' SVC.pickle'])
        with open(svc_file, 'wb') as f:
            pickle.dump(svc, f)
        
        gbc.fit(tuning_x, tuning_y, sample_weight=sample_weight)
        gbc_file = path_sep.join([path_ova, sp_unique[i] + ' GBC.pickle'])
        with open(gbc_file, 'wb') as f:
            pickle.dump(gbc, f)
            
        rfc.fit(tuning_x, tuning_y)
        rfc_file = path_sep.join([path_ova, sp_unique[i] + ' RFC.pickle'])
        with open(rfc_file, 'wb') as f:
            pickle.dump(rfc, f)
            
        if verbose:
            print("[ STATUS ]: ----------")
        
# ok, now that we have performed the cross validation for each model, plot out the cv results
if plot_ova_cv:
    plt.figure()
    plt_x = np.arange(len(sp_unique))
    plt.errorbar(plt_x-0.15, mean_svc, yerr = stdv_svc, fmt = 'o', label = 'SVC')
    plt.errorbar(plt_x, mean_gbc, yerr = stdv_gbc, fmt = 'o', label = 'GBC')
    plt.errorbar(plt_x+0.15, mean_rfc, yerr = stdv_rfc, fmt = 'o', label = 'RFC')
    plt.ylabel(scoring)
    plt.xticks(plt_x, sp_unique, rotation = 'vertical')
    plt.title('One-vs-all cross validation performance')
    plt.legend()
    plt.tight_layout()
    
    # save the figure
    path_fig = path_sep.join([path_plots, 'Cross validation results.png'])
    plt.savefig(path_fig, dpi=300)
    plt.close()

######
# 5b. One-vs-one model tuning

# create dictionaries to store the tuning parameters, and 2d arrays
#  to store the scores
best_params_svc, best_params_gbc, best_params_rfc = {}, {}, {}
best_score_svc = np.zeros((n_species, n_species), dtype = np.float)
best_score_gbc = copy.copy(best_score_svc)
best_score_rfc = copy.copy(best_score_svc)
mean_svc = copy.copy(best_score_svc)
stdv_svc = copy.copy(best_score_svc)
mean_gbc = copy.copy(best_score_svc)
stdv_gbc = copy.copy(best_score_svc)
mean_rfc = copy.copy(best_score_svc)
stdv_rfc = copy.copy(best_score_svc)

if verbose:
    print("[ STATUS ]: Beginning one-vs-one classification")

# loop through and classify as each species vs all    
for i in range(n_species):
    
    # since we don't need to repeat classifications (e.g., sp 1 vs 2, then 2 vs 1),
    #  remove all sp. of (i) that have already done the one vs one.
    #non_target_species = range(n_species)
    #for k in range(i+1):
    #    non_target_species.pop(k)
        
    non_target_species = np.arange(n_species - (i + 1)) + (i + 1)
        
    # ok, now loop through all of the novel one vs one combinations
    for j in non_target_species:
        if verbose:
            print("[ STATUS ]: Classifying species: {} vs. {}".format(sp_unique[i], sp_unique[j]))
        
        # we'll balance the data sets based on the lowest-sampled class
        ind_target = np.where(spid == i)
        ind_backgr = np.where(spid == j)
        
        ###
        # NEW SAMPLING - WEIGHTED CLASSES
        ###
        n_target = len(ind_target[0])
        n_backgr = len(ind_backgr[0])
        class_weight = float(n_backgr) / n_target
        sample_weight = np.ones(n_target+n_backgr)
        sample_weight[n_backgr:] = class_weight
        fit_params = {'sample_weight': sample_weight}
        
        tuning_y = np.zeros(n_backgr+n_target, dtype=np.uint8)
        tuning_y[n_backgr:] = 1
        tuning_x = np.zeros((n_backgr+n_target, target_data.shape[1]))
        tuning_x[0:n_backgr] = target_data[ind_backgr[0]]
        tuning_x[n_backgr:] = target_data[ind_target[0]]
        
        #n_per_class = np.min((len(ind_target[0]), len(ind_backgr[0])))
        
        # create the arrays to store the randomly sampled classes
        #tuning_x = np.zeros((n_per_class * 2, target_data.shape[1]))
        #tuning_y = np.zeros(n_per_class * 2, dtype = np.uint8)
        #tuning_y[n_per_class:] = 1
        
        # sample the background class first
        #rnd_backgr = ind_backgr[0][np.random.choice(len(ind_backgr[0]), n_per_class)]
        #tuning_x[0:n_per_class] = target_data[rnd_backgr]
        
        # then, the target class
        #rnd_target = ind_target[0][np.random.choice(len(ind_target[0]), n_per_class)]
        #tuning_x[n_per_class:] = target_data[rnd_target]
        
        # set the path for the output files
        sp_pair = "{}-{}".format(sp_unique[i], sp_unique[j])
        svc_pca_outputs = path_sep.join([path_ovo, "{} SVC params.pickle".format(sp_pair)])
        gbc_pca_outputs = path_sep.join([path_ovo, "{} GBC params.pickle".format(sp_pair)])
        rfc_pca_outputs = path_sep.join([path_ovo, "{} RFC params.pickle".format(sp_pair)])
        
        # now that we've got each class sampled, we'll tune a bunch of parameters
        if tune_ovo:
            if verbose:
                print("[ STATUS ]: Tuning model parameters")
                
            # first, SVC
            model_tuning = aei.model.tune(tuning_x, tuning_y)
            model_tuning.SVC(scoring='f1_weighted', class_weight={0: 1, 1: class_weight})
            model_tuning.SVC(scoring='f1_weighted')
            best_params_svc[sp_pair] = model_tuning.best_params
            best_score_svc[i, j] = model_tuning.best_score
            with open(svc_pca_outputs, 'wb') as f:
                pickle.dump(model_tuning.best_params, f)
            
            # second, gradient boosting
            model_tuning.param_grid = None
            model_tuning.GradientBoostClassifier(scoring='f1_weighted')
            best_params_gbc[sp_pair] = model_tuning.best_params
            best_score_gbc[i, j] = model_tuning.best_score
            with open(gbc_pca_outputs, 'wb') as f:
                pickle.dump(model_tuning.best_params, f)
            
            # finally, random forest
            model_tuning.param_grid = None
            model_tuning.RandomForestClassifier(scoring='f1_weighted', class_weight={0: 1, 1: class_weight})
            model_tuning.RandomForestClassifier(scoring='f1_weighted')
            best_params_rfc[sp_pair] = model_tuning.best_params
            best_score_rfc[i, j] = model_tuning.best_score
            with open(rfc_pca_outputs, 'wb') as f:
                pickle.dump(model_tuning.best_params, f)
                    
        else:
            with open(svc_pca_outputs, 'rb') as f:
                best_params_svc[sp_pair] = pickle.load(f)
            with open(gbc_pca_outputs, 'rb') as f:
                best_params_gbc[sp_pair] = pickle.load(f)
            with open(rfc_pca_outputs, 'rb') as f:
                best_params_rfc[sp_pair] = pickle.load(f)
                    
        # ok, now that we have the best parameters for each model type, let's perform some
        #  cross validation to assess model performance
        scoring = 'f1_weighted'
        
        if verbose:
            print("[ STATUS ]: Performing cross validation")
        
        # first, svm
        svc = svm.SVC(class_weight={0: 1, 1: class_weight},**best_params_svc[sp_pair])
        cv_score_svc = model_selection.cross_val_score(svc, tuning_x, tuning_y, scoring = scoring)
        mean_svc[i,j] = cv_score_svc.mean()
        stdv_svc[i,j] = cv_score_svc.std() * 2
        
        # next, gradient boosting
        gbc = ensemble.GradientBoostingClassifier(**best_params_gbc[sp_pair])
        cv_score_gbc = model_selection.cross_val_score(gbc, tuning_x, tuning_y, scoring = scoring)
        mean_gbc[i,j] = cv_score_gbc.mean()
        stdv_gbc[i,j] = cv_score_gbc.std() * 2
        
        # finally, random forest
        rfc = ensemble.RandomForestClassifier(class_weight={0: 1, 1: class_weight},**best_params_rfc[sp_pair])
        cv_score_rfc = model_selection.cross_val_score(rfc, tuning_x, tuning_y, scoring = scoring)
        mean_rfc[i,j] = cv_score_rfc.mean()
        stdv_rfc[i,j] = cv_score_rfc.std() * 2
        
        if verbose:
            print("[ STATUS ]: Fitting and saving final models")
        
        # train and save the final models using all data
        svc.fit(tuning_x, tuning_y)
        svc_file = path_sep.join([path_ovo, '{} SVC model.pickle'.format(sp_pair)])
        with open(svc_file, 'wb') as f:
            pickle.dump(svc, f)
        
        gbc.fit(tuning_x, tuning_y, sample_weight=sample_weight)
        gbc_file = path_sep.join([path_ovo, '{} GBC model.pickle'.format(sp_pair)])
        with open(gbc_file, 'wb') as f:
            pickle.dump(gbc, f)
            
        rfc.fit(tuning_x, tuning_y)
        rfc_file = path_sep.join([path_ovo, '{} RFC model.pickle'.format(sp_pair)])
        with open(rfc_file, 'wb') as f:
            pickle.dump(rfc, f)
            
        if verbose:
            print("[ STATUS ]: ----------")
                
# alright, next steps are going to be plotting out the results of the one vs one confusion matrix
if plot_ovo_cv:
    mask = mean_rfc == 0
    cv_plot_svc = sns.heatmap(mean_svc, xticklabels=sp_unique, 
        yticklabels=sp_unique, linewidths=0.3, annot=True, cmap="YlGnBu", 
        mask=mask, vmin = 0.75, vmax = 1.0)
    plt.title('SVC one-vs-one cross-val scores: {}'.format(scoring))
    plt.tight_layout()
    path_fig = path_sep.join([path_plots, 'Cross validation one-vs-one SVC.png'])
    plt.savefig(path_fig, dpi=300)
    plt.close()
    
    cv_plot_svc = sns.heatmap(mean_gbc, xticklabels=sp_unique, 
        yticklabels=sp_unique, linewidths=0.3, annot=True, cmap="YlGnBu", 
        mask=mask, vmin = 0.75, vmax = 1.0)
    plt.title('GBC one-vs-one cross-val scores: {}'.format(scoring))
    plt.tight_layout()
    path_fig = path_sep.join([path_plots, 'Cross validation one-vs-one GBC.png'])
    plt.savefig(path_fig, dpi=300)
    plt.close()
    
    cv_plot_svc = sns.heatmap(mean_rfc, xticklabels=sp_unique, 
        yticklabels=sp_unique, linewidths=0.3, annot=True, cmap="YlGnBu", 
        mask=mask, vmin = 0.75, vmax = 1.0)
    plt.title('RFC one-vs-one cross-val scores: {}'.format(scoring))
    plt.tight_layout()
    path_fig = path_sep.join([path_plots, 'Cross validation one-vs-one RFC.png'])
    plt.savefig(path_fig, dpi=300)
    plt.close()
    
# and save all the other weird paramaters to apply the model outputs
output_vars = {}
output_vars['gb'] = gb
output_vars['path_ovo'] = path_ovo
output_vars['path_ova'] = path_ova
output_vars['sp_unique'] = sp_unique
output_vars['ge_unique'] = ge_unique

output_var_file = path_sep.join([path_outputs, 'misc_variables.pickle'])
with open(output_var_file, 'wb') as f:
    pickle.dump(output_vars, f)
    