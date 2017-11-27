#
#
#####

# package imports
import aei
import copy
import pickle
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn import metrics
from sklearn import ensemble
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

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
plot_spectra = False
reduce_dims = True
tune_params = False
tune_ovo = True
use_transformed = True

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
    
    # save the number of spectra for each class for later
    n_spec.append(len(sp_ind[0]))
    
    if plot_spectra:
        # plot the mean spectra with the stdev highlighted
        plt.plot(wv, spec_mn, linewidth = 2, c = 'black')
        plt.fill_between(wv, spec_mn-spec_sd, spec_mn+spec_sd, facecolor = sp_colors[i], alpha = 0.8)
        plt.ylim(0, 0.6)
        
        # label the data
        plt.xlabel('Wavelength (um)')
        plt.ylabel('Reflectance (%)')
        plt.title('{}  |  n = {:d}'.format(sp_unique[i], len(sp_ind[0])))
        
        # save the figure for later
        path_fig = path_sep.join([path_plots, sp_unique[i]])
        plt.savefig(path_fig + '.png', dpi = 300)
        plt.close()
        
#####
# 4. Data dimensionality reduction and exploration

# Perform a full dimensionality reduction on the data with just the first 100 PCs
#  now that the outliers have been removed
if reduce_dims:
    n_components = 100
    reducer = PCA(n_components = n_components, whiten = True)
    transformed = reducer.fit_transform(refl[:,gb])
    
    # save the reducer to apply to data later
    reducer_output = path_sep.join([path_outputs, 'Data_transformer.pickle'])
    
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

# loop through and classify as each species vs all    
for i in range(n_species):
    if verbose:
        print("[ STATUS ]: Classifying species: {}".format(sp_unique[i]))
    ind_target = np.where(spid == i)
    ind_backgr = np.where(spid != i)
    
    n_per_class = np.int(np.ceil(len(ind_target[0]) / (n_species - 1)))
    
    # standardize the minimum number of samples per class to pull
    n_backgr = min(min_class, n_per_class) * (n_species - 1)
    n_target = min(len(ind_target[0]), (n_species - 1) * n_backgr)
    
    # create the arrays to store the randomly sampled classes
    tuning_x = np.zeros((n_backgr + n_target, n_components))
    tuning_y = np.zeros(n_backgr + n_target, dtype = np.uint8)
    tuning_y[n_backgr:] = 1
    
    # loop through each non-target species, randomly sample, and add to the tuning set
    non_target_species = range(n_species)
    non_target_species.pop(i)
    for j in range(len(non_target_species)):
        
        # get indices for the class, and pull a random subset of those
        ind_class = np.where(spid == non_target_species[j])
        ind_rnd = ind_class[0][np.random.choice(len(ind_class[0]), n_per_class)]
        tuning_x[j * n_per_class:(j + 1) * n_per_class] = target_data[ind_rnd]
    
    # then, randomly sample the target class
    tuning_x[n_backgr:] = target_data[ind_target[0][np.random.choice(len(ind_target[0]), n_target)]]

    # now we should have populated tuning x/y variables. now tune the models for each species
    #  using multiple classifiers
    
    if tune_params:
        print("[ STATUS ]: Tuning model parameters")
        # first, SVC
        model_tuning = aei.model.tune(tuning_x, tuning_y)
        model_tuning.SVC()
        best_params_svc.append(model_tuning.best_params)
        best_score_svc.append(model_tuning.best_score)
        
        # second, gradient boosting
        model_tuning.param_grid = None
        model_tuning.GradientBoostClassifier()
        best_params_gbc.append(model_tuning.best_params)
        best_score_gbc.append(model_tuning.best_score)
        
        # finally, random forest
        model_tuning.param_grid = None
        model_tuning.RandomForestClassifier()
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
            best_params_svc = pickle.load(f)
        with open(gbc_refl_outputs, 'rb') as f:
            best_params_gbc = pickle.load(f)
        with open(rfc_refl_outputs, 'rb') as f:
            best_params_rfc = pickle.load(f)
                
    # ok, now that we have the best parameters for each model type, let's perform some
    #  cross validation to assess model performance
    scoring = 'accuracy'
    
    if verbose:
        print("[ STATUS ]: Performing cross validation")
    
    # first, svm
    svc = svm.SVC(**best_params_svc[i])
    cv_score_svc = model_selection.cross_val_score(svc, tuning_x, tuning_y, scoring = scoring)
    mean_svc.append(cv_score_svc.mean())
    stdv_svc.append(cv_score_svc.std() * 2)
    
    # next, gradient boosting
    gbc = ensemble.GradientBoostingClassifier(**best_params_gbc[i])
    cv_score_gbc = model_selection.cross_val_score(gbc, tuning_x, tuning_y, scoring = scoring)
    mean_gbc.append(cv_score_gbc.mean())
    stdv_gbc.append(cv_score_gbc.std() * 2)
    
    # finally, random forest
    rfc = ensemble.RandomForestClassifier(**best_params_rfc[i])
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
    
    gbc.fit(tuning_x, tuning_y)
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
        
        n_per_class = np.min((len(ind_target[0]), len(ind_backgr[0])))
        
        # create the arrays to store the randomly sampled classes
        tuning_x = np.zeros((n_per_class * 2, target_data.shape[1]))
        tuning_y = np.zeros(n_per_class * 2, dtype = np.uint8)
        tuning_y[n_per_class:] = 1
        
        # sample the background class first
        rnd_backgr = ind_backgr[0][np.random.choice(len(ind_backgr[0]), n_per_class)]
        tuning_x[0:n_per_class] = target_data[rnd_backgr]
        
        # then, the target class
        rnd_target = ind_target[0][np.random.choice(len(ind_target[0]), n_per_class)]
        tuning_x[n_per_class:] = target_data[rnd_target]
        
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
            model_tuning.SVC()
            best_params_svc[sp_pair] = model_tuning.best_params
            best_score_svc[i, j] = model_tuning.best_score
            with open(svc_pca_outputs, 'wb') as f:
                pickle.dump(model_tuning.best_params, f)
            
            # second, gradient boosting
            model_tuning.param_grid = None
            model_tuning.GradientBoostClassifier()
            best_params_gbc[sp_pair] = model_tuning.best_params
            best_score_gbc[i, j] = model_tuning.best_score
            with open(gbc_pca_outputs, 'wb') as f:
                pickle.dump(model_tuning.best_params, f)
            
            # finally, random forest
            model_tuning.param_grid = None
            model_tuning.RandomForestClassifier()
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
        scoring = 'accuracy'
        
        if verbose:
            print("[ STATUS ]: Performing cross validation")
        
        # first, svm
        svc = svm.SVC(**best_params_svc[sp_pair])
        cv_score_svc = model_selection.cross_val_score(svc, tuning_x, tuning_y, scoring = scoring)
        mean_svc[i,j] = cv_score_svc.mean()
        stdv_svc[i,j] = cv_score_svc.std() * 2
        
        # next, gradient boosting
        gbc = ensemble.GradientBoostingClassifier(**best_params_gbc[sp_pair])
        cv_score_gbc = model_selection.cross_val_score(gbc, tuning_x, tuning_y, scoring = scoring)
        mean_gbc[i,j] = cv_score_gbc.mean()
        stdv_gbc[i,j] = cv_score_gbc.std() * 2
        
        # finally, random forest
        rfc = ensemble.RandomForestClassifier(**best_params_rfc[sp_pair])
        cv_score_rfc = model_selection.cross_val_score(rfc, tuning_x, tuning_y, scoring = scoring)
        mean_rfc[i,j] = cv_score_rfc.mean()
        stdv_rfc[i,j] = cv_score_rfc.std() * 2
        
        if verbose:
            print("[ STATUS ]: Fitting and saving final models")
        
        # train and save the final models using all data
        svc.fit(tuning_x, tuning_y)
        svc_file = path_sep.join([path_ova, '{} SVC model.pickle'.format(sp_pair)])
        with open(svc_file, 'wb') as f:
            pickle.dump(svc, f)
        
        gbc.fit(tuning_x, tuning_y)
        gbc_file = path_sep.join([path_ova, '{} GBC model.pickle'.format(sp_pair)])
        with open(gbc_file, 'wb') as f:
            pickle.dump(gbc, f)
            
        rfc.fit(tuning_x, tuning_y)
        rfc_file = path_sep.join([path_ova, '{} RFC model.pickle'.format(sp_pair)])
        with open(rfc_file, 'wb') as f:
            pickle.dump(rfc, f)
            
        if verbose:
            print("[ STATUS ]: ----------")
                
# alright, next steps are going to be plotting out the results of the one vs one confusion matrix