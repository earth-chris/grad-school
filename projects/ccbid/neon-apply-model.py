#
# script to apply species ID models from NEON reflectance data
#####

# package imports
import aei
import copy
import plotly
import pickle
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn import metrics
from sklearn import ensemble
from sklearn import preprocessing
from sklearn import model_selection
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA

# set the path to the ECODSEdataset files
path_sep = '/'
path_base = '/home/cba/Downloads/ECODSEdataset'
path_train_hyp = path_sep.join([path_base, 'Task3', 'GroundData', 'hyper_bands_train.csv'])
path_train_bands = path_sep.join([path_base, 'Task3', 'hyper_bands.csv'])
path_train_spp = path_sep.join([path_base, 'Task3', 'GroundData', 'species_id_train.csv'])
path_test_hyp = path_sep.join([path_base, 'Task3', 'GroundData', 'hyper_bands_test.csv'])

# select the classifier to use
classifier = 'RFC'
color_options = {'RFC': 'Orange', 'GBC': 'Blue', 'SVC': 'Green'}
color = color_options[classifier]

# set paths to new output files
path_plots = path_sep.join([path_base, 'plots'])
path_scratch = path_sep.join([path_base, 'scratch'])
path_outputs = path_sep.join([path_base, 'outputs'])
path_ova = path_sep.join([path_outputs, 'one_vs_all'])
path_ovo = path_sep.join([path_outputs, 'one_vs_one'])

# load other misc. variables from the classification script
misc_var_file = path_sep.join([path_outputs, 'misc_variables.pickle'])
with open(misc_var_file, 'r') as f:
    misc_vars = pickle.load(f)
    
sp_unique = misc_vars['sp_unique']
n_species = len(sp_unique)

# load the data reducer
reducer_file = path_sep.join([path_outputs, 'Data_transformer.pickle'])
with open(reducer_file, 'r') as f:
    reducer = pickle.load(f)
    
# load the one-vs-all models and the  one-vs-one models (again, just RF)
ova_models = []
ovo_models = {}
for i in range(n_species):
    
    # first, load the one-vs-all
    rfc_ova = path_sep.join([path_ova, '{} {}.pickle'.format(sp_unique[i], classifier)])
    with open(rfc_ova, 'r') as f:
        ova_models.append(pickle.load(f))
    
    # then loop through the other species for the one-vs-one
    non_target_species = np.arange(n_species - (i + 1)) + (i + 1)
    for j in non_target_species:
        # set the path for the input files
        sp_pair = "{}-{}".format(sp_unique[i], sp_unique[j])
        rfc_ovo = path_sep.join([path_ovo, "{} {} model.pickle".format(sp_pair, classifier)])
        with open(rfc_ovo, 'r') as f:
            ovo_models[sp_pair] = pickle.load(f)

# set various options for saving / printing outputs
verbose = True
remove_outliers = True
plot_hist = True
plot_ova_cv = False
plot_ovo_cv = True
reduce_dims = True
tune_params = False
tune_ovo = False
use_transformed = True

#####
# 
# Base functions
#
#####

def two_sp_id(xdata, s1, s2):
    sp_pair = "{}-{}".format(sp_unique[s1], sp_unique[s2])
    model = ovo_models[sp_pair]
    predicted = model.predict(xdata)
    choice = [s1, s2]
    return choice[predicted[0]]

#####
#
# Reading and preprocessing data
#
#####

test_hyp = pd.read_csv(path_test_hyp)
refl = np.array(test_hyp.loc[:,'band_1':])
crid = np.array(test_hyp.loc[:,'crown_id'])
hght = np.array(test_hyp.loc[:,'chm'])

# subset the reflectance data and transform it to the original data space
gb = misc_vars['gb']
refl = refl[:, gb]

if use_transformed:
    
    # transform the input data
    transformed = reducer.transform(refl)
    
    # then, remove outliers from these predictions
    if remove_outliers:
        
        if verbose:
            print("[ STATUS ]: Removing outliers from dataset")
        
        # we'll perform a PCA, just look through the first 20 or so components, and remove
        #  data points more than 3 stdv outside the mean
        mask = np.repeat(True, len(crid))
        n_components = transformed.shape[1]
        
        # since the data was whitened to unit variance, the stdv will be just 1
        std = 1
        
        # loop through each PC and maks out values > 4x std
        for i in range(n_components):
            outliers = np.where(abs(transformed[:, i]) > 4 * std)
            mask[outliers[0]] = False
            
        # finally, cull the data
        transformed = transformed[mask]
        crid = crid[mask]
        hght = hght[mask]
        
# ok, we should have (an outlier-removed set of) transformed spectra

# create an array to store the final species IDs
other_sp_val = 2
sp_final = np.zeros(transformed.shape[0]) + other_sp_val

# create a 2-d array to classify each species one-vs-all, then loop through and perform it
sp_pred = np.zeros((len(crid), n_species))

for i in range(n_species):
    sp_pred[:, i] = ova_models[i].predict(transformed)

# get the number of unique species classified
sp_sum = sp_pred.sum(axis=1)
    
# plot the results of each classification
if plot_hist:
    xbins = np.arange(sp_sum.max()+1) - 0.5
    plt.hist(sp_sum, bins=xbins, color=color, edgecolor='black', linewidth=1)
    plt.xlabel('# of species assigned')
    plt.ylabel('Count')
    plt.title('One-vs-one classification results: {}'.format(classifier))
    plt.tight_layout()
    output_hist = path_sep.join([path_plots, 'One-vs-one species assignments {}.png'.format(classifier)])
    plt.savefig(output_hist, dpi=300)
    plt.close()
    
# ok, we'll assign the species ID'd only once to the final array
sp1 = np.where(sp_sum == 1)
for i in range(len(sp1[0])):
    ind = np.where(sp_pred[sp1[0][i], :] == 1)
    sp_final[sp1[0][i]] = ind[0][0]
    
# then, run the one-vs-ones for the two-species cases
sp2 = np.where(sp_sum == 2)
for i in range(len(sp2[0])):
    ind = np.where(sp_pred[sp2[0][i], :] == 1)
    choice = two_sp_id(np.expand_dims(transformed[sp2[0][i], :], 0), ind[0][0], ind[0][1])
    sp_final[sp2[0][i]] = choice
    
# plot out the histograms for n pixels classified per species
xbins = np.arange(n_species+1)
plt.hist(sp_final, bins=xbins -0.5, color = color, edgecolor='black', linewidth=1, orientation='horizontal')
plt.yticks(xbins, sp_unique)#, rotation = 'vertical')
plt.ylabel('Count')
plt.title('Number of pixels identified per species: {}'.format(classifier))
plt.tight_layout()
path_count = path_sep.join([path_plots, 'One-vs-one species counts {}.png'.format(classifier)])
plt.savefig(path_count, dpi=300)
plt.close()

# plot out per-crown sp ID guesses
cr_unique = np.unique(crid)
for i in range(len(cr_unique)):
    cr_ind = np.where(crid == cr_unique[i])
    
    plt.hist(sp_final[cr_ind[0]], bins=xbins -0.5, color = color, edgecolor='black', linewidth=1, orientation='horizontal')
    plt.yticks(xbins, sp_unique)#, rotation = 'vertical')
    plt.ylabel('Count')
    plt.title('Species ID: {:03d}'.format(cr_unique[i]))
    plt.tight_layout()
    
    # save the output
    path_crowns = path_sep.join([path_plots, 'crowns', 'crown-{:03d}.png'.format(cr_unique[i])])
    plt.savefig(path_crowns, dpi=300)
    plt.close()