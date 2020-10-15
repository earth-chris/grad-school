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
from matplotlib import lines
%matplotlib tk

# set the paths to input files
path_sep = '/'
path_base = '/home/cba/Downloads/ECODSEdataset'
path_class_test = path_sep.join([path_base, 'class-metrics-testing.csv'])
path_class_train = path_sep.join([path_base, 'class-metrics-training.csv'])
path_class_prob = path_sep.join([path_base, 'class-metrics-testing-prop.csv'])

# set paths to new output files
path_plots = path_sep.join([path_base, 'plots'])

#####
# set some helper functions
# function to create legend proxies
def create_proxy(color, marker, linestyle='none'):
    line = lines.Line2D([0], [0], linestyle=linestyle, mfc=color,
        mec='black', marker=marker)
    return line

# function to label new colors
def label_color(row, colorby, unique_vals, colors):
    # get the index for the value passed and return the
    #  associated color
    ind = unique_vals.index(row[colorby])
    return colors[ind]
    
#####
# read the input data
class_test = pd.read_csv(path_class_test)
class_train = pd.read_csv(path_class_train)
class_prob = pd.read_csv(path_class_prob)

# set the color scheme
n_species = len(class_test)
sp_colors = aei.color.color_blind(n_species)
met_colors = aei.color.color_blind()
met_cols = [met_colors[2], met_colors[5], met_colors[3], met_colors[1]]
metrics = ['Accuracy', 'Specificity', 'Precision', 'Recall']

# plot a stacked bar chart of accuracy and specificity
plt.figure(figsize=(8,5), dpi=200)
width = 0.2
ind = np.arange(n_species)

plt.bar(ind - width/2 - width, class_test[metrics[0]], width, color=met_cols[0], alpha=0.9, edgecolor='black', label=metrics[0])
plt.bar(ind - width/2, class_test[metrics[1]], width, color=met_cols[1], alpha=0.9, edgecolor='black', label=metrics[1])
plt.bar(ind + width/2, class_test[metrics[2]], width, color=met_cols[2], alpha=0.9, edgecolor='black', label=metrics[2])
plt.bar(ind + width/2 + width, class_test[metrics[3]], width, color=met_cols[3], alpha=0.9, edgecolor='black', label=metrics[3])
#plt.bar(ind + width/2, class_test['Specificity'], width, color=sp_colors, alpha=0.5, edgecolor='black')

# set the labels
plt.xticks(ind, class_test['Species'], rotation='vertical', fontstyle='italic')
plt.ylabel('Score')
plt.title('Model performance on testing data')
plt.tight_layout()

# custom build the legend
#legend = ['Accuracy', 'Specificity']
#legend_colors = [(0.3, 0.3, 0.3, 0.9), (0.5, 0.5, 0.5, 0.5)]
#legend_marker = ['s', 's']

#proxies = []
#for i in range(len(legend)):
#    proxies.append(create_proxy(legend_colors[i], legend_marker[i]))#, legend_linestyle[i]))
    
plt.legend(ncol = 2, loc='upper center', bbox_to_anchor=(0.45, -0.45), 
    markerscale = 1.5)
    #-0.65
# save the figure
plt.savefig(path_plots + '/all-metrics-test-bin-no-ital.svg')
plt.savefig(path_plots + '/all-metrics-test-bin.png')
plt.savefig(path_plots + '/all-metrics-test-bin.pdf')
plt.close()

# do the same for the training data
plt.figure(figsize=(8,5), dpi=200)
width = 0.22
ind = np.arange(n_species)

plt.bar(ind - width/2 - width, class_train[metrics[0]], width, color=met_cols[0], alpha=0.9, edgecolor='black', 
    label='{}'.format(metrics[0], style='italics'))
plt.bar(ind - width/2, class_train[metrics[1]], width, color=met_cols[1], alpha=0.9, edgecolor='black', 
    label='{}'.format(metrics[1], style='italics'))
plt.bar(ind + width/2, class_train[metrics[2]], width, color=met_cols[2], alpha=0.9, edgecolor='black', 
    label='{}'.format(metrics[2], style='italics'))
plt.bar(ind + width/2 + width, class_train[metrics[3]], width, color=met_cols[3], alpha=0.9, edgecolor='black', 
    label='{}'.format(metrics[3], style='italics'))
#plt.bar(ind + width/2, class_test['Specificity'], width, color=sp_colors, alpha=0.5, edgecolor='black')

# set the labels
plt.xticks(ind, class_test['Species'], rotation='vertical', fontstyle='italic')
plt.ylabel('Score')
plt.title('Model performance on training data')
plt.tight_layout()
plt.legend(ncol = 2, loc='upper center', bbox_to_anchor=(0.45, -0.45),
    markerscale = 1.5)
    
# save the figure
plt.savefig(path_plots + '/all-metrics-train.svg')
plt.savefig(path_plots + '/all-metrics-train.png')
plt.savefig(path_plots + '/all-metrics-train.pdf')
plt.close()


# finally, for probability samples
plt.figure(figsize=(8,5), dpi=200)
width = 0.2
ind = np.arange(n_species)

plt.bar(ind - width/2 - width, class_prob[metrics[0]], width, color=met_cols[0], alpha=0.9, edgecolor='black', label=metrics[0])
plt.bar(ind - width/2, class_prob[metrics[1]], width, color=met_cols[1], alpha=0.9, edgecolor='black', label=metrics[1])
plt.bar(ind + width/2, class_prob[metrics[2]], width, color=met_cols[2], alpha=0.9, edgecolor='black', label=metrics[2])
plt.bar(ind + width/2 + width, class_prob[metrics[3]], width, color=met_cols[3], alpha=0.9, edgecolor='black', label=metrics[3])
#plt.bar(ind + width/2, class_test['Specificity'], width, color=sp_colors, alpha=0.5, edgecolor='black')

# set the labels
plt.xticks(ind, class_prob['Species'], rotation='vertical', fontstyle='italic')
plt.ylabel('Score')
plt.title('Model performance on test data using sample probabilities')
plt.tight_layout()

plt.legend(ncol = 2, loc='upper center', bbox_to_anchor=(0.45, -0.45), 
    markerscale = 1.5)
    #-0.65
# save the figure
plt.savefig(path_plots + '/all-metrics-test-prob.svg')
plt.savefig(path_plots + '/all-metrics-test-prob.png')
plt.savefig(path_plots + '/all-metrics-test-prob.pdf')
plt.close()

# do it again, individual plots per metric
plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_test[metrics[0]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[0])
plt.xticks(ind, class_test['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[0])
plt.tight_layout()
plt.savefig(path_plots + '/{}-test.svg'.format(metrics[0]))
plt.savefig(path_plots + '/{}-test.png'.format(metrics[0]))
plt.close()

plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_test[metrics[1]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[1])
plt.xticks(ind, class_test['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[1])
plt.tight_layout()
plt.savefig(path_plots + '/{}-test.svg'.format(metrics[1]))
plt.savefig(path_plots + '/{}-test.png'.format(metrics[1]))
plt.close()

plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_test[metrics[2]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[2])
plt.xticks(ind, class_test['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[2])
plt.tight_layout()
plt.savefig(path_plots + '/{}-test.svg'.format(metrics[2]))
plt.savefig(path_plots + '/{}-test.png'.format(metrics[2]))
plt.close()

plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_test[metrics[3]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[3])
plt.xticks(ind, class_test['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[3])
plt.tight_layout()
plt.savefig(path_plots + '/{}-test.svg'.format(metrics[3]))
plt.savefig(path_plots + '/{}-test.png'.format(metrics[3]))
plt.close()

# and duplicate for the training data
plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_train[metrics[0]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[0])
plt.xticks(ind, class_train['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[0])
plt.tight_layout()
plt.savefig(path_plots + '/{}-train.svg'.format(metrics[0]))
plt.savefig(path_plots + '/{}-train.png'.format(metrics[0]))
plt.close()

plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_train[metrics[1]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[1])
plt.xticks(ind, class_train['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[1])
plt.tight_layout()
plt.savefig(path_plots + '/{}-train.svg'.format(metrics[1]))
plt.savefig(path_plots + '/{}-train.png'.format(metrics[1]))
plt.close()

plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_train[metrics[2]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[2])
plt.xticks(ind, class_train['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[2])
plt.tight_layout()
plt.savefig(path_plots + '/{}-train.svg'.format(metrics[2]))
plt.savefig(path_plots + '/{}-train.png'.format(metrics[2]))
plt.close()

plt.figure(figsize=(4,4), dpi=200)
width = 0.85
ind = np.arange(n_species)
plt.bar(ind, class_train[metrics[3]], width, color=sp_colors, alpha=0.9, edgecolor='black', label=metrics[3])
plt.xticks(ind, class_train['Species'], rotation='vertical')
plt.ylabel('Score')
plt.title(metrics[3])
plt.tight_layout()
plt.savefig(path_plots + '/{}-train.svg'.format(metrics[3]))
plt.savefig(path_plots + '/{}-train.png'.format(metrics[3]))
plt.close()