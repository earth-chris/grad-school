# 
import aei
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn import metrics
from sklearn import ensemble
from sklearn import model_selection

%matplotlib tk

# set the paths to the input and output data
base = '/home/salo/Downloads/costa-rica/'
plots = base + 'plots/'
richness_file = base + 'sp-data/sites-covariates.csv'

# read the richness data into memory
rich = pd.read_csv(richness_file)

# find the unique taxa
taxa = rich['Taxon'].unique()
taxa.sort()
n_taxa = len(taxa)

# and the unique classes to plot
covar = rich.columns[11:]
n_covar = len(covar)

# set the polynomial orders
poly = [1,1,1,1,1,1,1,1,1,1,1,1,1]

# set the colors for each taxon
#cols = ['blue', 'brown', 'purple', 'orange', 'grey', 'red', 'green', 'magenta']
cols = aei.color.color_blind(n_taxa)

# set the response variable
yval = 'Observed_s'
xlabel = 'Observed species (count)'
ylabel = 'Frequency'
title = 'CCB plot network species richness'

# loop through each taxon and extract the richness data to use
for i in range(len(taxa)):
    
    # subset the taxon data
    subset = rich[rich['Taxon'] == taxa[i]]
    vals = np.array(subset[yval])
    
    # figure out the number of bins to plot
    ranges = vals.max() - vals.min()
    iterator = int(np.ceil(ranges / 10.))
    bins = np.arange(vals.min(), vals.max()+iterator, iterator) - (iterator/2.)
    
    # set up the plot device
    plt.figure(figsize = (5, 5), dpi=150)

    # plot the histogram
    plt.hist(vals, color=cols[i], bins=bins, edgecolor='black', 
        label='n = {}'.format(len(vals)))
    
    # label it up
    plt.xlabel(xlabel)
    plt.title('Species richness\n{}'.format(taxa[i]))
    plt.ylabel(ylabel)
    plt.legend()
    
    # and save the results
    plt.tight_layout()
    plt.savefig('{}Richness-{}.png'.format(plots, taxa[i]))
    plt.close()
    