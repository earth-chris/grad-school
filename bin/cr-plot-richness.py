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
base = '/home/cba/Downloads/costa-rica/'
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
cols = ['blue', 'brown', 'purple', 'orange', 'grey', 'red', 'green', 'magenta']

# set the response variable
yval = 'Observed_s'
ylabel = 'Observed species'

# set the units for each covariate
xunit = ['Mg C/ha', '%', '%', '%', 'dB', 'dB', 'unitless', 'unitless', 'unitless', 'C', 'C', 'C', '%']
yunit = 'count'

# loop through each covariate, then plot least squares fits for each taxon
for i in range(n_covar):
    for j in range(n_taxa):
        
        # extract the x and y data
        subset = rich[rich['Taxon'] == taxa[j]]
        x = subset[covar[i]]
        y = subset[yval]
        
        # calculate the fit 
        z = np.polyfit(x, y, poly[i])
        f = np.poly1d(z)
        
        # set the range for plotting fitted values
        x_new = np.linspace(x.min(), x.max(), 50)
        y_new = f(x_new)
        y_eval = f(x)
        
        # calculate some basic metrics
        rsq = metrics.r2_score(y, y_eval)
        rms = np.sqrt(metrics.mean_squared_error(y, y_eval))
        
        # plot the outputs
        plt.figure(figsize=(6,6), dpi=150)
        plt.scatter(x, y, color=cols[j], label='data', edgecolor='black')
        plt.plot(x_new, y_new, color='black', label='linear fit')
        
        # set the plot lables
        plt.xlabel('{} ({})'.format(covar[i], xunit[i]))
        plt.ylabel('{} ({})'.format(ylabel, yunit))
        plt.title('{} - {}\nr-squared: {:0.2f}, rmse: {:0.2f}'.format(taxa[j], covar[i], rsq, rms))
        
        # add the legend
        #plt.legend()
        
        # tighten it up
        plt.tight_layout()
        
        # then save the figure
        plt.savefig("{}{}-{}-linear.png".format(plots, taxa[j], covar[i]))
        plt.close()
        
# now we're going to do the silly thing of using a gradient boosting classifier to predict the observed # of species for each taxon
for j in range(n_taxa):

    # extract the x and y data
subset = rich[rich['Taxon'] == taxa[j]]
x = subset[covar]
y = subset[yval]

# find the best model params
tuner = aei.model.tune(x, y, n_splits=3)
tuner.GradientBoostRegressor(scoring='neg_mean_squared_error')

# clean up a deprecated param
del(tuner.best_params['min_impurity_split'])

# set up the model
gbr = ensemble.GradientBoostingRegressor(**tuner.best_params)

# run cross validation metrics
cv_score = model_selection.cross_validate(gbr, x, y, scoring=['r2', 'neg_mean_squared_error'])

# fit the model on all the data
gbr.fit(x, y)

# calculate the metrics
y_eval = gbr.predict(x)
rsq = metrics.r2_score(y, y_eval)
mse = metrics.mean_squared_error(y, y_eval)

# set the linear fit
z = np.polyfit(y_eval, y, 1)
f = np.poly1d(z)
x_new = np.linspace(y.min(), y.max(), 50)
y_new = f(x_new)

# plot the results
plt.figure(figsize=(6,6), dpi=150)
plt.scatter(y_eval, y, color=cols[j], edgecolor='black')
plt.plot(x_new, y_new, color='black', label='cross-val r-sq: {:0.2f}'.format(cv_score['test_r2'].mean()))

# set the limits to the bounds of the input data
plt.ylim([y.min()-1, y.max()+1])
plt.xlim([y.min()-1, y.max()+1])

# label the plots
plt.xlabel('Predicted species counts')
plt.ylabel('Observed species counts')
plt.title('GBM model\n{}'.format(taxa[j]))
plt.legend()
    
    # tight layout that shit
    plt.tight_layout()
    
    # save and close the figure
    plt.savefig('{}GBM-{}.png'.format(plots, taxa[j]))
    plt.close()