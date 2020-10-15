# 
import aei
import gdal
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

%matplotlib tk

# set the working directories
base = '/home/cba/Downloads/scale-conceptual/'
plots = base + 'plots/'
tch_file = base + 'marin_tch_byte.tif'
gd_file = base + 'marin_tch_mask.tif'
rnd_data = base + 'random-sampling-data.pck'

# set flag for whether to calculate or load the random sampling
calc_rs = False

# set the number of resolutions to assess
#res = [20, 30, 50, 100, 250, 500, 1000, 1250, 1500, 1750, 2000]
res = [20, 30, 50, 100, 250, 500, 1000]

# set the number of pixels to sample 
fullres = 2 * 2
testres = min(res) * min(res)
n_var = testres / fullres
n_loc = n_var
n_rep = 1000

if calc_rs:
    # open the suitability raster and read in good indices
    gdref = gdal.Open(gd_file)
    gd = gdref.ReadAsArray()
    #gd_ext1 = gd == 1
    gd_ext1 = np.where(gd == 1)
    bd = gd == 0
    #gd_ext2 = np.where(gd > 1)
    gd = None
    gdref = None
    
    # read the tch data into memory 
    tchref = gdal.Open(tch_file)
    #tchb1 = tchref.GetRasterBand(1)
    #ndval = tchb1.GetNoDataValue()
    tch = tchref.ReadAsArray().astype(np.float)
    tchb1 = None
    tchref = None
    
    # set the nodata values to nan
    tch[bd] = np.nan
    
    # create a set of arrays to store the data
    wg_var = np.zeros((len(res), n_rep, n_loc))
    bg_var = np.zeros((len(res), n_rep, n_loc))
    
    # loop through each resolution, and calculate within/between grain variance
    for i in range(len(res)):
        rdif = res[i] / 2
        for j in range(n_rep):
            
            # find the random pixel locations to select
            loc_ext1 = np.random.choice(len(gd_ext1[0]), n_loc)
            #loc_ext2 = np.random.choice(len(gd_ext2[0]), n_loc)
            
            # loop through each random location
            for k in range(n_loc):
                # sample the x/y space around this pixel
                xmin = gd_ext1[0][loc_ext1[k]] - rdif
                xmax = gd_ext1[0][loc_ext1[k]] + rdif
                ymin = gd_ext1[1][loc_ext1[k]] - rdif
                ymax = gd_ext1[1][loc_ext1[k]] + rdif
                subset = tch[xmin:xmax, ymin:ymax]
                avg_ext1 = np.nanmean(subset)
                var_ext1 = np.nanvar(subset)
                wg_var[i,j,k] = var_ext1
                bg_var[i,j,k] = avg_ext1
                
                # then do it for the second extent
                #xmin = gd[0][loc_ext2[k]] - rdif
                #xmax = gd[0][loc_ext2[k]] + rdif
                #ymin = gd[1][loc_ext2[k]] - rdif
                #ymax = gd[1][loc_ext2[k]] + rdif
                #subset = tch[xmin:xmax, ymin:ymax]
                #avg_ext2 = np.nanmean(subset)
                #var_ext2 = np.nanvar(subset)
                
    # save this loop to a data file
    with open(rnd_data, 'w+') as f:
        pickle.dump([wg_var, bg_var], f)
else:
    with open(rnd_data, 'r+') as f:
        wg_var, bg_var = pickle.load(f)
        
# now that we have the data, reduce it to the plot the within/between grain variance
wg_mean_loc = np.nanmean(wg_var, axis=2) # mean across all locations per rep
wg_mean_all = np.nanmean(wg_mean_loc, axis=1) # mean of all locs, all reps
wg_std = np.nanstd(wg_mean_loc, axis=1) # stdev of within grain variation across reps

bg_var_loc = np.nanvar(bg_var, axis=2), # variance of between-grain measurements
bg_mean = np.nanmean(bg_var_loc[0], axis=1) # mean of between-grain variance
bg_std = np.nanstd(bg_var_loc[0], axis=1) # stdev of between-grain variance

# plot these results
cols = aei.color.color_blind(5)
col = '#E3C16D'
#col=cols[4]
fill_alpha = 0.7
plt.figure(figsize=(4,3), dpi=200)

# plot the standard deviations for each
plt.fill_between(res, wg_mean_all-wg_std, wg_mean_all+wg_std, 
  color=col, alpha=fill_alpha, label='Bootstrapped\nstandard deviation')
plt.fill_between(res, bg_mean-bg_std, bg_mean+bg_std, 
  color=col, alpha=fill_alpha)

# set the line plots
plt.plot(res, wg_mean_all, color='black', linewidth=1.5, linestyle = '-',
    label='Within-grain\nvariance')
plt.plot(res, bg_mean, color='black', linewidth=1.5, linestyle = '--',
    label='Between-grain\nvariance')

# set the labels
plt.xlabel('Log grain size (m)')
plt.ylabel('Spatial\nvariance (m{})'.format(r'$^2$'))
plt.title('Scale-dependence in\nspatial tree height patterns')

# log transform the axes
#plt.yticks([], [])
#ax = plt.gca()
#ax.loglog()
plt.xscale('log')
#plt.yscale('log')

# replace the axis labels
#yticks=(20, 40, 60, 80, 100)
yticks=[25, 50, 75, 100]
#ylabels = ax.get_yticklabels()
plt.yticks(yticks, ['{}'.format(f) for f in yticks])
plt.xticks(res, res)

# set the custom legend
lgd = plt.legend(loc='right', bbox_to_anchor=(1.7, 0.5), fancybox=True)#, 
#    fancybox=True, shadow=True)

# save the figure
plt.tight_layout()
plt.savefig('{}{}.png'.format(plots, 'Variance-plots'), dpi=300,
    additional_artists = [lgd], bbox_inches="tight")
plt.savefig('{}{}.svg'.format(plots, 'Variance-plots'), 
    additional_artists = [lgd], bbox_inches="tight")