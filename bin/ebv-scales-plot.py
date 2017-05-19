# plots the satellite scale data
import aei
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from adjustText import adjust_text

# function to label new colors
def label_color(row, colorby, unique_vals, colors):
    # get the index for the value passed and return the
    #  associated color
    ind = unique_vals.index(row[colorby])
    return colors[ind]
    
# set path to the input csv
input_file = aei.params.environ['AEI_GS'] + '/data/EBV-Satellite-scales.csv'

# read the file using pandas
df = pd.read_csv(input_file)

# first work on a scatter plot of the points
xdata = 'Plot revisit time'
xlabel = "Frequency (days between revisit)"
ydata = 'Plot resolution'
ylabel = "Spatial resolution (m)"
title = "Spatiotemporal Scales of\nBiodiversity Measurements from Earth Observations"

# set variable to color by
colorby = 'EBV Class'
unique = list(df[colorby].unique())
ncolors = len(unique)
colors = aei.color.color_blind(ncolors)
colors = []
color_map = cm.Dark2
for val in np.arange(0 + 1./(ncolors+1), 1 + 1/(ncolors+1), 1./(ncolors+1)):
    colors.append(color_map(val))
#df['scatter_color'] = df.apply(lambda row: label_color(
#    row, colorby, unique, colors), axis = 1)

plt.figure()
# plot by unique color scheme
for i in range(ncolors):
    fl = df[(df[colorby] == unique[i])]
    plt.scatter(x = fl[xdata], y = fl[ydata], c = colors[i], 
        alpha = 0.9, label = unique[i], s = 175,
        edgecolor='black', linewidth='1')
    
# add titles/labels        
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.legend(loc = 'lower left', ncol = 2)

# handle axes
xmin = df[xdata].min() * 0.8
xmax = df[xdata].max() * 1.2
ymin = df[ydata].min() * 0.8
ymax = df[ydata].max() * 1.2
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
ax = plt.gca()
ax.invert_xaxis()
ax.invert_yaxis()
ax.loglog()

# annotate each point
annotate = 'Sensor Name'
for i in range(len(df)):
    ax.annotate(df[annotate][i], 
        #xy = (df[xdata][i], df[ydata][i]), xytext = (4, 4),
        xy = (df[xdata][i], df[ydata][i]), xytext = (df['xoff'][i]*2, df['yoff'][i]*2),
        arrowprops = dict(arrowstyle = '-', color = 'black'),
        #textcoords = 'offset points', ha = 'left', va = 'bottom',
        textcoords = 'offset points', ha = df['ha'][i], va = df['va'][i],
        bbox = dict(boxstyle = 'round, pad=0.2', fc = 'black', alpha = 0.1))

# full plot
plt.tight_layout()
ax.plot()

##########
# bar plot of timeline for missions

# set labels and data to plot
xlabel = "Years of Operation"
ydata = 'Launch Year'
wdata = 'Decomission Year'
ylabel = "Earth Observation Mission"
ytlabel = 'Sensor Name'
title = 'Timeline of Earth Observations of Biodiversity'

# sort the data frame by start year
df_sorted = df.sort_values(by = ['EBV Class', 'Launch Year'], ascending = [0,1]).reset_index(drop = True)

# set the number of locations to plot
w = 1.0
y = np.arange(len(df_sorted)) * w

# set the bar width
wdata = 'Decomission Year'
width = [i-j for i,j in zip(df_sorted[wdata], df_sorted[ydata])]

# set the minimum year for the plot base
min_year = 1970 # min(df_sorted[ydata]) - 1
max_year = 2018 # max(df_sorted[wdata]) + 1

# plot as a bar chart
fig, ax1 = plt.subplots()
for i in range(ncolors):
    ind = np.where(df_sorted[colorby] == unique[i])
    ax1.barh(np.array(y)[ind[0]], 
        np.array(width)[ind[0]],
        left = df_sorted[ydata][ind[0]],
        color = colors[i],
        alpha = 0.9,
        label = unique[i],
        edgecolor='black', linewidth='1')

# add some dashed lines to track each sensor
for i in range(len(df_sorted)):
    ax1.plot((min_year, df_sorted[ydata][i]),
        (y[i], y[i]), '--', c = 'black', alpha = 0.2)

# annotate with sensor information
ax = plt.gca()
annotate = 'Sensor Description'
for i in range(len(df_sorted)):
    ax.annotate(df_sorted[annotate][i],
        xy = (df_sorted[ydata][i]-0.5, y[i]), #+ (width[i]/2.), y[i]),
        ha = 'right', va = 'center',
        bbox = dict(boxstyle = 'round, pad=0.1', lw=0, fc = 'white'))
    
# label the axes
ax.invert_yaxis()
plt.yticks(y, df_sorted[ytlabel])
plt.xlim(min_year, max_year)
plt.ylabel(ylabel)
plt.xlabel(xlabel)
plt.title(title)
plt.legend(loc = 'lower left')
#ax2 = ax1.twinx()
#plt.ylim(ax1.get_ylim())
#ax2.yaxis.tick_right()
#plt.yticks(y, df_sorted['EBV Authors'])
#plt.ylabel("Reference")
plt.tight_layout()