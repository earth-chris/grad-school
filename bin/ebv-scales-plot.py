# plots the satellite scale data
import aei
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib import lines
from matplotlib import rc
from matplotlib import colors as cols

# activate latex text rendering
#rc('text', usetex=True)

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
    
# set path to the input csv
input_file = aei.params.environ['AEI_GS'] + '/data/EBV-Satellite-scales.csv'

# read the file using pandas
df = pd.read_csv(input_file)

# first work on a scatter plot of the points
xdata = 'Plot revisit time'
xlabel = "low  --------------->  high\nFrequency (days between revisit)"
ydata = 'Plot resolution'
ylabel = "Grain (m)\ncoarse  --------------->  fine"
title = "Spatiotemporal Scales of\nBiodiversity Measurements from Earth Observations"

# set variable to color by
colorby = 'EBV Class'
unique = list(df[colorby].unique())
unique.sort()
ncolors = len(unique)
#colors = aei.color.color_blind(ncolors)
#colors = aei.color.color_blind()
colors = aei.objects.color(palette = ['#E56C2D', '#00A583', '#F1A53A', '#0081B4', '#F5E369'], 
    n = ncolors).palette

#color_map = cm.Dark2
#for val in np.arange(0 + 1./(ncolors+1), 1 + 1/(ncolors+1), 1./(ncolors+1)):
#    colors.append(color_map(val))
#df['scatter_color'] = df.apply(lambda row: label_color(
#    row, colorby, unique, colors), axis = 1)

# add columns to the data frame for plotting in a single go
dfl = len(df)
markers = ['o', 'D']
marker_titles = ['Continuous', 'Discrete']
#linestyle = ['solid', 'dashed']
df['Shape'] = pd.Series(np.repeat(markers[0], dfl))
df['Shape'][df['Coverage'] == marker_titles[1]] = markers[1]
#alphas = [0.9, 0.9]

styles = ['normal', 'italic']
styles_titles = ['Open access', 'Restricted']
df['Style'] = pd.Series(np.repeat(styles[0], dfl))
df['Style'][df['Access'] == styles_titles[1]] = styles[1]

plt.figure()
# plot by unique color scheme
for i in range(ncolors):
    fl = df[(df[colorby] == unique[i])]
    #plt.scatter(x = fl[xdata], y = fl[ydata], c = colors[i], 
    #    alpha = 0.9, label = unique[i], s = 175, marker = fl['Shape'],
    #    edgecolor='black', linewidth='1')

    for j in range(len(markers)):
        #print("label: {label}".format(label=unique[i]))
        #print('color: {color}'.format(color=colors[i]))
        #print('marker: {marker}'.format(marker=markers[j]))
        fc = fl[(fl['Shape'] == markers[j])]
        c = cols.to_hex(colors[i])
        #if len(fc) == len(colors[i]):
        #    c = cols.to_hex(colors[i])
        #else:
        #    c = cols[i]
        plt.scatter(x = fc[xdata], y = fc[ydata], c = c, 
            alpha = 0.9, label = unique[i], s = 195, marker = markers[j],
            edgecolor='black')#, linewidth='1')
    
# add titles/labels        
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
#plt.legend(loc = 'lower left', ncol = 2)

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
#ax.arrow(xmin, 0, xmax-xmin, 0)
#ax.arrow(0, ymin, 0, ymax-ymin)

# set the ticks
plt.xticks([1e0, 1e1, 50], ['1', '10', '50'])
plt.yticks([1e0, 1e1, 1e2, 1e3], ['1', '10', '100', '1000'])

# annotate each point
annotate = 'Sensor Name'
for i in range(len(df)):
    ax.annotate(df[annotate][i], 
        #xy = (df[xdata][i], df[ydata][i]), xytext = (4, 4),
        xy = (df[xdata][i], df[ydata][i]), xytext = (df['xoff'][i]*2, df['yoff'][i]*2),
        arrowprops = dict(arrowstyle = '-', color = 'black'),
        #textcoords = 'offset points', ha = 'left', va = 'bottom',
        textcoords = 'offset points', ha = df['ha'][i], va = df['va'][i],
        bbox = dict(boxstyle = 'round, pad=0.2', fc = 'black', alpha = 0.15),
        fontstyle = df['Style'][i])

# custom build the legend
legend = list(unique)
legend.append(marker_titles[1])
legend.append(marker_titles[0])
legend.append('')
legend.append('{style}'.format(style=styles_titles[0]))
legend.append('$\it{style}$'.format(style=styles_titles[1]))
legend_colors = list(colors)
legend_colors.append((1,1,1))
legend_colors.append((1,1,1))
legend_colors.append((0,0,0))
legend_colors.append((0,0,0))
legend_colors.append((0,0,0))
legend_marker = list(np.repeat(markers[0], ncolors))
legend_marker.append(markers[1])
legend_marker.append(markers[0])
legend_marker.append('')
legend_marker.append('')
legend_marker.append('')
#legend_linestyle = list(np.repeat(linestyle[0], ncolors))
#legend_linestyle.append(linestyle[0])
#legend_linestyle.append(linestyle[1])
#legend_linestyle.append('')
#legend_linestyle.append('')
#legend_linestyle.append('')
proxies = []
for i in range(len(legend)):
    proxies.append(create_proxy(legend_colors[i], legend_marker[i]))#, legend_linestyle[i]))
    
plt.legend(proxies, legend, loc = 'lower left', ncol = 2,
    markerscale = 2)

# set the figure size explicitly to fit full page width
xwidth = 6.811 # 173 millimeters is full width limit
ywidth = xwidth * 0.75
scaler = 1.5
fig = plt.gcf()
fig.set_size_inches(xwidth*scaler, ywidth*scaler, forward=True)
plt.tight_layout()

# save the output file
plt.savefig('/home/cba/cba/aei-grad-school/figures/EO-Biodiv-spatiotemporal-by-group-hr-access.tif',
    dpi=600/scaler)
plt.savefig('/home/cba/cba/aei-grad-school/figures/EO-Biodiv-spatiotemporal-by-group-hr-access.png',
dpi=600/scaler)
    
plt.close()

# full plot
#plt.tight_layout()
#plt.show()

##########
# bar plot of timeline for missions

# set labels and data to plot
xlabel = "Years of Operation"
ydata = 'Launch Year'
wdata = 'Decomission Year'
ylabel = "Earth Observation Mission"
ytlabel = 'Sensor Name'
title = 'Timeline of Earth Observations of Biodiversity'

# set variable to color by
colorby = 'EBV Class'
unique = list(df[colorby].unique())
unique.sort()
ncolors = len(unique)
#colors = aei.color.color_blind(ncolors)
#colors = aei.color.color_blind()
    
#distinct colors
colors = aei.objects.color(palette = ['#e6194b', '#C05354', '#C28053', '#A89846', '#FEF866',
    '#98E55A', '#6CE4B3', '#69C6E2', '#71B4FF', '#617BE3',
    '#8F5EFF', '#F95FF6', '#A54396', '#F86791'], 
    n = ncolors).palette
    
#colors = aei.color.color_blind()
colors = aei.objects.color(palette = ['#E56C2D', '#00A583', '#F1A53A', '#0081B4', '#F5E369'], 
    n = ncolors).palette

# sort the data frame by start year
#df_sorted = df.sort_values(by = ['EBV Class', 'Launch Year'], ascending = [0,1]).reset_index(drop = True)
df_sorted = df.sort_values(by = ['Sensor order', 'Launch Year'], ascending = [1,1]).reset_index(drop = True)

# set the number of locations to plot
w = 1.0
y = np.arange(len(df_sorted)) * w

# set the bar width
wdata = 'Decomission Year'
width = [i-j for i,j in zip(df_sorted[wdata], df_sorted[ydata])]

# set the minimum year for the plot base
min_year = 1970 # min(df_sorted[ydata]) - 1
max_year = 2019 # max(df_sorted[wdata]) + 1

# find the y position to place a dividing line
ysplit = np.where(df_sorted['Sensor order'] == 2)[0].min()

# move the ydata for each one down a notch
y[ysplit:] += 1

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
        
# add a solid line to distinguish between active and passive sensors
ax1.plot((min_year-2, max_year+2), (ysplit, ysplit), '-', c = 'black', linewidth = 0.5)

# annotate with sensor information
ax = plt.gca()
annotate = 'Sensor Description'
for i in range(len(df_sorted)):
    ax.annotate(df_sorted[annotate][i],
        xy = (df_sorted[ydata][i]-0.5, y[i]), #+ (width[i]/2.), y[i]),
        ha = 'right', va = 'center',
        bbox = dict(boxstyle = 'round, pad=0.1', lw=0, fc = 'white'),
        fontstyle = df_sorted['Style'][i])
    
# set italics for ylabels
ytick_labels = []
for i in range(len(df_sorted)):
    if df_sorted['Style'][i] == 'italic':
        ytick_labels.append('$\it{label}$'.format(label=df_sorted[ytlabel][i]))
    else:
        ytick_labels.append('{label}'.format(label=df_sorted[ytlabel][i]))

# label the axes
ax.invert_yaxis()
#plt.yticks(y, df_sorted[ytlabel])
plt.yticks(y, ytick_labels)
plt.xlim(min_year, max_year)
plt.ylim(max(y)+1, min(y)-1)
#plt.ylabel(ylabel)
plt.xlabel(xlabel)
plt.title(title)

# custom build the legend
legend = list(unique)
legend.append('{style}'.format(style=styles_titles[0]))
legend.append('$\it{style}$'.format(style=styles_titles[1]))
legend_colors = list(colors)
legend_colors.append((0,0,0))
legend_colors.append((0,0,0))
legend_marker = list(np.repeat('s', ncolors))
legend_marker.append('')
legend_marker.append('')
proxies = []
for i in range(len(legend)):
    proxies.append(create_proxy(legend_colors[i], legend_marker[i]))#, legend_linestyle[i]))
    
plt.legend(proxies, legend, loc = 'lower left', ncol = 1,
    markerscale = 2)
#plt.legend(loc = 'lower left')

# set the labels on the right
ax2 = ax1.twinx()
plt.ylim(ax1.get_ylim())
ax2.yaxis.tick_right()
plt.yticks([y[:ysplit].mean(), y[ysplit:].mean()], ['Passive sensors', 'Active sensors'], rotation=270,
    horizontalalignment='left', verticalalignment='center')
ax2.tick_params(length=0)
#plt.yticks(y, df_sorted['EBV Authors'])
#plt.ylabel("Reference")

# set the figure size explicitly to fit full page width
xwidth = 6.811 # 173 millimeters is full width limit
ywidth = xwidth * 0.8
scaler = 1.8
fig = plt.gcf()
fig.set_size_inches(xwidth*scaler, ywidth*scaler, forward=True)
plt.tight_layout()

# save the output file
plt.savefig('/home/cba/cba/aei-grad-school/figures/EO-Biodiv-timeline-by-group-and-sensor.tif',
    dpi=600/scaler)
plt.savefig('/home/cba/cba/aei-grad-school/figures/EO-Biodiv-timeline-by-group-and-sensor.png',
dpi=600/scaler)