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
xdata = 'Min. revisit time (days)'
xlabel = "Temporal resolution (freq.)"
ydata = 'Min. spatial resolution (m)'
ylabel = "Spatial resolution (m)"

# set variable to color by
colorby = 'EBV Class'
unique = list(df[colorby].unique())
ncolors = len(unique)
colors = aei.color.color_blind(ncolors)
df['scatter_color'] = df.apply(lambda row: label_color(
    row, colorby, unique, colors), axis = 1)

plt.figure()
# plot by unique color scheme
for i in range(ncolors):
    fl = df[(df[colorby] == unique[i])]
    plt.scatter(x = fl[xdata], y = fl[ydata], c = colors[i], 
        alpha = 0.9, label = unique[i], s = 150)
    
# add titles/labels        
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()

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
annotate = 'Sensor/instrument'
for i in range(len(df)):
    ax.annotate(df[annotate][i], 
        xy = (df[xdata][i], df[ydata][i]), xytext = (8, 0),
        #arrowprops = dict(arrowstyle = '-', color = 'black'),
        textcoords = 'offset points', ha = 'left', va = 'bottom',
        bbox = dict(boxstyle = 'round, pad=0.2', fc = 'black', alpha = 0.1))

# full plot
ax.plot()