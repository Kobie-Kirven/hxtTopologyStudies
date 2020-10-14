### Plot Generation
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tmsr import readTMSR

##The input file
inputFile = readTMSR('hxt1.tmsr')

#The colors of the boxes
colors = ['green','yellow',
         'blue','orange','red','pink','brown']

tools = inputFile[0]
coordinates = inputFile[1]
fig = plt.figure()
ax = fig.add_subplot(111)

#Build the chart
height = len(tools)*3
height = height -1
for z in range(len(coordinates)):
    height = height - ((len(tools)-1.5)/2)
    for i in range(len(coordinates[z])):
        rect1 = matplotlib.patches.Rectangle((coordinates[z][i][0],height),
                                             (coordinates[z][i][1]-coordinates[z][i][0]),
                                             2, color =colors[z])
        ax.add_patch(rect1)

plt.xlim([0, 570])
plt.ylim([0, 20])

frame1 = plt.gca()
frame1.axes.get_yaxis().set_ticks([])

#Make the legend
handle = []
for g in range(len(tools)):
    handle.append(mpatches.Patch(color=colors[g], label=tools[g]))
    plt.legend(handles=handle,bbox_to_anchor=(1.05, 1), loc='upper left')

#Axis Labels
plt.xlabel('Residue Position')
plt.title('Transmembrane Region Predictions: hxt1')

#Output file
plt.tight_layout()
plt.savefig('topologyOutputs/hxt1_toplogy.pdf')
