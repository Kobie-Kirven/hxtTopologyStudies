### Plot Generation
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tmsr import readTMSR

coords = [[0,0],[0,1],[1,0],[1,1]]

#The colors of the boxes
colors = ['green','yellow','blue','orange','red','pink','brown']

fig, axs = plt.subplots(2, 2)

ax1, ax2, ax3, ax4 = axs.ravel()

axes = (ax1, ax2, ax3, ax4)

for r in range(len(coords)):

	file = 1
	##The input file
	inputFile = readTMSR('../tmsrFiles/hxt' + str(file) + '.tmsr')

	tools = inputFile[0]
	coordinates = inputFile[1]

	#Build the chart
	# height = len(tools)*3
	# height = height
	change = (20 / len(tools))
	height = 0
	for z in range(len(coordinates)):
	    for i in range(len(coordinates[z])):
	        rect1 = matplotlib.patches.Rectangle((coordinates[z][i][0],height),
	                                             (coordinates[z][i][1]-coordinates[z][i][0]),
	                                             change-0.1, color =colors[z])
	        axes[r].add_patch(rect1)

	    height += change

	plt.xlim([0, 570])
	plt.ylim([0, 20])

	plt.axes.get_yaxis().set_ticks([])

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
plt.savefig('../outputFiles/hxt1_toplogy.png')
