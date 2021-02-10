#Script for plotting the results from the rmsd analysis

#imports
import matplotlib
import matplotlib.pyplot as plt


#open the rmsd file and read it line by line
fn = open("hxt1_rmsd.dat")
lines = fn.readlines()
fn.close()


frame, rmsd = [],[]


#Get the frame and rmsd data
for line in lines:
	line = line.split('\t')
	frame.append(int(line[0]))
	rmsd.append(float(line[1]))

for i in range(len(frame)):
	frame[i] = ((frame[i]*250)*2)/1000000
#Create the plot
fig, ax = plt.subplots()
ax.plot(frame, rmsd)

#label the axies
ax.set(xlabel='Time (ns)', ylabel='RMSD',
       title='RMSD Over Time')
ax.grid()

#output the figure
fig.savefig("hxt1_rmsd_plot.png")
plt.show()
