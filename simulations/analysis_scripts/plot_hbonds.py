#Script for plotting the results from the rmsd analysis

#imports
import matplotlib
import matplotlib.pyplot as plt


#open the rmsd file and read it line by line
fn = open("hbonds.dat")
lines = fn.readlines()

frame, bonds = [],[]

#Get the frame and hbonds data
for line in lines:
	line = line.strip('\n').split(' ')
	frame.append(int(line[0]))
	bonds.append(int(line[1]))

for i in range(len(frame)):
	frame[i] = ((frame[i]*250)*2)/1000000

#Create the plot
fig, ax = plt.subplots()
ax.plot(frame, bonds)

#label the axies
ax.set(xlabel='Time (ns)', ylabel='H-Bonds',
       title='Hydrogen Bonds Between the N and C-Termani')
ax.grid()

#output the figure
fig.savefig("hxt1_hbonds_plot.png")
plt.show()
