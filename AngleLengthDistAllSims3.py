#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import csv
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy import misc
import matplotlib.gridspec as gridspec
plt.rcParams.update({'errorbar.capsize': 2})
import matplotlib.lines  as mlines
import os
from math import sqrt
import math
import statistics

outputwidth   = 3.7
outputheight  = 6.24
fontsize      = 8
lblfontsize   = 18
lblfontsize2  = 10
linewidth     = 2.5
lgndmarkersize= 3.5
markersize    = 5
linelength    = 5
lgndmarkersize = 5
elinewidth     = 1
fontweight     = 'bold'
xMax           = 80


#Just for debugging, it should be turned off
plt.ioff()


# Create the figure
fig = plt.figure()
spec1 = gridspec.GridSpec(3,2,width_ratios = [0.9,0.05], hspace=0, wspace=0.1)
ax1 = fig.add_subplot(spec1[0,0])
ax2 = fig.add_subplot(spec1[1,0])
ax3 = fig.add_subplot(spec1[2,0])

Strain    = []
Strain2   = []
Cid       = []
Cos       = []
Len       = []
Cen       = []
zMin      = []
zMax      = []
StrainAvg = []
CosAvg    = []
CosErr    = []
LenAvg    = []
LenErr    = []
CenAvg    = []
CenErr    = []
zMinAvg   = []
zMaxAvg   = []

UniqueStrain  = []
UniqueStrain2 = []
UniqueCosAvg  = []
UniqueLenAvg  = []
UniqueCenAvg  = []
UniquezMinAvg = []
UniquezMaxAvg = []
UniqueCosErr  = []
UniqueLenErr  = []
UniqueCenErr  = []


cidC  = []
cidCR = []
with open("1-Series1/dump/dump.equilibrate.lammpstrj", 'r') as eqFile :
	dump = eqFile.readlines()
tsFound = False
numAtoms = 0
xmin = xmax = ymin = ymax = zmin = zmax = 0.0
for i in range(len(dump)):
	line = dump[i]
	if "ITEM: TIMESTEP" in line:
		ts = int(dump[i+1])
		i = i + 1
		if (ts == 0):
			tsFound = True
	elif "ITEM: NUMBER OF ATOMS" in line:
		numAtoms = int(dump[i+1])
		i = i + 1
	elif "ITEM: BOX BOUNDS" in line:
		line = dump[i+1].split()
		i = i + 1
		xmin = float(line[0])
		xmax = float(line[1])
		line = dump[i+1].split()
		i = i + 1
		ymin = float(line[0])
		ymax = float(line[1])
		line = dump[i+1].split()
		i = i + 1
		zmin = float(line[0])
		zmax = float(line[1])
	elif "ITEM: ATOMS" in line:
		if tsFound == False :
			i = i + numAtoms
		else :
			for j in range(numAtoms):
				#print j
				i = i + 1
				line   = dump[i].split()
				aid    = int(line[0])
				pType  = int(line[1])
				x      = float(line[2])
				y      = float(line[3])
				z      = float(line[4])
				x      = xmin + (xmax-xmin)*x
				y      = ymin + (ymax-ymin)*y
				radius = sqrt((x+46.5)*(x+46.5)+(y-147.9)*(y-147.9))
	
				if (pType == 1):
					cid = int(aid/219+1)
					cidC.append(cid)
					cidCR.append(radius)
			break

Sims =["1-Series1","2-Series2","3-Series3","4-Series4","5-Series5"]

for sim in Sims:
	file = sim + "/dump/z-cId/strain-MFangle-MFlength_v02.csv"
	with open(file, 'r') as csvfile:
		plots = csv.reader(csvfile,delimiter=',')
		next(plots)
		for row in plots:
			Strain.append(float(row[0]))
			Cid.append(int(row[1]))
			Cos.append(float(row[2]))
			Len.append(float(row[3]))
			Cen.append(float(row[4]))
	file2 = sim + "/dump/z-cId/strain-zMin-zMax.csv"
	with open(file2, 'r') as csvfile:
		plots = csv.reader(csvfile,delimiter=',')
		next(plots)
		for row in plots:
			Strain2.append(float(row[0]))
			zMin.append(float(row[1]))
			zMax.append(float(row[2]))

for strain in Strain:
	if strain not in UniqueStrain:
		UniqueStrain.append(strain)

for strain in Strain2:
	if strain not in UniqueStrain2:
		UniqueStrain2.append(strain)

for strain in UniqueStrain:
	cidS = []
	cosS = []
	lenS = []
	cenS = []
	for i in range(len(Strain)):
		if (Strain[i] == strain):
			cidS.append(Cid[i])
			cosS.append(Cos[i])
			lenS.append(Len[i])
			cenS.append(Cen[i])
	UniqueCosAvg.append(statistics.mean(cosS))
	UniqueLenAvg.append(statistics.mean(lenS))
	UniqueCenAvg.append(statistics.mean(cenS))
	UniqueCosErr.append(statistics.stdev(cosS))
	UniqueLenErr.append(statistics.stdev(lenS))
	UniqueCenErr.append(statistics.stdev(cenS))

for strain in UniqueStrain2:
	zminS= []
	zmaxS= []
	for i in range(len(Strain2)):
		if (Strain2[i] == strain):
			zminS.append(zMin[i])
			zmaxS.append(zMax[i])
	UniquezMinAvg.append(statistics.mean(zminS))
	UniquezMaxAvg.append(statistics.mean(zmaxS))

for i in range(len(Cen)):
	zmn = UniquezMinAvg[UniqueStrain2.index(Strain[i])]
	zmx = UniquezMaxAvg[UniqueStrain2.index(Strain[i])]
	Cen[i] = (Cen[i]-zmn)/(zmx-zmn)



# colormap
ax = fig.add_subplot(spec1[:,1])
ax.set_title("Microfibril Radius (A)",fontsize=fontsize-1,fontweight='bold')
N = 5
cmap = plt.get_cmap('rainbow', N)
# Normalizer
norm = mpl.colors.Normalize(vmin=min(cidCR), vmax=max(cidCR))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

plt.colorbar(sm,ax)

colors = []
for cid in Cid:
	color = cidCR[cidC.index(cid)]
	colors.append(color)

color = 'blue'
ax1.scatter( Strain, Cos, s=0.01,alpha=1, marker='o', c = colors, cmap=cmap)
color = 'red'
#ax1.errorbar(UniqueStrain,UniqueCosAvg,yerr=UniqueCosErr, alpha= 0.8, fmt='rx', color=color , markersize=1.5,capthick=0.5,capsize=1.5,elinewidth=elinewidth-0.5)
ax1.set_xlim(0,xMax)
ax1.set_ylim(0,3)
ax1.set_xticks([])
ax1.set_xlabel("Strain (%)",color='k',fontsize=fontsize, fontweight=fontweight)
ax1.set_ylabel("Angle (deg)",fontsize=fontsize, fontweight=fontweight)

x1_line = mlines.Line2D([],[],color='blue',mfc='b', mec='b',marker='o',markersize=lgndmarkersize,linestyle='none',label='Individual MFs')
x2_line = mlines.Line2D([],[],color='red',mfc='r', mec='r',marker='x',markersize=lgndmarkersize,linestyle='none',label='Average all the MFs')
x3_line = mlines.Line2D([],[],color='green',mfc='g', mec='g',marker='_',markersize=lgndmarkersize,linestyle='none',label='External Strain')

plt.sca(ax1)
legend = plt.legend(handles=[x1_line,x2_line,x3_line], fontsize='x-small' ,loc='upper left')
legend.get_frame().set_facecolor((1., .8, 0.1))
legend.get_frame().set_edgecolor((1., 0, 0))
legend.get_frame().set_alpha(0.5)
plt.sca(ax2)

color = 'blue'
ax2.scatter( Strain, Len, s=0.01,alpha=1, marker='o', c = colors, cmap=cmap)
color = 'red'
#ax2.errorbar(UniqueStrain,UniqueLenAvg,yerr=UniqueLenErr, alpha= 0.8, fmt='rx', color=color , markersize=1.5,capthick=0.5,capsize=1.5,elinewidth=elinewidth-0.5)
t = np.arange(0., xMax, 1)
plt.plot(t, t/100.0, 'g--')
ax2.set_xlim(0,xMax)
ax2.set_ylim(0,1)
#ax2.set_yticks(np.arange(2000,6000,2000))
ax2.set_xticks([])
ax2.set_xlabel("Strain (%)",color='k',fontsize=fontsize, fontweight=fontweight)
ax2.set_ylabel("Microfibril Strain",fontsize=fontsize, fontweight=fontweight)

color = 'blue'
ax3.scatter( Strain, Cen, s=0.01,alpha=1, marker='o', c = colors, cmap=cmap)
color = 'red'
#ax3.errorbar(UniqueStrain,UniqueCenAvg,yerr=UniqueCenErr, alpha= 0.8, fmt='rx', color=color , markersize=1.5,capthick=0.5,capsize=1.5,elinewidth=elinewidth-0.5)
color = 'green'
#ax3.scatter(UniqueStrain2, UniquezMinAvg , s=0.1,alpha=1, marker='o', c = color)
#ax3.scatter(UniqueStrain2, UniquezMaxAvg , s=0.1,alpha=1, marker='o', c = color)
ax3.set_xlim(0,xMax)
ax3.set_ylim(0.2,0.6)
#ax3.set_yticks(np.arange(0,4000,2000))
#ax.set_yticks([])
ax3.set_xlabel("Strain (%)",color='k',fontsize=fontsize, fontweight=fontweight)
ax3.set_ylabel("Center of Mass (A)",fontsize=fontsize, fontweight=fontweight)
		
for axis in ['top','bottom','left','right']:
	ax1.spines[axis].set_linewidth(2.5)
	ax2.spines[axis].set_linewidth(2.5)
	ax3.spines[axis].set_linewidth(2.5)


ax1.xaxis.set_tick_params(width=2.5,length=5)
ax1.yaxis.set_tick_params(width=2.5,length=5)
ax2.xaxis.set_tick_params(width=2.5,length=5)
ax2.yaxis.set_tick_params(width=2.5,length=5)
ax3.xaxis.set_tick_params(width=2.5,length=5)
ax3.yaxis.set_tick_params(width=2.5,length=5)
plt.sca(ax1)
plt.xticks(fontsize=fontsize, fontweight=fontweight)
plt.yticks(fontsize=fontsize, fontweight=fontweight)
plt.sca(ax2)
plt.xticks(fontsize=fontsize, fontweight=fontweight)
plt.yticks(fontsize=fontsize, fontweight=fontweight)
plt.sca(ax3)
plt.xticks(fontsize=fontsize, fontweight=fontweight)
plt.yticks(fontsize=fontsize, fontweight=fontweight)





fig.tight_layout()

plt.subplots_adjust(top=0.95, bottom=0.12, left=0.18, right=0.90, hspace=0.,wspace=0.05)


#Setting the output resolution
fileName = "MFLength-MFAngles-AllSims3"
fig.set_dpi(300)
fig.set_size_inches(outputwidth, outputheight)
fig.savefig(fileName+"-test.tif")

# Compressing the image
command = "tiffcp -c lzw " + fileName + "-test.tif " +fileName + ".tif"
os.popen(command)
command = "rm " + fileName + "-test.tif"
os.popen(command)

plt.close()

