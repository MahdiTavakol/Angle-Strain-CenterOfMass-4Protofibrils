#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import csv
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy import misc
import matplotlib.gridspec as gridspec
plt.rcParams.update({'errorbar.capsize': 2})
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
fig, (ax1,ax2,ax3) = plt.subplots(3,1)

Strain    = []
Cos       = []
Len       = []
Cen       = []
Len2      = []
StrainAvg = []
CosAvg    = []
CosErr    = []
LenAvg    = []
LenErr    = []
Len2Avg   = []
Len2Err   = []
CenAvg    = []
CenErr    = []
zMax      = []
zMin      = []
CID       = []


Len0      = []
cid0      = []       


cid       = []
cidA0     = []

cidStrain = []
cidDelA   = []


with open("../dump.defo.lammpstrj", 'r') as eqFile :
	dump = eqFile.readlines()
tsFound = False
numAtoms = 0
xmin = xmax = ymin = ymax = zmin = zmax = 0.0
for i in range(len(dump)):
	line = dump[i]
	if "ITEM: TIMESTEP" in line:
		ts = int(dump[i+1])
		i = i + 1
		if ((ts%400000)==0):
			tsFound = True
			print "working on frame " , ts
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
			#if (ts == 15600000):
			#	break
			tsFound = False
			X1   = []
			Y1   = []
			Z1   = []
			cid1 = []
			X2   = []
			Y2   = []
			Z2   = []
			cid2 = []
			Xcen = 0.0
			Ycen = 0.0
			num  = 0
			
			StrainCurr = []
			CosCurr    = []
			LenCurr    = []
			Len2Curr   = []
			CenCurr    = []

			for j in range(numAtoms):
				#print j
				i = i + 1
				line   = dump[i].split()
				aid    = int(line[1])-1 # The Lammps has 1-based indexing system
				pType  = int(line[2])
				if (pType == 1):
					if (True):
					#if ((aid%219 == 0) or (aid%219 == 218)):
						x      = float(line[3])
						y      = float(line[4])
						z      = float(line[5])
						x      = xmin + (xmax-xmin)*x
						y      = ymin + (ymax-ymin)*y
						z      = zmin + (zmax-zmin)*z
						Xcen = x + Xcen
						Ycen = y + Ycen
						num = num + 1
					if (aid%219 ==0):
							X1.append(x)
							Y1.append(y)
							Z1.append(z)
							cid1.append(int(aid/219))
					elif (aid%219 == 218):
							X2.append(x)
							Y2.append(y)
							Z2.append(z)
							cid2.append(int(aid/219))
			Xcen = Xcen/num
			Ycen = Ycen/num
			for i in range(len(cid1)):
				x1 = X1[i]
				y1 = Y1[i]
				z1 = Z1[i]
				if (cid1[i] not in cid2):
					continue
				x2 = X2[cid2.index(cid1[i])]
				y2 = Y2[cid2.index(cid1[i])]
				z2 = Z2[cid2.index(cid1[i])]
				r1 = sqrt((x1-Xcen)*(x1-Xcen)+(y1-Ycen)*(y1-Ycen))
				r2 = sqrt((x2-Xcen)*(x2-Xcen)+(y2-Ycen)*(y2-Ycen))
				if (z2 < z1):
					if ((z2-zmin)/(zmax-zmin) < 0.5):
						z2 = z2 + zmax - zmin
					else:
						z1 = z1 - zmax + zmin
				a = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
				cos = math.acos((z2-z1)/a)*180/math.pi
				lenn = z2 - z1
				cen  = (z1+z2)/2.0
				index = int((cen - zmin)/((zmax-zmin)/5))
				cindex = float(index*(zmax-zmin)/5.0) + zmin
				cen = cen - cindex
				lenn2 = lenn
				if (ts == 0):
					Len0.append(lenn)
					cid0.append(cid1[i])
				lenn = (lenn-Len0[cid0.index(cid1[i])])/Len0[cid0.index(cid1[i])]
				if (True): #cos < 10): #True): #i > 66 and i < 85) : #(r1 > 70 and r2 > 70):
					Strain.append(float(ts/400000))
					Cos.append(cos)
					Len.append(lenn)
					Len2.append(lenn2)
					Cen.append(cen)
					StrainCurr.append(float(ts/400000))
					CosCurr.append(cos)
					LenCurr.append(lenn)
					Len2Curr.append(lenn2)
					CenCurr.append(cen)
					CID.append(i+1)
			StrainAvg.append(sum(StrainCurr)/len(StrainCurr))
			CosAvg.append(sum(CosCurr)/len(CosCurr))
			CosErr.append(statistics.stdev(CosCurr))
			LenAvg.append(sum(LenCurr)/len(LenCurr))
			LenErr.append(statistics.stdev(LenCurr))
			Len2Avg.append(sum(Len2Curr)/len(Len2Curr))
			Len2Err.append(statistics.stdev(Len2Curr))
			CenAvg.append(sum(CenCurr)/len(CenCurr))
			CenErr.append(statistics.stdev(CenCurr))
			zMax.append(zmax)
			zMin.append(zmin)

#for i in range(1,len(Len)):
	#Len[i] = (Len[i]-LenAvg[0])/LenAvg[0]
#for i in range(1,len(LenAvg)):
	#LenAvg[i] = (LenAvg[i]-LenAvg[0])/LenAvg[0]

color = 'blue'
ax1.scatter( Strain, Cos, s=0.01,alpha=1, marker='o', c = color)
color = 'red'
ax1.errorbar(StrainAvg,CosAvg,yerr=CosErr, alpha= 0.8, fmt='ro', color=color , markersize=1,capthick=0.5,capsize=1.5,elinewidth=elinewidth-0.5)
ax1.set_xlim(0,xMax)
ax1.set_ylim(0,5)
ax1.set_xticks([])
ax1.set_xlabel("Strain (%)",color='k',fontsize=fontsize, fontweight=fontweight)
ax1.set_ylabel("Angle (deg)",fontsize=fontsize, fontweight=fontweight)


color = 'blue'
ax2.scatter( Strain, Len, s=0.01,alpha=1, marker='o', c = color)
color = 'red'
ax2.scatter( StrainAvg, LenAvg, s=1,alpha=0.8, marker='o', c = color)
#ax2.errorbar(StrainAvg,LenAvg,yerr=LenErr, alpha= 0.8, fmt='ro', color=color , markersize=1,capthick=0.5,capsize=1.5,elinewidth=elinewidth-0.5)
ax2.set_xlim(0,xMax)
ax2.set_ylim(0,1)
ax2.set_xticks([])
ax2.set_xlabel("Strain (%)",color='k',fontsize=fontsize, fontweight=fontweight)
ax2.set_ylabel("Microfibril Strain",fontsize=fontsize, fontweight=fontweight)

color = 'blue'
ax3.scatter( Strain, Cen, s=0.01,alpha=1, marker='o', c = color)
color = 'red'
ax3.errorbar(StrainAvg,CenAvg,yerr=LenErr, alpha= 0.8, fmt='ro', color=color , markersize=1,capthick=0.5,capsize=1.5,elinewidth=elinewidth-0.5)
color = 'green'
ax3.scatter( StrainAvg, zMin, s=0.1,alpha=1, marker='o', c = color)
ax3.scatter( StrainAvg, zMax, s=0.1,alpha=1, marker='o', c = color)
ax3.set_xlim(0,xMax)
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

plt.subplots_adjust(top=0.95, bottom=0.12, left=0.18, right=0.95, hspace=0.,wspace=0.05)


with open('strain-MFangle-MFlength.csv', 'w') as csvfilew:
	plotsw = csv.writer(csvfilew,delimiter=',')
	plotsw.writerow(["Strain(%)","MF-Angle","MF-Length","MF-CenterOfMass"])
	plotsw.writerows(zip(Strain,Cos,Len,Cen))

with open('strain-MFangle-MFlength_v02.csv', 'w') as csvfilew:
	plotsw = csv.writer(csvfilew,delimiter=',')
	plotsw.writerow(["Strain(%)","cid","MF-Angle","MF-Strain","MF-CenterOfMass"])
	plotsw.writerows(zip(Strain,CID,Cos,Len,Cen))

with open('strain-MFangle-MFlength_v03.csv', 'w') as csvfilew:
	plotsw = csv.writer(csvfilew,delimiter=',')
	plotsw.writerow(["Strain(%)","cid","MF-Length"])
	plotsw.writerows(zip(Strain,CID,Len2))

with open('strain-zMin-zMax.csv', 'w') as csvfilew:
	plotsw = csv.writer(csvfilew,delimiter=',')
	plotsw.writerow(["Strain(%)","zmin","zmax"])
	plotsw.writerows(zip(StrainAvg,zMin,zMax))


#Setting the output resolution
fileName = "MFLength-MFAngles"
fig.set_dpi(300)
fig.set_size_inches(outputwidth, outputheight)
fig.savefig(fileName+"-test.tif")

# Compressing the image
command = "tiffcp -c lzw " + fileName + "-test.tif " +fileName + ".tif"
os.popen(command)
command = "rm " + fileName + "-test.tif"
os.popen(command)


