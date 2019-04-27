import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import sys
import os
import shutil
Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
dt=0.01             #Size of timestep, in Ga


def rel_shells(file,startline):
	f=open(file,'r')	
	lines=f.readlines()[startline:]
	f.close()
	radii=[]
	surfa=[]
	vol_cumul=[]
	weight_total=[]
	weight_local=[]
	count=0
	for line in lines:
		temp=line.split(',')
		R=float(temp[1])*1000 #converts to meters
		SA=4*np.pi*R**2
		V=(4./3.)*np.pi*R**3
		radii.append(R)
		surfa.append(SA)
		vol_cumul.append(V)
	Vtot=max(vol_cumul)-min(vol_cumul)
	for i in range(len(radii)-1):
		SA1=surfa[i]/(surfa[i]+surfa[i+1])
		SA2=surfa[i+1]/(surfa[i]+surfa[i+1])
		Vsh=(vol_cumul[i+1]-vol_cumul[i])
		weight_local.append([SA1,SA2])
		weight_total.append(Vsh/Vtot)
	return radii, weight_local, weight_total


def plot_rel_contributions(radii,wt_local):
	fig = figure(1)
	ax = fig.add_subplot(111, autoscale_on=True)
	weight_local=np.asarray(wt_local)
	ax.plot(radii[:-1],weight_local[:,0],'k',linewidth=2)
	maxx=max(radii)
	maxy=max(weight_local[:,0])
	minx=min(radii)
	miny=min(weight_local[:,0])

	plt.xlim(minx,maxx)
	plt.ylim(miny,maxy)
	#title='depth_rel_contribution.pdf'
	#plt.savefig(title)
	#return title
	return ' '

def weightedavg1(file,startline,wt_local,wt_total):
	f=open(file,'r')	
	lines=f.readlines()[startline+1:]
	f.close()
	layernumber=0
	sumtotals=[]
	for line in lines:
		temp=line.split(',')
		temp2=[]
		i=0
		for item in temp:
			temp2.append(float(item)*wt_total[layernumber])
		layernumber=layernumber+1
		sumtotals.append(temp2)
	a=np.asarray(sumtotals)
	return a

def read_cols(file):
	f=open(file,'r')
	lines=f.readlines()[0].split(',')
	lines[0]=lines[0][1:] #removes comment character from header
	lines[-1]=lines[-1][0:-1] #removes newline character 
	return lines

