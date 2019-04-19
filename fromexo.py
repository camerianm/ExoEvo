import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import sys
import os
import shutil
Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
Qe = 34.0e-12       #Early Earth heat production in W/kg mantle
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
	for line in lines:
		temp=line.split(', ')
		R=float(temp[1])*1000 #converts to meters
		SA=4*np.pi*R**2
		V=(4./3.)*np.pi*R**3
		radii.append(R)
		surfa.append(SA)
		vol_cumul.append(V)
	Vtot=max(vol_cumul)
	for i in range(len(radii)-1):
		SA1=surfa[i]/(surfa[i]+surfa[i+1])
		SA2=surfa[i+1]/(surfa[i]+surfa[i+1])
		Vsh=abs(vol_cumul[i+1]-vol_cumul[i])
		weight_local.append([SA1,SA2])
		weight_total.append(Vsh/Vtot)
	return radii, weight_local, weight_total

#def calc_averages()

#def heat_production(Mp,Mc,Qpl,t):
#	Q0=(Qe*Qpl)*(Mp-Mc)       #initial radionuclide abundance total - w/kg versus Earth, times #kg mantle
#	t=0.0
#   Ht=[0.0,0.0,0.0,0.0]
#	for i in range(3):
#		Ht[i]=radio['wtpercent'][i]*np.exp(radio['lambda'][i]*t)
#	print('H')

#       
#    heat=sum(Ht)
#    return heat
def plot_rel_contributions(radii,wt_local):
	fig = figure(1)
	ax = fig.add_subplot(111, autoscale_on=True)
	weight_local=np.asarray(wt_local)
	ax.plot(radii[:-1],weight_local[:,0],'k',linewidth=2)
	maxx=max(radii)
	maxy=max(weight_local[:,0])
	minx=min(radii)
	miny=min(weight_local[:,0])

	#plt.xlim(float(sys.argv[7]),float(sys.argv[8]))
	#plt.ylim(float(sys.argv[9]),float(sys.argv[10]))
	plt.xlim(minx,maxx)
	plt.ylim(miny,maxy)
	title='depth_rel_contribution.pdf'
	plt.savefig(title)
	return title

def weightedavg1(file,startline,wt_local,wt_total):
	f=open(file,'r')	
	lines=f.readlines()[startline+1:]
	f.close()
	layernumber=0
	sumtotals=[]
	for line in lines:
		temp=line.split(', ')
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
	lines=f.readlines()[0].split(', ')
	return lines

