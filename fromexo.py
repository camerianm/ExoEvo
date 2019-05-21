import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show

#import scipy as sp
#import sys
#import os
#import shutil


Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters

def weights_by_volume(file,startline): #returns: radius of each shell, relative contribution of each depth
	f=open(file,'r')					# vs the one following it, relative contribution of each shell to volume 
	lines=f.readlines()[startline:]
	f.close()
	radii=[]
	surfa=[]
	vol_cumul=[]
	weight_local=[]
	fraction_of_volume=[]
	count=0
	for line in lines:
		temp=line.split(separator)
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
		fraction_of_volume.append(Vsh/Vtot)
	vol_weights=np.asarray(weight_local)
	return radii, vol_weights, fraction_of_volume

def weights_by_mass(file,startline): #returns: mass of each shell, relative contribution of each depth to mass, relative contribution of each shell to mass
	radii, weight_local, weight_total = weights_by_volume(file,startline)
	Vtot=(4./3.)*np.pi*(max(radii)**3 - min(radii)**3)
	f=open(file,'r')	
	lines=f.readlines()[startline:]
	f.close()
	densi=[]
	shellmasses=[]
	weight_masses=[]
	fraction_of_mass=[]
	count=0
	for line in lines:
		temp=line.split(separator)
		rho=float(temp[2])*1000 #converts to kg/m3
		densi.append(rho)
	for i in range(len(weight_local[:,0])-1):
		rho1=densi[i]/(densi[i]+densi[i+1]) # radius 1 contains this much (more/less) mass per m3 than radius 2
		rho2=densi[i+1]/(densi[i]+densi[i+1]) # note that weight_local indicates this much (more/less) volume.

		mass_sum=rho1*weight_local[i,0] + rho2*weight_local[i,1]
		M1=rho1*weight_local[i,0] / mass_sum  # normalized contribution to total mass of shell, for calculating
		M2=rho2*weight_local[i,1] / mass_sum # how to weigh each mass percent
		weight_masses.append([M1,M2])

		rho_avg=(weight_local[i,0]*densi[i])+(weight_local[i,1]*densi[i+1]) # volumetric average for density itself
		shellmasses.append(rho_avg*(weight_total[i]*Vtot))
	totalmass=sum(shellmasses)
	for i in range(len(shellmasses)):
		fraction_of_mass.append(shellmasses[i]/totalmass)
	return shellmasses, weight_masses, fraction_of_mass

def find_average(file,startline,relative_fractions,fraction_of_total):
	relevant = np.genfromtxt(file, delimiter=separator, skip_header=startline+1)
	section_average = []
	rel_fractions=np.asarray(relative_fractions)
	runningtotal=np.zeros(len(list(relevant[0,:])))
	#layer=0
	for i in range(len(list(relevant[:,0]))-1):
		part1=relevant[i,:]*rel_fractions[i,0] #first weight * first item
		part2=relevant[i+1,:]*rel_fractions[i,1] #second weight * second item
		average=part1+part2 #element-wise addition - average between first and second
		contribution=average*fraction_of_total[i] #each section's contribution to total
		runningtotal=runningtotal+contribution
		section_average.append([list(average)])
	return list(runningtotal), section_average

def read_cols(file):
	f=open(file,'r')
	lines=f.readlines()[0].split(separator)
	lines[0]=lines[0][1:] #removes comment character from header
	lines[-1]=lines[-1][0:-1] #removes newline character 
	return lines

def bulk_mass_fraction(file,startline):
	shellmasses, weight_masses, fraction_of_mass = weights_by_mass(file,startline)
	M_runningtotal, M_section_average = find_average(file,startline,weight_masses,fraction_of_mass)
	columns=read_cols(file)[10:]
	bulkmassfraction=dict(zip(columns, 0.01*np.asarray(M_runningtotal[10:])))
	return bulkmassfraction

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

#To test these functions independent of ExoEvo, uncomment the below:
#file='test_exoplex_file_Ca0.05_Si0.954_Al0.06_Fe1.0.csv'
#startline=1000  # This is the line where the mantle begins in the sample file. All NANs except Fe above it.
#separator=','  # Default ExoPlex output is tab-delimited. The included sample output is in .csv format.
#print('\nTesting calculation of mass fractions:\n',bulk_mass_fraction(file,startline))

#Column 0 in sample ExoPlex output is depth. 1 is radius, 2 is density.
#Columns 5 and 6 are alpha and Cp assuming Tp=1600K. Useful if investigating static Cp and a in ExoEvo.
#Cp likely (?) needs mass-based averaging scheme. Alpha likely (?) needs volume-based averaging scheme.
#Columns 10 onward are minerals. They are in variable order, and there are a variable number of them.



