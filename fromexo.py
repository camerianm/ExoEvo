import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import getall as get
#import scipy as sp
#import sys
#import os
#import shutil


Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
separator=','  # Default ExoPlex output is tab-delimited. The included sample output is in .csv format.

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

def bulk_mass_fraction(file,startline,frac_Fe_phases):
	from mineralDB import exonames as e

	shellmasses, weight_masses, fraction_of_mass = weights_by_mass(file,startline)
	M_runningtotal, M_section_average = find_average(file,startline,weight_masses,fraction_of_mass)
	phasefrac = 0.01*np.asarray(M_runningtotal[11:-2])
	columns=read_cols(file)[11:-2] #Excludes iron phase and mass
	fFe = frac_Fe_phases
	emfrac = []
	alts=[]

	for colno, col in enumerate(columns): #converts exoplex columns to sample pure phases
		for i in range(len(e[col])):
			alts.append(e[col][i])
		#alts.append(e[col])
		if len(e[col])>1:
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		else:
			emfrac.append(phasefrac[colno])
		'''
		if col == 'C2/c':
			alts.append('hpcEn')
			emfrac.append(phasefrac[colno])
		elif col == 'Wus':
			alts.append('Per')
			alts.append('Wus')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'Pv':
			alts.append('MgPrv')
			alts.append('FePrv')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'O':
			alts.append('Fo')
			alts.append('Fa')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'Wad':
			alts.append('MgWds')
			alts.append('FeWds')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'Ring': 
			alts.append('MgRwd')
			alts.append('FeRwd')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'Opx':
			alts.append('En')
			alts.append('Fs')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'Cpx':
			alts.append('cEn')
			emfrac.append(phasefrac[colno])
		elif col == 'Aki':
			alts.append('MgAki')
			alts.append('FeAki')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'Gt_maj':
			alts.append('Maj')
			emfrac.append(phasefrac[colno])
		elif col == 'Ppv':
			alts.append('MgPpv')
			alts.append('FePpv')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'CF':
			alts.append('MgCf')
			alts.append('FeCf')
			emfrac.append(phasefrac[colno]*(1-fFe))
			emfrac.append(phasefrac[colno]*(fFe))
		elif col == 'ca-pv':
			alts.append('CaPrv')
			emfrac.append(phasefrac[colno])
		elif col == 'cfs': 
			alts.append('hpcFs')
			emfrac.append(phasefrac[colno])
		elif col == 'Sp':
			alts.append('Spl')
			emfrac.append(phasefrac[colno])
		elif col == 'st':
			alts.append('Sti')
			emfrac.append(phasefrac[colno])
		elif col == 'q':
			alts.append('Qz')
			emfrac.append(phasefrac[colno])
		else:
			alts.append(col.capitalize())
			emfrac.append(phasefrac[colno])
		'''
		# elif col == 'an': alts.append(col)
		# elif col == 'coe': alts.append(col)
		# elif col == 'ky': alts.append(col)
		# elif col == 'seif': alts.append(col)
	bulkmassfraction=dict(zip(alts, emfrac))
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
#file='earth_nomantleFe_FeMg0.9_0.07_0.9_0.09_0.9.csv'
#startline=1000  # This is the line where the mantle begins in the sample file. All NANs except Fe above it.
#separator=','  # Default ExoPlex output is tab-delimited. The included sample output is in .csv format.
#print('\nTesting calculation of mass fractions:\n',bulk_mass_fraction(file,startline))

#Column 0 in sample ExoPlex output is depth. 1 is radius, 2 is density.
#Columns 5 and 6 are alpha and Cp assuming Tp=1600K. Useful if investigating static Cp and a in ExoEvo.
#Cp likely (?) needs mass-based averaging scheme. Alpha likely (?) needs volume-based averaging scheme.
#Columns 10 onward are minerals. They are in variable order, and there are a variable number of them.
'''
def build(Mpl,file,Tp0):
    startline=1000
    with open(file) as f:
        for i, line in enumerate(f):
            if i == 0:
                minerals=line.split(separator)[10:-1]
            if i == 1:
                Rpl=np.float(line.split(separator)[0])*1000/Re
                Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb=get.build(Mpl=Mpl,Rpl=Rpl,Tp0=Tp0)
            if i == startline+1:
                Rc=np.float(line.split(separator)[1])*1000
                d=Rp-Rc
                Vm=(4./3.)*np.pi*(Rp**3 - Rc**3)
                pm=(Mp-Mc)/Vm
                Pcmb=np.float(line.split(separator)[3])
    return Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb

'''

def build(file,Tp0):
    startline=1000
    with open(file) as f:
        for i, line in enumerate(f):
            if i == 0:
                pass
                #minerals=line.split(separator)[11:-2]
            elif i == 1:
                Rp=np.float(line.split(separator)[0])*1000
                Rpl=Rp/Re
                Tcmb=get.CMB_T(Rp,Tp0)
                Sa=4*np.pi*Rp**2
                #Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb=get.build(Mpl=Mpl,Rpl=Rpl,Tp0=Tp0)
            elif i == startline+1:
                Rc=np.float(line.split(separator)[1])*1000
                d=Rp-Rc
                Vm=(4./3.)*np.pi*(Rp**3 - Rc**3)
                Pcmb=np.float(line.split(separator)[3])
                Mc=np.float(line.split(separator)[-1])
            elif i == 3000:
                Mp=np.float(line.split(separator)[-1])
                pm=(Mp-Mc)/Vm
                g=Grav*Mp/(Rp**2)
    return Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb
