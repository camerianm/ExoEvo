import numpy as np
import getall as get
from mineralDB import minerals
import pandas as pd
from constants import *

Pe = lambda n: format(n, '.4e')
Pf = lambda n: format(n, '.4f')

def planets_from_summary():
	# Purpose: Imports a tidy-format file containing all necessary (and any optional 
    #    e.g. mineral) planetary parameters. Enables users to run planets in batch.
    # Inputs: CSV file, without headers. Each line has format of: PlanetID,Parameter,Value
    # Required parameters: alpha, CMBP, CMF, Cp, CRF, Mantle_depth, Mantle_mass, Mantle_rho, Mantle_vol,
        # Mass_kg, Mass_Me, Radius_m, Radius_Re ... minerals are optional, any # of mins works for a given ID
        # but the mineral abbreviation must match those in mineralDB.py, and end in '_percent' (e.g. Cpx_percent)
    # Outputs: a nested dictionary whose keys are the planet IDs; each planet is a dict with parameters from file.
    # Limitations: Formatting is specific, needs all params. See default cumulative.csv for sample formatting.
    # Calls: none
    # Tasks: Couple with get.build function to populate parameters for planets with missing parameters.
    # Refs: n/a
	file = 'cumulative.csv'
	f = open(file,'r')
	lines = f.readlines()
	f.close()
	files = {}
	for line in lines:
		temp = line.split(',')
		if temp[0] not in files:
			files[temp[0]] = {}
		files[temp[0]][temp[1]] = float(temp[2])

	for file in files:
		files[file]['composition'] = {}
		for mineral in minerals.keys():
			if (mineral in files[file]) and (files[file][mineral]>0):
					files[file]['composition'][mineral] = files[file][mineral]
			try:
				del files[file][mineral]
			except KeyError:
				pass
		files[file]['composition'] = get.adds_up(files[file]['composition'])
		# print("\n" + file + "\n" + "\n".join("{}: {}".format(k, v) for k, v in files[file]['composition'].items()) + "}")
	return(files)

def read_cols(file):
	f=open(file,'r')
	lines=f.readlines()[0].split(separator)
	lines[0]=lines[0][1:] #removes comment character from header
	lines[-1]=lines[-1][0:-1] #removes newline character 
	return lines

def weights_by_volume(file,startline):
	#returns: radius of each shell, relative contribution of each depth
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

def weights_by_mass(file,startline):
	#returns: mass of each shell, relative contribution of each depth to mass, relative contribution of each shell to mass
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

def bulk_mass_fraction(file,startline):
	shellmasses, weight_masses, fraction_of_mass = weights_by_mass(file,startline)
	M_runningtotal, M_section_average = find_average(file,startline,weight_masses,fraction_of_mass)
	columns=read_cols(file)[11:-2] #Excludes iron phase and mass
	bulkmassfraction=dict(zip(columns, 0.01*np.asarray(M_runningtotal[11:-2])))
	#print('Done importing composition from ExoPlex.')
	return bulkmassfraction

def build(planet,file):
    if not('Tp0' in planet.keys()):
        planet['Tp0'] = DEFAULT['Tp0']
    startline=1000
    with open(file) as f:
        for i, line in enumerate(f):
            if i == 0:
                pass
                #minerals=line.split(separator)[11:-2]
            elif i == 1:
                planet['Rp']=np.float(line.split(separator)[0])*1000
                planet['Rpl']=planet['Rp']/Re
                planet['Tcmb']=get.CMB_T(planet['Rp'],planet['Tp0'])
                planet['Sa']=4*np.pi*planet['Rp']**2
                #Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb=get.build(Mpl=Mpl,Rpl=Rpl,Tp0=Tp0)
            elif i == startline+1:
                planet['Rc']=np.float(line.split(separator)[1])*1000
                planet['d']=planet['Rp']-planet['Rc']
                planet['Vm']=(4./3.)*np.pi*(planet['Rp']**3 - planet['Rc']**3)
                planet['Pcmb']=np.float(line.split(separator)[3])
                planet['Mc']=np.float(line.split(separator)[-1])
            elif i == 3000:
                planet['Mp']=np.float(line.split(separator)[-1])
                planet['pm']=(planet['Mp']-planet['Mc'])/planet['Vm']
                planet['g']=Grav*planet['Mp']/(planet['Rp']**2)
                planet['CMF']=planet['Mc']/planet['Mp']
                planet['CRF']=planet['Rc']/planet['Rp']
    planet=thermals_from_file(planet, file, startline)
    print('Done importing ExoPlex output file.')
    return planet

def thermals_from_file(planet, file, startline):
	shellmasses, weight_masses, fraction_of_mass = weights_by_mass(file,startline)
	M_runningtotal, M_section_average = find_average(file,startline,weight_masses,fraction_of_mass)
	
	radii, vol_weights, fraction_of_volume = weights_by_volume(file,startline)
	V_runningtotal, V_section_average = find_average(file,startline,vol_weights,fraction_of_volume)

	#Get thermal conductivity from uppermost mantle
	UMradii, UMvol_weights, UMfraction_of_volume = weights_by_volume(file, 2900)
	UMV_runningtotal, UMV_section_average = find_average(file,2900,UMvol_weights,UMfraction_of_volume)
	planet['k'] = 0.0681*np.exp(0.0006 * UMV_runningtotal[8] * 1000.0) 
	# planet['k'] = get.average_property(get.adds_up(bulk_mass_fraction(file, startline)), 'k', 5.0)
	planet['alpha'] = V_runningtotal[5]
	planet['Cp'] = M_runningtotal[6]
	lith = get.adds_up(bulk_mass_fraction(file, 2900)) #lithosphere composition - last few lines of exoplex file
	# This is where viscosity params could be added, from lith proportions
	print('Done importing thermal parameters from ExoPlex.')
	return planet

def lith_rheology(file, startline):
	shellmasses, weight_masses, fraction_of_mass = weights_by_mass(file,startline)
	M_runningtotal, M_section_average = find_average(file,startline,weight_masses,fraction_of_mass)
	columns=read_cols(file)[11:-2] #Excludes iron phase and mass
	bulkmassfraction=dict(zip(columns, 0.01*np.asarray(M_runningtotal[11:-2])))
	#print('Done importing composition from ExoPlex.')
	return bulkmassfraction


#thermals = pd.read_csv(file, sep=', header=0, skiprows=startline)
#note that pandas takes much longer to load the csv for MC cases.