
#Import external packages
import numpy as np
import time
import matplotlib.pyplot as plt


#Import internal packages
import evolve
import fromexo
import getall as get
import printall as prnt
from printall import Pe as Pe #print scientific notation, 4 decimal
from printall import Pf as Pf #print float, 4 decimal
import read_cumulative
'''
# Although you can design your own composition, let's start with a sample planet from ExoPlex.
file='earth_nomantleFe_FeMg0.9_0.07_0.9_0.09_0.9.csv'
startline=1000 #This is where the core stops and the mantle begins, in that file.

# User input values:
# Composition is in weight percent. All solid solutions are represented by their Mg endmembers.
exo_composition=fromexo.bulk_mass_fraction(file,startline)
composition = exo_composition


# Build-your-own-planet mode:
composition = {'C2/c':5.605938492, 'Wus':0.196424301, 'Pv':58.03824705, 'an':0.00, \
               'O':0.249338793, 'Wad':0.072264906, 'Ring':0.028707673, 'Opx':14.88882685, \
               'Cpx':1.099284717, 'Aki':0.0000703828958617, 'Gt_maj':9.763623743, 'Ppv':6.440039009, \
               'CF':0.00, 'st':0.00, 'q':0.00, 'ca-pv':3.617239801, \
               'cfs':0.00, 'coe':0.00, 'ky':0.00, 'seif':0.00}

# Make sure your planet is self-consistent
composition=get.adds_up(composition)
'''
files=read_cumulative.planets_from_summary()
Hts=[]

print("FINAL VALUES:")
print('file\t\t\t alpha \t\tCp \ttemp(K) \tRayleigh \tHeatLoss(W) \tUreyRatio')

for file in files:

    method='dynamic'    #static, dynamic, or benchmark thermal parameters
    Mpl=files[file]['Mass_Me']            #Planet mass in Me - usually between 0.5 and 5     Earth = 1.0
    Rpl=files[file]['Radius_Re']             #Relative heat production per kg mantle, vs Earth  Earth = 1.0
    Tp0=2000          #starting mantle potential temperature in K        Earth = 2000.0 (initial), 1600 (present)
    tmax=4.55           #ending time, in Ga                                Earth=4.55
    Qpl=1.0
    Pref=12.0

    planet={'Mpl':Mpl, 'Rpl':Rpl, 'Qpl':Qpl, 'Tp0':Tp0}     # Build your mantle and acquire its unchanging material properties.
    #Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb=get.build(Mpl=Mpl,Rpl=Rpl,Tp0=Tp0)
    planet['Mp']=files[file]['Mass_kg']
    planet['Mc']=planet['Mp']*files[file]['CMF']
    planet['Rp']=files[file]['Radius_m']
    planet['Rc']=planet['Rp']*files[file]['CRF']
    planet['d']=files[file]['Mantle_depth']
    planet['Vm']=files[file]['Mantle_vol']
    planet['Sa']=4*np.pi*planet['Rp']**2
    planet['pm']=files[file]['Mantle_rho']
    planet['g']=get.Grav*planet['Mp']/(planet['Rp']**2)
    planet['Pcmb']=files[file]['CMBP']
    planet['Tcmb']=get.CMB_T(planet['Rp'],planet['Tp0'])
    planet['Pref']=Pref

    composition=files[file]['composition']

    '''
    params={"Mp":Mp,"Mc":Mc,"CMF":Mc/Mp,\
            "Rp":Rp,"Rc":Rc,"CRF":Rc/Rp,\
            "d":d,"Vm":Vm,"Sa":Sa,\
            "pm":pm,"g":g,"Pcmb":Pcmb,\
            "c1":c1,"Ev":Ev,"visc0":visc0}

    prnt.unchanging(params)
    '''

    # Compare differing bulk compositions.
    Tp=Tp0
    Ts=300.0
    dt=0.01
    #Hts=[]             #A list of lists; column names are in get.keys['columns']
    t=0.0              #Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.


    planet['c1'],planet['Ev'],planet['visc0']=get.TdepVisc(composition)
    thermals=get.thermals_at_P_ave(composition,planet['Pref'])

    start=time.time()
    while t <= tmax:

        if method=='dynamic': alpha,cp,k=get.Tdep_thermals(thermals,Tp)
        if method=='static': alpha,cp,k=get.Tdep_thermals(thermals,1625)

        viscT=get.viscosity(planet,Tp)
        Ra=get.rayleigh(planet,Tp,Ts,viscT,alpha,cp,k)
        
        production=evolve.produce_heat(planet,t)
        loss=evolve.flux_heat(planet,k,Tp,Ra)
        dTp=(dt*get.seconds*(production-loss))/(cp*planet['pm']*planet['Vm']) #Potentially change to (cp*Mp)?
        Hts.append([t,Tp,Ra,production,loss,production/loss])
        
        Tp=Tp+dTp
        t=t+dt

    print(file, '\t', Pe(alpha), Pf(cp), '\t', Pf(Tp), '\t', Pe(Ra), '\t',  Pe(loss),  '\t', Pf(production/loss))

    #Evolution=np.asarray(Hts)

    #Temps=evolve.plot_heat(Evolution[:,(0,1)],"Temperature (K) vs Time (Ga)")
    #Rayleighs=evolve.plot_heat(Evolution[:,(0,2)],"Rayleigh number vs Time (Ga)")
    #Heat_loss=evolve.plot_heat(Evolution[:,(0,4)],"Heat loss (W) vs Time (Ga)")
    #Heat_production=evolve.plot_heat(Evolution[:,(0,3)],"Heat production (W) vs Time (Ga)")
    #plt.savefig('radiogenic.pdf')#Rayleighs=evolve.plot_heat(Evolution[:,(0,2)],"Rayleigh number vs Time (Ga)")

Evolution=np.asarray(Hts)
Temps=evolve.plot_heat(Evolution[:,(0,1)],"Temperature (K) vs Time (Ga)")
