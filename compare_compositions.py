
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
file='test_exoplex_file_Ca0.05_Si0.954_Al0.06_Fe1.0.csv'
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
    Qp=0.8


    # Build your mantle and acquire its unchanging material properties.
    #Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb=get.build(Mpl=Mpl,Rpl=Rpl,Tp0=Tp0)
    Mp=files[file]['Mass_kg']
    Mc=Mp*files[file]['CMF']
    Rp=files[file]['Radius_m']
    Rc=Rp*files[file]['CRF']
    d=files[file]['Mantle_depth']
    Vm=files[file]['Mantle_vol']
    Sa=4*np.pi*Rp**2
    pm=files[file]['Mantle_rho']
    g=get.Grav*Mp/(Rp**2)
    Pcmb=files[file]['CMBP']
    Tcmb=get.CMB_T(Rp,Tp0)

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


    c1,Ev,visc0=get.TdepVisc(composition)
    thermals=get.thermals_at_P_ave(composition,0.5*Pcmb)

    start=time.time()
    while t <= tmax:

        if method=='dynamic': alpha,cp,k=get.Tdep_thermals(thermals,Tp)
        if method=='static': alpha,cp,k=get.Tdep_thermals(thermals,1625)

        viscT=get.viscosity(Ev,visc0,Tp)
        Ra=get.rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k)

        production=evolve.produce_heat(Mp,Mc,Qp,t)
        loss=evolve.flux_heat(Sa,c1,k,Tp,d,Ra,Ev)
        dTp=(dt*get.seconds*(production-loss))/(cp*pm*Vm) #Potentially change to (cp*Mp)?
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
