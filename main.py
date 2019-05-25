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
from mineralDB import minerals as minerals 

# Although you can design your own composition, let's start with a sample planet from ExoPlex.
file='test_exoplex_file_Ca0.05_Si0.954_Al0.06_Fe1.0.csv'
startline=1000 #This is where the core stops and the mantle begins, in that file.

# User input values:
# Composition is in weight percent. All solid solutions are represented by their Mg endmembers.
exo_composition=fromexo.bulk_mass_fraction(file,startline)
composition = exo_composition

''' # Build-your-own-planet mode:
composition = {'C2/c':5.605938492, 'Wus':0.196424301, 'Pv':58.03824705, 'an':0.00, \
               'O':0.249338793, 'Wad':0.072264906, 'Ring':0.028707673, 'Opx':14.88882685, \
               'Cpx':1.099284717, 'Aki':0.0000703828958617, 'Gt_maj':9.763623743, 'Ppv':6.440039009, \
               'CF':0.00, 'st':0.00, 'q':0.00, 'ca-pv':3.617239801, \
               'cfs':0.00, 'coe':0.00, 'ky':0.00, 'seif':0.00}
'''

method='dynamic'    #static, dynamic, or benchmark thermal parameters
Mpl=1.0            #Planet mass in Me - usually between 0.5 and 5     Earth = 1.0
Rpl=1.0             #Relative heat production per kg mantle, vs Earth  Earth = 1.0
Tp0=2000          #starting mantle potential temperature in K        Earth = 2000.0 (initial), 1600 (present)
tmax=4.55           #ending time, in Ga                                Earth=4.55
Qp=1.0

composition=get.adds_up(composition)
CMF = get.CMF_estimate(Mpl,Rpl)
try_CRF = get.CRF_estimate(Mpl,CMF,Rpl)


# Build your mantle and acquire its unchanging material properties.
Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb=get.build(Mpl=Mpl,Rpl=Rpl,Tp0=Tp0)
c1,Ev,visc0=get.TdepVisc(composition)

params={"Mp":Mp,"Mc":Mc,"CMF":Mc/Mp,\
        "Rp":Rp,"Rc":Rc,"CRF":Rc/Rp,\
        "d":d,"Vm":Vm,"Sa":Sa,\
        "pm":pm,"g":g,"Pcmb":Pcmb,\
        "c1":c1,"Ev":Ev,"visc0":visc0}

prnt.unchanging(params)



#Evolve your planet over time
Tp=Tp0
Ts=300.0
dt=0.01
Hts=[]             #A list of lists; column names are in get.keys['columns']
t=0.0              #Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.

print("FINAL VALUES:")
print('time\t temp(C) \tRayleigh \tRadioHeat(W)\tHeatLoss(W) \tUreyRatio')

while t <= tmax:

    if method=='dynamic': alpha,cp,k=get.thermals(composition,Tp)
    if method=='static': alpha,cp,k=get.thermals(composition,1625)
    if method=='benchmark': alpha,cp,k,pm=3.7e-5,1250.,5.0,3340. #common benchmark values
    
    viscT=get.viscosity(Ev,visc0,Tp)
    Ra=get.rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k)
    
    production=evolve.produce_heat(Mp,Mc,Qp,t)
    loss=evolve.flux_heat(Sa,c1,k,Tp,d,Ra,Ev)
    dTp=(dt*get.seconds*(production-loss))/(cp*pm*Vm) #Potentially change to (cp*Mp)?
    Hts.append([t,Tp-273.15,Ra,production,loss,production/loss])
    
    Tp=Tp+dTp
    t=t+dt

Evolution=np.asarray(Hts)

print(Pf(t), '\t', Pf(Tp), '\t', Pe(Ra), '\t',  Pe(production), '\t', Pe(loss),  '\t', Pf(production/loss))
print()

Temps=evolve.plot_heat(Evolution[:,(0,1)],"Temperature (C) vs Time (Ga)")


#Compare approaches to calculating Cp and alpha
methods=('dynamic','static','benchmark')    #static or dynamic thermal parameters
Tp=Tp0
Ts=300.0
dt=0.01
Hts=[]             #A list of lists; column names are in get.keys['columns']
t=0.0              #Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.

print("FINAL VALUES:")
print('method\t\t temp(C) \tRayleigh \tHeatLoss(W) \tUreyRatio')

for method in methods:
    
    Tp=Tp0
    Ts=300.0
    dt=0.01
    t=0.0  
    
    while t <= tmax:
        alpha,cp,k=get.thermals(composition,Tp)

        if method=='dynamic': alpha,cp,k=get.thermals(composition,Tp)
        if method=='static': alpha,cp,k=get.thermals(composition,1625)
        if method=='benchmark': alpha,cp,k=3.7e-5,1250.,5.0 #common benchmark values

        viscT=get.viscosity(Ev,visc0,Tp)
        Ra=get.rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k)

        production=evolve.produce_heat(Mp,Mc,Qp,t)
        loss=evolve.flux_heat(Sa,c1,k,Tp,d,Ra,Ev)
        dTp=(dt*get.seconds*(production-loss))/(cp*pm*Vm) #Potentially change to (cp*Mp)?
        Hts.append([t,Tp-273.15,Ra,production,loss,production/loss])

        Tp=Tp+dTp
        t=t+dt
        
    print(method, '  \t', Pf(Tp), '\t', Pe(Ra), '\t',  Pe(loss),  '\t', Pf(production/loss))

Evolution=np.asarray(Hts)

Temps=evolve.plot_heat(Evolution[:,(0,1)],"Temperature (C) vs Time (Ga)")




# Compare differing bulk compositions.
Tp=Tp0
Ts=300.0
dt=0.01
Hts=[]             #A list of lists; column names are in get.keys['columns']
t=0.0              #Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.

print("FINAL VALUES:")
print('mineral\t temp(C) \tRayleigh \tHeatLoss(W) \tUreyRatio')

for m in minerals:
    onemin = {m:1.0}
    composition = onemin

    Tp=Tp0
    Ts=300.0
    dt=0.01
    t=0.0  
    
    start=time.time()
    while t <= tmax:
        alpha,cp,k=get.thermals(composition,Tp)

        if method=='dynamic': alpha,cp,k=get.thermals(composition,Tp)
        if method=='static': alpha,cp,k=get.thermals(composition,1625)

        viscT=get.viscosity(Ev,visc0,Tp)
        Ra=get.rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k)

        production=evolve.produce_heat(Mp,Mc,Qp,t)
        loss=evolve.flux_heat(Sa,c1,k,Tp,d,Ra,Ev)
        dTp=(dt*get.seconds*(production-loss))/(cp*pm*Vm) #Potentially change to (cp*Mp)?
        Hts.append([t,Tp-273.15,Ra,production,loss,production/loss])

        Tp=Tp+dTp
        t=t+dt
    
    print(m, '\t', Pf(Tp), '\t', Pe(Ra), '\t',  Pe(loss),  '\t', Pf(production/loss))

Evolution=np.asarray(Hts)

Temps=evolve.plot_heat(Evolution[:,(0,1)],"Temperature (C) vs Time (Ga)")
Rayleighs=evolve.plot_heat(Evolution[:,(0,2)],"Rayleigh number vs Time (Ga)")
Heat_loss=evolve.plot_heat(Evolution[:,(0,4)],"Heat loss (W) vs Time (Ga)")
#Heat_production=evolve.plot_heat(Evolution[:,(0,3)],"Heat production (W) vs Time (Ga)")


