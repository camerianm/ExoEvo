# Import external packages
import numpy as np
import time
import matplotlib.pyplot as plt

# Import internal packages
import evolve
import fromexo
import getall as get
import printall as prnt
from printall import Pe as Pe #print scientific notation, 4 decimal
from printall import Pf as Pf #print float, 4 decimal


# User input values:
minerals=('forsterite','fayalite','orthoenstatite','clinoenstatite','periclase', \
          'corundum','spinel','diopside','diamond','ca-al pyroxene')
fractions=[0.60, 0.10, 0.05, 0.00, 0.249, \
           0.00, 0.00, 0.00, 0.001, 0.00]
method='dynamic'      #static or dynamic thermal parameters
Mpl=1.0             #Planet mass in Me                                 Earth = 1.0
Rpl=1.0             #Planet radius in Re                               Earth = 1.0
CMF=0.3333          #Core mass fraction                                Earth = 0.3333
CRF=0.5437          #Core radius fraction                              Earth = 0.5437
Qp=1.0              #Relative heat production per kg mantle, vs Earth  Earth=1.0
Ts=300.0            #Surface temperature in K                          Earth=300.0
Tp0=2000.0           #starting mantle potential temperature in K       Earth=2000.0
dt=0.1              #timestep size, in Ga
tmax=4.55           #ending time, in Ga                                Earth=4.55


# Build your mantle and acquire its unchanging material properties.
Mp,Mc,Rp,Rc=get.SIunits(Mpl,CMF,Rpl,CRF)
d,Vm,Sa,pm,g=get.build(Mp=Mp,Mc=Mc,Rp=Rp,Rc=Rc)
c1,Ev,visc0=get.TdepVisc(minerals=minerals, fractions=fractions)

params={"Mp":Mp,"Mc":Mc,"Rp":Rp,"Rc":Rc,"d":d,"Vm":Vm,"Sa":Sa,"pm":pm,"g":g,"c1":c1,"Ev":Ev,"visc0":visc0}
prnt.unchanging(params)


# Evolve your planet over time.
Hts=[]             #A list of lists; column names are in get.keys['columns']
t=0.0              #Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.
Tp=Tp0
dt=0.01
#start=time.time()
while t <= tmax:
    alpha,cp,k=get.thermals(minerals,fractions,Tp)
    
    if method=='static':
        alpha,cp,k=get.thermals(minerals,fractions,1625)
        #alpha,cp,k,pm=3.7e-5,1250.,5.0,3340. #Uncomment for common benchmark values
        
    viscT=get.viscosity(Ev,visc0,Tp)
    Ra=get.rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k)
    
    production=evolve.produce_heat(Mp,Mc,Qp,t)
    loss=evolve.flux_heat(Sa,c1,k,Tp,d,Ra,Ev)
    dTp=(dt*get.seconds*(production-loss))/(cp*pm*Vm) #Potentially change to (cp*Mp)?
    Hts.append([t,Tp-273.15,Ra,production,loss,production/loss])
    
    Tp=Tp+dTp
    t=t+dt
#end=time.time()
#print("Program running time: ", Pf(end-start), " seconds")

Evolution=np.asarray(Hts)

# Summarize final state of planet and progression over time.
print(get.keys['columns'])
print([Pf(t), Pf(Tp), Pe(Ra), Pe(production), Pe(loss), Pf(production/loss)])

Temps=evolve.plot_heat(Evolution[:,(0,1)],"Temperature (C) vs Time (Ga)")
#Rayleighs=evolve.plot_heat(Evolution[:,(0,2)],"Rayleigh number vs Time (Ga)")
#Heat_production=evolve.plot_heat(Evolution[:,(0,3)],"Heat production (W) vs Time (Ga)")
#Heat_loss=evolve.plot_heat(Evolution[:,(0,4)],"Heat loss (W) vs Time (Ga)")



