# Import external packages
import numpy as np

# Import internal packages
import evolve
import fromexo
import getall as get
import printall as prnt
from printall import Pe #print scientific notation, 4 decimal
from printall import Pf #print float, 4 decimal
from mineralDB import minerals

########################################################################
# User input values:
########################################################################

# Should I import planetary params from ExoPlex? If not, my_composition is assumed.
ExoPlex='TRUE'

# What exoplex file should I import, and where does the mantle start?
file='earth_nomantleFe_FeMg0.9_0.07_0.9_0.09_0.9.csv'
startline=1000 #This is where the core stops and the mantle begins, in that file.

# Info about your planet - note that Mpl and Rpl are *overwritten* if ExoPlex='TRUE'
Mpl=1.0             # Planet mass in Me - usually between 0.5 and 5. ignored if ExoPlex = 'TRUE' Earth = 1.0
Rpl=1.0             # Relative heat production per kg mantle, vs Earth  Earth = 1.0
Qpl=1.0             # Planet's starting radiogenic abundance, per kg mantle
Tp0=1800.           # starting mantle potential temperature in K        Earth = 1800.0 (initial), 1650 (present)
Pref=5.0           # reference pressure for thermal calculations, in GPa. If Pref is less than 4, Pref is set to half the CMB pressure.
tmax=4.55           # ending time, in Ga - how long to cool the planet         Earth = 4.55
method='dynamic'    # 'static', 'dynamic', or 'benchmark' thermal parameters
my_composition = {'C2/c':5.605938492, 'Wus':0.196424301, 'Pv':58.03824705, 'an':0.00, \
                  'O':0.249338793, 'Wad':0.072264906, 'Ring':0.028707673, 'Opx':14.88882685, \
                  'Cpx':1.099284717, 'Aki':0.0000703828958617, 'Gt_maj':9.763623743, 'Ppv':6.440039009, \
                  'CF':0.00, 'st':0.00, 'q':0.00, 'ca-pv':3.617239801, \
                  'cfs':0.00, 'coe':0.00, 'ky':0.00, 'seif':0.00}

########################################################################
# End of user input values
########################################################################

columnkeys = ['time', 'temp', 'rayleigh', 'production', 'loss', 'urey', 'alpha', 'cp']
planet={'Mpl':Mpl, 'Rpl':Rpl, 'Qpl':Qpl, 'Tp0':Tp0, 'Pref':Pref}

# Composition is in weight percent. All solid solutions are represented by their Mg endmembers.
if ExoPlex == 'TRUE':
  composition=fromexo.bulk_mass_fraction(file,startline)
  planet=fromexo.build(planet=planet,file=file)

else: # custom composition
  composition = my_composition
  planet=get.build(planet)

if Pref<4.001:
    planet['Pref']=0.5*planet['Pcmb']

# Make sure your planet is self-consistent
composition=get.adds_up(composition)

# Build your mantle and acquire its unchanging material properties.
planet['c1'],planet['Ev'],planet['visc0']=get.TdepVisc(composition)
thermals=get.thermals_at_P_ave(composition, Pref)
prnt.unchanging(planet, composition)

#Evolve your planet over time
Tp=Tp0
Ts=300.0
dt=0.01
Hts=[]             # A list of lists; column names are in get.keys['columns']
t=0.0              # Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.

while t <= tmax:

    if method=='dynamic': alpha,cp,k=get.Tdep_thermals(thermals,Tp)
    if method=='static': alpha,cp,k=get.Tdep_thermals(thermals,1600.)
    if method=='benchmark': alpha,cp,k,planet['pm']=3.7e-5,1250.,5.0,3340. #common benchmark values
    
    viscT=get.viscosity(planet,Tp)
    Ra=get.rayleigh(planet,Tp,Ts,viscT,alpha,cp,k)
    
    production=evolve.produce_heat(planet,t)
    loss=evolve.flux_heat(planet,k,Tp,Ra)
    dTp=(dt*get.seconds*(production-loss))/(cp*planet['pm']*planet['Vm']) #Potentially change to (cp*Mp)?
    Hts.append([t,Tp,Ra,production,loss,production/loss,alpha,cp])
    
    Tp=Tp+dTp
    t=t+dt

Evolution=np.asarray(Hts)
evolve.evolution_colorcoded(Evolution, columnkeys, 'cp', 'continuous')
