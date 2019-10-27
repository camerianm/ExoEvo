# Import external packages
import numpy as np
import pandas as pd

# Import internal packages
import evolve
import fromexo
import getall as get
import printall as prnt
from printall import Pe #print scientific notation, 4 decimal
from printall import Pf #print float, 4 decimal
from mineralDB import minerals
from constants import STO

#import plotly.express as px

########################################################################
# User input values:
########################################################################

# Should I import planetary params from ExoPlex? If not, my_composition is assumed.
ExoPlex='TRUE'
method='static'    # 'static', 'dynamic', or 'benchmark' thermal parameters

# Where should I send outputs from this run?
outfolder = 'OUTPUT/'
outfile = 'earth_nomantleFe.csv'

# What exoplex file should I import, and where does the mantle start?
file='earth_nomantleFe_FeMg0.9_0.07_0.9_0.09_0.9.csv'
startline=1000 #This is where the core stops and the mantle begins, in that file.

# Info about your planet - note that Mpl and Rpl are *overwritten* if ExoPlex='TRUE'
Mpl=1.0             # Planet mass in Me - usually between 0.5 and 5. ignored if ExoPlex = 'TRUE' Earth = 1.0
Rpl=1.0             # Relative heat production per kg mantle, vs Earth  Earth = 1.0
Qpl=1.0             # Planet's starting radiogenic abundance, per kg mantle
Tp0=2000.           # starting mantle potential temperature in K        Earth = 1800.0 (initial), 1650 (present)
Pref=5.0           # reference pressure for thermal calculations, in GPa. If Pref is less than 4, Pref is set to half the CMB pressure.
tmax=4.5           # ending time, in Ga - how long to cool the planet         Earth = 4.55
my_composition = {'O': 1.0}
'''
my_composition = {'C2/c':5.605938492, 'Wus':0.196424301, 'Pv':58.03824705, 'an':0.00, \
                  'O':0.249338793, 'Wad':0.072264906, 'Ring':0.028707673, 'Opx':14.88882685, \
                  'Cpx':1.099284717, 'Aki':0.0000703828958617, 'Gt_maj':9.763623743, 'Ppv':6.440039009, \
                  'CF':0.00, 'st':0.00, 'q':0.00, 'ca-pv':3.617239801, \
                  'cfs':0.00, 'coe':0.00, 'ky':0.00, 'seif':0.00}
'''
########################################################################
# End of user input values
########################################################################
columnkeys = ['time', 'temp', 'rayleigh', 'production', 'loss', 'urey', 'alpha', 'cp','viscT', 'k']
planet={'Mpl':Mpl, 'Rpl':Rpl, 'Qpl':Qpl, 'Tp0':Tp0, 'Pref':Pref}

# Composition is in weight percent. All solid solutions are represented by their Mg endmembers.
# Build your mantle and acquire its unchanging material properties.
if ExoPlex == 'TRUE':
    composition = fromexo.bulk_mass_fraction(file,startline)
    composition = get.adds_up(composition)
    planet = fromexo.build(planet=planet,file=file)
    thermals = get.thermals_at_P_ave(composition, Pref) #makes sure to get thermal conductivity
    thermals[:,1] = planet['alpha'] #replaces values in thermals with those in file
    thermals[:,2] = planet['Cp']
    thermals[:,3] = planet['k']

else: # custom composition
    planet=get.build(planet)
    composition = get.adds_up(my_composition)
    thermals=get.thermals_at_P_ave(composition, Pref)

planet['c1'],planet['Ev'],planet['visc0']=get.TdepVisc(composition)

if method == 'benchmark':
        planet.update(STO)

if Pref<4.001:
    planet['Pref'] = 0.5*planet['Pcmb']

# Evolve your planet over time.
Evolution = evolve.ThermEv(planet, thermals, method, planet['Tp0'], tmax)
#evolve.plot_heat(Evolution,method)
evolve.evolution_colorcoded(Evolution, columnkeys, 'k', 'discrete')
ev = pd.DataFrame(data=Evolution, columns=columnkeys)
ev.to_csv(outfolder+outfile)

for i in sorted(planet.keys()):
    print(i, '\t', planet[i])
print()
for i in sorted(composition.keys()):
      print(i, '\t', composition[i])
print()
print('Done! See output file: ', outfile)

