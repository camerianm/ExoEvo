# Import external packages
import numpy as np
import pandas as pd
import plotly.io as pio

# Import internal packages
import evolve
import fromexo
import getall as get
from mineralDB import minerals
from constants import *
import plot
Pe = lambda n: format(n, '.4e')
Pf = lambda n: format(n, '.4f')

########################################################################
# User input values:
########################################################################

# Should I import planetary params from ExoPlex? If not, my_composition is assumed.
ID = 'test'+method
method='dynamic'    # dynamic, benchmark, MC, or DEFAULT parameters if not provided explicitly
ExoPlex='TRUE'

# Where should I send outputs from this run?
outfolder = 'OUTPUT/'
outfile = 'earth_nomantleFe.csv'

# What exoplex file should I import, and where does the mantle start?
file='earth_nomantleFe_FeMg0.9_0.07_0.9_0.09_0.9.csv'
startline=1000 #This is where the core stops and the mantle begins, in that file.

# Info about your planet - note that Mpl and Rpl are *overwritten* if ExoPlex='TRUE'
Tp0=2000.           # starting mantle potential temperature in K        Earth = 1800.0 (initial), 1650 (present)
Mpl=1.0             # Planet mass in Me - usually between 0.5 and 5. ignored if ExoPlex = 'TRUE' Earth = 1.0
Rpl=1.0             # Relative heat production per kg mantle, vs Earth  Earth = 1.0
Qpl=1.0             # Planet's starting radiogenic abundance, per kg mantle
Pref=6.0            # reference pressure for thermal calculations, in GPa. If Pref is less than 4, Pref is set to half the CMB pressure.
tmax=4.55           # ending time, in Ga - how long to cool the planet         Earth = 4.55
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
planet={'ID': ID, 'Mpl':Mpl, 'Rpl':Rpl, 'Qpl':Qpl, 'Tp0':Tp0, 'Pref':Pref, 
     'outcols': ['ID', 'time', 'temp', 'Ra', 'H', 'Q', 'Urey', 'viscT', 
     'visc0', 'Ev', 'log10visc', 'beta']}
if Pref<4.001:
    planet['Pref'] = 0.5*planet['Pcmb']
# Composition is in weight percent. All solid solutions are represented by their Mg endmembers.
# Build your mantle and acquire its unchanging material properties.
if ExoPlex == 'TRUE':
    planet['composition'] = get.adds_up(fromexo.bulk_mass_fraction(file,startline))
    planet = fromexo.build(planet=planet,file=file)
    thermals = {'alpha': planet['alpha'], 'Cp': planet['Cp'], 'k': planet['k']} #get.thermals_at_P_ave(composition, Pref)

else: # custom composition
    planet['composition'] = get.adds_up(my_composition)
    planet=get.build(planet)
    composition = planet['X'] 
    thermals=get.thermals_at_P_ave(composition, Pref)

# Evolve your planet over time.
Evolution = evolve.ThermEv(planet, thermals, method, planet['Tp0'], tmax)
Evolution.to_csv(outfolder+planet['ID']+outfile)
p = plot.evolution_colorcoded(Evolution, 'viscT', 'continuous')
pio.write_html(p, outfolder+planet['ID']+'.html')
print('Done! See '+outfolder+planet['ID']+'.html for plot.')


exit()

for i in sorted(planet.keys()):
    print(i, '\t', planet[i])
print()
for i in sorted(composition.keys()):
      print(i, '\t', composition[i])
