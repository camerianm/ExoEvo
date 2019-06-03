cpgrid_README.txt

The attached grid of heat capacities (by weight) is organized as follows (where P=pressure in GPa, T=temperature in K):

NaN, P1, P2,... P_N
T1, Cp(T1,P1), Cp(T1,P2),... Cp(T1,P_N)
T2, Cp(T2,P1), Cp(T2,P2),... Cp(T2,P_N)
...
T_M, Cp(T_M,P1), Cp(T200,P2), Cp(T_M,P_N)

That the grids are spaced evenly, from 10K to 2990K (rows) and from 1.0GPa to 649.0GPa (columns), with differences of 1.0eN (N=0, for pressure in GPa) and 1.0eM (M=1, for temperature in K) allows for minimal array size and use, with easy access and indexing on an as-needed basis. The function Cps_from_grid shows more about the trimming process.


Values were obtained using ENKIPortal's ThermoEngine tool, which serves as an interface for various thermodynamic databases. ThermoEngine is not yet available to the public (that's in the works!), but ENKIPortal users may access it on the ENKI computing server, at https://enki.ofm-research.org/. Documentation is available at https://enki-portal.gitlab.io/ThermoEngine/.

The database used here is Stixrude & Lithgow-Bertelloni (2011).



# --------------------------------------------------------
#   Code used to generate the grid files, in ENKIportal:
# --------------------------------------------------------

import numpy as np
from thermoengine import phases
from thermoengine import model

#Designate database as Stixrude & Lithgow-Bertelloni (2011)
modelDBStix = model.Database(database='Stixrude')

#Uncomment to see what pure phases are available.
#print(modelDBStix.phase_info[15:]) 

#Uncomment to show solid solutions.
#print(modelDBStix.phase_info)

'''
Note on solid solutions: The below dictionary designates a single - and usually Mg-bearing - end member for each solid solution in ExoPlex output. This is because the primary science case for which this was developed investigated planets without Fe in their mantles - and most (simpler) exoplanet models assume an Mg-silicate mantle. However, ExoPlex output files include Ca and Al bearing phases, if given a non-zero Al/Mg and Ca/Mg ratio as input. Those phases aren't included here, but you could theoretically edit this to add end-members, so the format of an entry is 'Endmember': 'DB_abbrev_for_Endmember'. NOTE that Wus =/= Wus in this dictionary, but periclase (Mg endmember of ferropericlase solid solution series).

'''

mins={'C2/c':'hpcEn', 'Wus':'Per', 'Pv':'MgPrv', 'an':'An', 'O':'Fo', 'Wad':'MgWds', 'Ring':'MgRwd', \
           'Opx':'En', 'Cpx':'cEn', 'Aki':'MgAki', 'Gt_maj':'Maj', 'Ppv':'MgPpv', 'CF':'hpcEn', \
           'st':'Sti', 'q':'Qz', 'ca-pv':'CaPrv', 'cfs':'MgCf', 'coe':'Coe', 'ky':'Ky', 'seif':'Seif'}

# More non-Fe silicates are available to assume as present in mantles w/o Fe, but they are not included in the above list.
#      Aki phases (MgAki, AlAki)
#      Opx phases (4 - En, oDi, MgTs?)
#      Cpx phases (5 - cEn, Di, CaTs)
#      Maj phases (both Maj and NaMaj have Al in them)

#Change the below '200' to a larger number for a finer-resolution grid
T=np.arange(10,3000,10) 
P=np.arange(10000,6500000,10000) 

#Generate a grid for each mineral.
for item in mins:
    minstring=str(item)
    minstring=''.join(char for char in minstring if char.isalnum())
    fname=minstring+'_cpgrid.csv'
    f=open(fname,'a+')
    mineral=modelDBStix.get_phase(mins[item])
    molwt=mineral.props['molwt'][0]/1000 #grams to kg
    P_Gpa=P/(10000)
    f.write('NaN,'+(','.join(map(str, P_Gpa)))+'\n')
    for Ti in T:
        heatcap=(mineral.heat_capacity(Ti, P, const='P'))/molwt #J/K*kg
        heatcap=','.join(map(str, heatcap))
        f.write(str(Ti)+','+heatcap+'\n')
    f.close()

