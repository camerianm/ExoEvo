# compare_compositions.py

# Import external packages
import numpy as np
import matplotlib.pyplot as plt
import time

# Import internal packages
import evolve
import fromexo
import getall as get
import read_cumulative
import printall as prnt
from printall import Pe  # print scientific notation, 4 decimal
from printall import Pf  # print float, 4 decimal

# For plotting purposes:
columnkeys = ['time', 'temp', 'rayleigh', 'production', 'loss', 'urey', 'Mg', 'Si', 'Ca', 'Al', 'MgSi', 'alpha', 'cp', 'Water']


# USER INPUT VALUES
method = 'dynamic'  # static or dynamic calculation of thermal parameters
Pref = 70         # Reference pressure for Cp and alpha, GPa
Ts = 300.0          # surface temperature, K
Tp0 = 2000          # starting mantle potential temperature in K        Earth = 2000.0 (initial), 1600 (present)
Qpl = 1.0           # Relative heat production per kg mantle, vs Earth  Earth = 1.0
tmax = 4.55         # ending time, in Ga                                Earth = 4.55
output_file='compare_compositions_results.csv'
output_file_2='bulk_elemental_fraction.csv'

# Import planet compositions from summary file.
files = read_cumulative.planets_from_summary()

#Establish files which will hold final values
out = open(output_file, 'w+')
out.write('file,Mg/Si,Ca/Si,Al/Si,alpha,Cp,k,temp(K),Rayleigh,HeatLoss(W),UreyRatio,WaterFrac\n')
out2 = open(output_file_2, 'w+')
out2.write('f_Mg,f_Si,f_Ca,f_Al,tempK,Water\n')

Hts = []  # A list of lists; column names are in get.keys['columns']
sep = ','
nplanets = 0

for file in files:

    Mpl = files[file]['Mass_Me']            # Planet mass in Me - usually between 0.5 and 5     Earth = 1.0
    Rpl = files[file]['Radius_Re']
    planet = {'Mpl': Mpl, 'Rpl': Rpl, 'Qpl': Qpl, 'Tp0': Tp0}
    planet['Mp'] = files[file]['Mass_kg']
    planet['Mc'] = planet['Mp'] * files[file]['CMF']
    planet['Rp'] = files[file]['Radius_m']
    planet['Rc'] = planet['Rp'] * files[file]['CRF']
    planet['d'] = files[file]['Mantle_depth']
    planet['Vm'] = files[file]['Mantle_vol']
    planet['Sa'] = 4 * np.pi * planet['Rp']**2
    planet['pm'] = files[file]['Mantle_rho']
    planet['g'] = get.Grav * planet['Mp']/(planet['Rp']**2)
    planet['Pcmb'] = files[file]['CMBP']
    planet['Tcmb'] = get.CMB_T(planet['Rp'], planet['Tp0'])
    planet['Pref'] = Pref

    # Get elemental abundances for easy plotting
    MgSi = 1 / files[file]['SiMg']
    CaSi = files[file]['CaMg'] * MgSi
    AlSi = files[file]['AlMg'] * MgSi
    totmol = MgSi + CaSi + AlSi + 1
    Mg = MgSi / totmol
    Ca = CaSi / totmol
    Al = AlSi / totmol
    Si = 1 / totmol

    # Given what your planet's made of, find out its water storage capacity and rheological constraints.
    composition = files[file]['composition']
    composition = get.adds_up(composition)
    planet['c1'], planet['Ev'], planet['visc0'] = get.TdepVisc(composition)
    thermals = get.thermals_at_P_ave(composition, planet['Pref'])
    Water = get.average_property(composition, 'water', 0.0)

    # Now start the evolutionary clock...
    # Be sure to keep Hts, Tp, and t=0.0 here, so at the end of each run, values are reset
    # and the clock starts again.
    Tp = Tp0
    dt = 0.01
    t = 0.0

    while t <= tmax:

        if method == 'dynamic': alpha, cp, k = get.Tdep_thermals(thermals, Tp)
        if method == 'static': alpha, cp, k = get.Tdep_thermals(thermals, 1625)

        viscT = get.viscosity(planet, Tp)
        Ra = get.rayleigh(planet, Tp, Ts, viscT, alpha, cp, k)

        production = evolve.produce_heat(planet, t)
        loss = evolve.flux_heat(planet, k, Tp, Ra)
        dTp = (dt * get.seconds * (production - loss)) / (cp * planet['pm'] * planet['Vm'])

        # This is what columnkey applies to
        Hts.append([t, Tp, Ra, production, loss, production / loss, Mg, Si, Ca, Al, MgSi, alpha, cp, Water])

        Tp = Tp + dTp
        t = t + dt

    nplanets = nplanets + 1
    print('|')

    line = (file, Pf(MgSi), Pf(CaSi), Pf(AlSi), Pe(alpha), Pf(cp), Pf(k), Pf(Tp), Pe(Ra), Pe(loss), Pf(production/loss), Pe(Water))
    out.write(sep.join(line))
    out.write('\n')

    line2 = (Pf(Mg), Pf(Si), Pf(Ca), Pf(Al), Pf(Tp), Pe(Water))
    out2.write(sep.join(line2))
    out2.write('\n')

print()

print('Output file names: ' + output_file + '\t' + output_file_2)
Evolution = np.asarray(Hts)
Temps = evolve.plot_heat(Evolution[:, (0, 1)], "Temperature (K) vs Time (Ga)")
# plt.savefig('radiogenic.pdf')
evolve.evolution_colorcoded(nparray=Hts, columnkeys=columnkeys, colorcolumn='Water', colortype='continuous')
evolve.evolution_colorcoded(nparray=Hts, columnkeys=columnkeys, colorcolumn='MgSi', colortype='discrete')
evolve.evolution_colorcoded(nparray=Hts, columnkeys=columnkeys, colorcolumn='Mg', colortype='continuous')
evolve.evolution_colorcoded(nparray=Hts, columnkeys=columnkeys, colorcolumn='Si', colortype='continuous')
evolve.evolution_colorcoded(nparray=Hts, columnkeys=columnkeys, colorcolumn='alpha', colortype='continuous')
evolve.evolution_colorcoded(nparray=Hts, columnkeys=columnkeys, colorcolumn='cp', colortype='continuous')

