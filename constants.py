# constants.py
# Description: Contains constants used in all calculations, plus benchmark values

# Common to all runs:
pi = 3.1415926535897932384626433
e = 2.71828182845904523536028747135266249775724709369995
R = 8.3144598  # Ideal gas constant
Re = 6.371e6  # Earth radius in m
Grav = 6.67408e-11  # Newtonian constant of gravitation
Me = 5.97236473e24  # Earth mass in kg
seconds = 3.1536e16  # number of s in 1 Gyr

# Customize structure and outputs:
sep = ','  # outputs are comma-separated
separator = sep
tmax = 4.55  # in Gyr
dt = 0.01  # in Gyr
verbose = "false"  # if "true": all print statements activated
error_tolerance = 1.0e-6

# -----------BEGIN BENCHMARK PARAMETERS-----------
# Schubert, Turcotte, and Olson, 2001 ; doi.org/10.1017/CBO9780511612879.014
STO = {
    'method': 'benchmark',
    'scaletemp': 1623.15,
    'urey': 0.75,
    'Ra_cr': 1100.0,
    'beta': 0.3,
    'Ts': 273.0,
    'A0': 7.0e4,
    'Ev': R*7.0e4,
    'visc0': 1.65e2, #3.084e4, #5.815e4, #
    'Tp0': 3273.0,
    'pm': 3400.0,
    'k': 4.18,
    'kappa': 1.0e-6,
    'Cp': 4.18/(1.0e-6*3400.0),
    'g': 10.0,
    'alpha': 3.0e-5,
    'd': 2.8e6,
    'Tf': 1950.0,
    'decay': 1.42e-17 * seconds,
    'Mp': Me,
    'Mc': Me-4.06e24,    #'Vm': (4./3. * pi * ((Re**3 - (Re - 2.8e6)**3)))
    }

# Group 1:
STO['Qp'] = 4.317e-14 * (STO['Mp']-STO['Mc'])*STO['Cp'] 
STO['Sa'] = STO['Cp'] * (STO['Mp']-STO['Mc']) * 1.377e-13 #4*pi*(Re)**2,
STO['Rp'] = (STO['Sa']/(4*pi))**(0.5)
STO['Rc'] = STO['Rp'] - 2.8e6
STO['Vm'] = (4./3. * pi * ((STO['Rp']**3-(STO['Rc'])**3)))

# Gets closer to 1950:
# STO['Qe'] = 34.5e-12 * (STO['Mp'] - STO['Mc'])
# STO['Sa'] = 4*pi*(Re)**2
# STO['Rc'] = Re - 2.8e6 #3485.0e3,
# STO['Rp'] = Re
# STO['Vm'] = (4./3. * pi * (STO['Rp']**3-(STO['Rc'])**3))


# For a value based solely on the Urey Ratio and present-day radiogenic abundances:
# STO['Qe'] = ((STO['urey'] * 60.0e-3 * (4*pi*Re*Re))/(e**(-1 * STO['decay'] * 4.5)))/(Me-STO['Mc'])  #4.06e24, #(4.86398 * (0.75*36.5e12))/(4.06e24)



# Seales and Lenardic 2019 - parameter space exploration
MC = {
    'method': 'MC',
    'Ra_cr': 1100.0,
    'Ts': 273.0,
    'gs': 10.0e6,
    'pm': 4000.0,
    'k': 5.0,
    'kappa': 1.0e-6,
    'Cp': 1250.,
    'c1': 1.0,
    'g': 9.8,
    'alpha': 3.0e-5,
    'd': 2.8e6,
    'Tf': 1950.0,
    'Mp': 5.97219e24,
    'Mc': 1.883e24,    #'Vm': (4./3. * pi * ((Re**3 - (Re - 2.8e6)**3)))
    'Rc': 3.4881e6    #'decay': 1.42e-17 * seconds,
    }
MC['Rp'] = MC['Rc'] + MC['d']
MC['Sa'] = 4*pi*(MC['Rp'])**2
MC['Vm'] = (4./3. * pi * ((MC['Rp']**3-(MC['Rc'])**3)))
MC['DensityDerived'] = (MC['Mp']-MC['Mc'])/MC['Vm']

# Karato and Korenaga 2008 - DOI: 10.1029/2007JB005100.
# Consider: adding 10.1029/2018JB016558
KK = {
    'method': 'KK', #Karato and Korenaga olivine rheology
    'urey': 0.75,
    'Ra_cr': 1100.0,
    'Ts': 273.0,
    'A0': 10**4.32,
    'p2': 2.56,
    'COH': 1000./1.0e6, #upper mantle 0.1% OH by weight
    'rCOH': 1.93,
    'gs': 10.0e6,
    'Ev': 387.0e3,
    'pm': 4000.0,
    'k': 5.0,
    'kappa': 1.0e-6,
    'Cp': 1250.,
    'c1': 1.0,
    'g': 9.8,
    'alpha': 3.0e-5,
    'd': 2.89e6,
    'Tf': 1950.0,
    'Mp': 5.97219e24,
    'Mc': 1.883e24,    #'Vm': (4./3. * pi * ((Re**3 - (Re - 2.8e6)**3)))
    'decay': 1.42e-17 * seconds
    }
MC['Rp'] = MC['Rc'] + MC['d']
MC['Sa'] = 4*pi*(MC['Rp'])**2
MC['Vm'] = (4./3. * pi * ((MC['Rp']**3-(MC['Rc'])**3)))
MC['DensityDerived'] = (MC['Mp']-MC['Mc'])/MC['Vm']
KK['visc0'] = 1./(KK['A0'] * KK['gs'] ** (-1 * KK['p2']) * KK['COH'] ** KK['rCOH'])

# -----------END BENCHMARK PARAMETERS-----------
