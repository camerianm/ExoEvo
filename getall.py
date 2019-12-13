# A module for calculating thermodynamic parameters.
# Recommended import style (to be distinguished from "get" function for dictionaries):
#   import getall as get
from constants import *
import numpy as np
from mineralDB import minerals as mins
Pe = lambda n: format(n, '.4e')
Pf = lambda n: format(n, '.4f')
from scipy import stats


def adds_up(composition):
    # Purpose: weighted averaging schemes are vital in this code. adds_up protects against user error
    #    by checking that the sum of all weights is 1. This eases iterating over compositional spaces.
    # Inputs: dictionary containing fractional or percentage weight abundances of mineral phases
    # Outputs: the same dictionary, but normalized such that sum of phase proportions is 1
    # Limitations: may be picky w.r.t. error tolerance; % amts are recast as fractions
    # Calls: none
    # Tasks: n/a
    # Refs: n/a
    subtotal = 0
    for i in composition:
        subtotal = subtotal + composition[i]
    if abs(subtotal-1.0) > error_tolerance:
        for i in composition:
            new = composition[i] / subtotal
            composition[i] = new
            # if verbose == "true":
            #     print(i,'\t\t',Pf(composition[i]))
    return(composition)


def TWdepVisc(planet):
    # Purpose: Establish constants to be used in Arrhenius expressions for temperature-dependent viscosity.
    # Draws these from 2 sigma (95%CI) about a normal distribution about Jain, et al. 2019 for OL-WA1.
    # Not used for most cases.
    # Inputs: dictionary containing fractional or percentage weight abundances of mineral phases
    # Outputs: 3 floats: a constant, mass-averaged activation energy for diffusion creep, viscosity baseline
    # Limitations: If Ev isn't populated in minDB.py (common), a default value is assumed, possibly skewing results.
    #   Viscosity prefactors and activation energies are subject to considerable experimental uncertainty,
    #   and translating from stress/strain expt results to (strictly) temp-dependent viscosity form can present issues.
    #  Accommodating grain size may ease this, but hasn't been implemented. So for now, olivine is largely hardcoded.
    #   C1=0.5 in Foley & Smye 2017, considered flux at base of a modeled stagnant lid.
    #   Uncertain how should play out for implied flux from surface w/o considering lid.
    #   Ev is populated for some pyroxene phases, but diffusion creep literature is lacking. Dislocation abundant.
    #   Viscosity prefactors required for Arrhenius form (A*exp(Ev/RT)), but numerical studies are normalized to present viscosity.
    #   "Wet" and "dry" diffusion creep activation energies differ. "Dry" values are in minDB.py.
    # Calls: average_property
    # Tasks: figure out incorporation of prefactors for Cpx-dominated planets
    #   consider implementation of grain size and/or pressure, to accommodate rheological data whose prefactors depend on this
    #       add wet values to minDB.py and/or use balance of wet/dry to calculate weighted Ev_tot.
    # Refs: DOI: 10.1089/ast.2017.1695
    # Ev_tot = average_property(composition, 'Ev', Ev_default)
    # visc0_tot = average_property(composition, 'visc0', visc0_default)
    gs = 100.0 #grain size in microns - arbitrary right now
    c1 = 0.5
    Ev = 1000.0 * stats.truncnorm.rvs(-2, 2, loc=364., scale=61., size=1)
    visc0 = (1.00e6 * (10.0**stats.truncnorm.rvs(-2, 2, loc=5.49, scale=0.76, size=1)) * #prefactor
             (planet['Water']*1.0e2*1.0e4)**stats.truncnorm.rvs(-2,2,loc=0.85, scale=0.25, size=1) * #wt fraction water to ppm H/Si
             gs ** stats.truncnorm.rvs(-2,2,loc=1.74, scale=0.13)) #grain size
    flowparams = {'c1': c1, 'Ev': Ev[0], 'visc0': visc0[0]}
    planet.update(flowparams)
    return planet

def viscosity(planet,Tp):
    # Purpose: Arrhenius expression for strictly temperature-dependent viscosity given diffusion creep
    # Inputs: planet dictionary, potential temperature
    # Outputs: 1 float: viscosity, in Pa * s
    # Limitations: all those listed in TdepVisc; note that visc also affected by pressure, grain size, and water content
    # Calls: n/a
    # Tasks: consider implementation of grain size and/or pressure, to accommodate rheological data whose prefactors depend on this
    # Refs: DOI: 10.1089/ast.2017.1695
    for i in ['Ev', 'visc0', 'c1']: #if user hasn't provided viscosity law data, assume default case - i.e., visc0=scaled viscosity.
        if not(i in planet.keys()):
            planet[i]=DEFAULT[i]
    if not('scaletemp' in planet.keys()): #if no scaling temp is provided...
        if (planet['visc0']<1.0e13): #assume prefactor isn't normalized to a temp
            visc = planet['visc0'] * np.exp(planet['Ev']/(R * Tp)) #and get viscosity
        else: #but if the prefactor was too high for that to be reasonable, assume it's a viscosity scaled to default scaletemp
            planet['scaletemp'] = DEFAULT['scaletemp']
            visc = planet['visc0'] * np.exp((planet['Ev']/(R * Tp)) - (planet['Ev']/(R*planet['scaletemp'])))
    else:
        visc = planet['visc0'] * np.exp((planet['Ev']/(R * Tp)) - (planet['Ev']/(R*planet['scaletemp'])))
    return visc

def rayleigh(planet,Tp,Ts,viscT):
    # Purpose: calculate Rayleigh number, representing vigor of convection relative to conductive heat transport
    # Inputs: planet dictionary, plus supporting potentially time-variant float parameters
    # Outputs: 1 float: Rayleigh number, non-dimensional
    # Limitations: uses surface temp and potential temp; does not explicitly treat internal heating rate
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.1029/164GM03
    Ra = (planet['pm'] ** 2) * planet['g'] * planet['alpha'] * (Tp-Ts) * (planet['d'] ** 3) * planet['Cp']/(planet['k'] * viscT)
    return Ra


def CMF_estimate(Mpl,Rpl):
    # Purpose: estimate core mass fraction from planet's given mass and radius; empirical fit of hydrostatic model
    # Inputs: floats: planetary mass and radius, in Me and Re (Earth mass and Earth radius)
    # Outputs: 1 float: core mass fraction, non-dimensional
    # Limitations: assumes solid material throughout, i.e. no volatile envelope; assumes rocky, i.e. no mantle stripping
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.3847/0004-637X/819/2/127, close to DOI: 10.1029/2018JE005844
    CMF = (1/0.21) * (1.07 - (Rpl/(Mpl ** (1.0/3.7))))  
    return CMF


def CRF_estimate(Mpl,CMF,Rpl):
    # Purpose: estimate core radius fraction based on the total mass of core and amount of overlying material
    # Inputs: floats: planetary mass, core mass fraction, and planetary radius; mass and radius in Me and Re
    # Outputs: 1 float: core radius fraction, non-dimensional
    # Limitations: empirical fit to core masses and radii from 1-D hydrostatic planet model, from 2 CMFs across many planet masses.
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.1029/2018JE005844, data from which fits are obtained at https://tinyurl.com/y7wlvamw
    CM_me = CMF * Mpl   # error margin under ~1% (of actual core radius) for most cases 
    cr_scaling_exponent = -0.0156950415578063 * (CMF-0.45547) + 0.259709441585
    cr_scaling_coefficient = 0.116283877635892 * (CMF-0.45547) + 0.704856846608
    CR_re = cr_scaling_coefficient * (CM_me ** cr_scaling_exponent)
    # Note that the below /mean/ scaling of core mass to core radius, across all CMF, has (on average) 5x % error w.r.t. modeled core radii.
    # CR_re = 0.694571363775 * (CM_me ** 0.268828669083)
    CRF = CR_re/Rpl
    return CRF


def CMB_T(Rp,Tp):
    # Purpose: estimates (adiabatic) core-mantle boundary temperature given a potential temperature and radius
    # Inputs: 2 floats: planet radius (in m, then converted to Re) and planet potential temperature (in K)
    # Outputs: 1 float: core-mantle boundary temperature, in K
    # Limitations: adiabatic! reported fit to results of 1-D hydrostatic planet model
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.1029/2018JE005844
    R = Rp/Re
    Tcmb = 4180 * R - 2764 * R ** 2 + 1219 * R ** 3 + ((Tp - 1600) * (0.82 + R ** 1.81))
    return Tcmb


def CMB_P(Rp):
    # Purpose: estimates (adiabatic) core-mantle boundary pressure given a planet radius
    # Inputs: 1 float: planet radius (in m) which is then converted to Re
    # Outputs: 1 float: core-mantle boundary pressure, in GPa
    # Limitations: adiabatic! reported fit to results of 1-D hydrostatic planet model
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.1029/2018JE005844
    Rpl = Rp/Re
    Pcmb = 262 * Rpl - 550 * Rpl ** 2 + 432 * Rpl ** 3
    return Pcmb


def reasonable_mass(R,CMF):
    # Purpose: suggests a reasonable (starting-point) mass for a rocky planet of a given internal structure
    # Inputs: 2 floats: planet radius (in Re) and core mass fraction (non-dimensional)
    # Outputs: 1 float: planet mass, in Me
    # Limitations: 2-layer model fit, based on PREM model of Earth structure; does not demand hydrostatic equilibrium
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.3847/0004-637X/819/2/127
    M = (R/(1.07-0.21 * CMF)) ** 3.7
    return M


def reasonable_radius(M,CMF):
    # Purpose: suggests a (starting-point) radius for a rocky planet of known mass, given an internal structure
    # Inputs: 2 floats: planet mass (in Me) and core mass fraction (non-dimensional)
    # Outputs: 1 float: planet radius, in Re
    # Limitations: 2-layer model fit, based on PREM model of Earth structure; does not demand hydrostatic equilibrium
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.3847/0004-637X/819/2/127
    R = (1.07-0.21 * CMF) * (M) ** (1/3.7)
    return R


def sample_masses(R):
    # Purpose: suggests a range of masses for a plausibly rocky planet, given an observed radius (e.g. transiting planets)
    # Inputs: 1 float: planet radius, in Re
    # Outputs: a print statement for the user, to guide parameter selection
    # Limitations: 2-layer model fit, based on PREM model of Earth structure; does not demand hydrostatic equilibrium
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.3847/0004-637X/819/2/127
    CMF_list = [0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a radius of ',R,'Re, here\'s a reasonable range of masses:')
    print('CMF\tMass(Me)')
    for m in CMF_list:
        print(Pf(m),'\t',Pf(reasonable_mass(R,m)))
    return


def sample_radii(M):
    # Purpose: suggests a (starting-point) radius for a planet of constrained mass (e.g. RV planets)
    # Inputs: 1 float: planet mass (in Me)
    # Outputs: a print statement for the user, to guide parameter selection
    # Limitations: 2-layer model fit, based on PREM model of Earth structure; does not demand hydrostatic equilibrium
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.3847/0004-637X/819/2/127
    CMF_list = [0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a mass of ', M,'Me, here\'s a reasonable range of radii:')
    print('CMF\tRadius(Re)')
    for cmf in CMF_list:
        print(Pf(cmf),'\t',Pf(reasonable_radius(M,cmf)))
    return


def build(planet): 
    # Purpose: estimates planet properties given sparse inputs (mass in Me, radius in Re)
    # Inputs: dictionary of fixed planet parameters
    # Outputs: the same dictionary, but with additional derived and/or dimensional properties
    # Limitations: constrains user inputs to those planets one could consider "rocky"
    # Calls: n/a
    # Tasks: n/a
    # Refs: DOI: 10.1029/2018JE005844, with conclusions drawn from (e.g.) DOI: 10.1088/0004-637X/801/1/41
    if planet['Rpl']>1.50001:
        print('Oops! At',Rpl,'earth radii, your planet isn\'t likely to be rocky.')
        print('I\'ll run the largest plausibly "rocky" radius, 1.5Re.')
        planet['Rpl'] = 1.50
    
    if planet['Mpl']>reasonable_mass(planet['Rpl'],0.8):
        print('The mass provided is too large to be plausible for this radius.')
        print('I\'ll run the largest plausible mass I can, at CMF = 0.8, given the provided radius.')
        planet['Mpl'] = reasonable_mass(planet['Rpl'],0.8)
        print('New mass:\t',str(Pf(planet['Mpl'])))
        sample_masses(planet['Rpl'])
    
    if planet['Mpl']<reasonable_mass(planet['Rpl'],0.1):
        print('This planet isn\'t dense enough to be rocky!')
        print('I\'ll run the least massive plausibly "rocky" planet I can, at CMF = 0.1, given the provided radius.')
        planet['Mpl'] = reasonable_mass(planet['Rpl'],0.1)
        print('New mass:\t',str(Pf(planet['Mpl'])))
        sample_masses(planet['Rpl'])

    planet['CMF'] = CMF_estimate(planet['Mpl'],planet['Rpl'])  # Estimate core mass fraction from mass-radius relationships
    planet['CRF'] = CRF_estimate(planet['Mpl'],planet['CMF'],planet['Rpl'])  # Estimate core radius fraction from mass-radius relationships
    planet = SIunits(planet)
    
    # Calculate derived values
    planet['d'] = planet['Rp']-planet['Rc']                # mantle depth/thickness in m
    planet['Vm'] = (4./3.) * np.pi * ((planet['Rp'] ** 3)-(planet['Rc'] ** 3))    # mantle volume in m ** 3
    planet['Sa'] = 4 * np.pi * (planet['Rp'] ** 2)    # surface area in m ** 2
    planet['pm'] = (planet['Mp']-planet['Mc'])/planet['Vm']         # mantle density in kg/m ** 3
    planet['g'] = Grav * planet['Mp']/(planet['Rp'] ** 2)    # surface gravity
    
    # Find core-mantle boundary adiabatic temperature and pressure, per relationships in arxiv.org/abs/1905.06530
    planet['Tcmb'] = CMB_T(planet['Rp'],planet['Tp0'])
    planet['Pcmb'] = CMB_P(planet['Rp'])

    planet['k'] = average_property(adds_up(planet['X']), 'k', 5.0)
    return planet


def SIunits(planet): #Mpl,CMF,Rpl,CRF):
    # Converts from Earth-equivalents to SI units
    planet['Mp'] = planet['Mpl'] * Me       # planet mass in kg
    planet['Mc'] = planet['Mp'] * planet['CMF']       # core mass in kg
    planet['Rp'] = planet['Rpl'] * Re       # planet radius in m
    planet['Rc'] = planet['Rp'] * planet['CRF']       # core radius in m
    return planet


def sample_planets(Rpl,Tp0):
    # Generates a grid of planets of varying mass which are plausible for a given radius and potential temperature.
    filename = 'Rpl_' + str(Pf(Rpl)) + '_Tp_' + str(Pf(Tp0)) + '_sample_planets.csv'
    CMF_range = np.linspace(0.1,0.8,20)
    f = open(filename,'w+')
    f.write('Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb\n')
    for CMF in CMF_range:
        Mpl = reasonable_mass(Rpl,CMF)
        p = {'Rpl': Rpl,'Mpl': Mpl, 'Tp0': Tp0}
        p = build(p)
        line = (Pe(p['Mp']), Pe(p['Mc']), Pe(p['Rp']), Pe(p['Rc']), Pf(p['d']), Pe(p['Vm']), Pe(p['Sa']), Pf(p['pm']), Pf(p['g']), Pf(p['Pcmb']), Pf(p['Tcmb']), '\n')
        f.write(sep.join(line))
        f.write('\n')
    print('Sample planets available at: ', filename)
    return


def thermals_at_P_ave(composition,P):
    # Returns an array, whose columns are: T_P, alpha_P, Cp_P, and k_P
    # i.e., a temperature range, with corresponding compositionally-averaged
    # alpha, Cp, and k for a (predetermined) "average" mantle pressure.
    # Note that this smears out the thermal consequences of pressure-dependent
    # compositional heterogeneity (i.e. it treats mantle as isobaric).

    # The below can be changed if indexing scheme changes from 1 Gpa resolution

    gridstandard = 'alphagrid/O_alphagrid.csv'  # will not be changed when composition method is converted to end-member structure

    lP_index = int(np.floor(P))
    uP_index = lP_index + 2  # Only two values print due to string read-in mode, one GPa apart.
    T_P = np.genfromtxt(gridstandard,delimiter=',', usecols=0, skip_header=1)
    nTs = len(T_P[:])
    thermals = np.zeros((nTs,4))
    thermals[:,0] = T_P
    thermals[:,3] = average_property(composition, 'k', 5.0)
    for i in range(len(thermals[:,3])):
        krad = (8.5*thermals[i,0]**3)/(1.0e11) # radiative portion from 10.1126/science.283.5408.1699
        thermals[i,3] = thermals[i,3] + krad
    fgridstandard = open(gridstandard,'r')
    Ps = fgridstandard.readline().split(',')[lP_index:uP_index]
    loP,hiP = np.float(Ps[0]),np.float(Ps[1])
    hi_wt = (P-loP)/(hiP-loP)
    lo_wt = (hiP-P)/(hiP-loP)

    for i in composition:
        istring = str(i)
        istring = ''.join(char for char in istring if char.isalnum())
        cpfile = 'CPgrid/' + istring + '_cpgrid.csv'
        alphafile = 'alphagrid/' + istring + '_alphagrid.csv'
        try:
            cp_arr = np.genfromtxt(cpfile, delimiter=',', skip_header=1, usecols=(lP_index,lP_index + 1))
            alpha_arr = np.genfromtxt(alphafile, delimiter=',', skip_header=1, usecols=(lP_index,lP_index + 1))
            thermals[:,1] =  thermals[:,1] + (composition[i] * (alpha_arr[:,0] * lo_wt + alpha_arr[:,1] * hi_wt))
            thermals[:,2] = thermals[:,2] +  (composition[i] * ((cp_arr[:,0] * lo_wt) + (cp_arr[:,1] * hi_wt)))
        except:
            cp_arr = np.empty_like(T_P)
            alpha_arr = np.empty_like(T_P)
            j = 0
            for Tp in T_P:
                alpha_arr[j], cp_arr[j] = alphadefault, Cpdefault
                j = j + 1

            thermals[:,1] = thermals[:,1] + composition[i] * alpha_arr
            thermals[:,2] = thermals[:,2] + composition[i] * cp_arr

    return thermals

def Tdep_thermals(thermals,Tp):
    if len(thermals) < 4: #if it's just a dictionary of only 3 values...
        ack=thermals
    else:
        gap = thermals[1,0]-thermals[0,0]
        lT_index = int(np.floor(Tp/gap))-1
        loT = thermals[lT_index,0]
        hiT = thermals[lT_index + 1,0]
        hi_wt = (Tp-loT)/(hiT-loT)
        lo_wt = (hiT-Tp)/(hiT-loT)
        thermals[lT_index,0]
        alpha = thermals[lT_index,1] * lo_wt + thermals[lT_index + 1,1] * hi_wt
        Cp = thermals[lT_index,2] * lo_wt + thermals[lT_index + 1,2] * hi_wt
        k = thermals[lT_index,3]  # right now, all are same. changes are artifacts of slight heterogeneity.
        ack = {'alpha': alpha, 'Cp': Cp, 'k': k}
    return ack #alpha,Cp,k

def representative_mantle(Rp,Rc):
    # for calculating mantle's volume-averaged properties
    nR = 100
    SAtot = 0.0
    nSA = 0
    shell_radii = np.linspace(Rp,Rc,nR)
    for r in shell_radii:
        SA = 4 * np.pi * (r ** 2)
        SAtot = SAtot + SA
        nSA = nSA + 1
    SArep = SAtot/nSA
    rrep = (SArep/(4 * np.pi)) ** (0.5)
    fracdepth = (Rp - rrep) / (Rp - Rc)
    frac_height = 1 - fracdepth
    return frac_height

def average_property(composition, property, defaultval):
    # Linear average of species properties by weight in composition.
    # If a species doesn't have [property] in minDB.py,
    # this assumes it has defaultval.
    subtotal = 0.0
    wtused = 0.0

    for i in composition:
        mineral = i # for clarity
        wt = composition[i]
        try:
            part = mins[mineral][property]
            # wtused = wtused + wt
        except:
            part = defaultval
        subtotal = subtotal + (part * wt)
    # subtotal = Ev_tot / wtused
    return subtotal

