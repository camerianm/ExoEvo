# A module for calculating thermodynamic parameters.
# Recommended import style (to be distinguished from "get" function for dictionaries):
#   import getall as get

import numpy as np
from mineralDB import minerals as mins
from printall import Pe as Pe #print scientific notation, 4 decimal
from printall import Pf as Pf #print float, 4 decimal

#Constants
Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
Ts=300.0
seconds=3.1536e16    #billion years to seconds conversion
Qe = 3.611610290257257e-11
error_tolerance = 1.0e-10
Plithbase=13

old_minlist=['forsterite','fayalite', 'orthoenstatite','clinoenstatite','periclase',
	'corundum', 'spinel', 'diopside', 'diamond', 'ca-al pyroxene']

keys={'columns':['time', 'temp', 'Ra', 'production','loss','urey']}


#composition = {'C2/c':0.00, 'Wus':0.00, 'Pv':0.20, 'an':0.00, 'O':0.03, 'Wad':0.40, 'Ring':0.01, 'Opx':0.20, 'Cpx':0.10, 'Aki':0.00, 'Gt_maj':0.05, 'Ppv':0.04, 'CF':0.00, 'st':0.00, 'q':0.00, 'ca-pv':0.00, 'cfs':0.00, 'coe':0.00, 'ky':0.00, 'seif':0.00}

#print(thermals(composition,1600))

def adds_up(composition):
    subtotal=0
    for i in composition:
        subtotal=subtotal+composition[i]
    if abs(subtotal-1.0)>error_tolerance:
        print('Fractional abundances don\'t add up; normalizing for total.\nNew fractional abundances:')
        for i in composition:
            new=composition[i]/subtotal
            composition[i]=new
            print(i,'\t\t',str(Pf(composition[i])))
    return(composition)

def TdepVisc(composition):
    # Look-up table for baseline constants used in equations for temperature-dependent viscosity.
    # Currently doesn't implement different minerals, only olivine. Currently seeking rheological data.
    # Will then add to mineral dictionary file.
    
    # Olivine values
    c1_default=0.5
    Ev_default=300.0e3
    visc0_default=4.0e10
    
    c1_tot=0.0
    Ev_tot=0.0
    visc0_tot=0.0
    
    for i in composition:
        mineral=i
        wt=composition[i]
        
        c1=0.5
        Ev=300.0e3
        visc0=4.0e10
        
        c1_tot=c1_tot+(c1*wt)
        Ev_tot=Ev_tot+(Ev*wt)
        visc0_tot=visc0_tot+(visc0*wt)

    return c1_tot,Ev_tot,visc0_tot

#Calculates thermal parameters
def thermals(composition,Tp):

    def berman(Tp,k):
        molcp = k[0] + k[1]*Tp**(-0.5) + k[2]*Tp**(-2) + k[3]*Tp**(-3) + k[4]*Tp**(-1) + k[5]*Tp + k[6]*Tp**2
        return molcp

    def alpha_coeffs(Tp,k): 
        alpha = k[0] + k[1]*Tp**(1/3) + k[2]*Tp**(-1/3) + k[3]*Tp**(-1) + k[4]*Tp**(2) + k[5]*Tp**(-3) - k[6]*Tp**(1/2)
        return alpha

	#Not yet obtained for non-olivine minerals
    #alpha=3.7e-5
    k=5.0
    k_default=5.0

    cp_tot=0.0
    alpha_tot=0.0
    k_tot=0.0

    for i in composition:
        mineral=i
        wt=composition[i]

        alpha_tot = alpha_tot + wt * (alpha_coeffs(Tp,mins[mineral]['alpha']))
        cp_tot = cp_tot + wt * (berman(Tp,mins[mineral]['Cp'])) * 1000./mins[mineral]['MW']
        k_tot = k_default

        #k0 through k6 are derived from extended form of eq.4 in Berman 1988, doi:10.1093/petrology/29.2.445
        #k0 through k6 have units J mol-1 K-1. commented values are from 1988; uncommented are from doi:10.4095/223425
        #final 'cp' variable converts Cp from J mol-1 K-1 to J kg-1 K-1
        #alpha has units K-1
        '''
        if mineral == 'forsterite':  # Foley & Smye 2018; doi:10.1089/ast.2017.1695
            MW=140.69 #	FORSTERITE molar weight, in grams
            #k0,k1,k2,k3=238.6400,-2001.300,0.0,-116240000.0
            k0,k1,k2,k3=233.18030,-1801.580,0.000,-267937600. #fit to richet data
            k4,k5,k6=0,0,0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            #molcp = k0 + k1*Tp**(-0.5) + k2*Tp**(-2) + k3*Tp**(-3) + k4*Tp**(-1) + k5*Tp + k6*Tp**2
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

            k=k_default

        elif mineral == 'fayalite':
            MW=203.78
            #k0,k1,k2,k3=248.93,-1923.9,0.0,-139100000.
            k0,k1,k2,k3=251.99620,-2013.697,0.000,-62189100.
            k4,k5,k6=0.0,0.0,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        elif mineral == 'orthoenstatite':
            MW=200.78
            #k0,k1,k2,k3=166.58,-1200.6,-2270600.,279150000.
            k0,k1,k2,k3=1332.63600,-9604.704,-18164480.000,2233202400.
            k4,k5,k6=0.0,0.0,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        elif mineral == 'clinoenstatite':
            MW=200.78
            #k0,k1,k2,k3=139.96,-497.,-4400200.,535710000.
            k0,k1,k2,k3=139.95824,-497.034,-4400237.000,535708928.
            k4,k5,k6=0.0,0.0,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        elif mineral == 'periclase':
            MW=40.3
            #k0,k1,k2,k3=61.11,-296.2,-621200.,5840000.
            k0,k1,k2,k3=61.10965,-296.199,-621154.000,5844612.
            k4,k5,k6=0.0,0.0,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        elif mineral == 'corundum':
            MW=101.96
            #k0,k1,k2,k3=155.02,-828.4,-3861400.,409080000.
            k0,k1,k2,k3=155.01888,-828.387,-3861363.000,409083648.
            k4,k5,k6=0.0,0.0,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        elif mineral == 'spinel':
            MW=142.27
            k0,k1,k2,k3=235.9,-1766.6,-1710400.,40620000.
            k4,k5,k6=0.0,0.0,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        elif mineral == 'diopside':
            MW=216.55
            k0,k1,k2,k3=305.41333,-1604.931,-7165973.000,921837568.
            k4,k5,k6=0.0,0.0,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        elif mineral == 'diamond':
            MW=12.01
            k0,k1,k2,k3=24.30000,-273.400,-377400.000,0.0
            k4,k5,k6=0.0,0.006272,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=2200.

        elif mineral == 'ca-al pyroxene': # CA(1)AL(2)SI(1)O(6)
            MW=218.12
            k0,k1,k2,k3=310.69775,-1671.627,-7455263.000,948781568.
            k4,k5,k6=0.0,0.006272,0.0
            molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
            cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1
            k=k_default

        else:
            #print('COULD NOT OBTAIN THERMAL BASELINE. Assuming static values from DOI:10.1089/ast.2017.1695.')
            alpha,cp,k=3.0e-5,1250.0,5.0
        
        cp_tot=cp_tot+(cp*wt)
        '''
        #k_tot=k_tot+(k*wt)

	#molcp = k0 + k1*Tp**(-0.5) + k2*Tp**(-2) + k3*Tp**(-3) + k4*Tp**(-1) + k5*Tp + k6*Tp**2
	#cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

    #print(alpha_tot,cp_tot,k_tot)
    return alpha_tot,cp_tot,k_tot


def viscosity(Ev,visc0,Tp):
	visc=visc0*np.exp(Ev/(R*Tp))
	return visc

def rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k):
	Ra=(pm**2)*g*alpha*(Tp-Ts)*(d**3)*cp/(k*viscT)
	return Ra


def CMF_estimate(Mpl,Rpl):     # Estimate core mass fraction from mass-radius relationships
    CMF=(1/0.21)*(1.07 - (Rpl/(Mpl**(1.0/3.7))))  #per doi:10.3847/0004-637X/819/2/127, close fit to arxiv.org/abs/1905.06530
    return CMF

def CRF_estimate(Mpl,CMF,Rpl):     #Estimate core radius fraction from model results in https://tinyurl.com/y7wlvamw, from arxiv.org/abs/1905.06530
    CM_me=CMF*Mpl   # error margin under ~1% (of actual core radius) for most cases 
    cr_scaling_exponent=0.266858062163334 - 0.0156950415578063*CM_me
    cr_scaling_coefficient=0.65189302886118 + 0.116283877635892*CM_me
    CR_re=cr_scaling_coefficient*(CM_me**cr_scaling_exponent)
    #Note that the below /mean/ scaling of core mass to core radius, across all CMF, has (on average) 5x % error with respect to modeled core radii.
    #CR_re = 0.694571363775 * (CM_me**0.268828669083)
    CRF=CR_re/Rpl
    return CRF

def CMB_T(Rp,Tp): #Estimate core-mantle boundary temperature, per arxiv.org/abs/1905.06530
    R=Rp/Re #below fits are in Re ; results are in K
    Tcmb = 4180*R - 2764*R**2 + 1219*R**3 + ((Tp - 1600)*(0.82 + R**1.81))
    return Tcmb

def CMB_P(Rp): #Estimate core-mantle boundary pressure, per arxiv.org/abs/1905.06530
    Rpl=Rp/Re #below fits are in Re ; results are in GPa
    Pcmb = 262 * Rpl - 550*Rpl**2 + 432*Rpl**3
    return Pcmb

def representative_mantle(Rp,Rc): #for calculating mantle's volume-averaged properties
    nR=100
    SAtot=0.0
    nSA=0
    shell_radii=np.linspace(Rp,Rc,nR)
    for r in shell_radii:
        SA=4*np.pi*(r**2)
        SAtot=SAtot+SA
        nSA=nSA+1
    SArep=SAtot/nSA
    rrep=(SArep/(4*np.pi))**(0.5)
    fracdepth=(Rp-rrep)/(Rp-Rc)
    frac_height=1-fracdepth
    return frac_height

def reasonable_mass(R,CMF):
    M=(R/(1.07-0.21*CMF))**3.7
    return M

def reasonable_radius(M,CMF):
    R=(1.07-0.21*CMF)*(M)**(1/3.7)
    return R

def sample_masses(R):
    CMF_list=[0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a radius of ',R,'Re, here\'s a reasonable range of masses:')
    print('CMF\tMass(Me)')
    for m in CMF_list:
        print(str(Pf(m)),'\t',str(Pf(reasonable_mass(R,m))))
    return

def sample_radii(M):
    CMF_list=[0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a mass of ', M,'Me, here\'s a reasonable range of radii:')
    print('CMF\tRadius(Re)')
    for cmf in CMF_list:
        print(str(Pf(cmf)),'\t',str(Pf(reasonable_radius(M,cmf))))
    return

def build(Mpl,Rpl,Tp0):
    
    print()
    print('Now building your planet. One moment...')
    print()
    
    if Mpl>6.0:
        sample_radii(Mpl)
    
    if Rpl>1.50001:
        print('\nOops! You\'ve exceeded the range of radii where planets are likely to be rocky (<1.5Re).')
        print('I\'ll run the largest plausibly "rocky" radius I can - 1.5Re.')
        Rpl=1.50
    
    if Mpl>reasonable_mass(Rpl,0.9):
        print('\nSeems you\'ve exceeded the range of masses plausible for this radius.')
        print('I\'ll run the most massive plausibly "rocky" mass I can, at CMF=0.8, given the provided radius.')
        Mpl=reasonable_mass(Rpl,0.7)
        print('New mass:\t',str(Pf(Mpl)))
        sample_masses(Rpl)
    
    CMF = CMF_estimate(Mpl,Rpl)  # Estimate core mass fraction from mass-radius relationships
    CRF = CRF_estimate(Mpl,CMF,Rpl) # Estimate core radius fraction from mass-radius relationships

    if (CMF<0) or (CRF<0):
        print('Either increase your mass, or decrease your radius, for a reasonable structure.')
    if (CMF>1) or (CRF>1):
        print('Either decrease your mass, or increase your radius, for a reasonable structure.')
    if (abs(np.floor(CMF))>0) or (abs(np.floor(CRF))>0):
        CMF=0.33
        print('\nLet\'s first try changing the mass and keeping the current radius of',Rpl,' with an Earth-like CMF of 0.33.')
        print('If results seem weird, try changing your planet\'s mass or radius! Keep in mind that "rocky" worlds have radii <1.5Re.')        
        Mpl=reasonable_mass(Rpl,0.33)
    
    CMF = CMF_estimate(Mpl,Rpl)  # Estimate core mass fraction from mass-radius relationships - https://tinyurl.com/y7wlvamw, from arxiv.org/abs/1905.06530
    CRF = CRF_estimate(Mpl,CMF,Rpl) # Estimate core radius fraction from mass-radius relationships
    Mpl=reasonable_mass(Rpl,CMF)

    
    #Convert to SI units for further calculations
    Mp,Mc,Rp,Rc=SIunits(Mpl,CMF,Rpl,CRF)
    
    #Calculate derived values
    d=Rp-Rc				#mantle depth/thickness in m
    Vm=(4./3.)*np.pi*((Rp**3)-(Rc**3))	#mantle volume in m**3
    Sa=4*np.pi*(Rp**2)	#surface area in m**2
    pm=(Mp-Mc)/Vm 		#mantle density in kg/m**3
    g=Grav*Mp/(Rp**2)	#surface gravity
    
    #Find core-mantle boundary adiabatic temperature and pressure, per relationships in arxiv.org/abs/1905.06530
    Tcmb=CMB_T(Rp,Tp0)
    Pcmb=CMB_P(Rp)

    #Determine volume-averaged pressure - not currently implemented - helpful for calibrating against P implicit in Cp calculation
    frac_height=representative_mantle(Rp,Rc)
    Prep=Pcmb-frac_height*(Pcmb-Plithbase)
    #print('Representative mantle pressure: ', str(Pf(Prep)))
    
    print('\nDone!')
    
    #Report-out
    return Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb

def SIunits(Mpl,CMF,Rpl,CRF):
    Mp=Mpl*Me	#planet mass in kg
    Mc=Mp*CMF   #core mass in kg
    Rp=Rpl*Re		#planet radius in m
    Rc=Rp*CRF			#core radius in m
    return Mp,Mc,Rp,Rc

