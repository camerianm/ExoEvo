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
        
        c1=c1_default
        Ev=Ev_default
        visc0=visc0_default
        
        c1_tot=c1_tot+(c1*wt)
        Ev_tot=Ev_tot+(Ev*wt)
        visc0_tot=visc0_tot+(visc0*wt)

    return c1_tot,Ev_tot,visc0_tot

#Calculates thermal parameters: Heat capacity Cp, thermal expansivity alpha, thermal conductivity k_tot
#Units for Cp: J kg-1 K-1   Units for alpha: K-1   Units for k: W m-1 K-1
def thermals(composition,Tp):

    def berman(Tp,k):
        #Returns: Heat capacity by weight, at constant pressure of 75 GPa - units of J kg-1 K-1
        #The form of this expression (coefficients and exponents) is derived from extended form of eq.4 in Berman 1988, doi:10.1093/petrology/29.2.445
        #Coefficient values, found in mineralDB.py, are obtained using scipy.optimize curvefit function on data from Stixrude & Lithgow-Bertelloni 2011
        wtcp = k[0] + k[1]*Tp**(-0.5) + k[2]*Tp**(-2) + k[3]*Tp**(-3) + k[4]*Tp**(-1) + k[5]*Tp + k[6]*Tp**2
        return wtcp

    def alpha_coeffs(Tp,k): 
        #Returns: Coefficient of thermal expansivity at constant pressure of 75 GPa - units of K-1
        #Coefficients in this expression, found in mineralDB.py, are obtained using scipy.optimize curvefit on values obtained from Stixrude &
        #Lithgow-Bertelloni 2011. Process: V(P,T) * dV/dT
        alpha = k[0] + k[1]*Tp**(1/3) + k[2]*Tp**(-1/3) + k[3]*Tp**(-1) + k[4]*Tp**(2) + k[5]*Tp**(-3) - k[6]*Tp**(1/2)
        return alpha

	#Not yet obtained for non-olivine minerals
    k=5.0
    k_default=5.0

    cp_tot=0.0
    alpha_tot=0.0
    k_tot=0.0

    for i in composition:
        mineral=i
        wt=composition[i]

        alpha_tot = alpha_tot + wt * (alpha_coeffs(Tp,mins[mineral]['alpha']))
        cp_tot = cp_tot + wt * (berman(Tp,mins[mineral]['Cp']))# * 1000./mins[mineral]['MW']
        k_tot = k_default

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
    cr_scaling_exponent=-0.0156950415578063*(CMF-0.45547)+0.259709441585
    cr_scaling_coefficient=0.116283877635892*(CMF-0.45547)+0.704856846608
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
    # Estimates a mass for a planet whose radius and internal structure are only vaguely known
    M=(R/(1.07-0.21*CMF))**3.7
    return M

def reasonable_radius(M,CMF):
    # Estimates a radius for a planet whose mass and internal structure are only vaguely known
    R=(1.07-0.21*CMF)*(M)**(1/3.7)
    return R

def sample_masses(R):
    # Suggsts plausible masses for a given radius
    CMF_list=[0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a radius of ',R,'Re, here\'s a reasonable range of masses:')
    print('CMF\tMass(Me)')
    for m in CMF_list:
        print(str(Pf(m)),'\t',str(Pf(reasonable_mass(R,m))))
    return

def sample_radii(M):
    # Suggsts plausible radii for a given mass
    CMF_list=[0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a mass of ', M,'Me, here\'s a reasonable range of radii:')
    print('CMF\tRadius(Re)')
    for cmf in CMF_list:
        print(str(Pf(cmf)),'\t',str(Pf(reasonable_radius(M,cmf))))
    return

def build(Mpl,Rpl,Tp0):
    # derives planet properties: Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb
    print('Now building your planet. One moment...')

    if Rpl>1.50001:
        print('Oops! At',Rpl,'earth masses, your planet isn\'t likely to be rocky.')
        print('I\'ll run the largest plausibly "rocky" radius we expect - 1.5Re.')
        Rpl=1.50
    
    if Mpl>reasonable_mass(Rpl,0.95):
        print('Seems you\'ve exceeded the range of masses plausible for this radius.')
        print('I\'ll run the most massive plausibly "rocky" mass I can, at CMF=0.8, given the provided radius.')
        Mpl=reasonable_mass(Rpl,0.8)
        print('New mass:\t',str(Pf(Mpl)))
        sample_masses(Rpl)
    
    if Mpl<reasonable_mass(Rpl,0.05):
        print('This planet isn\'t dense enough to be rocky!')
        print('I\'ll run the least massive plausibly "rocky" planet I can, at CMF=0.1, given the provided radius.')
        Mpl=reasonable_mass(Rpl,0.1)
        print('New mass:\t',str(Pf(Mpl)))
        sample_masses(Rpl)

    CMF = CMF_estimate(Mpl,Rpl)  # Estimate core mass fraction from mass-radius relationships
    CRF = CRF_estimate(Mpl,CMF,Rpl) # Estimate core radius fraction from mass-radius relationships
    
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
    
    print('Done!')
    
    #Report-out
    return Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb

def SIunits(Mpl,CMF,Rpl,CRF):
    #Converts from Earth-equivalents to SI units
    Mp=Mpl*Me	#planet mass in kg
    Mc=Mp*CMF   #core mass in kg
    Rp=Rpl*Re		#planet radius in m
    Rc=Rp*CRF			#core radius in m
    return Mp,Mc,Rp,Rc

def sample_planets(Rpl,Tp0):
     # Generates a grid of planets of varying mass which are plausible for a given radius and potential temperature.
     CMF_range=np.linspace(0.05,0.95,20)
     filename='Rpl_'+str(Pf(Rpl))+'_Tp_'+str(Pf(Tp0))+'_sample_planets.csv'
     f=open(filename,'a+')
     f.write('Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb\n')
     for CMF in CMF_range:
        Mpl=reasonable_mass(Rpl,CMF)
        Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb=build(Mpl,Rpl,Tp0)
        f.write(str(Pe(Mp))+','+str(Pe(Mc))+','+str(Pe(Rp))+','+str(Pe(Rc))+','+str(Pf(d))+','+str(Pe(Vm))+','+str(Pe(Sa))+','+str(Pf(pm))+','+str(Pf(g))+','+str(Pf(Pcmb))+','+str(Pf(Tcmb))+'\n')
     print('Sample planets available at: ', filename)
     return 0