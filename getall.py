# A module for calculating thermodynamic parameters.
# Recommended import style (to be distinguished from "get" function for dictionaries):
#   import getall as get

import numpy as np
from mineralDB import minerals as mins
from printall import Pe #print scientific notation, 4 decimal
from printall import Pf #print float, 4 decimal

#Constants
Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
Ts=300.0
seconds=3.1536e16    #billion years to seconds conversion
error_tolerance = 1.0e-6
Plithbase=13

keys={'columns':['time', 'temp', 'Ra', 'production','loss','urey']}

def adds_up(composition):
    subtotal=0
    for i in composition:
        subtotal=subtotal+composition[i]
    if abs(subtotal-1.0)>error_tolerance:
        print('Normalized fractional abundances: \n')
        for i in composition:
            new=composition[i]/subtotal
            composition[i]=new
            print(i,'\t\t',str(Pf(composition[i])))
    return(composition)

def TdepVisc(composition):
    # Look-up table for baseline constants used in equations for temperature-dependent viscosity.
    # Linear average of species whose activation energies are currently in dictionary. 
    # Not all species have Ev populated.
    
    # Olivine values
    c1_default=1.0
    Ev_default=300.0e3
    visc0_default=4.0e10
    
    Ev_tot=average_property(composition, 'Ev', Ev_default)
    visc0_tot=visc0_default
    c1_tot=c1_default
    return c1_tot,Ev_tot,visc0_tot

def viscosity(planet,Tp):
	visc=planet['visc0']*np.exp(planet['Ev']/(R*Tp))
	return visc

def rayleigh(planet,Tp,Ts,viscT,alpha,cp,k):
	Ra=(planet['pm']**2)*planet['g']*alpha*(Tp-Ts)*(planet['d']**3)*cp/(k*viscT)
	return Ra

def CMF_estimate(Mpl,Rpl):
    # Estimate core mass fraction from mass-radius relationships
    CMF=(1/0.21)*(1.07 - (Rpl/(Mpl**(1.0/3.7))))  #per doi:10.3847/0004-637X/819/2/127, close fit to arxiv.org/abs/1905.06530
    return CMF

def CRF_estimate(Mpl,CMF,Rpl):
    #Estimate core radius fraction from model results in https://tinyurl.com/y7wlvamw, from arxiv.org/abs/1905.06530
    CM_me=CMF*Mpl   # error margin under ~1% (of actual core radius) for most cases 
    cr_scaling_exponent=-0.0156950415578063*(CMF-0.45547)+0.259709441585
    cr_scaling_coefficient=0.116283877635892*(CMF-0.45547)+0.704856846608
    CR_re=cr_scaling_coefficient*(CM_me**cr_scaling_exponent)
    #Note that the below /mean/ scaling of core mass to core radius, across all CMF, has (on average) 5x % error with respect to modeled core radii.
    #CR_re = 0.694571363775 * (CM_me**0.268828669083)
    CRF=CR_re/Rpl
    return CRF

def CMB_T(Rp,Tp):
    #Estimate core-mantle boundary temperature, per arxiv.org/abs/1905.06530
    R=Rp/Re #below fits are in Re ; results are in K
    Tcmb = 4180*R - 2764*R**2 + 1219*R**3 + ((Tp - 1600)*(0.82 + R**1.81))
    return Tcmb

def CMB_P(Rp):
    #Estimate core-mantle boundary pressure, per arxiv.org/abs/1905.06530
    Rpl=Rp/Re #below fits are in Re ; results are in GPa
    Pcmb = 262 * Rpl - 550*Rpl**2 + 432*Rpl**3
    return Pcmb

def reasonable_mass(R,CMF):
    # Estimates a mass for a planet whose radius and internal structure are only vaguely known
    M=(R/(1.07-0.21*CMF))**3.7
    return M

def reasonable_radius(M,CMF):
    # Estimates a radius for a planet whose mass and internal structure are only vaguely known
    R=(1.07-0.21*CMF)*(M)**(1/3.7)
    return R

def sample_masses(R):
    # Suggests plausible masses for a given radius
    CMF_list=[0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a radius of ',R,'Re, here\'s a reasonable range of masses:')
    print('CMF\tMass(Me)')
    for m in CMF_list:
        print(str(Pf(m)),'\t',str(Pf(reasonable_mass(R,m))))
    return

def sample_radii(M):
    # Suggests plausible radii for a given mass
    CMF_list=[0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    print('\nFor a mass of ', M,'Me, here\'s a reasonable range of radii:')
    print('CMF\tRadius(Re)')
    for cmf in CMF_list:
        print(str(Pf(cmf)),'\t',str(Pf(reasonable_radius(M,cmf))))
    return

def build(planet): 
    # derives planet properties: Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb
    if planet['Rpl']>1.50001:
        print('Oops! At',Rpl,'earth radii, your planet isn\'t likely to be rocky.')
        print('I\'ll run the largest plausibly "rocky" radius, 1.5Re.')
        planet['Rpl']=1.50
    
    if planet['Mpl']>reasonable_mass(planet['Rpl'],0.8):
        print('The mass provided is too large to be plausible for this radius.')
        print('I\'ll run the largest plausible mass I can, at CMF=0.8, given the provided radius.')
        planet['Mpl']=reasonable_mass(planet['Rpl'],0.8)
        print('New mass:\t',str(Pf(planet['Mpl'])))
        sample_masses(planet['Rpl'])
    
    if planet['Mpl']<reasonable_mass(planet['Rpl'],0.1):
        print('This planet isn\'t dense enough to be rocky!')
        print('I\'ll run the least massive plausibly "rocky" planet I can, at CMF=0.1, given the provided radius.')
        planet['Mpl']=reasonable_mass(planet['Rpl'],0.1)
        print('New mass:\t',str(Pf(planet['Mpl'])))
        sample_masses(planet['Rpl'])

    planet['CMF'] = CMF_estimate(planet['Mpl'],planet['Rpl'])  # Estimate core mass fraction from mass-radius relationships
    planet['CRF'] = CRF_estimate(planet['Mpl'],planet['CMF'],planet['Rpl']) # Estimate core radius fraction from mass-radius relationships
    
    #Convert to SI units for further calculations
    Mp,Mc,Rp,Rc=SIunits(planet['Mpl'],planet['CMF'],planet['Rpl'],planet['CRF'])
    planet['Mp']=Mp
    planet['Mc']=Mc
    planet['Rp']=Rp
    planet['Rc']=Rc
    
    #Calculate derived values
    planet['d']=planet['Rp']-planet['Rc']				#mantle depth/thickness in m
    planet['Vm']=(4./3.)*np.pi*((planet['Rp']**3)-(planet['Rc']**3))	#mantle volume in m**3
    planet['Sa']=4*np.pi*(planet['Rp']**2)	#surface area in m**2
    planet['pm']=(planet['Mp']-planet['Mc'])/planet['Vm'] 		#mantle density in kg/m**3
    planet['g']=Grav*planet['Mp']/(planet['Rp']**2)	#surface gravity
    
    #Find core-mantle boundary adiabatic temperature and pressure, per relationships in arxiv.org/abs/1905.06530
    planet['Tcmb']=CMB_T(planet['Rp'],planet['Tp0'])
    planet['Pcmb']=CMB_P(planet['Rp'])

    return planet #Mp,Mc,Rp,Rc,d,Vm,Sa,pm,g,Pcmb,Tcmb

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
     return

def thermals_at_P_ave(composition,P):
    # Returns an array, whose columns are: T_P, alpha_P, Cp_P, and k_P
    # i.e., a temperature range, with corresponding compositionally-averaged
    # alpha, Cp, and k for a (predetermined) "average" mantle pressure.
    # Note that this smears out the thermal consequences of pressure-dependent
    # compositional heterogeneity (i.e. it treats mantle as isobaric).

    #The below can be changed if indexing scheme changes from 1 Gpa resolution

    gridstandard='alphagrid/q_alphagrid.csv' #will not be changed when composition method is converted to end-member structure

    lP_index=int(np.floor(P))
    uP_index=lP_index+2 #Only two values print due to string read-in mode, one Gpa apart.
    T_P=np.genfromtxt(gridstandard,delimiter=',',usecols=0,skip_header=1)
    nTs=len(T_P[:])
    thermals=np.zeros((nTs,4))
    thermals[:,0]=T_P
    thermals[:,3]=average_property(composition, 'k', 5.0)

    fgridstandard=open(gridstandard,'r')
    Ps=fgridstandard.readline().split(',')[lP_index:uP_index]
    loP,hiP=np.float(Ps[0]),np.float(Ps[1])
    hi_wt=(P-loP)/(hiP-loP)
    lo_wt=(hiP-P)/(hiP-loP)

    for i in composition:
        istring=str(i)
        istring=''.join(char for char in istring if char.isalnum())
        cpfile='CPgrid/'+istring+'_cpgrid.csv'
        alphafile='alphagrid/'+istring+'_alphagrid.csv'
        try:
            cp_arr=np.genfromtxt(cpfile,delimiter=',',skip_header=1,usecols=(lP_index,lP_index+1))
            alpha_arr=np.genfromtxt(alphafile,delimiter=',',skip_header=1,usecols=(lP_index,lP_index+1))
            thermals[:,1]= thermals[:,1] + (composition[i] * (alpha_arr[:,0]*lo_wt + alpha_arr[:,1]*hi_wt))
            thermals[:,2] = thermals[:,2]+ (composition[i] * ((cp_arr[:,0]*lo_wt) + (cp_arr[:,1]*hi_wt)))
        except:
            cp_arr=np.empty_like(T_P)
            alpha_arr=np.empty_like(T_P)
            idict={istring:composition[i]}
            j=0
            for Tp in T_P:
                alpha_arr[j], cp_arr[j], k_tot = thermals_75GPa(idict,Tp,0.0)
                j=j+1

            thermals[:,1] = thermals[:,1] + composition[i] * alpha_arr
            thermals[:,2] = thermals[:,2] + composition[i] * cp_arr

    return thermals

def Tdep_thermals(thermals,Tp):

    gap=thermals[1,0]-thermals[0,0]
    lT_index=int(np.floor(Tp/gap))-1
    loT=thermals[lT_index,0]
    hiT=thermals[lT_index+1,0]
    hi_wt=(Tp-loT)/(hiT-loT)
    lo_wt=(hiT-Tp)/(hiT-loT)
    thermals[lT_index,0]
    alpha=thermals[lT_index,1]*lo_wt + thermals[lT_index+1,1]*hi_wt
    Cp=thermals[lT_index,2]*lo_wt + thermals[lT_index+1,2]*hi_wt
    k=thermals[lT_index,3]*lo_wt + thermals[lT_index+1,3]*hi_wt
    return alpha,Cp,k

def representative_mantle(Rp,Rc):
    #for calculating mantle's volume-averaged properties
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

def average_property(composition, property, defaultval):
    # Linear average of species properties by weight in composition.
    # If a species doesn't have [property] in minDB.py,
    # this assumes it has defaultval.
    subtotal = 0.0
    wtused = 0.0

    for i in composition:
        mineral=i #for clarity
        wt=composition[i]
        try:
            part=mins[mineral][property]
            #wtused=wtused+wt
        except:
            part=defaultval
        subtotal=subtotal+(part*wt)
    #subtotal=Ev_tot/wtused
    return subtotal