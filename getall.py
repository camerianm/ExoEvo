# A module for calculating thermodynamic parameters.
# Recommended import style (to be distinguished from "get" function for dictionaries):
#   import getall as get

import numpy as np

#Constants
Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
Ts=300.0
seconds=3.1536e16    #billion years to seconds conversion
Qe = 3.611610290257257e-11


keys={
	'minerals':['forsterite','fayalite', 'orthoenstatite','clinoenstatite','periclase',
	'corundum', 'spinel', 'diopside', 'diamond', 'ca-al pyroxene'],
	'columns':['time', 'temp', 'Ra', 'production','loss','urey']
}

#Look-up table for baseline constants used in equations for temperature-dependent viscosity.
def TdepVisc(mineral):
	if mineral == 'forsterite':   # Foley & Smye 2018; doi:10.1089/ast.2017.1695
		c1=0.5
		Ev=300.0e3
		visc0=4.0e10
	else:
		print('Could not obtain viscosity baseline values. Assuming olivine values.')
		c1,Ev,visc0=0.5,300.0e3,4.0e10
	return c1,Ev,visc0

#Calculates thermal parameters
def thermals(mineral,Tp):

	def berman(Tp,k0,k1,k2,k3,k4,k5,k6):
		molcp = k0 + k1*Tp**(-0.5) + k2*Tp**(-2) + k3*Tp**(-3) + k4*Tp**(-1) + k5*Tp + k6*Tp**2
		return molcp

	#Not yet obtained for non-olivine minerals
	alpha=3.7e-5
	k=5.0

	#k0 through k6 are derived from extended form of eq.4 in Berman 1988, doi:10.1093/petrology/29.2.445
	#k0 through k6 have units J mol-1 K-1. commented values are from 1988; uncommented are from doi:10.4095/223425
	#final 'cp' variable converts Cp from J mol-1 K-1 to J kg-1 K-1
	#alpha has units K-1
	if mineral == 'forsterite':  # Foley & Smye 2018; doi:10.1089/ast.2017.1695
		MW=140.69 #	FORSTERITE molar weight, in grams
		#k0,k1,k2,k3=238.6400,-2001.300,0.0,-116240000.0
		k0,k1,k2,k3=233.18030,-1801.580,0.000,-267937600. #fit to richet data
		k4,k5,k6=0,0,0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		#molcp = k0 + k1*Tp**(-0.5) + k2*Tp**(-2) + k3*Tp**(-3) + k4*Tp**(-1) + k5*Tp + k6*Tp**2
		cp=(1000./MW)*molcp 
		
	elif mineral == 'fayalite':
		MW=203.78
		#k0,k1,k2,k3=248.93,-1923.9,0.0,-139100000.
		k0,k1,k2,k3=251.99620,-2013.697,0.000,-62189100.
		k4,k5,k6=0.0,0.0,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'orthoenstatite':
		MW=200.78
		#k0,k1,k2,k3=166.58,-1200.6,-2270600.,279150000.
		k0,k1,k2,k3=1332.63600,-9604.704,-18164480.000,2233202400.
		k4,k5,k6=0.0,0.0,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'clinoenstatite':
		MW=200.78
		#k0,k1,k2,k3=139.96,-497.,-4400200.,535710000.
		k0,k1,k2,k3=139.95824,-497.034,-4400237.000,535708928.
		k4,k5,k6=0.0,0.0,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'periclase':
		MW=40.3
		#k0,k1,k2,k3=61.11,-296.2,-621200.,5840000.
		k0,k1,k2,k3=61.10965,-296.199,-621154.000,5844612.
		k4,k5,k6=0.0,0.0,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'corundum':
		MW=101.96
		#k0,k1,k2,k3=155.02,-828.4,-3861400.,409080000.
		k0,k1,k2,k3=155.01888,-828.387,-3861363.000,409083648.
		k4,k5,k6=0.0,0.0,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'spinel':
		MW=142.27
		k0,k1,k2,k3=235.9,-1766.6,-1710400.,40620000.
		k4,k5,k6=0.0,0.0,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'diopside':
		MW=216.55
		k0,k1,k2,k3=305.41333,-1604.931,-7165973.000,921837568.
		k4,k5,k6=0.0,0.0,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'diamond':
		MW=12.01
		k0,k1,k2,k3=24.30000,-273.400,-377400.000,0.0
		k4,k5,k6=0.0,0.006272,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	elif mineral == 'ca-al pyroxene': # CA(1)AL(2)SI(1)O(6)
		MW=218.12
		k0,k1,k2,k3=310.69775,-1671.627,-7455263.000,948781568.
		k4,k5,k6=0.0,0.006272,0.0
		molcp=berman(Tp,k0,k1,k2,k3,k4,k5,k6)
		cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	else:
		#print('COULD NOT OBTAIN THERMAL BASELINE. Assuming static values from DOI:10.1089/ast.2017.1695.')
		alpha,cp,k=3.0e-5,1250.0,5.0

	#molcp = k0 + k1*Tp**(-0.5) + k2*Tp**(-2) + k3*Tp**(-3) + k4*Tp**(-1) + k5*Tp + k6*Tp**2
	#cp=(1000./MW)*molcp #converts Cp from J mol-1 K-1 to J kg-1 K-1

	return alpha,cp,k

def SIunits(Mpl,CMF,Rpl,CRF):
	Mp=Mpl*Me	#planet mass in kg
	Mc=Mp*CMF			#core mass in kg
	Rp=Rpl*Re		#planet radius in m
	Rc=Rp*CRF			#core radius in m
	return Mp,Mc,Rp,Rc

def build(Mp,Mc,Rp,Rc):
	d=Rp-Rc				#mantle depth in m
	Vm=(4./3.)*np.pi*((Rp**3)-(Rc**3))	#mantle volume in m**3
	SA=4*np.pi*(Rp**2)	#surface area in m**2
	pm=(Mp-Mc)/Vm 		#mantle density in kg/m**3
	g=Grav*Mp/(Rp**2)	#surface gravity
	return d,Vm,SA,pm,g

def viscosity(Ev,visc0,Tp):
	visc=visc0*np.exp(Ev/(R*Tp))
	return visc

def rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k):
	Ra=(pm**2)*g*alpha*(Tp-Ts)*(d**3)*cp/(k*viscT)
	return Ra