# get.py

#import sys
#import os
#import time
import numpy as np
import scipy as sp

#Constants
Grav = 6.67408e-11  #Gravitational constant
Me = 5.97e24        #Earth mass in kg
Qe = 34.0e-12       #Early Earth heat production in W/kg mantle
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
dt=0.01             #Size of timestep, in Ga

radio={'species':['238U','235U','232Th','40K'],
	'rel_amt':[1.0,1.0,1.0,1.0],	#relative to Earth's values
	'wtpercent':[0.14715,0.29649,0.10464,0.45172], #early earth values - normalized to total U
	'lambda':[0.155,0.985,0.0495,0.555]}

#Baseline returns, in order: c1, Ev, visc0
def TdepVisc(mineral):
	if mineral == 'olivine':
		c1=0.5
		Ev=300.0
		visc0=4.0e10
	else:
		print('COULD NOT OBTAIN VISCOSITY BASELINE. Assuming olivine values.')
		c1,Ev,visc0=0.5,300.0,4.0e10

	return c1,Ev,visc0

def thermals(mineral,Tp):
	if mineral == 'olivine':
		alpha=3.0e-5
		cp=1250.0
		k=5.0
	else:
		print('COULD NOT OBTAIN THERMAL BASELINE. Assuming olivine values.')
		alpha,cp,k=3.0e-5,1250.0,5.0
		
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

#def heat_production(Mp,Mc,Qpl,t):
#	Q0=(Qe*Qpl)*(Mp-Mc)       #initial radionuclide abundance total - w/kg versus Earth, times #kg mantle
#	t=0.0
#   Ht=[0.0,0.0,0.0,0.0]
#	for i in range(3):
#		Ht[i]=radio['wtpercent'][i]*np.exp(radio['lambda'][i]*t)
#	print('H')

#       
#    heat=sum(Ht)
#    return heat
