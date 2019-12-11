import numpy as np
import matplotlib.pyplot as plt
import getall as get
import plotly.express as px
import pandas as pd
from constants import *
Pe = lambda n: format(n, '.4e')
Pf = lambda n: format(n, '.4f')

Ts=273.0
UreyRatio = 0.23
Qe = 1.055e-11 		# If bulk silicate earth Urey ratio (BSE) is assumed, crust is included. 16TW currently; ~80TW past; core mass excluded; Qe=2.00e-11
				 	# If not, Qe=1.055e-11. Estimated from convective Urey Ratio for present day, multiplied by 5 to simulate starting Earth values
					# i.e., Qe = ~5.00 * (0.23 * 36.5e12 W)/(M(earth)-M(core)-M(crust)). 
					# circa 4.55Ga. This would imply a closed reservoir of radionuclides in the mantle, with the total mass
					# of continental crust remaining constant.
					# If crustal growth module is added, should use this present day starting point, then add partitioning
					# coefficient for HPE in melt vs. solid.
					# Current HPE estimates from: doi:10.1029/2007RG000241
					#    Mcrust=0.006Mmantle ; Mcore=0.33*Mearth.

radio = np.array([
	# '238U','235U','232Th','40K'; 
	[1.00, 1.00, 1.00, 1.00],  # 'rel_amt':relative to Earth's values
	#[0.15363, 0.2837, 0.11046, 0.45221], # at 4.50Ga
	#[0.15053, 0.28976, 0.10767, 0.45204], # 'wtpercent' early earth values, i.e. 4.55Ga - normalized to total U
	[0.372, 0.0164, 0.430, 0.181],  # heat produced by a given isotope
	[0.155, 0.985, 0.0495, 0.555]])  #decay constants in 1/Ga

def produce_heat(planet,t):
	if 'decay' in planet.keys():
		heat = planet['Qp'] * np.exp(-1*planet['decay']*t)
	else:
		Ht=[0.0,0.0,0.0,0.0]
		for i in range(4):
			Ht[i]=radio[0,i] * radio[1,i]*np.exp(radio[2,i] * (-1*t)) #CHANGED TO REFLECT GOING BACK FROM PRESENT DAY
		heat=planet['Qp']*sum(Ht)
	return heat

def frank_kamenetskii(Ev,Tp):
	theta=Ev*(Tp-Ts)/(R*(Tp**2))
	return theta

def flux_heat(planet,Tp,Ra):
	if 'Q0' in planet:
		Fman = (planet['Q0'] * (Tp / planet['scaletemp'])**(1 + planet['beta']) * 
			(get.viscosity(planet,planet['scaletemp']) / get.viscosity(planet,Tp))**(planet['beta']))
	else:  # naive method - no initial flux required.
		theta=frank_kamenetskii(planet['Ev'],Tp)
		Fman=(planet['Sa']*(planet['c1']*planet['k']*(Tp-Ts)/planet['d'])*
			(theta**(-(1+planet['beta'])))*(Ra**(planet['beta'])))

	return Fman

def ThermEv(planet, thermals, method, Tp0, tmax):
	if not('beta' in planet.keys()):
		planet['beta'] = DEFAULT['beta']
    Tp = Tp0
    #dt = 0.01
    t = 0.0              # Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.
    Hts = []             # A list of lists; column names are in get.keys['columns']
    
    if method=='MC':
        planet.update(MC)
    if method=='benchmark':
        planet.update(STO)
    #planet['Qp'] = (7.38e-12 * UreyRatio/0.75862069) * np.exp(1.42e-17 * seconds * tmax) * (planet['Mp'] - planet['Mc']) * planet['Qpl']
    while t <= tmax:
        #if method =='dynamic':
        #    planet.update(get.Tdep_thermals(thermals,Tp))
        #if method=='static':
        #    planet.update(get.Tdep_thermals(thermals,Tp))

        viscT=get.viscosity(planet,Tp)
        Ra=get.rayleigh(planet,Tp,Ts,viscT)

        production=produce_heat(planet,t)
        loss=flux_heat(planet,Tp,Ra)
        dTp=(dt*seconds*(production-loss))/(planet['Cp']*(planet['Mp']-planet['Mc']))
        Hts.append([str({k:v for k,v in planet.items() if k in ['scaletemp', 'beta', 'Q0', 'Qp', 'Ev', 'visc0']}), t,Tp-273.15,Ra,production,loss,production/loss, viscT, np.log10(planet['visc0']),planet['Ev']/1000.,np.log10(viscT), planet['beta']])  #planet['k']])

        Tp=Tp+dTp
        t=t+dt
    
    Evolution=pd.DataFrame(Hts, columns = planet['outcols'])
    Evolution['passfail'] = len(Evolution) * [50>=abs(Evolution.at[int(np.floor((tmax)/0.01)), 'temp']-1350.0)]
    return Evolution


# low_mg = {'C2/c':0.02553006, 'Wus':0.000000, 'Pv':0.40009268, 'an':0.00, \
#                'O':0.00000, 'Wad':0.00000, 'Ring':0.00000, 'Opx':0.06821423, \
#                'Cpx':0.00715906, 'Aki':4.31e-05, 'Gt_maj':0.06264212, 'Ppv':0.05175689, \
#                'CF':0.00, 'st':0.01990724, 'q':0.00125053, 'ca-pv':0.02987515, \
#                'cfs':0.00, 'coe':0.00165384, 'ky':0.00, 'seif':0.00154011}
