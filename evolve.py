import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import sys
import os
import shutil
import getall as get

Me = 5.97e24        #Earth mass in kg
Qe = 3.611610290257257e-11 #Calculated so Qe= ~5.016 * 28.8e12/(Mp-Mc) - i.e. 5x present day, as primordial Earth.
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
Ts=300.0
#Tp=1923.0

radio = np.array([
	#'238U','235U','232Th','40K'; 
	[1.00, 1.00, 1.00, 1.00], #'rel_amt':relative to Earth's values
	[0.14715, 0.29649, 0.10464, 0.45172], #'wtpercent' early earth values - normalized to total U
	[0.155, 0.985, 0.0495, 0.555]]) #decay constants in 1/Ga

def produce_heat(Mp,Mc,Qpl,t):
	Q0=(Qe*Qpl)*(Mp-Mc)       #initial radionuclide abundance total - w/kg versus Earth, times #kg mantle
	Ht=[0.0,0.0,0.0,0.0]
	for i in range(4):
		Ht[i]=radio[0,i]*radio[1,i]*np.exp(radio[2,i]*-t)
	heat=Q0*sum(Ht)
	return heat

def plot_heat(source,title):
	fig = figure(1)
	ax = fig.add_subplot(111, autoscale_on=True)
	production=source
	ax.plot(production[:,0],production[:,1],'red',linewidth=2)
	minx=min(production[:,0])
	maxx=max(production[:,0])
	miny=min(production[:,1])
	maxy=max(production[:,1])+((max(production[:,1])-miny)*0.1)
	plt.xlim(minx,maxx)
	plt.ylim(miny,maxy)
	plt.title(str(title))
	fname=str(str(title).split(' '))+'_temp_evolution.png'
	#title=str(mineral)+'_temp_evolution.pdf'
	#plt.savefig(title)
	plt.show()

def frank_kamenetskii(Ev,Tp):
	theta=Ev*(Tp-Ts)/(R*(Tp**2))
	return theta

def lose_heat(Mpl,CMF,Rpl,CRF,mineral,Tp,Ra):
	Mp,Mc,Rp,Rc=get.SIunits(Mpl,CMF,Rpl,CRF)
	d,Vm,Sa,pm,g=get.build(Mp,Mc,Rp,Rc)
	alpha,cp,k=get.thermals(mineral,Tp)
	c1,Ev,visc0=get.TdepVisc(mineral)
	theta=frank_kamenetskii(Ev,Tp)
	Fman=Sa*(c1*k*(Tp-Ts)/d)*(theta**(-4./3.))*(Ra**(1./3.))
	return Fman

def flux_heat(Sa,c1,k,Tp,d,Ra,Ev):
	theta=frank_kamenetskii(Ev,Tp)
	Fman=Sa*(c1*k*(Tp-Ts)/d)*(theta**(-4./3.))*(Ra**(1./3.))
	return Fman