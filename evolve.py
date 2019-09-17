import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import getall as get
import plotly.express as px
import pandas as pd

Me = 5.97e24        # Earth mass in kg
Qe = 1.055e-11 		# If bulk silicate earth Urey ratio (BSE) is assumed, crust is included. 16TW currently; ~80TW past; core mass excluded; Qe=2.00e-11
				 	# If not, Qe=1.055e-11. Estimated from convective Urey Ratio for present day, multiplied by 5 to simulate starting Earth values
					# i.e., Qe = ~5.00 * (0.23 * 36.5e12 W)/(M(earth)-M(core)-M(crust)). 
					# circa 4.55Ga. This would imply a closed reservoir of radionuclides in the mantle, with the total mass
					# of continental crust remaining constant.
					# If crustal growth module is added, should use this present day starting point, then add partitioning
					# coefficient for HPE in melt vs. solid.
					# Current HPE estimates from: doi:10.1029/2007RG000241
					#    Mcrust=0.006Mmantle ; Mcore=0.33*Mearth.
R = 8.3145          # Ideal gas constant
Re = 6.371e6        # Earth radius in meters
Ts = 300.0


radio = np.array([
	# '238U','235U','232Th','40K'; 
	[1.00, 1.00, 1.00, 1.00],  # 'rel_amt':relative to Earth's values
	[0.15053, 0.28976, 0.10767, 0.45204], # 'wtpercent' early earth values, i.e. 4.55Ga - normalized to total U
	[0.155, 0.985, 0.0495, 0.555]])  #decay constants in 1/Ga

def produce_heat(planet,t):
	Q0=(Qe*planet['Qpl']) * (planet['Mp'] - planet['Mc'])       # initial radionuclide abundance total - w/kg versus Earth, times #kg mantle
	Ht=[0.0,0.0,0.0,0.0]
	for i in range(4):
		Ht[i]=radio[0,i] * radio[1,i]*np.exp(radio[2,i] * -t)
	heat=Q0*sum(Ht)
	return heat

def plot_heat(source,title):
	plt.rcParams['figure.dpi'] = 150
	fig = figure(1)
	ax = fig.add_subplot(211, autoscale_on=True)
	production=source
	ax.scatter(production[:,0],production[:,1],c='blue',alpha=0.2,s=0.5)
	minx=min(production[:,0])
	maxx=max(production[:,0])
	miny=min(production[:,1])
	maxy=max(production[:,1])+((max(production[:,1])-miny)*0.1)
	plt.xlim(minx,maxx)
	plt.ylim(miny,maxy)
	plt.xlabel('Time, Ga')
	plt.ylabel('Temperature, K')
	plt.grid(which='both', linestyle='--')
	plt.title(str(title))
	fname = str(str(title).split(' ')) + '_temp_evolution.pdf'
	plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
	plt.savefig(fname)
	plt.show()
	return 0

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

def flux_heat(planet,k,Tp,Ra):
	theta=frank_kamenetskii(planet['Ev'],Tp)
	Fman=planet['Sa']*(planet['c1']*k*(Tp-Ts)/planet['d'])*(theta**(-4./3.))*(Ra**(1./3.))
	return Fman

def evolution_colorcoded(nparray, columnkeys, colorcolumn, colortype):
	# Input: nparray = Numpy array; columnkeys = list of strings; colorcolumn = col name; colortype = continuous or discrete
	df = pd.DataFrame(data=nparray, columns=columnkeys)
	if colortype == "continuous":
		plot = px.scatter(df, x="time", y="temp", color=colorcolumn,
			color_continuous_scale=px.colors.diverging.Spectral[::-1])
		plot.show()
	else:
		plot = px.scatter(df, x="time", y="temp", color=colorcolumn,
            colors=px.colors.qualitative.Safe[::-1])
		plot.show()
	return df


def ThermEv(planet, thermals, method, Tp0, tmax):
    Tp=Tp0
    Ts=300.0
    dt=0.01
    Hts=[]             # A list of lists; column names are in get.keys['columns']
    t=0.0              # Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.

    while t <= tmax:

        if method=='dynamic': alpha,cp,k=get.Tdep_thermals(thermals,Tp)
        if method=='static': alpha,cp,k=get.Tdep_thermals(thermals,1600.)
        if method=='benchmark': alpha,cp,k,planet['pm']=3.7e-5,1250.,5.0,3340. #common benchmark values
    
        viscT=get.viscosity(planet,Tp)
        Ra=get.rayleigh(planet,Tp,Ts,viscT,alpha,cp,k)

        production=produce_heat(planet,t)
        loss=flux_heat(planet,k,Tp,Ra)
        dTp=(dt*get.seconds*(production-loss))/(cp*planet['pm']*planet['Vm']) #Potentially change to (cp*Mp)?
        Hts.append([t,Tp,Ra,production,loss,production/loss,alpha,cp])

        Tp=Tp+dTp
        t=t+dt

    Evolution=np.asarray(Hts)
    return Evolution
