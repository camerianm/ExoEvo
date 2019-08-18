import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import sys
import os
import shutil
import getall as get

Me = 5.97e24        #Earth mass in kg
Qe = 1.055e-11 		# If bulk silicate earth Urey ratio (BSE) is assumed, crust is included. 16TW currently; ~80TW past; core mass excluded; Qe=2.00e-11
				 	# If not, Qe=1.055e-11. Estimated from convective Urey Ratio for present day, multiplied by 5 to simulate starting Earth values
					# i.e., Qe = ~5.00 * (0.23 * 36.5e12 W)/(M(earth)-M(core)-M(crust)). 
					# circa 4.55Ga. This would imply a closed reservoir of radionuclides in the mantle, with the total mass
					# of continental crust remaining constant.
					# If crustal growth module is added, should use this present day starting point, then add partitioning
					# coefficient for HPE in melt vs. solid.
					# Current HPE estimates from: doi:10.1029/2007RG000241
					#    Mcrust=0.006Mmantle ; Mcore=0.33*Mearth.
R = 8.3145          #Ideal gas constant
Re = 6.371e6        #Earth radius in meters
Ts = 300.0

radio = np.array([
	#'238U','235U','232Th','40K'; 
	[1.00, 1.00, 1.00, 1.00], #'rel_amt':relative to Earth's values
	[0.15053, 0.28976, 0.10767, 0.45204], #'wtpercent' early earth values, i.e. 4.55Ga - normalized to total U
	[0.155, 0.985, 0.0495, 0.555]]) #decay constants in 1/Ga

def produce_heat(planet,t): #Mp,Mc,Qpl,t):
	Q0=(Qe*planet['Qpl'])*(planet['Mp']-planet['Mc'])       #initial radionuclide abundance total - w/kg versus Earth, times #kg mantle
	Ht=[0.0,0.0,0.0,0.0]
	for i in range(4):
		Ht[i]=radio[0,i]*radio[1,i]*np.exp(radio[2,i]*-t)
	heat=Q0*sum(Ht)
	return heat

def dorn_heat(Mp,Mc,Qpl,t):
    Q0=(2.42e-11*Qpl)*(Mp-Mc)
    heat=Q0*np.exp((-1*t)/2.85)
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
	plt.grid(which='both',linestyle='--')
	plt.title(str(title))
	fname=str(str(title).split(' '))+'_temp_evolution.pdf'
	plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
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

def flux_heat(planet,k,Tp,Ra):# Sa,c1,k,Tp,d,Ra,Ev):
	theta=frank_kamenetskii(planet['Ev'],Tp)
	Fman=planet['Sa']*(planet['c1']*k*(Tp-Ts)/planet['d'])*(theta**(-4./3.))*(Ra**(1./3.))
	return Fman

def show_structure(Rp,CRF):   #beta
    from mpl_toolkits.mplot3d import Axes3D
    plt.rcParams['figure.dpi'] = 100
    fig = plt.figure()
    ax1 = fig.add_subplot(111)#, projection='3d')

    u, v = np.mgrid[0:1.5*np.pi:100j, 0:1*np.pi:100j]
    x = Rp * np.cos(u)*np.sin(v)
    y = Rp * np.sin(u)*np.sin(v)
    z = Rp * np.cos(v)
    ax1.plot_area(x, y, color='red', alpha=0.5) #, z,
    ax1.set_aspect('equal', 'box')
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
    x = Rp*CRF * np.cos(u)*np.sin(v)
    y = Rp*CRF * np.sin(u)*np.sin(v)
    z = Rp*CRF * np.cos(v)
    ax1.plot_area(x, y, color='black', alpha=1.0) # z, 
    plt.axis('off')
    plt.show()
    
    
# Hts=[]             #A list of lists; column names are in get.keys['columns']
# t=0.0              #Keep Hts=[], Tp=Tp0, and t=0.0 here, so we can reset values and run again.
# Tp=Tp0
# dt=0.01

# start=time.time()
# while t <= tmax:
#     alpha,cp,k=get.thermals(mineral,Tp)
    
#     if method=='static':
#         alpha,cp,k=get.thermals(mineral,1625)
#         #alpha,cp,k,pm=3.7e-5,1250.,5.0,3340. #Uncomment for common benchmark values
        
#     viscT=get.viscosity(Ev,visc0,Tp)
#     Ra=get.rayleigh(d,g,pm,Tp,Ts,viscT,alpha,cp,k)
    
#     production=evolve.produce_heat(Mp,Mc,Qp,t)
#     loss=evolve.flux_heat(Sa,c1,k,Tp,d,Ra,Ev)
#     dTp=(dt*get.seconds*(production-loss))/(cp*pm*Vm) #Potentially change to (cp*Mp)?
#     Hts.append([t,Tp-273.15,Ra,production,loss,production/loss])
    
#     Tp=Tp+dTp
#     t=t+dt
# end=time.time()
# Evolution=np.asarray(Hts)

# print("Program running time: ", Pf(end-start), " seconds")
# print("Final values:")
# print(get.keys['columns'])
# print([Pf(t), Pf(Tp), Pe(Ra), Pe(production), Pe(loss), Pf(production/loss)])
