import numpy as np
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
	fname=str(str(title).split(' '))+'_temp_evolution.png'
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

def flux_heat(Sa,c1,k,Tp,d,Ra,Ev):
	theta=frank_kamenetskii(Ev,Tp)
	Fman=Sa*(c1*k*(Tp-Ts)/d)*(theta**(-4./3.))*(Ra**(1./3.))
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