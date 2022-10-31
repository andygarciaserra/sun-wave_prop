
##PACKAGES & GENERAL VARIABLES##    
import numpy as np
import matplotlib.pyplot as plt
filename = 'model_jcd.dat'
figdir = 'figures/'



##PLOT-FORMAT##
plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=20)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels


##IMPORT##

data = np.loadtxt(filename,skiprows=1)
Z = data[:,0]
P = data[:,1]
Rho = data[:,2]
T = data[:,3]


##PLOTTING##
    #Pressure
plt.figure(figsize=(12,8))
plt.plot(Z,P,lw=2,color='green',label='P(Z)')
plt.xlabel('Z [Km]')
plt.ylabel('P [dyn/cm^2]')
plt.legend()
plt.tight_layout()
plt.savefig(figdir+'model_P(Z).png',dpi=100)
plt.show()

    #Density
plt.figure(figsize=(12,8))
plt.plot(Z,Rho,lw=2,color='orange',label='Rho(Z)')
plt.xlabel('Z [Km]')
plt.ylabel('Rho [g/cm^2]')
plt.legend()
plt.tight_layout()
plt.savefig(figdir+'model_Rho(Z).png',dpi=100)
plt.show()

    #Temperature
plt.figure(figsize=(12,8))
plt.plot(Z,T,lw=2,color='red',label='T(Z)')
plt.xlabel('Z [Km]')
plt.ylabel('T [K]')
plt.legend()
plt.tight_layout()
plt.savefig(figdir+'model_T(Z).png',dpi=100)
plt.show()

























