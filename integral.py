#----------------------------------------------------------------------------------------#
#
# This script uses a solar interior model script holding the variation of P, T and Rho vs
# Z to calculate the propagation of acoustic-gravitational waves in the solar interior at
# different frequencies (2, 3, 3.5, 5 mHz) from the "Lower Turning Point (LTP)" (k_z = 0).
# 
# A 2.5mHz wave trajectory is also presented as an example.
#
#----------------------------------------------------------------------------------------#





##PACKAGES, DIRECTORIES & CONTANTS##
import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
model = 'model_jcd.dat'
figdir = 'figures/'
EPS = 0.000025
G = 6.674*1e-8          # G in cgs
RSUN = 696340           # Solar radius in km 
MSUN = 1.98e33          # Solar mass in g
GAMMA =5./3             # Adiabatic index
DELTAZ = 23.3           # Model provided has a delta z of 23.3
DELTAX = 23.3           # We set the delta x equal to delta z
FRECUENCIES = np.array([0.002, 0.003, 0.0035, 0.005])   # Frecuencies in Hz
Z_LTP = 500             # Z value where we calculate the LTP





##FITTING SOLAR MODEL##
    # Importing
data = np.loadtxt(model,skiprows=1)
Z = data[:,0]
P = data[:,1]
Rho = data[:,2]
T = data[:,3]
cs = np.sqrt(GAMMA*P/Rho)
r = (RSUN + Z)*1e5
g = G*MSUN/r**2
H = P/Rho*g
wc = cs/(2*H)
N = np.sqrt((g/H)*(GAMMA-1)/GAMMA)
dN_dz = np.gradient(N)
dcs_dz = np.gradient(cs)
dwc_dz = np.gradient(wc)

    # Fitting
p_z = np.polyval(np.polyfit(Z,P,4),Z)
rho_z = np.polyval(np.polyfit(Z,Rho,4),Z)
T_z = np.polyval(np.polyfit(Z,T,4),Z)

    # Plotting
for i in [[p_z,P],[rho_z,Rho],[T_z,T]]:
    plt.figure()
    plt.plot(Z,i[0],label='fit')
    plt.plot(Z,i[1],label='model')
    plt.legend()
    plt.show()





##DEFINING DERIVATIVES##
kx, kz, w, cs, wc, N = sym.symbols('kx kz w cs wc N')
cs_z, wc_z, N_z = sym.symbols('csz wcz Nz')
                                       #! This minus sign is what changes between fast P-mode (+) and slow G-mode (-)
F = (cs**2 * (kx**2 + kz**2) + wc**2)/2 - (1/2) * \
    sym.sqrt((cs**2 * (kx**2 + kz**2) + wc**2)**2 - 4*cs**2 * N**2 * kx**2) \
    - w**2

dkx_ds = 0
dkz_ds = -(sym.diff(F, cs)*cs_z + sym.diff(F, wc)*wc_z + sym.diff(F, N)*N_z)
dx_ds = sym.diff(F, kx)
dz_ds = sym.diff(F, kz)

print(dkx_ds)
print(dkz_ds)
print(dx_ds)
print(dz_ds)

dkx = dkx_ds
dkz = sym.lambdify([N,N_z,cs,kx,kz,wc,cs_z,wc_z],dkz_ds,'numpy')
dx = sym.lambdify([cs,kx,N,kz,wc],dx_ds,'numpy')
dz = sym.lambdify([cs,kz,kx,wc],dz_ds,'numpy')

