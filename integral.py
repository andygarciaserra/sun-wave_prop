##PACKAGES##
import sympy as sym
import numpy as np



##DEFINING DERIVATIVES##
kx, kz, w, cs, wc, N = sym.symbols('kx kz w cs wc N')
cs_z, wc_z, N_z = sym.symbols('csz wcz Nz')

F = (cs**2 * (kx**2 + kz**2) + wc**2)/2 + (1/2) * \
    sym.sqrt((cs**2 * (kx**2 + kz**2) + wc**2)**2 - 4*cs**2 * N**2 * kx**2) \
    - w**2

dkx_ds = 0
dkz_ds = -(sym.diff(F, cs)*cs_z + sym.diff(F, wc)*wc_z + sym.diff(F, N)*N_z)
dx_ds = sym.diff(F, kx)
dz_ds = sym.diff(F, kz)

dkx = dkx_ds
dkz = sym.lambdify([N,N_z,cs,kx,kz,wc,cs_z,wc_z],dkz_ds,'numpy')
dx = sym.lambdify([cs,kx,N,kz,wc],dx_ds,'numpy')
dz = sym.lambdify([cs,kz,kx,wc],dz_ds,'numpy')



##INTEGRATING VIA RUNGE-KUTTA##
    #dkx_ds





