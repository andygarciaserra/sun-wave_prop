# Task 1 - Acoustic-gravitational wave propagation in the solar interior
# ULL - 19/10/2022
# Sergio Guerra Arencibia

from sympy import *
import matplotlib.pyplot as plt
import numpy as np

# Definition of constants
EPS = 0.000025
RADIUS_SUN = 696340  # Solar radius in km 
GAMMA = 5/3          # Adiabatic index
DELTAZ = 23.3        # Model provided has a delta z of 23.3
DELTAX = 23.3        # We set the delta x equal to delta z
FRECUENCIES = np.array([0.002, 0.003, 0.0035, 0.005])   # Frecuencies in Hz
ZLTP = 500           # Z value where we calculate the LTP

# Function for plotting 2 magnitudes
def plotMagnitudes(xValues, yValues, points = []):
    fig, ax = plt.subplots()
    fig.set_figheight(10)
    fig.set_figwidth(15)
    ax.plot(xValues, yValues)
    if (len(points) != 0):
        for i in range(len(points)):
            ax.scatter(points[i], FRECUENCIES[i], s = 20)
    plt.show()

# Function to change from Z values provided in the solar model to Radius values
def modelZtoRadius(zValues):
    return (RADIUS_SUN + zValues) / RADIUS_SUN

# RK4 generic method for a system of 4 equations
def rk4(function, initialValues, t0, tf, steps):
    t = np.linspace(t0, tf, steps)
    deltaT = t[1] - t[0]
    values = np.zeros(shape=(4, steps))
    values[0, 0] = initialValues[0]
    values[1, 0] = initialValues[1]
    values[2, 0] = initialValues[2]
    values[3, 0] = initialValues[3]

    for i in range(1, steps):
        k1 = function(*values[:, i-1])
        k2 = function(*values[:, i-1] + 0.5*k1*deltaT)
        k3 = function(*values[:, i-1] + 0.5*k2*deltaT)
        k4 = function(*values[:, i-1] + k3*deltaT)
        values[:, i] = values[:, i-1] + (1/6)*deltaT*(k1 + 2*k2 + 2*k3 + k4)
    return (values)

# Function to integrate
def fFunction(x, z, kx, kz):
    index = np.argmin(np.abs(z - zArray))

    dkx_ds = 0

    # Using negative sign
    # dx_ds = (kx * (cs[index]**2)) - (0.5*(2*kx*(((kx**2) + (kz**2))*(cs[index]**2) +( wc[index]**2)) * (cs[index]**2) - (4*kx*(cs[index]**2)*(n[index]**2)))) \
    #     /(np.sqrt(-4*(kx**2)*(cs[index]**2)*(n[index]**2) + ((cs[index]**2)*(kx**2 + kz**2) +( wc[index]**2))**2 ))


    # dz_ds = (kz * (cs[index]**2)) - ( kz * (cs[index]**2) * ((kx**2 + kz**2)*(cs[index]**2) + (wc[index]**2)) ) \
    #     / (np.sqrt( -4*(kx**2)*(cs[index]**2)*(n[index]**2) + ((kx**2 + kz**2) * (cs[index]**2) + (wc[index]**2) )**2) )


    # dkz_ds = (kx**2 + kz**2) * cs[index] * gradcs[index] + (wc[index] * gradwc[index]) - \
    #     (0.5 *( ( -4*(kx**2)*(cs[index]**2)*n[index]*gradn[index] ) - (4*(kx**2)*cs[index]*(n[index]**2)*gradcs[index] ) + \
    #     ( ( ((kx**2 + kz**2)*(cs[index]**2) + (wc[index]**2)) * (4*(kx**2 + kz**2)*cs[index]*gradcs[index] + 4*wc[index]*gradwc[index]) ) / 2)        ) ) / \
    #     ( np.sqrt( -4*(kx**2)*(cs[index]**2)*(n[index]**2) + ((kx**2 + kz**2)*(cs[index]**2) + (wc[index]**2))**2 ) )

    # Using positive sign
    dx_ds = (kx * (cs[index]**2)) + (0.5*(2*kx*(((kx**2) + (kz**2))*(cs[index]**2) +( wc[index]**2)) * (cs[index]**2) - (4*kx*(cs[index]**2)*(n[index]**2)))) \
        /(np.sqrt(-4*(kx**2)*(cs[index]**2)*(n[index]**2) + ((cs[index]**2)*(kx**2 + kz**2) +( wc[index]**2))**2 ))


    dz_ds = (kz * (cs[index]**2)) + ( kz * (cs[index]**2) * ((kx**2 + kz**2)*(cs[index]**2) + (wc[index]**2)) ) \
        / (np.sqrt( -4*(kx**2)*(cs[index]**2)*(n[index]**2) + ((kx**2 + kz**2) * (cs[index]**2) + (wc[index]**2) )**2) )


    dkz_ds = (kx**2 + kz**2) * cs[index] * gradcs[index] + (wc[index] * gradwc[index]) + \
        (0.5 *( ( -4*(kx**2)*(cs[index]**2)*n[index]*gradn[index] ) - (4*(kx**2)*cs[index]*(n[index]**2)*gradcs[index] ) + \
        ( ( ((kx**2 + kz**2)*(cs[index]**2) + (wc[index]**2)) * (4*(kx**2 + kz**2)*cs[index]*gradcs[index] + 4*wc[index]*gradwc[index]) ) / 2)        ) ) / \
        ( np.sqrt( -4*(kx**2)*(cs[index]**2)*(n[index]**2) + ((kx**2 + kz**2)*(cs[index]**2) + (wc[index]**2))**2 ) )

    dz_ds = -dz_ds
    return np.array([dx_ds, dz_ds, dkx_ds, dkz_ds])

# results = rk4(fFunction, [0, zArray[ZLTP], k_xArray[0], 0], 0, 1000, 3)

# -------------------------------------------------------------------------
# Reading the model data
zArray = np.array([])
pArray = np.array([])
rhoArray = np.array([])
tArray = np.array([])

with open("model_jcd.dat") as openFileObject:
    next(openFileObject)
    for line in openFileObject:
        currentValues = line.split()
        zArray = np.append(zArray, float(currentValues[0]))
        pArray = np.append(pArray, float(currentValues[1]))
        rhoArray = np.append(rhoArray, float(currentValues[2]))
        tArray = np.append(tArray, float(currentValues[3]))

# # Plotting preassure, density and temperature vs Z
# fig, ax = plt.subplots(3, 1)
# fig.set_figheight(10)
# fig.set_figwidth(15)
# ax[0].plot(zArray, tArray)
# ax[1].plot(zArray, pArray)
# ax[2].plot(zArray, rhoArray)
# plt.show()

# -------------------------------------------------------------------------
# Sound velocity
cs = np.array([])

def soundVelocityFormula(preassure, density):
    cs = (GAMMA * preassure) / density
    return(np.sqrt(cs))
cs = soundVelocityFormula(pArray, rhoArray)
# plotMagnitudes(modelZtoRadius(zArray), cs)

# Height scale
sunGravity = 274 # m/s**2
sunGravity = sunGravity * 100 # to c.g.s

def heightScaleFormula(preassure, density):
    h = preassure / (density * sunGravity)
    return h
scaleHeight = heightScaleFormula(pArray, rhoArray)
# plotMagnitudes(modelZtoRadius(zArray), scaleHeight)

# w_c
def wcFormula(cs, h):
    return (cs / (2 * h))
wc = wcFormula(cs, scaleHeight)
# plotMagnitudes(modelZtoRadius(zArray), wc)

# N
def bruntVaisala (h):
    firstMember = sunGravity / h
    return np.sqrt(firstMember * ((GAMMA - 1) / GAMMA))
n = bruntVaisala(scaleHeight)
# plotMagnitudes(modelZtoRadius(zArray), n)

# ---------------------------------------------------------

# Get the external Z at which the wave reflects
externalZ = np.array([])
frecuencyw_c = wc / (2*np.pi)
for i in range(len(FRECUENCIES)):
    for j in range(len(wc)):
        if (np.abs(FRECUENCIES[i] - frecuencyw_c[j]) < EPS):
            externalZ = np.append(externalZ, zArray[j])
            break
# plotMagnitudes(modelZtoRadius(zArray), frecuencyw_c, modelZtoRadius(externalZ))

# Get the k_x of the wave (it's constant since in this axis the medium does not vary)
# We get it from the LTP expression, where we know that k_z = 0
k_xArray = np.array([])  # Units will be cm^-1
for i in range(len(FRECUENCIES)):
    k_xArray = np.append(k_xArray, (2*np.pi*FRECUENCIES[i]) / cs[ZLTP])


# -------------------------------------------------------

# Derivation
# x, z, kx, kz, w = symbols('x z kx kz w')
# csSym = Function('cs')(z)
# wcSym = Function('wc')(z)
# nSym = Function('n')(z)

# firstMember = (csSym**2 * (kx**2 + kz**2) + wcSym**2) / 2
# insideRoot = (csSym**2 * (kx**2 + kz**2) + wcSym**2)**2 - (4 * csSym**2 * nSym**2 * kx**2)
# secondMember = 0.5*sqrt(insideRoot) - w**2
# f = firstMember + secondMember

# This way we can see the derivatives and create the function to integrate
# print(f.diff(x))
# print(f.diff(z))
# print(f.diff(kx))
# print(f.diff(kz))

# Numeric integration

gradcs = np.gradient(cs, -DELTAZ)
gradwc = np.gradient(wc, -DELTAZ)
gradn = np.gradient(n, -DELTAZ)

lines = []
results1 = rk4(fFunction, [0, zArray[ZLTP], k_xArray[0], 0], 0, 2, 1500)
results2 = rk4(fFunction, [0, zArray[ZLTP], k_xArray[1], 0], 0, 2, 1500)
results3 = rk4(fFunction, [0, zArray[ZLTP], k_xArray[2], 0], 0, 2, 1500)
results4 = rk4(fFunction, [0, zArray[ZLTP], k_xArray[3], 0], 0, 2, 1500)



fig, ax = plt.subplots()
fig.set_figheight(10)
fig.set_figwidth(15)
ax.set_xlabel('Distancia [km]')
ax.set_ylabel('Altura [km]', rotation = 90)
ax.yaxis.label.set_size(15)
ax.xaxis.label.set_size(15)
plt.rcParams.update({'font.size': 17})
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)


lines += ax.plot(results1[0], results1[1], color='blueviolet')
lines += ax.plot(results2[0], results2[1], color='cadetblue')
lines += ax.plot(results3[0], results3[1], color='gold')
lines += ax.plot(results4[0], results4[1], color='peru')

ax.legend(lines, ["2 mHz", "3 mHz", "3.5 mHz", "5 mHz"], fontsize = "medium", loc = "upper right")
# ax.legend(lines, ["LTP 11km", "LTP 8.8km", "LTP 6.5km", "LTP 4.1km"], fontsize = "medium", loc = "upper right")

plt.show()