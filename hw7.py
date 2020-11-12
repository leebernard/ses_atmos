"""
SES 598 - Atmospheres
Lee Bernard
Homework 7 code
"""

import numpy as np
import matplotlib.pyplot as plt

'''
Question 1 for all of these
'''

'''
Part c)
'''
# from part b)

def o_depth(pressure, n, p0=1, k0=(1e-3*1e5)):
    # assumes pressure in bars
    # k0: m^2/kg
    # k0 * 1e5 converts to units of 1/barr
    g = 10  # m/s^2

    return k0/g * 1/(n + 1) * pressure**(n + 1) / p0**n


pressure = np.linspace(1e-7, 1, num=500)
n_space = [0, 1, 2]

plt.figure('pressure dependence', figsize=(8, 6))
for n in n_space:
    plt.plot(pressure, o_depth(pressure, n))

plt.legend(('n=%d' %n_space[0],
            'n=%d' %n_space[1],
            'n=%d' %n_space[2]))
plt.xlabel('Pressure (barrs)')
plt.ylabel('Opacity τ(p)')
plt.yscale('log')
plt.xscale('log')
plt.title('Optical depth with pressure-dependant k')


'''
Part d)
'''

def transmittance(pressure, n, p0=1):
    return np.exp(-1 * o_depth(pressure, n, p0=p0))


plt.figure('transmittance', figsize=(8, 6))
for n in n_space:
    plt.plot(pressure, transmittance(pressure, n))

plt.legend(('n=%d' %n_space[0],
            'n=%d' %n_space[1],
            'n=%d' %n_space[2]))
plt.xlabel('Pressure (barrs)')
plt.ylabel('Transmittance')
plt.title('Transmittance t as a function of pressure')


'''
Part g)
'''

def weighting(pressure, n, H=8e3, p0=1, k0=1e-3, angle=0):
    # expect angle in degrees
    mu = np.cos(np.radians(angle))
    g = 10
    # convert k0 to barr units
    k0 = k0*1e5
    return k0/g * 1/mu * (pressure/p0)**n * pressure/H * np.exp(-1/mu * o_depth(pressure, n, p0=p0, k0=k0))


pressure = np.linspace(10e-7, 1, num=500)
k0_values = [10**-1, 10**-2, 10**-3, 10**-4]

plt.figure('weighting opacities', figsize=(8, 6))

for k0 in k0_values:
    plt.plot(pressure, weighting(pressure, n=0, H=8e3, k0=k0))

plt.legend(('k0=%.1e' %k0_values[0],
            'k0=%.1e' %k0_values[1],
            'k0=%.1e' %k0_values[2],
            'k0=%.1e' %k0_values[3]))
plt.xlabel('Pressure (barrs)')
plt.ylabel('Weighting')
plt.title('Weighting function for various k0 values')


angles = [0, 30, 60, 85]

plt.figure('angles incidence', figsize=(8, 6))
for angle in angles:
    plt.plot(pressure, weighting(pressure, n=0, H=8e3, k0=1e-3, angle=angle))

plt.legend(('angle=%d' %angles[0],
            'angle=%d' %angles[1],
            'angle=%d' %angles[2],
            'angle=%d' %angles[3]))
plt.xlabel('Pressure (barrs)')
plt.ylabel('Weighting')
plt.title('Weighting function at different angles of incidence, k0=10^-3')


pressure = np.linspace(1e-7, 0.1, num=10000)
plt.figure('optical depth weighting', figsize=(8, 6))
for k0 in k0_values:
    plt.plot(o_depth(pressure, n=0, k0=(k0*1e5)), weighting(pressure, n=0, H=8e3, k0=k0))

plt.legend(('k0=%.1e' %k0_values[0],
            'k0=%.1e' %k0_values[1],
            'k0=%.1e' %k0_values[2],
            'k0=%.1e' %k0_values[3]))
plt.xlabel('optical depth')
plt.ylabel('Weighting')
plt.title('Comparison of optical depth to weighting')
plt.xscale('log')



