"""
Stuff about equilibrium chemistry and law of mass action
"""

import numpy as np
import matplotlib.pyplot as plt

'''
Question 5
'''
def K_eq(T):
    R = 8.314 * 0.2389
    return np.exp(-G(T)/(R*T))


def G(T):
    aa = 1.106131e6
    bb = -5.6895e4
    cc = 62.565
    dd = -5.81396E-4
    ee = 2.3465E-8
    return aa/T + bb + cc*T + dd*(T**2) + ee*(T**3)


def f_co(T, P=1):
    f_h2 = .85
    f_c = 5.91e-4
    f_o = 1.04e-3

    return 1/2*(P**2 * f_h2**3 * K_eq(T) + f_c + f_o) - \
           1/2*np.sqrt((P**2 * f_h2**3 * K_eq(T) + f_c + f_o)**2 - 4*f_o*f_c)


T = 1000
f_c = 5.91e-4
f_o = 1.04e-3

f_co(T)
f_h20 = f_o - f_co(T)
f_ch4 = f_c - f_co(T)


'''
Question 6
'''
temps = np.linspace(500, 2100, num=500)

p10_f_co = f_co(temps, P=10)
p10_f_h20 = f_o - p10_f_co
p10_f_ch4 = f_c - p10_f_co

p1_f_co = f_co(temps, P=1)
p1_f_h20 = f_o - p1_f_co
p1_f_ch4 = f_c - p1_f_co

p001_f_co = f_co(temps, P=0.01)
p001_f_h20 = f_o - p001_f_co
p001_f_ch4 = f_c - p001_f_co

plt.figure('hw8q6', figsize=(8,6))

plt.plot(temps, p1_f_co, c='tab:blue')
plt.plot(temps, p10_f_co, c='tab:blue', ls='--')
plt.plot(temps, p001_f_co, c='tab:blue', ls=':')
plt.plot(temps, p1_f_h20, c='tab:orange')
plt.plot(temps, p1_f_ch4, c='tab:olive')

plt.plot(temps, p10_f_h20, c='tab:orange', ls='--')
plt.plot(temps, p10_f_ch4, c='tab:olive', ls='--')

plt.plot(temps, p001_f_h20, c='tab:orange', ls=':')
plt.plot(temps, p001_f_ch4, c='tab:olive', ls=':')

plt.xlabel('Temperature (K)')
plt.ylabel('Mixing ratio')
plt.legend(('CO, 1 barr', '10 barr', '0.01 bar', 'H2O', 'CH4'))


'''
Question 7
'''
f_h2 = .85
P = 0.01  # barr
T = 1000


def C_O_tracer(ratio, P=1):
    # ratio = f_c/f_o
    # returns f_co

    # calculate new factions of Oxygen and Carbon
    f_total = 5.91e-4 + 1.04e-3
    f_o = f_total/(1 + ratio)
    f_c = f_total/(1 + 1/ratio)

    f_co = 1/2*(P**2 * f_h2**3 * K_eq(T) + f_c + f_o) - \
           1/2*np.sqrt((P**2 * f_h2**3 * K_eq(T) + f_c + f_o)**2 - 4*f_o*f_c)
    return f_co, f_o, f_c


ratios = np.geomspace(.01, 100, num=1000)

f_co, f_o, f_c = C_O_tracer(ratios, P=.01)
f_h20 = f_o - f_co
f_ch4 = f_c - f_co


plt.figure('hw8q7', figsize=(8, 6))
plt.plot(ratios, f_co)
plt.plot(ratios, f_h20)
plt.plot(ratios, f_ch4)

f_co57, _, _ = C_O_tracer(.57, P=.01)
plt.plot(.57, f_co57, 'ro')


plt.title('Temp = 1000 K, Pressure = .01 bar')
plt.xlabel('C/O ratio')
plt.ylabel('Mixing ratio')
plt.legend(('CO', 'H2O', 'CH4', 'Ratio=.57'))
plt.xscale('log')


