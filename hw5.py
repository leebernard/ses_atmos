"""
SES 598 atmos HW 5
"""

import numpy as np
import matplotlib.pyplot as plt

'''
Problem 1
'''
.454 * 4180 * 25 / (3.336e5 + 2108*10)
# 0.133 kg

'''
Problem 3
'''
es = 611.2 * np.exp((18.729*25 - 25**2/227.3) / (25 + 257.87))

# 3168 Pa

8e3 * np.log(es/101.3e3)
# -27719

'''
Problem 4
'''
# part a solved in wolfram alpha
# part b)
R = 5.67e-8 * (148)**4 / (1000*5.9e5)
# 4.61e-8
R * 687 * 86400

# part c)

# mass of ice
2 * np.pi * (3390e3)**2 * (1 - np.sin(np.radians(70))) * 2.7 * 1000
# 1.176e+16

# scale height
k = 1.38e-23  # boltzmann constant k_b in J/K
amu_kg = 1.66e-27  # kg/amu

H = k*148.4/(44*amu_kg * 3.76)

'''
Problem 5
'''
L = 2264.76  # kJ/kg
ul = 18  # amu
ua = 29  # amu
amu_kg = 1.66e-27  # kg/amu
cp = 1000  # J/kg/K

T0 = 296  # K
P0 = 1  # bar

k = 1.38e-23  # boltzmann constant k_b in J/K


def dry_adiabat(P):
    kappa = k/(ua*amu_kg * cp)

    # return np.log(T0) + kappa*np.log(P/P0)

    # new version
    lnP = np.log(P)
    # T_out = [T0]
    lnT_out = [np.log(T0)]
    for i, curr_p in enumerate(lnP):
        if i + 1 < len(lnP):
            # return the current dlnT/dlnP
            dlnt_dlnp = kappa
            # calculate dlnP as simple delta
            dlnp = lnP[i + 1] - lnP[i]
            # calculate the next temperature
            next_lnT = dlnt_dlnp * dlnp + lnT_out[i]
            # and save it
            lnT_out.append(next_lnT)

    return np.array(lnT_out)


# wet adiabat

def e_s(temp):
    # dew point
    # this is only valid for C, convert K to C
    temp = temp - 273.15
    # Buck equation
    # returns in barrs
    return 6.1121 * np.exp((18.729 * temp - temp ** 2 / 227.3) / (temp + 257.87))
    # return 6.1121e-3 * np.exp((18.678 * temp - temp ** 2 / 234.5) / (temp + 257.14))


def w_s(P, T):
    # expects pressure in bars and temperature in kelvin
    return -(ul/ua * e_s(T)) / P


def wet_lapse(P, T):
    kappa = k / (ua * amu_kg * cp)
    specific_lapse = \
        (1 + (ua*amu_kg)*L*w_s(P, T) / (k*T)) / (1 + (ul*amu_kg)*(L**2)*w_s(P, T) / (cp*k*(T**2))) * kappa
    return specific_lapse


def wet_adiabat(P):
    lnP = np.log(P)
    # T_out = [T0]
    es_out = [e_s(T0)]
    lnT_out = [np.log(T0)]
    for i, curr_p in enumerate(P):
        if i+1 < len(P):
            # return the current dlnT/dlnP
            dlnt_dlnp = wet_lapse(curr_p, np.exp(lnT_out[i]))
            # calculate dlnP as simple delta
            dlnp = lnP[i+1] - lnP[i]
            # calculate the next temperature
            next_lnT = dlnt_dlnp * dlnp + lnT_out[i]
            # and save it
            lnT_out.append(next_lnT)

            # save the current es for debugging
            es_out.append(e_s(np.exp(next_lnT)))

    return np.array(lnT_out), np.array(es_out)


# T = np.linspace(75, 325, num=100)
# lnT = np.log(T)
P = np.linspace(1, .001, num=1000)
lnP = np.log(P)


lnT = dry_adiabat(P)
T = np.exp(lnT)

wet_lnT, es = wet_adiabat(P)

plt.figure('Dew_point')
plt.plot(es, P)
plt.yscale('log')
plt.xscale('log')
plt.xlim(left=10e-24)
plt.gca().invert_yaxis()

plt.figure('Adiabatic_lapse', figsize=(8, 6))
plt.plot(lnT, lnP)
plt.plot(wet_lnT, lnP)
plt.xlabel('lnT (K)')
# plt.xlim(left=0)
plt.ylabel('lnP (barr)')
plt.legend(('Dry', 'Wet'))
plt.xlim(4.2, 5.8)
plt.ylim(-6, 0)
plt.gca().invert_yaxis()

plt.savefig('hw5_q5')


