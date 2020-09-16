"""
SES 598 - Atmopheres
Homework 1 code
"""
import numpy as np
import matplotlib.pyplot as plt

import yaml

from tabulate import tabulate

"""
Question 1
"""
# equation for solar constant
# s_new = (d_0/d_new)**2 * s_0 = 1/d_new**2 * 1365 au^2 W/m^2
solar_fun = lambda d_new: 1365/(d_new**2)

# equation for equilibrium temperature
'''
T_new = 1/sqrt(d) * T_earth works, but does not account for albedo
Besides, I would need to calculate T_earth. 
Instead, start from earth's solar flux.
'''
# T_new = (s_new * (1 - A)/4*sigma)**(1/4)
# sigma = stefan-boltzmann constant = 5.67e10−8 W⋅m−2⋅K−4
sigma = 5.67e-8
equil_temp = lambda d_new, A: (solar_fun(d_new) * (1 - A)/(4*sigma))**(1/4)

with open('hw1_data/data.yaml') as file:
    planet_data = yaml.load(file, Loader=yaml.FullLoader)

solar_constants = [solar_fun(planet_data[planet]['distance']) for planet in planet_data.keys()]

equiv_temps = [equil_temp(planet_data[planet]['distance'], planet_data[planet]['albedo'])
               for planet in planet_data.keys()]


table_entries = [list(planet_data.keys()), solar_constants, equiv_temps]
transposed = list(map(list, zip(*table_entries)))

colLabels = ("Planet", "Solar Constant", "Equivalent Temp")
print(tabulate(transposed, headers=colLabels))

"""
Question 2:
Equation use for plotting equilibrium temperature
(0.7 / 4)^(1 / 4) * (RSTAR[m] / SEP[m])^(1 / 2) * TEFF[k]
"""

"""
Question 3
"""

s_1 = 1365 * (3600/5400)**4 * (.5)**2
# output: 67.40
s_01 = s_1 + 1365
# output: 1432
T_eq = 245*(1432/1365)**(1/4)
# output: 247

# fraction of solar constant due to m dwarf
frac = 67/s_01
# output : 0.047

"""
Question 4
"""


def q_day(d, phi, delta, s_0=1365):
    output = []
    # in degrees, need radians
    phi = np.radians(phi)
    delta = np.radians(delta)

    for subphi in phi:
        data_sanitizer = -np.tan(subphi) * np.tan(delta)
        # if -tan(phi)tan(delta) < -1, sun never sets
        data_sanitizer[data_sanitizer < -1] = -1
        # if -tan(phi)tan(delta) > 1, sun never rises
        data_sanitizer[data_sanitizer > 1] = 1

        # now that the data is sanitized, calculate sunrise/sunset angle
        # cos(-1) = pi, sun never sets
        # cos(1) = 0, sun never rises
        h_0 = np.arccos(data_sanitizer)

        output.append(s_0/np.pi * (1/d)**2 * (h_0 * np.sin(subphi) * np.sin(delta) +
                                              np.cos(subphi) * np.cos(delta) * np.sin(h_0)))

    return np.array(output)

"""
with open('ephemeris.txt', 'r') as file:
    raw_data = file.readlines()
    earth_data = []
    for x in raw_data:
        earth_data.append(x.split())
    # transpose it from row (lines) to columns
    earth_data = list(map(list, zip(*earth_data)))
    # print(earth_data)

declination = np.array([float(number) for number in earth_data[3]])
distance = np.array([float(number) for number in earth_data[4]])

polar_angles = np.linspace(-90, 90, num=360)
insolation = q_day(d=distance, phi=polar_angles, delta=declination)

plot_name = 'Earth'
plt.figure('earth_inso', figsize=(8, 6))
conplot = plt.contour(np.arange(declination.size), polar_angles, insolation)
plt.clabel(conplot)
plt.ylabel('Latitude')
plt.xlabel('Days')
plt.title(r'Daily Averaged Insolation for ' + plot_name + r' $(W/m^2)$')
plt.savefig(plot_name + '_inso')
"""


"""
Question 5
"""


def generate_isolation_contour(data_file='ephemeris.txt', output_name='Earth'):
    with open(data_file, 'r') as file:
        raw_data = file.readlines()
        earth_data = []
        for x in raw_data:
            earth_data.append(x.split())
        # transpose it from row (lines) to columns
        earth_data = list(map(list, zip(*earth_data)))
        # print(earth_data)

    declination = np.array([float(number) for number in earth_data[3]])
    distance = np.array([float(number) for number in earth_data[4]])

    polar_angles = np.linspace(-90, 90, num=360)
    insolation = q_day(d=distance, phi=polar_angles, delta=declination)

    plt.figure(output_name, figsize=(8, 6))
    conplot = plt.contour(np.arange(declination.size), polar_angles, insolation)
    plt.clabel(conplot)
    plt.ylabel('Latitude')
    plt.xlabel('Days')
    plt.title(r'Daily Averaged Insolation for ' + output_name + r' $(W/m^2)$')
    plt.savefig(output_name + '_inso')


# do mars!
generate_isolation_contour(data_file='hw1_data/mars_ephemeris.txt', output_name='Mars')

# now Jupiter!
generate_isolation_contour(data_file='hw1_data/jupiter_emphemeris.txt', output_name='Jupiter')

# now with Titan. This is gonna get weird
generate_isolation_contour(data_file='hw1_data/titan_emphemeris.txt', output_name='Titan')

# and last and definitely the strangest: Uranus
generate_isolation_contour(data_file='hw1_data/uranus_ephemeris.txt', output_name='Uranus')


