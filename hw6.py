"""
SES 598 - Atmospheres
Lee Bernard
Homework 6 code
"""

import numpy as np
import matplotlib.pyplot as plt
import yaml


def solar_fun(d_new, baseline_constant=1365):
    # equation for solar constant
    # s_new = (d_0/d_new)**2 * s_0 = 1/d_new**2 * 1365 au^2 W/m^2
    return baseline_constant/(d_new**2)


def equil_temp(d_new, A):
    # equation for equilibrium temperature
    '''
    T_new = 1/sqrt(d) * T_earth works, but does not account for albedo
    Besides, I would need to calculate T_earth.
    Instead, start from earth's solar flux.
    '''

    # T_new = (s_new * (1 - A)/4*sigma)**(1/4)
    # sigma = stefan-boltzmann constant = 5.67e10−8 W⋅m−2⋅K−4
    sigma = 5.67e-8
    return (solar_fun(d_new) * (1 - A)/(4*sigma))**(1/4)


def blackbody(wavelength, temp=300):
    # expects wavelength in meters

    c = 3e8  # um/s
    h = 6.626e-34  # Js
    k = 1.38e-23  # J/K

    # returns W m^-2 m^-1 sr^-1
    return (2*h*(c**2) / (wavelength**5)) * 1/(np.exp(h*c / (wavelength*temp*k)) - 1)


def incident_radiance(wavelength, distance, temp=5800, r_star=.00465):
    # units of solar radius and distance must match

    # solar radiance per sq area is attunated by (r_star/distance)^2
    return (r_star/distance)**2 * blackbody(wavelength, temp=temp)


def reflected_rad(wavelength, distance, albedo, temp=5800, r_star=.00465):
    # default temp and radius are for the Sun
    # return the reflected light, which is incident simply scaled with albedo
    # note this assumes albedo is grey
    return incident_radiance(wavelength, distance, temp=temp, r_star=r_star) * albedo


def find_intersection(set1, set2):
    # size of the datasets must match
    diff = np.abs(set1 - set2)
    return np.argmin(diff)

c = 3e8  # m/s
h = 6.626e-34  # Js
k = 1.38e-23  # J/K
'''
Problem 1
'''

# refer back to the data from HW1
with open('hw1_data/data.yaml') as file:
    planet_data = yaml.load(file, Loader=yaml.FullLoader)


def prob1_plot_data(raw_wavelengths, distance, albedo, name):
    # expects wavelengths in meters, outputs in um
    earth_temp = equil_temp(distance, albedo)

    earth_spectrum = blackbody(raw_wavelengths, temp=earth_temp)
    reflected_spectrum = reflected_rad(raw_wavelengths, distance, albedo)

    # convert from per meter to per micrometer
    earth_spectrum = earth_spectrum * 1e-6
    reflected_spectrum = reflected_spectrum * 1e-6
    wavelengths = raw_wavelengths * 1e6

    idx = find_intersection(earth_spectrum, reflected_spectrum)
    intersection_wavelength = wavelengths[idx]

    plt.figure(name, figsize=(8, 6))
    plt.plot(wavelengths, earth_spectrum)
    plt.plot(wavelengths, reflected_spectrum)
    plt.plot(wavelengths[idx], earth_spectrum[idx], 'ro')
    plt.legend((name + ' emission',
                'Reflected light',
                'Intersection at %.2f μm' %intersection_wavelength))
    plt.title(name)
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('W m^-2 μm^-1 sr^-1')
    # plt.yscale('log')

    return intersection_wavelength


raw_wavelengths = np.linspace(1e-6, 30e-6, num=500)
name_list = ['Earth', 'Venus', 'Mars', 'Jupiter']

intersection_list = []

for curr_name in name_list:
    curr_data = planet_data[curr_name]
    intersection_list.append(prob1_plot_data(raw_wavelengths,
                                             curr_data['distance'],
                                             curr_data['albedo'],
                                             curr_name))

# do planet 9 with special bounds
curr_data = planet_data['Planet9']
intersection_list.append(prob1_plot_data(np.linspace(2e-6, 100e-6, num=500),
                                         curr_data['distance'],
                                         curr_data['albedo'],
                                         'Planet9'))




'''
Problem 5
'''
wave = 10

power = blackbody(wave, temp=300)
power * wave / (h*c)
# 4.97e26
