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


def equil_temp(d_new, A, baseline_constant=1365):
    # equation for equilibrium temperature
    '''
    T_new = 1/sqrt(d) * T_earth works, but does not account for albedo
    Besides, I would need to calculate T_earth.
    Instead, start from earth's solar flux.
    '''

    # T_new = (s_new * (1 - A)/4*sigma)**(1/4)
    # sigma = stefan-boltzmann constant = 5.67e10−8 W⋅m−2⋅K−4
    sigma = 5.67e-8
    return (solar_fun(d_new, baseline_constant=baseline_constant) * (1 - A)/(4*sigma))**(1/4)


def blackbody(wavelength, temp=300):
    # expects wavelength in meters

    c = 3e8  # um/s
    h = 6.626e-34  # J
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


def prob1_plot_data(raw_wavelengths, distance, albedo, name,
                    r_star=.00465, star_temp=5800, baseline_constant=1365):
    # expects wavelengths in meters, outputs in um
    planet_temp = equil_temp(distance, albedo, baseline_constant=baseline_constant)
    print('Planet temp:', planet_temp)

    planet_spectrum = blackbody(raw_wavelengths, temp=planet_temp)
    reflected_spectrum = reflected_rad(raw_wavelengths, distance, albedo,
                                       r_star=r_star, temp=star_temp)

    # convert from per meter to per micrometer
    planet_spectrum = planet_spectrum * 1e-6
    reflected_spectrum = reflected_spectrum * 1e-6
    wavelengths = raw_wavelengths * 1e6

    idx = find_intersection(planet_spectrum, reflected_spectrum)
    intersection_wavelength = wavelengths[idx]

    plt.figure(name, figsize=(8, 6))
    plt.plot(wavelengths, planet_spectrum)
    plt.plot(wavelengths, reflected_spectrum)
    plt.plot(wavelengths[idx], planet_spectrum[idx], 'ro')
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

# for exoplanet, rescale the solar constant
new_insolation = 1365 * (1.15)**2 * (6065/5800)**4
intersection_list.append(prob1_plot_data(raw_wavelengths=np.linspace(.5e-6, 5e-6, num=500),
                         distance=.047,
                         albedo=.01,
                         name='Hot Jup',
                         r_star=.00465*1.15, star_temp=6065, baseline_constant=new_insolation))



'''
Problem 4
'''

def rayleigh_jeans(wavelength, temp=300):
    # expects wavelength in meters

    c = 3e8  # um/s
    h = 6.626e-34  # Js
    k = 1.38e-23  # J/K

    # returns W m^-2 m^-1 sr^-1
    return 2*c*k/(wavelength**4) * temp


raw_wavelengths = np.linspace(1e-6, 2500e-6, num=1000)
wavelengths = raw_wavelengths * 1e6

bb_emission = blackbody(raw_wavelengths, temp=300) * 1e-6  # in um
approx_emission = rayleigh_jeans(raw_wavelengths, temp=300) * 1e-6

percent_difference = (approx_emission - bb_emission)/bb_emission * 100
percent_difference = np.around(percent_difference, decimals=3)
idx = np.nonzero(percent_difference == 1.0)[0][0]
one_percent_wavelength = wavelengths[idx]

plt.figure('RJ approx', figsize=(8, 6))
plt.plot(wavelengths, bb_emission)
plt.plot(wavelengths, approx_emission)
plt.plot(wavelengths[idx], approx_emission[idx], 'ro')

plt.legend(('Black body emission',
            'Rayleigh Jeans approx',
            '1%% difference at %d μm' %int(one_percent_wavelength)))
plt.xlabel('Wavelength (μm)')
plt.ylabel('W m^-2 μm^-1 sr^-1')
plt.yscale('log')
plt.title('Blackbody at 300K')




'''
Problem 5
'''
wave = 10

power = blackbody(wave, temp=300)
power * wave / (h*c)
# 4.97e26


'''
Problem 6
'''

def freq_blackbody(freq, temp=300):
    # expects wavelength in meters

    c = 3e8  # um/s
    h = 6.626e-34  # Js
    k = 1.38e-23  # J/K

    # returns W m^-2 m^-1 sr^-1
    return (2*h*freq**3 / c**2) * 1/(np.exp(h*freq / (temp*k)) - 1)


raw_wavelengths = np.linspace(1e-6, 20e-6, num=1000)

wave_spectrum = blackbody(raw_wavelengths, temp=700)

plt.figure('test', figsize=(8, 6))
plt.plot(raw_wavelengths*1e6, wave_spectrum*1e-6)

plt.xlabel('Wavelength (μm)')
plt.ylabel('W m^-2 μm^-1 sr^-1')

# convert to freq
freq_spectrum = wave_spectrum * (raw_wavelengths**2)/c

freqs = c/raw_wavelengths

theo_freq_spectrum = freq_blackbody(freqs, temp=700)

tera_freqs = freqs * 1e-12  # gigahertz
plt.figure('Freqencies', figsize=(8, 6))
plt.plot(tera_freqs, freq_spectrum*1e12, 'o')
plt.plot(tera_freqs, theo_freq_spectrum*1e12)
plt.xlabel('Frequency (THz)')
plt.ylabel('W m^-2 GHz^-1 sr^-1')
plt.legend(('Converted',
            'Actual'))
plt.title('Frequency converted from wavelength for a 700K blackbody')





