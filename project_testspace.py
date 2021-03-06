import numpy as np
import matplotlib.pyplot as plt


def delta_z(sigma, scale_h, wavelengths):
    # need wavelengths to properly scale the derivative
    return scale_h * np.gradient(np.log(sigma), wavelengths)


def delta_alpha(sigma, scale_h, wavelengths, star_radius, planet_radius):
    return (2 * planet_radius * delta_z(sigma, scale_h, wavelengths))/(star_radius**2)


'''generate derivative of absorption profile'''
with open('project_data/1H2-16O_6250-12500_300K_20.000000.sigma') as file:
    raw_data = file.readlines()
    wave_numbers = []
    cross_sections = []
    for x in raw_data:
        wave_string, cross_string = x.split()
        wave_numbers.append(float(wave_string))
        cross_sections.append(float(cross_string))
    wave_numbers = np.array(wave_numbers)
    cross_sections = np.array(cross_sections)

# convert wavenumber to wavelength
cross_wavelengths = 1e7/wave_numbers

# scale height
k = 1.38e-23  # boltzmann constant k_b in J/K
amu_kg = 1.66e-27  # kg/amu
g = 10  # m/s^2
T = 290  # K
mass = 18  # amu
H = k*T/(mass*amu_kg * g)

first_order_depth = delta_alpha(sigma=cross_sections,
                                scale_h=H,
                                wavelengths=cross_wavelengths,
                                star_radius=6.957e8,
                                planet_radius=6.3710e6)

# mean spectral resolution
spec_res = np.mean(np.diff(np.flip(cross_wavelengths)))

plt.figure('transit depth %.2f' %spec_res, figsize=(8, 8))
plt.subplot(212)
plt.plot(cross_wavelengths, cross_sections)
plt.title('Cross section of H2O')
plt.xlabel('Wavelength (nm)')
plt.ylabel('cm^2/molecule')

plt.subplot(211)
plt.plot(cross_wavelengths, first_order_depth)
plt.title('1st order transit depth, resolution %.2f nm' %spec_res)
plt.axhline(c='k')
# plt.xlabel('Wavelength (nm)')




plt.figure('test')
plt.plot(cross_wavelengths, cross_sections)
plt.plot(cross_wavelengths, delta_z(sigma=cross_sections,
                                    scale_h=H,
                                    wavelengths=cross_wavelengths))



