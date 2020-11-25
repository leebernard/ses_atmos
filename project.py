"""
model a spectrum being eclipsed by a planet

coding up Line et al 2016
"""

import numpy as np
import matplotlib.pyplot as plt
import yaml

from astropy.io import fits
from scipy.ndimage import gaussian_filter


def spectrum_slicer(start_angstrom, end_angstrom, angstrom_data, spectrum_data):
    start_index = np.where(angstrom_data == start_angstrom)[0][0]
    end_index = np.where(angstrom_data == end_angstrom)[0][0]
    spectrum_slice = spectrum_data[start_index:end_index]
    angstrom_slice = angstrom_data[start_index:end_index]

    return angstrom_slice, spectrum_slice

with fits.open('sun.fits') as hdul:
    # pull the data
    # hdul = fits.open('sun.fits')
    sun = hdul[0]
    sun_data = sun.data

    angstrom_per_pix = sun.header['CDELT1']
    intial_angstrom = sun.header['CRVAL1']


# generate an x-data set
angstrom = np.arange(0, sun_data.size, step=1) * angstrom_per_pix + intial_angstrom

# find a slice of data
start_ang = 7300
end_ang = angstrom[-1]
angstrom_slice, sun_slice = spectrum_slicer(start_ang, end_ang, angstrom, sun_data)

# filter the spectrum slice
# this should produce a resolution of 1 angstrom
resolution = 500  # 1/.002  # the resolution of the spectrum in angstroms. This corresponds to FWHM
sigma = resolution / (2.0 * np.sqrt(2.0 * np.log(2.0)))
filtered_sun = gaussian_filter(sun_slice, sigma)

'''down sample the data to a lower resolution'''
sim_angstroms_per_pixel = .1  # resolution of the simlulated pixel grid
bin_factor = int(sim_angstroms_per_pixel/angstrom_per_pix)

excess_data_index = int(sun_slice.size % bin_factor)
if excess_data_index:
    pared_angstrom_slice = angstrom_slice[:-excess_data_index]
    pared_sun_slice = filtered_sun[:-excess_data_index]

binned_angstroms = np.reshape(pared_angstrom_slice, (int(pared_angstrom_slice.size/bin_factor), bin_factor))
binned_spectrum = np.reshape(pared_sun_slice, (int(pared_sun_slice.size/bin_factor), bin_factor))

binned_angstroms = np.mean(binned_angstroms, axis=1)
binned_spectrum = np.mean(binned_spectrum, axis=1)  # * sim_angstroms_per_pixel

plt.figure('binned spectrum')
plt.scatter(angstrom_slice, filtered_sun, s=1)
plt.scatter(binned_angstroms, binned_spectrum, s=1)
plt.title('Normalized Slice of Solar spectrum')
plt.xlabel('Angstroms')
plt.legend(['.002 angstroms resolution', '%.3f angstroms resolution' %sim_angstroms_per_pixel])

with open('project_data/1H2-16O_13513-13698_300K_0.185000.sigma') as file:
    raw_data = file.readlines()
    cross_sections = []
    for x in raw_data:
        cross_sections.append(float(x))
    cross_sections = np.array(cross_sections)


def alpha_lambda(z, star_radius, planet_radius):
    return (planet_radius / star_radius)**2 + (2 * planet_radius * z)/star_radius**2


def z_lambda():
    pass



