"""
SES 598 - Atmospheres
HW 2 code
"""

import numpy as np
import matplotlib.pyplot as plt


def open_tabular(data_file, headerlines=0):
    """
    Opens a txt file with tabulated data, and returns it as column lists

    Parameters
    ----------
    data_file

    Returns
    -------

    """
    with open(data_file, 'r') as file:
        while headerlines > 0:
            file.readline()
            print('File line skipped')
            headerlines -= 1

        raw_data = file.readlines()
        # print(raw_data)
        out_data = []
        for x in raw_data:
            out_data.append(x.split())
        # transpose it from row (lines) to columns
        out_data = list(map(list, zip(*out_data)))
        # print(out_data)

    return out_data


def string_to_array_1d(string_list):
    """
    Converts a 1D list of strings into an array of floats

    Parameters
    ----------
    string_list

    Returns
    -------
    1D numpy array of floats
    """
    return np.array([float(number) for number in string_list])


def q_day(e, obliquity, biglamda, omega, phi, s_0=1365):
    """
    q_day = s0/pi (d0/d)**2 (h0*sin(phi)*sin(delta) + cos(phi)cos(delta)(sinh0))

    where h0= arccos(tan(phi) * tan(delta)), with special cases for the sun never setting or
    never rising at the poles.

    This calculation is done for every phi.

    Need to replace (d0/d)**2 with a more complex model, one that reflects actual orbital
    mechanics.

    distance model:
    ( (1 + e*cos(omega - biglambda)/(1 - e**2) )**2

    Parameters
    ----------
    e: eccentricity
    obliquity: radians
        used to calculate delta, the angle of subsolar point
    biglamda: radians
        angle between spring equinox and perihelion
    omega: radians
        current angle of orbit, measured from spring equinox (obliquity=0)
    phi: degrees
        latitude

    s_0: average solar irradiance

    Returns
    -------
    Average daily insolation at the latitudes given
    """
    output = []
    # in degrees, need radians
    # delta = np.radians(delta)
    biglamda = biglamda + np.pi  # pi radians phase shift
    # omega = np.radians(omega)
    phi = np.radians(phi)

    for subphi in phi:
        # calculate the relative tilt to the sun
        delta = obliquity * np.sin(omega)

        data_sanitizer = -np.tan(subphi) * np.tan(delta)
        # if -tan(phi)tan(delta) < -1, sun never sets
        data_sanitizer[data_sanitizer < -1] = -1
        # if -tan(phi)tan(delta) > 1, sun never rises
        data_sanitizer[data_sanitizer > 1] = 1

        # now that the data is sanitized, calculate sunrise/sunset angle
        # cos(-1) = pi, sun never sets
        # cos(1) = 0, sun never rises
        h_0 = np.arccos(data_sanitizer)

        # calculate the fractional difference in distance from the sun
        d = ((1 + e*np.cos(omega - biglamda)) / (1 - e**2))**2

        output.append(s_0/np.pi * d * (h_0 * np.sin(subphi) * np.sin(delta) +
                                       np.cos(subphi) * np.cos(delta) * np.sin(h_0)))

    return np.array(output)

"""
Question 1
"""
filename = 'hw2_data/Earth_Milankovitch_Laskar.txt'

data = open_tabular(filename, headerlines=1)
year = string_to_array_1d(data[0])
eccentricity = string_to_array_1d(data[1])
obliquity = string_to_array_1d(data[2])
arg_perhi = string_to_array_1d(data[3])

omega = np.pi/2  # keep omega in radians
phi = [33.4]  # keep latitude in degrees

insolation = q_day(e=eccentricity,
                   obliquity=obliquity,
                   biglamda=arg_perhi,
                   omega=omega,
                   phi=phi)

plt.figure('phoenix_inso', figsize=(8, 6))
plt.plot(year, insolation[0])
# limit to 1000 Kyr
plt.xlim(0, 1000)
plt.title('Insolation at 33.4° N latitude (Pheonix)')
plt.ylabel('Daily average insolation (W/m^2)')
plt.xlabel('Time in the future (Kyr)')
plt.savefig('phoenix_inso')





