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


def q_day(e, delta, biglamda, omega, phi, s_0=1365):
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
    delta: angle of subsolar point, or obliquity
    phi: latitude

    s_0: average solar irradiance

    Returns
    -------

    """
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

        # calculate the fractional difference in distance from the sun
        d = ((1 + e*np.cos(omega - biglamda)) / (1 - e**2))**2

        output.append(s_0/np.pi * d * (h_0 * np.sin(subphi) * np.sin(delta) + np.cos(subphi) * np.cos(delta) * np.sin(h_0)))

    return np.array(output)


filename = 'hw2_data/Earth_Milankovitch_Laskar.txt'

data = open_tabular(filename, headerlines=1)



