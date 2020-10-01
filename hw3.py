"""
homework 3

"""
import numpy as np

"""
Question 2
"""
.85*2 + .15*4
# 2.3

(1.38e-23 * 120) / (2.3 * 1.66e-27 * 20)
# 21687

1.7/8.5
# 2.55


"""
Question 3
"""
'%.2E' % (3/(4*np.pi) * 6e24/(6372e3**3))
# '5.54E+3'

'%.2E' % (2/3 * np.pi * 6.67e-11 * 5.54e3**2 * 6372e3**2)
# '3.05E+24'


"""
Question 4
"""
k = 1.38e-23  # boltzmann constant k_b
amu_kg = 1.66e-27  # kg/amu
G = 6.67e-11

(k * 90 * 2600e3**2) / (28*amu_kg * 1.34e23 * G)
# 20210 m


