"""
hw 4 ses 598

calculations
"""


'''
problem 1
'''
k = 1.38e-23  # boltzmann constant k_b in J/K

amu_kg = 1.66e-27  # kg/amu

60 * k * 300 / (10e5 * 10e-3) * 29/(4 * (29-4) * amu_kg)


218 * (750/200)**0.286
# 318.15


'''
Question 3
'''

# Venus
8.87 / (4 * k/(44.0 * amu_kg))
# 0.0117

# Mars
3.71/8.87 * .01174
3.71 / (4 * k/(44.0 * amu_kg))

# Jupiter
c_bar = .85*(7/2)/2.02 + .15*(5/2)/4.00
21 / (c_bar/amu_kg * k)
# 0.00161

# water exoplanet
10 / (4 * k/(18.02 * amu_kg))


'''
Question 4
'''
amu = .21*32.0 + .79*28.0
# 28.8

dT = 5 / (10 * 60**2)

watts_kg = 5/2 * (k * dT)/(amu * amu_kg)

watts_kg * 1.25 * 10
# 1.25

