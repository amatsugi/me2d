#! /usr/bin/env python3

"""
ILT (inverse laplace transform) calculation of k(E) for H2O2 => OH + OH
output file: h2o2_iltE_dE10.dat
"""

from me2d import RoVib
from me2d import ilt

nsym = 2
rotA = 10.3560
rotB2D = 0.84680
freq = [877, 1266, 1402, 3599, 3608]
# HO-OH rotational energy levels
levels = [0.0, 13.28417, 251.96034, 375.2889, 576.38451, 786.65549, 1014.27698,
          1251.27288, 1492.82778, 1738.07658, 1965.91838, 2224.64158, 2354.85548,
          2728.62328, 2750.95608, 3300.27678, 3301.84558, 3960.96818, 3961.04278,
          4709.42418, 4709.42688, 5542.65208, 5542.65218, 6458.95218, 6458.95218,
          7457.38088, 7457.38088, 8537.37598, 8537.37598, 9698.58218, 9698.58218,
          10940.76388, 10940.76388, 12263.76088, 12263.76088, 13667.45788, 13667.45788,
          15151.77188, 15151.77188, 16716.64188, 16716.64188, 18362.02188, 18362.02188,
          20087.87688, 20087.87688, 21894.17788, 21894.17788, 23780.90488, 23780.90488,
          25748.03888, 25748.03888, 27795.56688, 27795.56688, 29923.47788, 29923.47788,
          32131.76188, 32131.76188, 34420.41188, 34420.41188, 36789.42088, 36789.42088,
          39238.78388, 39238.78388, 41768.49588, 41768.49588, 44378.55388, 44378.55388,
          47068.95488, 47068.95488, 49839.69488, 49839.69488, 52690.77188, 52690.77188,
          55622.18388, 55622.18388, 58633.92988, 58633.92988, 61726.00688, 61726.00688,
          64898.41488, 64898.41488]
states = [(1, x) for x in levels]
rovib = RoVib(nsym, rotA, rotB2D, freq, states=states)
    
# ILT [k(T) = A*exp(-E/RT)]
ilt_A = 6.7e14    # s^-1
ilt_E = 17236.884 # cm^-1

maxE = 65000
dE = 10
ilt(maxE, dE, rovib, ilt_A, ilt_E, outf="h2o2_iltE_dE%d.dat" % (dE))
