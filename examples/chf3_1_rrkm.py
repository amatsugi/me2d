#! /usr/bin/env python3

"""
RRKM calculation for CHF3 => CF2 + HF
[Ref: A. Matsugi, J. Phys. Chem. A 124 (2020) 6645.]

This code generates chf3_1_rrkmE_dE10.dat and chf3_1_rrkmEJ_dE50.dat files
for 1D and 2D ME calculations, respectively.
"""

from me2d import RoVib
from me2d import rrkmE, rrkmEJ

# reactant
reactant = RoVib(\
    nsym = 3,
    rotA = 0.46925, # A [cm^-1]
    rotB2D = 0.22, # Beff [cm^-1]
    freq = [510.7, 510.7, 707.6, 1158.9, 1189.0, 1189.0,
            1408.1, 1408.1, 3158.1], # cm^-1
    fscale = 0.94
    )

# ts
ts = RoVib(\
    nsym = 1,
    rotA = 0.37266,
    rotB2D = 0.17679,
    freq = [213.0, 284.6, 654.6, 843.7, 1190.0, 1252.0, 1349.9, 2317.0],
    freqimg = 1347. # for Eckart tunneling correction [cm^-1]
    )

E0 = 25471.  # 0K (ZPE-corrected) barrier height [cm^-1]

# product (for Eckart tunneling correction)
deltaH0 = 18441. # 0K (ZPE-corrected) reaction energy [cm^-1]
rotB2Dprod = 0   # B2D of product (zero for bimolecular products) [cm^-1]


# rho(E) and k(E)
maxE = 55000   # max. value of energy [cm^-1]
dE = 10   # energy grain size [cm^-1]
rrkmE(maxE, dE, reactant, [ts], [E0], [deltaH0], outf="chf3_1_rrkmE_dE%d.dat" % (dE))

# rho(E,J) and k(E,J)
dE = 50
maxJ = 260  # max. value of J [hbar]
rrkmEJ(maxE, dE, maxJ, reactant, [ts], [E0], [deltaH0], [rotB2Dprod],
       outf="chf3_1_rrkmEJ_dE%d.dat" % dE, verbose=True)

