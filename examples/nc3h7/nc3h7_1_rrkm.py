#! /usr/bin/env python3

"""
RRKM calculation for multichannel decomposition of n-C3H7 radical
nC3H7 => C2H4 + CH3
nC3H7 => C3H6 + H

generates three files:
nc3h7_1_rrkmE_dE10.dat    --- for 1D ME as a function of E
nc3h7_1_rrkmEeps_dE10.dat --- for 1D ME as a function of eps with centrifugal correction
nc3h7_1_rrkmEJ_dE50.dat   --- for 2D ME as a function of E and J
"""

from me2d import RoVib
from me2d import rrkmE, rrkmEJ

# reactant
nc3h7 = RoVib(\
    nsym = 1,
    rotA = 1.10189,
    rotB2D = 0.279377,
    freq = [250.196, 379.029, 467.993, 762.122, 890.026, 932.685, 1052.69, 1100.62,
            1174.23, 1272.97, 1366.72, 1418.08, 1466.2, 1471.32, 1502.51, 1508.29,
            2965.42, 3031.31, 3041.89, 3115.12, 3122.37, 3146.24, 3253.11],
    fscale = 0.95,
    introt = [(10.5666, 2, 0, 0)],  # free rotor
    states = [(1, 0)]
    )


# C-C fission ts
ts_cc = RoVib(\
    nsym = 1,
    rotA = 0.921745,
    rotB2D = 0.202788,
    freq = [227.001, 357.09, 515.069, 537.441, 812.099, 829.171, 882.218, 961.145,
            1025.23, 1240.15, 1309.98, 1418.73, 1424.42, 1464.61, 1590.76, 3094.79,
            3143.73, 3155.49, 3226.41, 3251.87, 3258.36, 3268.56],
    fscale = 0.95,
    introt = [(5.9042, 3, 0, 0)],  # free rotor
    states = [(1, 0)],
    freqimg = 471.661
    )
E0_cc = 10934.13
deltaH0_cc = 7740.52 # C2H4 + CH3
rotB2Dprod_cc = 0

# C-H fission ts
ts_ch = RoVib(\
    nsym = 1,
    rotA = 1.08246,
    rotB2D = 0.278345,
    freq = [181.989, 372.051, 429.281, 439.386, 596.959, 925.403, 930.241, 952.843,
            1031.16, 1073.08, 1194.57, 1313.9, 1409.77, 1449.54, 1485.65, 1499.23,
            1654.67, 3043.62, 3111.73, 3133.05, 3146.45, 3160.58, 3242.67],
    fscale = 0.95,
    states = [(2, 0)], # optical isomer
    freqimg = 758.838
    )
E0_ch = 12565.05
deltaH0_ch = 11138.38  # C3H6 + H
rotB2Dprod_ch = 0


# rho(E) and k(E)
maxE = 55000
dE = 10
rrkmE(maxE, dE, nc3h7, [ts_cc, ts_ch], [E0_cc, E0_ch], [deltaH0_cc, deltaH0_ch],
      outf="nc3h7_1_rrkmE_dE%d.dat" % (dE))

# rho(eps) and k(eps) (with centrifugal correction)
rrkmE(maxE, dE, nc3h7, [ts_cc, ts_ch], [E0_cc, E0_ch], [deltaH0_cc, deltaH0_ch],
      convJ=False, centcorr=True, outf="nc3h7_1_rrkmeps_dE%d.dat" % (dE))

# rho(E,J) and k(E,J)
dE = 50
maxJ = 260
rrkmEJ(maxE, dE, maxJ, nc3h7, [ts_cc, ts_ch], [E0_cc, E0_ch], [deltaH0_cc, deltaH0_ch],
       [rotB2Dprod_cc, rotB2Dprod_ch],
       outf="nc3h7_1_rrkmEJ_dE%d.dat" % dE, verbose=True)

