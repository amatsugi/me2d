#! /usr/bin/env python3

"""
RRKM calculation for n- and s-C5H11
nC5H11 <=> sC5H11
nC5H11 => C2H4 + C3H7
sC5H11 => C3H6 + C2H5

generates six files
c5h11_1_rrkmE_[ns]c5h11_dE10.dat    --- for 1D ME as a function of E
c5h11_1_rrkmEeps_[ns]c5h11_dE10.dat --- for 1D ME as a function of eps with centrifugal correction
c5h11_1_rrkmEJ_[ns]c5h11_dE50.dat   --- for 2D ME as a function of E and J
"""

from me2d import RoVib
from me2d import rrkmE, rrkmEJ

# reactants
nc5h11 = RoVib(\
    nsym = 1,
    rotA = 0.609391,
    rotB2D = 0.0654364,
    freq = [113.598, 120.817, 184.873, 247.846, 401.834, 410.665, 469.376, 731.748,
            780.563, 872.445, 922.568, 947.387, 1043.24, 1079.07, 1098.31, 1105.21,
            1168.56, 1235.09, 1289.54, 1312.72, 1337.06, 1371.45, 1416.27, 1420.43,
            1466.27, 1468.73, 1493.98, 1500.6, 1503.72, 1512.2, 2952.34, 3021.08,
            3024.77, 3032.49, 3036.95, 3054.91, 3077.94, 3109.68, 3116.83, 3146.86,
            3253.57],
    fscale = 0.95,
    introt = [(10.1083, 2, 0, 0)],  # free rotor
    states = [(1, 0), (2, 9.01), (2, 89.1), (2, 143.23), (2, 361.86)] # conformers
    )

sc5h11 = RoVib(\
    nsym = 1,
    rotA = 0.64564,
    rotB2D = 0.0633999,
    freq = [117.076, 181.245, 241.578, 338.637, 405.582, 409.762, 740.463, 862.016,
            887.444, 947.624, 982.493, 1052.52, 1079.03, 1093.84, 1146.64, 1180.52,
            1253.02, 1286.53, 1326.41, 1370.7, 1409.52, 1420.08, 1435.14, 1470.15,
            1477.71, 1488.83, 1498.07, 1502.87, 1510.26, 2940.31, 2972.28, 2998.95,
            3036.22, 3038.24, 3041.47, 3074.74, 3109.06, 3115.43, 3117.51, 3173.02],
    fscale = 0.95,
    introt = [(1.56492, 3, 0, 0),  # free rotor
              (5.94229, 3, 0, 0)],  # free rotor
    states = [(1, 0), (2, 1.99), (2, 5.86), (2, 38.73), (2, 658.54)] # conformers
    )

# nc5h11 <=> sc5h11 isomerization ts
ts_isom = RoVib(\
    nsym = 1,
    rotA = 0.3031,
    rotB2D = 0.0954303,
    freq = [133.483, 187.022, 219.143, 345.646, 480.448, 545.692, 572.832, 810.097,
            855.332, 879.916, 905.965, 926.784, 990.452, 1006.91, 1044.83, 1101.9,
            1137.56, 1182.43, 1229.42, 1249.92, 1307.2, 1341.15, 1374.26, 1403.64,
            1423.6, 1458.3, 1488.63, 1491.28, 1495.01, 1506.2, 1658.78, 3011.58,
            3035.28, 3048.17, 3075.45, 3078.43, 3088.38, 3095.01, 3100.43, 3115.94,
            3185.69],
    fscale = 0.95,
    states = [(2, 0), (2, 32.46)], # conformers
    freqimg = 1878.06
    )

# nc5h11, C-C fission ts
ts_ncc = RoVib(\
    nsym = 1,
    rotA = 0.229743,
    rotB2D = 0.0854228,
    freq = [119.12, 166.293, 250.424, 295.145, 361.419, 413.49, 571.327, 745.144,
            808.473, 826.663, 861.237, 892.759, 930.289, 954.743, 1018.81, 1083.88,
            1103.39, 1190.6, 1238.7, 1286.39, 1305.51, 1372.36, 1418.59, 1463.33,
            1465.82, 1477.05, 1504.31, 1508.79, 1586.45, 2992.47, 3040.75, 3052.41,
            3112.97, 3115.87, 3121.05, 3141.69, 3153.36, 3213.03, 3223.87, 3248.19],
    fscale = 0.95,
    introt = [(1.05143, 3, 0, 0)],  # free rotor
    states = [(2, 0), (2, 235.81), (2, 247.77), (1, 328.21), (2, 390.09)], # conformers
    freqimg = 474.82
    )


# sc5h11, C-C fission ts
ts_scc = RoVib(\
    nsym = 1,
    rotA = 0.259238,
    rotB2D = 0.0758349,
    freq = [107.052, 136.519, 175.27, 190.914, 316.005, 432.051, 545.709, 663.579,
            771.005, 835.18, 892.786, 928.455, 939.933, 995.211, 1036, 1048.61,
            1081.85, 1194.39, 1217.57, 1294.18, 1402.86, 1404.58, 1439.14, 1473.22,
            1481.02, 1488.19, 1490.39, 1495.77, 1600.04, 3008.6, 3018.98, 3073.54,
            3076.19, 3108.26, 3113.13, 3122.64, 3137.03, 3158.63, 3216.14, 3229.78],
    fscale = 0.95,
    introt = [(0.939439, 3, 0, 0)],  # free rotor
    states = [(2, 0), (2, 0.09), (2, 165.67)], # conformers
    freqimg = 494.809
    )


# relative energies [cm^-1]
E0_nc5h11 = 948.38
E0_sc5h11 = 0.00
E0_ts_isom = 9033.04
E0_ts_ncc = 11391.70
H0_ncc =  8625.83
E0_ts_scc = 10392.83
H0_scc = 7591.92


# rho(E) and k(E)
maxE = 70000
dE = 10
rrkmE(maxE, dE, nc5h11, [ts_isom, ts_ncc], [E0_ts_isom - E0_nc5h11, E0_ts_ncc - E0_nc5h11],
      [E0_sc5h11 - E0_nc5h11, H0_ncc - E0_nc5h11],
      outf="c5h11_1_rrkmE_nc5h11_dE%d.dat" % (dE))
rrkmE(maxE, dE, sc5h11, [ts_isom, ts_scc], [E0_ts_isom - E0_sc5h11, E0_ts_scc - E0_sc5h11],
      [E0_nc5h11 - E0_sc5h11, H0_scc - E0_sc5h11],
      outf="c5h11_1_rrkmE_sc5h11_dE%d.dat" % (dE))

# rho(eps) and k(eps) (with centrifugal correction)
rrkmE(maxE, dE, nc5h11, [ts_isom, ts_ncc], [E0_ts_isom - E0_nc5h11, E0_ts_ncc - E0_nc5h11],
      [E0_sc5h11 - E0_nc5h11, H0_ncc - E0_nc5h11], convJ=False, centcorr=True, 
      outf="c5h11_1_rrkmeps_nc5h11_dE%d.dat" % (dE))
rrkmE(maxE, dE, sc5h11, [ts_isom, ts_scc], [E0_ts_isom - E0_sc5h11, E0_ts_scc - E0_sc5h11],
      [E0_nc5h11 - E0_sc5h11, H0_scc - E0_sc5h11], convJ=False, centcorr=True, 
      outf="c5h11_1_rrkmeps_sc5h11_dE%d.dat" % (dE))


# rho(E,J) and k(E,J)
dE = 50
maxJ = 450
rrkmEJ(maxE, dE, maxJ, nc5h11, [ts_isom, ts_ncc], [E0_ts_isom - E0_nc5h11, E0_ts_ncc - E0_nc5h11],
       [E0_sc5h11 - E0_nc5h11, H0_ncc - E0_nc5h11], [sc5h11.rotB2D, 0.],
       outf="c5h11_1_rrkmEJ_nc5h11_dE%d.dat" % dE, verbose=True)
rrkmEJ(maxE, dE, maxJ, sc5h11, [ts_isom, ts_scc], [E0_ts_isom - E0_sc5h11, E0_ts_scc - E0_sc5h11],
       [E0_nc5h11 - E0_sc5h11, H0_scc - E0_sc5h11], [nc5h11.rotB2D, 0.],
       outf="c5h11_1_rrkmEJ_sc5h11_dE%d.dat" % dE, verbose=True)

