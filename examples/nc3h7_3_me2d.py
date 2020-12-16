#! /usr/bin/env python3

"""
Thermal decomposition of n-C3H7
(single-well, two-channel, 2D ME)

sample output (nc3h7_3_me2d.dat):
   T[K]    p[bar]        k1[s-1]      k2[s-1]    ktot[s-1]
  1000.0    HPL       6.6249e+06   3.3274e+05   6.9576e+06
  1000.0  1.000e+02   5.1990e+06   2.0915e+05   5.4081e+06
  1000.0  1.000e+01   2.8637e+06   7.5828e+04   2.9396e+06
  1000.0  1.000e+00   9.8014e+05   1.3255e+04   9.9340e+05
  1000.0  1.000e-01   2.2732e+05   1.1960e+03   2.2852e+05
"""

from me2d import ME2D

maxE = 50000  # cm^-1
maxJ = 200
dE = 50
dJ = 10
rrkmEJfn = "nc3h7_1_rrkmEJ_dE%d.dat" % (dE)
outfn = "nc3h7_3_me2d.dat"

solver = "InvIter,cho" # inverse iteration with Cholesky decomposition
bandpcrit = 5e-7       # truncation threshold for banded matrix
nthreads = 2
maxmemGB = 2
verbose = True

T = 1000.  # K
pl = [100., 10., 1., 0.1]  # bar

y_e = 0.5
y_J = 1.0
Z = 1.07e-09  # cm^3 molecule^-1 s^-1
a_e = lambda E0,J0: 22.8 + J0 * 0.254  # cm^-1
a_J = lambda E0,J0: 9.91

me = ME2D.read_from(rrkmEJfn, dJ=dJ, maxE=maxE, maxJ=maxJ)
me.set_params(Z, y_e, y_J, a_e, a_J)
outfp = open(outfn, "w")
outfp.write("   T[K]    p[bar]   %s\n" % 
            (" ".join("%12s" % x for x in me.get_channel_strings()+["ktot[s-1]"])))

ktot, kchl, ga = me.hpl(T)
outfp.write("%8.1f    HPL     %s\n" % (T, " ".join("%12.4e" % x for x in kchl+[ktot])))
outfp.flush()

for p in pl:
    ktot, kchl, ga, vals, vec = me.solve(T, p, gguess=ga, solver=solver,
                                         verbose=verbose, bandpcrit=bandpcrit,
                                         nthreads=nthreads, maxmemGB=maxmemGB)
    outfp.write("%8.1f  %.3e %s\n" % (T, p, " ".join("%12.4e" % x for x in kchl+[ktot])))
    outfp.flush()

