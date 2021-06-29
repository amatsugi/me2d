#! /usr/bin/env python3

"""
Thermal decomposition of n-C3H7
(single-well, two-channel, 1D ME as a fuction of eps with centrifugal correction)

sample output (nc3h7_2b_me1d_eps.dat):
   T[K]    p[bar]        k1[s-1]      k2[s-1]    ktot[s-1]
  1000.0    HPL       6.6850e+06   3.3825e+05   7.0232e+06
  1000.0  1.000e+02   4.9303e+06   1.9731e+05   5.1276e+06
  1000.0  1.000e+01   2.5653e+06   6.9820e+04   2.6351e+06
  1000.0  1.000e+00   8.4877e+05   1.2749e+04   8.6151e+05
  1000.0  1.000e-01   1.9506e+05   1.2732e+03   1.9633e+05
"""

from me2d import ME1D


maxE = 50000  # cm^-1
rrkmEfn = "nc3h7_1_rrkmeps_dE10.dat"
outfn = "nc3h7_2b_me1d_eps.dat"

solver = "InvIter,cho" # inverse iteration with Cholesky decomposition
bandpcrit = None       # use dense matrix
nthreads = 2
maxmemGB = 1
verbose = True

T = 1000. # K
pl = [100., 10., 1., 0.1]  # bar

y = 0.5
alpha = 37.1  # cm^-1
Z = 1.07e-09  # cm^3 molecule^-1 s^-1

me = ME1D.read_from(rrkmEfn, maxE=maxE)
me.set_params(Z, y, alpha)
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

