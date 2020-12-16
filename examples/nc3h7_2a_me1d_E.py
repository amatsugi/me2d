#! /usr/bin/env python3

"""
Thermal decomposition of n-C3H7
(single-well, two-channel, 1D ME as a function of E)

sample output (nc3h7_2a_me1d_E.dat):
   T[K]    p[bar]        k1[s-1]      k2[s-1]    ktot[s-1]
  1000.0    HPL       6.6850e+06   3.4891e+05   7.0339e+06
  1000.0  1.000e+02   5.4397e+06   2.4148e+05   5.6811e+06
  1000.0  1.000e+01   3.2042e+06   1.0383e+05   3.3080e+06
  1000.0  1.000e+00   1.1987e+06   2.3415e+04   1.2221e+06
  1000.0  1.000e-01   3.0430e+05   2.9522e+03   3.0726e+05
"""

from me2d import ME1D


maxE = 50000  # cm^-1
rrkmEfn = "nc3h7_1_rrkmE_dE10.dat"
outfn = "nc3h7_2a_me1d_E.dat"

solver = "InvIter,cho" # inverse iteration with Cholesky decomposition
bandpcrit = None       # use dense matrix
nthreads = 2
maxmemGB = 1
verbose = True

T = 1000. # K
pl = [100., 10., 1., 0.1]  # bar

y = 0.5
alpha = 48.9  # cm^-1
Z = 1.09e-09  # cm^3 molecule^-1 s^-1

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


