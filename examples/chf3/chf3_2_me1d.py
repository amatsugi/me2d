#! /usr/bin/env python3

"""
Thermal decomposition of CHF3 (single-well, single-channel, 1D ME)

sample output (chf3_2_me1d.dat):
   T[K]    p[bar]        k1[s-1]    ktot[s-1]
  2000.0    HPL       4.5901e+06   4.5901e+06
  2000.0  1.000e+02   1.1342e+06   1.1342e+06
  2000.0  1.000e+01   3.3541e+05   3.3541e+05
  2000.0  1.000e+00   7.3630e+04   7.3630e+04
  2000.0  1.000e-01   1.2989e+04   1.2989e+04
  1500.0    HPL       8.5745e+03   8.5745e+03
  1500.0  1.000e+02   4.6348e+03   4.6348e+03
  1500.0  1.000e+01   2.0554e+03   2.0554e+03
  1500.0  1.000e+00   6.3833e+02   6.3833e+02
  1500.0  1.000e-01   1.4836e+02   1.4836e+02
"""

from me2d import ME1D
from me2d import collfreq_lj, collfreq_capture


maxE = 50000  # cm^-1
rrkmEfn = "chf3_1_rrkmE_dE10.dat"
outfn = "chf3_2_me1d.dat"


# solver (see "import me2d; me2d.show_solvers()" for detail)
solver = "InvIter,cho"   # inverse iteration with Cholesky decomposition
#solver = "InvIter,lu"   # inverse iteration with LU decomposition
#solver = "InvIter,ldlt" # inverse iteration with LDLT decomposition
#solver = "InvIter,cg"   # inverse iteration with conjugate gradient
#solver = "Eigen"        # direct eigen solver
#solver = "EigIter,cho"  # ARPACK (invert mode) with Cholesky decomposition
#solver = "EigIter,lu"   # ARPACK (invert mode) with LU decomposition
#solver = "EigIter,ldlt" # ARPACK (invert mode) with LDLT decomposition
#solver = "EigIter,cg"   # ARPACK (invert mode) with conjugate gradient

bandpcrit = None   # use dense matrix
#bandpcrit = 5e-7  # truncation threshold for banded matrix

nthreads = 2
maxmemGB = 1
verbose = True

Tl = [2000., 1500.]  # K
pl = [100., 10., 1., 0.1]  # bar

sigmaA, epsK = 3.74, 216. # LJ parameters

me = ME1D.read_from(rrkmEfn, maxE=maxE)
outfp = open(outfn, "w")
outfp.write("   T[K]    p[bar]   %s\n" % 
            (" ".join("%12s" % x for x in me.get_channel_strings()+["ktot[s-1]"])))

for T in Tl:
    y = 1.0
    alpha = 300. * (T/1000.)  # cm^-1
    # Z [cm^3 molecule^-1 s^-1]
    #Z = collfreq_lj("CHF3", "Ar", sigmaA=sigmaA, epsK=epsK, T=T)  # use LJ collision frequency
    Z = collfreq_capture("CHF3", "Ar", sigmaA=sigmaA, epsK=epsK, T=T)  # use capture rate constant
    
    # P(E,E') = C exp[-(|E-E'|/alpha)^y] for downward transition
    me.set_params(Z, y, alpha)
    
    ktot, kchl, ga = me.hpl(T)
    outfp.write("%8.1f    HPL     %s\n" % (T, " ".join("%12.4e" % x for x in kchl+[ktot])))
    outfp.flush()
    for p in pl:
        ktot, kchl, ga, vals, vec = me.solve(T, p, gguess=ga, solver=solver,
                                             verbose=verbose, bandpcrit=bandpcrit,
                                             nthreads=nthreads, maxmemGB=maxmemGB)
        outfp.write("%8.1f  %.3e %s\n" % (T, p, " ".join("%12.4e" % x for x in kchl+[ktot])))
        outfp.flush()

