#! /usr/bin/env python3

"""
Thermal decomposition of CHF3 (single-well, single-channel, 2D ME)

sample output (chf3_3_me2d.dat):
   T[K]    p[bar]        k1[s-1]    ktot[s-1]
  2000.0    HPL       4.5976e+06   4.5976e+06
  2000.0  1.000e+02   6.4722e+05   6.4722e+05
  2000.0  1.000e+01   1.5991e+05   1.5991e+05
  2000.0  1.000e+00   3.0433e+04   3.0433e+04
  2000.0  1.000e-01   4.7967e+03   4.7967e+03
  1500.0    HPL       8.5780e+03   8.5780e+03
  1500.0  1.000e+02   3.0698e+03   3.0698e+03
  1500.0  1.000e+01   1.0815e+03   1.0815e+03
  1500.0  1.000e+00   2.7406e+02   2.7406e+02
  1500.0  1.000e-01   5.3815e+01   5.3815e+01
"""

from me2d import ME2D

maxE = 50000  # cm^-1
maxJ = 250
dJ = 10
rrkmEJfn = "chf3_1_rrkmEJ_dE50.dat"
outfn = "chf3_3_me2d.dat"


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

#bandpcrit = None  # use dense matrix
bandpcrit = 5e-7   # truncation threshold for banded matrix

nthreads = 2
maxmemGB = 2
verbose = True

Tl = [2000., 1500.]  # K
pl = [100., 10., 1., 0.1]  # bar


me = ME2D.read_from(rrkmEJfn, dJ=dJ, maxE=maxE, maxJ=maxJ)
outfp = open(outfn, "w")
outfp.write("   T[K]    p[bar]   %s\n" % 
            (" ".join("%12s" % x for x in me.get_channel_strings()+["ktot[s-1]"])))

for T in Tl:
    y_e = 0.5
    y_J = 1.0
    Z = 7.94e-10 * (T/1000.)**0.18  # cm^3 molecule^-1 s^-1
    a_e = lambda E0,J0: 8.03*(T/1000.)**1.0 + J0 * 0.212*(T/1000.)**0.6  # cm^-1
    a_J = lambda E0,J0: 8.90*(T/1000.)**0.5
    
    # P(E,J;E',J') = C [S(E,J)/S(E',J')] exp[-(|eps-eps'|/a_e)^y_e] exp[-(|J-J'|/a_J)^y_J]
    me.set_params(Z, y_e, y_J, a_e, a_J)
    
    ktot, kchl, ga = me.hpl(T)
    outfp.write("%8.1f    HPL     %s\n" % (T, " ".join("%12.4e" % x for x in kchl+[ktot])))
    outfp.flush()
    for p in pl:
        ktot, kchl, ga, vals, vec = me.solve(T, p, gguess=ga, solver=solver,
                                             verbose=verbose, bandpcrit=bandpcrit,
                                             nthreads=nthreads, maxmemGB=maxmemGB)
        outfp.write("%8.1f  %.3e %s\n" % (T, p, " ".join("%12.4e" % x for x in kchl+[ktot])))
        outfp.flush()

