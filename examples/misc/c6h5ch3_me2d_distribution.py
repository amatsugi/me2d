#! /usr/bin/env python3

"""
The thermal decomposition of toluene
(single-well, two-channel, 2D ME)

output files:
c6h5ch3_me2d_distribution_E.dat  --- distribution over E
c6h5ch3_me2d_distribution_J.dat  --- distribution over J
"""

import numpy as np
from me2d import ME2D

maxE = 90000  # cm^-1
maxJ = 500
dE = 50
dJ = 20
rrkmEJfn = "c6h5ch3_vrrkmEJ_dE%d.dat" % (dE)
outfn_pE = "c6h5ch3_me2d_distribution_E.dat"
outfn_pJ = "c6h5ch3_me2d_distribution_J.dat"

solver = "InvIter,cho" # inverse iteration with Cholesky decomposition
bandpcrit = 5e-7       # truncation threshold for banded matrix
nthreads = 2
maxmemGB = 4
verbose = True

T = 2000.  # K
p = 1.  # bar

y_e = 0.5
y_J = 1.0
Z = 1.53e-09 # cm^3 molecule^-1 s^-1
a_e = lambda E0,J0: 40.6 + J0 * 0.221  # cm^-1
a_J = 27.7

me = ME2D.read_from(rrkmEJfn, dJ=dJ, maxE=maxE, maxJ=maxJ)
me.set_params(Z, y_e, y_J, a_e, a_J)
outfp_pE = open(outfn_pE, "w")
outfp_pJ = open(outfn_pJ, "w")
outfp_pE.write(" T[K] = %s, p[bar] = %s\n" % (T, p))
outfp_pJ.write(" T[K] = %s, p[bar] = %s\n" % (T, p))
outfp_pE.flush()
outfp_pJ.flush()

ktot, kchl, ga = me.hpl(T)
ga_therm = np.copy(ga)

ktot, kchl, ga, vals, vec = me.solve(T, p, gguess=ga, solver=solver,
                                     verbose=verbose, bandpcrit=bandpcrit,
                                     nthreads=nthreads, maxmemGB=maxmemGB)

outfp_pE.write(" %s [s^-1]\n" % 
               (", ".join("kch%d = %.3e" % (ich+1, kchl[ich]) for ich in range(len(kchl)))))
outfp_pJ.write(" %s [s^-1]\n" % 
               (", ".join("kch%d = %.3e" % (ich+1, kchl[ich]) for ich in range(len(kchl)))))
outfp_pE.flush()
outfp_pJ.flush()

Jl = [int(x) for x in me.Jl]
fE = [0. for i in range(me.sizE)]
fJ = [0. for i in range(len(Jl))]
gE = [0. for i in range(me.sizE)]
gJ = [0. for i in range(len(Jl))]
fluxE = [[0. for i in range(me.sizE)] for ich in range(len(kchl))]
fluxJ = [[0. for i in range(len(Jl))] for ich in range(len(kchl))]
for i in range(me.nsiz):
    iE = int(me.Ea[i]/me.dE)
    iJ = Jl.index(int(me.Ja[i]))
    fE[iE] += ga_therm[i]
    fJ[iJ] += ga_therm[i]
    gE[iE] += ga[i]
    gJ[iJ] += ga[i]
    for ich in range(len(kchl)):
        fluxE[ich][iE] += ga[i] * me.kchl[ich][i] / kchl[ich]
        fluxJ[ich][iJ] += ga[i] * me.kchl[ich][i] / kchl[ich]

outfp_pE.write(" f = thermal population distribution\n")
outfp_pE.write(" g = steady-state population distribution\n")
outfp_pE.write(" ch# = normalized dissociation flux distribution for channel #\n")
outfp_pE.write(" E[cm-1]   f(E)         g(E)         %s\n" %
               ("          ".join("ch%d" % (ich+1) for ich in range(len(kchl)))))
for i in range(me.sizE):
    outfp_pE.write("%6d %12.4e %12.4e %s\n" %
                   (me.dE*i, fE[i], gE[i], " ".join("%12.4e" % (x[i]) for x in fluxE)))
outfp_pE.close()

outfp_pJ.write(" f = thermal population distribution\n")
outfp_pJ.write(" g = steady-state population distribution\n")
outfp_pJ.write(" ch# = normalized dissociation flux distribution for channel #\n")
outfp_pJ.write(" J[hbar]   f(J)         g(J)         %s\n" % 
               ("          ".join("ch%d" % (ich+1) for ich in range(len(kchl)))))
for i in range(len(Jl)):
    outfp_pJ.write("%6d %12.4e %12.4e %s\n" %
                   (Jl[i], fJ[i], gJ[i], " ".join("%12.4e" % (x[i]) for x in fluxJ)))
outfp_pJ.close()

