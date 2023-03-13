#! /usr/bin/env python3

"""
Thermal decomposition in the n-,s-C5H11 system
(two-well, two-channel, 1D ME as a function of E)

Steady-state decomposition of n-C5H11
sample output (c5h11_2a_me1d_E_w1.dat):
   T[K]    p[bar]         w1-ch2       w2-ch2    ktot[s-1]        x(w1)        x(w2)
  1000.0  1.000e+02   4.7452e+06   3.2981e+06   8.0433e+06   6.9810e-01   3.0190e-01
  1000.0  1.000e+01   3.5927e+06   2.8640e+06   6.4566e+06   6.9099e-01   3.0901e-01
  1000.0  1.000e+00   1.6175e+06   1.8428e+06   3.4603e+06   6.7063e-01   3.2937e-01
  1000.0  1.000e-01   3.8500e+05   7.9742e+05   1.1824e+06   6.2665e-01   3.7335e-01
"""

from me2d import ME1DMW


maxE = 60000  # cm^-1

# list of (name, rrkm_filename, relative_energy)
well_list = [("w1", "c5h11_1_rrkmE_nc5h11_dE10.dat", 948.38),
             ("w2", "c5h11_1_rrkmE_sc5h11_dE10.dat", 0.)]
# list of ((name, ch), (name, ch))
connections = [(("w1", 1), ("w2", 1))]

outfn = "c5h11_2a_me1d_E_w1.dat"

solver = "InvIter,cg" # inverse iteration with conjugate gradient
reactant = "w1"       # solve steady-state decomposition of w1
neig = 1              # neig = 1 for InvIter solver
bandpcrit = 1e-9      # truncation threshold for banded matrix
nthreads = 2
maxmemGB = 1
verbose = True

T = 1000. # K
pl = [100., 10., 1., 0.1]  # bar

y = 0.5
alpha_w1 = 54.8   # cm^-1
Z_w1 = 1.36e-09   # cm^3 molecule^-1 s^-1
alpha_w2 = 54.4   # cm^-1
Z_w2 = 1.43e-09   # cm^3 molecule^-1 s^-1

memw = ME1DMW.read_from(well_list, connections, maxE=maxE)
memw["w1"].set_params(Z_w1, y, alpha_w1)
memw["w2"].set_params(Z_w2, y, alpha_w2)

outfp = open(outfn, "w")
kstrl = memw.get_channel_strings()
xstrl = ["x(%s)" % name for name in memw.names]
outfp.write("   T[K]    p[bar]   %s\n" % 
            (" ".join("%12s" % x for x in kstrl+["ktot[s-1]"]+xstrl)))

ga = None
for p in pl:
    ktot, kchl, ga, popl, vals, vec = \
          memw.solve(T, p, gguess=ga, solver=solver, reactant=reactant,
                     neig=neig, verbose=verbose, bandpcrit=bandpcrit,
                     nthreads=nthreads, maxmemGB=maxmemGB)
    outfp.write("%8.1f  %.3e %s\n" % (T, p, " ".join("%12.4e" % x for x in kchl+[ktot]+popl)))
    outfp.flush()

