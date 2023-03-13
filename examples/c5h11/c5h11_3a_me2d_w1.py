#! /usr/bin/env python3

"""
Thermal decomposition in the n-,s-C5H11 system
(two-well, two-channel, 2D ME)

Steady-state decomposition of n-C5H11
sample output (c5h11_3a_me2d_w1.dat):
   T[K]    p[bar]         w1-ch2       w2-ch2    ktot[s-1]        x(w1)        x(w2)
  1000.0  1.000e+02   4.7328e+06   3.2296e+06   7.9623e+06   7.0538e-01   2.9462e-01
  1000.0  1.000e+01   3.4178e+06   2.7322e+06   6.1500e+06   7.0220e-01   2.9780e-01
  1000.0  1.000e+00   1.3940e+06   1.6514e+06   3.0454e+06   6.8860e-01   3.1140e-01
  1000.0  1.000e-01   2.9046e+05   6.5476e+05   9.4522e+05   6.5329e-01   3.4671e-01
"""

from me2d import ME2DMW


maxE = 60000  # cm^-1
maxJ = 400
dE = 50
dJ = 20

# list of (name, rrkm_filename, relative_energy)
well_list = [("w1", "c5h11_1_rrkmEJ_nc5h11_dE%d.dat" %(dE), 948.38),
             ("w2", "c5h11_1_rrkmEJ_sc5h11_dE%d.dat" %(dE), 0.)]
# list of ((name, ch), (name, ch))
connections = [(("w1", 1), ("w2", 1))]

outfn = "c5h11_3a_me2d_w1.dat"

solver = "InvIter,cg" # inverse iteration with conjugate gradient
reactant = "w1"       # solve steady-state decomposition of w1
neig = 1              # neig = 1 for InvIter solver
bandpcrit = 5e-7      # truncation threshold for banded matrix
nthreads = 2
maxmemGB = 4
verbose = True

T = 1000. # K
pl = [100., 10., 1., 0.1]  # bar

y_e = 0.5
y_J = 1.0
Z_w1 = 1.35e-09  # cm^3 molecule^-1 s^-1
a_e_w1 = lambda E0,J0: 32.2 + J0 * 0.119  # cm^-1
a_J_w1 = lambda E0,J0: 22.4
Z_w2 = 1.35e-09  # cm^3 molecule^-1 s^-1
a_e_w2 = lambda E0,J0: 36.6 + J0 * 0.104  # cm^-1
a_J_w2 = lambda E0,J0: 22.1

memw = ME2DMW.read_from(well_list, connections, dJ=dJ, maxE=maxE, maxJ=maxJ)
memw["w1"].set_params(Z_w1, y_e, y_J, a_e_w1, a_J_w1)
memw["w2"].set_params(Z_w2, y_e, y_J, a_e_w2, a_J_w2)

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

