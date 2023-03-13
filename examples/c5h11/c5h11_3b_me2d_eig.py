#! /usr/bin/env python3

"""
Thermal decomposition in the n-,s-C5H11 system
(two-well, two-channel, 2D ME)

Steady-state decomposition of equilibrated mixture of n- and s-C5H11 and eigenvalues
sample output (c5h11_3b_me2d_eig.dat):
   T[K]    p[bar]         w1-ch2       w2-ch2    ktot[s-1]        x(w1)        x(w2)          ev1          ev2          ev3          ev4          ev5
  1000.0  1.000e+02   1.6402e+06   7.9684e+06   9.6086e+06   2.4322e-01   7.5678e-01  -9.6131e+06  -1.3071e+07  -6.8163e+09  -7.9500e+09  -1.2926e+10
  1000.0  1.000e+01   1.0036e+06   6.1290e+06   7.1326e+06   2.0007e-01   7.9993e-01  -7.1355e+06  -9.8239e+06  -7.6854e+08  -8.8914e+08  -1.5797e+09
  1000.0  1.000e+00   3.3975e+05   2.9581e+06   3.2979e+06   1.4968e-01   8.5032e-01  -3.2990e+06  -4.8171e+06  -9.9798e+07  -1.1159e+08  -2.1325e+08
  1000.0  1.000e-01   6.9802e+04   8.8642e+05   9.5622e+05   1.1595e-01   8.8405e-01  -9.5657e+05  -1.5561e+06  -1.3857e+07  -1.5192e+07  -2.9512e+07
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

outfn = "c5h11_3b_me2d_eig.dat"

solver = "EigIter,cg" # ARPACK (invert mode) with conjugate gradient
reactant = None       # no reactant specification for EigIter
neig = 5              # find 5 least-negative eigenvalues
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
if neig > 1: addstrl = ["ev%s" % (i+1) for i in range(neig)]
else: addstrl = []
outfp.write("   T[K]    p[bar]   %s\n" % 
            (" ".join("%12s" % x for x in kstrl+["ktot[s-1]"]+xstrl+addstrl)))

ga = None
for p in pl:
    ktot, kchl, ga, popl, vals, vec = \
          memw.solve(T, p, gguess=ga, solver=solver, reactant=reactant,
                     neig=neig, verbose=verbose, bandpcrit=bandpcrit,
                     nthreads=nthreads, maxmemGB=maxmemGB)
    if neig > 1: addl = list(vals)
    else: addl = []
    outfp.write("%8.1f  %.3e %s\n" % (T, p, " ".join("%12.4e" % x for x in kchl+[ktot]+popl+addl)))
    outfp.flush()

