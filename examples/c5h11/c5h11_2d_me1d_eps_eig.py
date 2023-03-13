#! /usr/bin/env python3

"""
Thermal decomposition in the n-,s-C5H11 system
(two-well, two-channel, 1D ME as a function of eps with centrifugal correction)

Steady-state decomposition of equilibrated mixture of n- and s-C5H11 and eigenvalues
sample output (c5h11_2d_me1d_eps_eig.dat):
   T[K]    p[bar]         w1-ch2       w2-ch2    ktot[s-1]        x(w1)        x(w2)          ev1          ev2          ev3          ev4          ev5
  1000.0  1.000e+02   1.5586e+06   7.7691e+06   9.3276e+06   2.3660e-01   7.6340e-01  -9.3276e+06  -1.2711e+07  -5.1937e+09  -6.1847e+09  -9.5693e+09
  1000.0  1.000e+01   8.8089e+05   5.6689e+06   6.5498e+06   1.8964e-01   8.1036e-01  -6.5498e+06  -9.1457e+06  -6.2688e+08  -7.3827e+08  -1.2983e+09
  1000.0  1.000e+00   2.8634e+05   2.5993e+06   2.8857e+06   1.4310e-01   8.5690e-01  -2.8857e+06  -4.3597e+06  -8.8340e+07  -9.9899e+07  -1.9309e+08
  1000.0  1.000e-01   5.9480e+04   7.6938e+05   8.2886e+05   1.1367e-01   8.8633e-01  -8.2886e+05  -1.4282e+06  -1.3146e+07  -1.4759e+07  -2.9102e+07
"""

from me2d import ME1DMW


maxE = 60000  # cm^-1

# list of (name, rrkm_filename, relative_energy)
well_list = [("w1", "c5h11_1_rrkmeps_nc5h11_dE10.dat", 948.38),
             ("w2", "c5h11_1_rrkmeps_sc5h11_dE10.dat", 0.)]
# list of ((name, ch), (name, ch))
connections = [(("w1", 1), ("w2", 1))]

outfn = "c5h11_2d_me1d_eps_eig.dat"

solver = "EigIter,cg" # ARPACK (invert mode) with conjugate gradient
reactant = None       # no reactant specification for EigIter
neig = 5              # find 5 least-negative eigenvalues
bandpcrit = 1e-9      # truncation threshold for banded matrix
nthreads = 2
maxmemGB = 1
verbose = True

T = 1000. # K
pl = [100., 10., 1., 0.1]  # bar

y = 0.5
alpha_w1 = 40.9  # cm^-1
Z_w1 = 1.27e-09  # cm^3 molecule^-1 s^-1
alpha_w2 = 43.3  # cm^-1
Z_w2 = 1.34e-09  # cm^3 molecule^-1 s^-1

memw = ME1DMW.read_from(well_list, connections, maxE=maxE)
memw["w1"].set_params(Z_w1, y, alpha_w1)
memw["w2"].set_params(Z_w2, y, alpha_w2)

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


