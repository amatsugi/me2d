#! /usr/bin/env python3

"""
Thermal decomposition in the n-,s-C5H11 system
(two-well, two-channel, 1D ME as a function of E)

Steady-state decomposition of equilibrated mixture of n- and s-C5H11 and eigenvalues
sample output (c5h11_2b_me1d_E_eig.dat):
   T[K]    p[bar]     w1-k2(dis)   w2-k2(dis)    ktot[s-1]        x(w1)        x(w2)          ev1          ev2          ev3          ev4          ev5
  1000.0  1.000e+02   1.6080e+06   8.0638e+06   9.6718e+06   2.3558e-01   7.6442e-01  -9.6718e+06  -1.3293e+07  -8.9940e+09  -9.5226e+09  -1.5845e+10
  1000.0  1.000e+01   1.0101e+06   6.4007e+06   7.4107e+06   1.8928e-01   8.1072e-01  -7.4107e+06  -1.0432e+07  -1.0473e+09  -1.1134e+09  -2.0894e+09
  1000.0  1.000e+00   3.7006e+05   3.2826e+06   3.6527e+06   1.3847e-01   8.6153e-01  -3.6527e+06  -5.5958e+06  -1.4219e+08  -1.5109e+08  -3.0697e+08
  1000.0  1.000e-01   8.6884e+04   1.0617e+06   1.1486e+06   1.0787e-01   8.9213e-01  -1.1486e+06  -2.0244e+06  -2.0391e+07  -2.3291e+07  -4.5239e+07
"""

from me2d import ME1DMW


maxE = 60000  # cm^-1

# list of (name, rrkm_filename, relative_energy)
well_list = [("w1", "c5h11_1_rrkmE_nc5h11_dE10.dat", 948.38),
             ("w2", "c5h11_1_rrkmE_sc5h11_dE10.dat", 0.)]
# list of ((name, ch), (name, ch))
connections = [(("w1", 1), ("w2", 1))]

outfn = "c5h11_2b_me1d_E_eig.dat"

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
alpha_w1 = 54.8   # cm^-1
Z_w1 = 1.36e-09   # cm^3 molecule^-1 s^-1
alpha_w2 = 54.4   # cm^-1
Z_w2 = 1.43e-09   # cm^3 molecule^-1 s^-1

memw = ME1DMW.read_from(well_list, connections, maxE=maxE)
memw["w1"].set_params(Z_w1, y, alpha_w1)
memw["w2"].set_params(Z_w2, y, alpha_w2)

outfp = open(outfn, "w")
kstrl, xstrl = memw.get_channel_strings()
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


