#! /usr/bin/env python3

"""
Thermal decomposition in the n-,s-C5H11 system
(two-well, two-channel, 1D ME as a function of E)

Phenomenological rate constants from steady-state MEs
(see Supporting Information of [A. Matsugi, J. Phys. Chem. A 125 (2021) 2532])

sample output (c5h11_4a_me1d_E_phenomenological.dat):
   T[K]    p[bar]       w1->w1-ch2     w2->w1-ch2     w1->w2-ch2     w2->w2-ch2         w2->w1         w1->w2
  1000.0  1.000e+02     6.7929e+06     1.0100e+04     1.8738e+05     1.0491e+07     6.5773e+05     4.8258e+06
  1000.0  1.000e+01     5.1834e+06     3.5558e+04     6.8599e+05     7.7342e+06     5.0427e+05     3.7001e+06
  1000.0  1.000e+00     2.3900e+06     4.4524e+04     9.5277e+05     3.6551e+06     2.6458e+05     1.9470e+06
  1000.0  1.000e-01     6.0012e+05     2.3930e+04     6.0857e+05     1.1144e+06     9.8458e+04     7.3686e+05
"""

from me2d import ME1DMW


maxE = 60000  # cm^-1

# list of (name, rrkm_filename, relative_energy)
well_list = [("w1", "c5h11_1_rrkmE_nc5h11_dE10.dat", 948.38),
             ("w2", "c5h11_1_rrkmE_sc5h11_dE10.dat", 0.)]
# list of ((name, ch), (name, ch))
connections = [(("w1", 1), ("w2", 1))]

outfn = "c5h11_4a_me1d_E_phenomenological.dat"

solver = "InvIter,cg" # inverse iteration with conjugate gradient
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
kdstrl, kwstrl = memw.get_channel_strings_phnm()
kstrl = []
for x in kdstrl+kwstrl: kstrl.extend(x)
outfp.write("   T[K]    p[bar]   %s\n" % (" ".join("%14s" % x for x in kstrl)))

gal = [None for i in range(memw.nwell)]
for p in pl:
    kll, popll = [], []
    for iwell in range(memw.nwell):
        reactant = memw.names[iwell]
        ktot, kchl, ga, popl, vals, vec = \
              memw.solve(T, p, gguess=gal[iwell], solver=solver, reactant=reactant,
                         neig=neig, verbose=verbose, bandpcrit=bandpcrit,
                         nthreads=nthreads, maxmemGB=maxmemGB)
        kll.append(kchl)
        popll.append(popl)
        gal[iwell] = ga
    kdl, kwl = memw.kphnm_from_ss(kll, popll)
    kl = []
    for x in kdl+kwl: kl.extend(x)
    outfp.write("%8.1f  %.3e %s\n" % (T, p, " ".join("%14.4e" % x for x in kl)))
    outfp.flush()

