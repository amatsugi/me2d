#! /usr/bin/env python3

"""
Thermal decomposition / chemical activation of toluene
(single-well, two-channel, 1D MEs as a function of E)

sample output (c6h5ch3_me1d_E_chemact.dat):
   T[K]    p[bar]         kd1[s-1]       kd2[s-1]    kca1_1[s-1]    kca1_2[s-1]  kb1_1[cm3s-1]  kb1_2[cm3s-1]    kr1[cm3s-1]    kca2_1[s-1]    kca2_2[s-1]  kb2_1[cm3s-1]  kb2_2[cm3s-1]    kr2[cm3s-1]
  2000.0    HPL         3.7576e+05     2.9843e+05
  2000.0  1.000e+02     3.4707e+05     2.6785e+05     3.7177e+05     2.9399e+05     6.1043e-12     6.4571e-12     1.5193e-10     3.8093e+05     3.0419e+05     1.8815e-12     2.0191e-12     3.4168e-11
  2000.0  1.000e+01     2.6473e+05     1.8721e+05     3.5992e+05     2.8159e+05     2.4407e-11     2.4200e-11     1.1589e-10     3.9817e+05     3.2225e+05     7.0509e-12     7.1356e-12     2.3882e-11
  2000.0  1.000e+00     1.4512e+05     8.6066e+04     3.3797e+05     2.6074e+05     5.2987e-11     4.7990e-11     6.3518e-11     4.3945e+05     3.6197e+05     1.3981e-11     1.3106e-11     1.0981e-11
  2000.0  1.000e-01     5.6225e+04     2.5542e+04     3.1017e+05     2.3668e+05     7.6388e-11     6.3512e-11     2.4596e-11     5.1981e+05     4.3405e+05     1.8502e-11     1.6303e-11     3.2633e-12
"""

from me2d import ME1D


maxE = 90000  # cm^-1
rrkmEfn = "c6h5ch3_vrrkmE_dE10.dat"
outfn = "c6h5ch3_me1d_E_chemact.dat"

solver = "InvIter,cho" # inverse iteration with Cholesky decomposition
bandpcrit = None       # use dense matrix
nthreads = 2
maxmemGB = 1
verbose = True


# for chemical activation
solver_ca = "LinEq,cho"
khpl_list = [1.644954e-10, 3.806832e-11] # CVTST rate constants

T = 2000. # K
pl = [100., 10., 1., 0.1]  # bar

y = 0.5
alpha = 97.7  # cm^-1
Z = 1.66e-09  # cm^3 molecule^-1 s^-1

me = ME1D.read_from(rrkmEfn, maxE=maxE)
me.set_params(Z, y, alpha)
outfp = open(outfn, "w")
strl = me.get_channel_strings("kd")
for ich in range(me.nchan):
    strl += me.get_channel_strings("kca%d_" % (ich+1), "s-1")
    strl += me.get_channel_strings("kb%d_" % (ich+1), "cm3s-1")
    strl += ["kr%d[cm3s-1]" % (ich+1)]
outfp.write("   T[K]    p[bar]   %s\n" % (" ".join("%14s" % x for x in strl )))

ktot, kchl, ga = me.hpl(T)
outfp.write("%8.1f    HPL     %s\n" % (T, " ".join("%14.4e" % x for x in kchl)))
outfp.flush()

for p in pl:
    ktot, kchl, ga, vals, vec = me.solve(T, p, gguess=ga, solver=solver,
                                         verbose=verbose, bandpcrit=bandpcrit,
                                         nthreads=nthreads, maxmemGB=maxmemGB)
    kdl = kchl
    koutl = list(kdl)
    for ich in range(me.nchan):
        ktot, kchl, ga, vals, vec = me.solve(T, p, gguess=None, solver=solver_ca,
                                             chemact_ch=ich+1,
                                             verbose=verbose, bandpcrit=bandpcrit,
                                             nthreads=nthreads, maxmemGB=maxmemGB)
        kr, kbl = me.k_chemical_activation(khpl_list[ich], kchl, kdl)
        koutl += list(kchl) + list(kbl) + [kr]

    outfp.write("%8.1f  %.3e %s\n" % (T, p, " ".join("%14.4e" % x for x in koutl)))
    outfp.flush()


