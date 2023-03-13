#! /usr/bin/env python3

"""
RRKM program
  - calculates density of states and microscopic rate constants
    for 1D and 2D master equation calculations
"""

import os, sys, time
import numpy as np

from . import constants
from .rovib import RoVib
from .utils import findmin
from .tst import tstrates
from .tst import cvtrates


def tunnel_eck(sums, dE, E0, deltaH0, freqimg):
    """ One-dimensional correction for tunneling through Eckart potential
    [Refs: C. Eckart, Phys. Rev. 35 (1930) 1303;
           W. H. Miller, J. Am. Chem. Soc. 101 (1979) 6810.]
    arguments:
      sums: array, sum of states
      dE: energy train size
      E0, deltaH0, freqimg: potential parameters
    """
    Bv = (np.sqrt(E0) + np.sqrt(E0-deltaH0))**2
    alpha_hb = freqimg * np.sqrt(Bv / (2.*E0*(E0-deltaH0)))
    tmpabd = 2. * np.pi * np.sqrt(2.) / alpha_hb
    d = tmpabd * np.sqrt(Bv - 0.5 * (alpha_hb/2.)**2)
    
    nbin = len(sums)
    nE0 = int(float(E0)/dE)
    naboveE0 = nbin - nE0
    Ecrit = max(0., float(deltaH0))
    ncritDH0 = int(Ecrit/dE)
    naboveDH0 = nbin - ncritDH0
    ntun = nE0 - ncritDH0
    
    rhode = np.zeros(naboveE0)
    rhode[0] = sums[nE0]
    rhode[1:] = sums[nE0+1:] - sums[nE0:nbin-1]
    
    # use energy bins relative to E0 (compatible with rrkmth)
    Ep = E0 - ntun*dE + np.arange(naboveDH0+1, dtype=float) * dE
    # correct for negative values
    Ep[Ep < deltaH0] = deltaH0
    nc_add = int(float(E0-Ecrit)/dE) - ntun
    
    a = tmpabd * np.sqrt(Ep)
    b = tmpabd * np.sqrt(Ep-deltaH0)
    ex2a, ex2b = np.exp(-2.*a), np.exp(-2.*b)
    ex2ab = ex2a*ex2b
    ptran = (1. + ex2ab - ex2a - ex2b) \
            / (1. + ex2ab + np.exp(d-(a+b)) + np.exp(-d-(a+b)))
    
    for i in range(naboveDH0):
        nc = i + 1 + nc_add
        if nc > naboveE0: nc = naboveE0
        pc = ptran[i-nc+1:i+1][::-1]
        sums[ncritDH0+i] = (pc * rhode[:nc]).sum()
    
    return sums


def numrates(T, rho, kEl, Ea):
    """ thermal rate constants from rho(E) and l(E) """
    T = np.atleast_1d(T)
    nchan = len(kEl)
    kl = [np.zeros(len(T)) for ich in range(nchan)]
    for i in range(len(T)):
        ga = rho * np.exp(- Ea * constants.cm2k / T[i])
        ga /= ga.sum()
        for ich in range(nchan): kl[ich][i] = (kEl[ich]*ga).sum()
    return kl


def rrkmE(maxE, dE, rovibm, rovibcl, E0l, deltaH0l, convK=True, convJ=True,
          centcorr=False, dEint=1, klow=2e-99, rho_sums=None, Tl=None, outf=sys.stdout):
    """ 1D RRKM calculation """
    if outf is None: ofp = None
    elif hasattr(outf, "write"): ofp = outf
    else: ofp = open(outf, "w")
    nbin = int(1. + maxE/dE)
    nchan = len(rovibcl)
    
    if ofp is not None:
        ofp.write("# reactant:\n")
        rovibm.write_to(ofp, prefix="#   ")
        for ich in range(nchan):
            ofp.write("# channel-%d:\n" % (ich+1))
            rovibcl[ich].write_to(ofp, prefix="#   ", 
                                  addinfo=(("E0", E0l[ich]), ("deltaH0", deltaH0l[ich])))
        ofp.write("# RRKM calculation:\n")
        if rho_sums is not None:
            ofp.write("#   use given density and sum of states\n")
        else:
            ofp.write("#   convK, convJ = %s, %s\n" % (convK, convJ))
            if (not convJ) and centcorr: 
                ofp.write("#   centrifugal correction applied to k\n")
                ofp.write("#   rho devided by B**(dim/2)*sig\n")
        ofp.write("#   nbin, dE = %d, %g\n" % (nbin, dE))
        ofp.write("#   nchan = %d\n" % (nchan))
        if hasattr(ofp, "flush"): ofp.flush()

    if rho_sums is None:
        rho = rovibm.dens(nbin, dE, convK=convK, convJ=convJ, dEint=dEint)
    else:
        if len(rho_sums[0]) < nbin:
            raise ValueError("insufficient len(rho): given %d, required %d"
                             % (len(rho_sums[0]), nbin))
        rho = np.copy(rho_sums[0])[:nbin]
    kl = []
    for ich in range(nchan):
        rovibc, E0, deltaH0 = rovibcl[ich], E0l[ich], deltaH0l[ich]
        nE0 = int(E0/float(dE))
        naboveE0 = nbin - nE0
        if ((rovibc.freqimg is not None) and (rovibc.freqimg > 0.) and (deltaH0 is not None)
            and (E0 > deltaH0 + dE) and (nE0 > 0)):
            tunnel = True
            # add some bins to ensure convergence in tunneling calcn.
            nadd = max(int(nbin/100), 10)
        else:
            tunnel = False
            nadd = 0
        ntmps = naboveE0 + nadd
        if rho_sums is None:
            tmps = rovibc.sums(ntmps, dE, convK=convK, convJ=convJ, dEint=dEint)
        else:
            if len(rho_sums[1][ich]) < ntmps:
                raise ValueError("insufficient len(sums): given %d, required %d"
                                 % (len(rho_sums[1][ich]), ntmps))
            tmps = np.copy(rho_sums[1][ich])[:ntmps]
        sums = np.zeros(nbin+nadd)
        if nE0 >= 0: sums[nE0:] += tmps
        else: sums += tmps[-nE0:]
        if tunnel: sums = tunnel_eck(sums, dE, E0, deltaH0, rovibc.freqimg)
        
        k = np.zeros(nbin)
        for i in range(nbin):
            if rho[i] <= 0.: continue
            if (deltaH0 is not None) and (dE*i < deltaH0): continue
            k[i] = sums[i] / (constants.h/constants.cm2j) / rho[i]
        if (not convJ) and centcorr: # centrifugal correction
            if centcorr is True: pwr = 1. # default is 2D
            else: pwr = 0.5*centcorr # centcorr is dimension
            k *= (rovibm.rotB2D/rovibc.rotB2D)**pwr  # nsym is already accounted for in mod_bs
        if klow is not None:
            for i in range(nbin):
                if k[i] < klow: k[i] = 0.
        kl.append(k)

    if (not convJ) and centcorr:  # devide rho by B**(dim/2)
        if centcorr is True: pwr = 1. # default is 2D
        else: pwr = 0.5*centcorr # centcorr is dimension
        rho /= rovibm.rotB2D**pwr  # nsym is already accounted for in mod_bs
    
    if ofp is not None:
        # thermal rates
        if Tl is None: Tl = [300., 500., 750., 1000., 1500., 2000.]
        ktstl = tstrates(Tl, rovibm, rovibcl, E0l, deltaH0l,
                         convK=convK, convJ=(convJ or centcorr))
        knuml = numrates(Tl, rho, kl, np.arange(nbin, dtype=float)*dE)
        ofp.write("# HPL rate constants [s-1]:\n")
        tmpl = ["k%02d(num)   k%02d(tst)" %
                (ich+1, ich+1) for ich in range(nchan)]
        ofp.write("#    T[K]  %s\n" % ("   ".join(tmpl)))
        for i in range(len(Tl)):
            T = Tl[i]
            tmpl = ["%10.4e %10.4e" % 
                    (knuml[ich][i], ktstl[ich][i]) for ich in range(nchan)]
            ofp.write("#%8.1f %s\n" % (T, " ".join(tmpl)))
            
        # rho(E) and k(E)
        ofp.write("# rho(E) and k(E):\n")
        ofp.write("# E[cm-1] rho[cm] %s\n" % 
                  (" ".join("    k%02d[s-1]" % (ich+1) for ich in range(nchan))))
        for i in range(nbin):
            ofp.write(" %6d %12.6e %s\n" % 
                      (i*dE, rho[i], " ".join("%12.6e" % k[i] for k in kl)))
        if hasattr(ofp, "flush"): ofp.flush()
    return rho, kl


def rrkmEJ(maxE, dE, maxJ, rovibm, rovibcl, E0l, deltaH0l, rotB2Dprodl,
           dEint=1, Tl=None, outf=sys.stdout, verbose=False):
    """ 2D RRKM calculation: rho(E,J) and k(E,J) """
    t0 = time.time()
    if outf is None: ofp = None
    elif hasattr(outf, "write"): ofp = outf
    else: ofp = open(outf, "w")
    nbin = int(1. + maxE/dE)
    nchan = len(rovibcl)
    
    if ofp is not None:
        ofp.write("# reactant:\n")
        rovibm.write_to(ofp, prefix="#   ")
        for ich in range(nchan):
            ofp.write("# channel-%d:\n" % (ich+1))
            rovibcl[ich].write_to(ofp, prefix="#   ",
                                  addinfo=(("E0", E0l[ich]), ("deltaH0", deltaH0l[ich]),
                                           ("rotB2Dprod", rotB2Dprodl[ich])))
        ofp.write("# RRKM(E,J) calculation:\n")
        ofp.write("#   nbin = %d\n" % (nbin))
        ofp.write("#   dE = %d\n" % (dE))
        ofp.write("#   maxJ = %d\n" % (maxJ))
        ofp.write("#   nchan = %d\n" % (nchan))
        ofp.write("#   B2D = %g\n" % (rovibm.rotB2D))
        ofp.write("# E[cm-1] rho[cm] %s\n" % 
                  (" ".join("    k%02d[s-1]" % (ich+1) for ich in range(nchan))))
        if hasattr(ofp, "flush"): ofp.flush()

    rho = rovibm.dens(nbin, dE, convK=True, convJ=False, dEint=dEint)
    sumsl = []
    nsums = nbin + max(int(nbin/100), 10)
    for rovibc in rovibcl:
        sumsl.append(rovibc.sums(nsums, dE, convK=True, convJ=False, dEint=dEint))
    rho_sums = (rho, sumsl)
    
    J = 0
    totcount = 0
    offsetl, rhol, kll = [], [], [[] for ich in range(nchan)]
    while J <= maxJ:
        Erot = rovibm.rotB2D * J * (J + 1.)
        offset = int(Erot/dE)
        if offset >= nbin: break
        
        E0Jl, deltaH0Jl = [], []
        for ich in range(nchan):
            Erotts = rovibcl[ich].rotB2D * J * (J + 1.)
            if rotB2Dprodl[ich] is None: Erotprod = 0.
            else: Erotprod = rotB2Dprodl[ich] * J * (J + 1.)
            E0J = E0l[ich] + (Erotts - Erot)
            if E0J < 0.: E0J = 0.
            E0Jl.append(E0J)
            deltaH0Jl.append(deltaH0l[ich] + (Erotprod - Erot))
        
        rho, kl = rrkmE(maxE, dE, rovibm, rovibcl, E0Jl, deltaH0Jl,
                        rho_sums=rho_sums, outf=None)
        rho = rho * (2.*J + 1.)  # 2J+1 factor is applied to rho
        offsetl.append(offset)
        rhol.append(rho)
        for ich in range(nchan): kll[ich].append(kl[ich])
        
        if ofp is not None:
            ofp.write("# J = %d\n" % J)
            for ig in range(nbin-offset):
                ofp.write(" %6d %12.6e %s\n" % 
                          (ig*dE, rho[ig], " ".join("%12.6e" % k[ig] for k in kl)))
            ofp.write("\n\n")
            if hasattr(ofp, "flush"): ofp.flush()
        totcount += nbin-offset
        if verbose and (J % 20) == 0:
            print("# rrkmEJ: J = %d, %.1f s elapsed" % (J, time.time()-t0))
            sys.stdout.flush()
        J += 1
    
    if ofp is not None:
        ofp.write("# total count = %d\n" % totcount)

        # thermal rates
        if Tl is None: Tl = [300., 500., 750., 1000., 1500., 2000.]
        ktstl = tstrates(Tl, rovibm, rovibcl, E0l, deltaH0l)
        El = []
        for i in range(len(offsetl)):
            El.append((np.arange(nbin, dtype=float)+offsetl[i])*dE)
        Ea = np.ravel(El)
        rhoa = np.ravel(rhol)
        kal = [np.ravel(kl) for kl in kll]
        knuml = numrates(Tl, rhoa, kal, Ea)
        ofp.write("# HPL rate constants [s-1]:\n")
        tmpl = ["k%02d(num)   k%02d(tst)" %
                (ich+1, ich+1) for ich in range(nchan)]
        ofp.write("#    T[K]  %s\n" % ("   ".join(tmpl)))
        for i in range(len(Tl)):
            T = Tl[i]
            tmpl = ["%10.4e %10.4e" % 
                    (knuml[ich][i], ktstl[ich][i]) for ich in range(nchan)]
            ofp.write("#%8.1f %s\n" % (T, " ".join(tmpl)))
        if hasattr(ofp, "flush"): ofp.flush()
    if verbose:
        print("# rrkmEJ: done, %.1f s elapsed" % (time.time()-t0))
        sys.stdout.flush()

    return offsetl, rhol, kll



def vrrkmE(maxE, dE, rovibm, rovibcl, E0l, deltaH0l, rcoordl, convK=True, convJ=True,
           centcorr=False, dEint=1, klow=2e-99, rho_sums=None,
           Tl=None, outf=sys.stdout, full_output=False):
    """ 1D variational RRKM calculation """
    if outf is None: ofp = None
    elif hasattr(outf, "write"): ofp = outf
    else: ofp = open(outf, "w")
    nbin = int(1. + maxE/dE)
    nchan = 1
    if ofp is not None:
        ofp.write("# reactant:\n")
        rovibm.write_to(ofp, prefix="#   ")
        ofp.write("# variational RRKM calculation (%d points):\n" % len(rovibcl))
        if rho_sums is not None:
            ofp.write("#   use given density and sum of states\n")
        else:
            ofp.write("#   convK, convJ = %s, %s\n" % (convK, convJ))
            if (not convJ) and centcorr: 
                ofp.write("#   centrifugal correction applied to k\n")
                ofp.write("#   rho devided by B**(dim/2)*sig\n")
        ofp.write("#   nbin, dE = %d, %g\n" % (nbin, dE))
        ofp.write("#   nchan = %d\n" % (nchan))
    
    rho, kl = rrkmE(maxE, dE, rovibm, rovibcl, E0l, deltaH0l, convK=convK, convJ=convJ, 
                    centcorr=centcorr, dEint=dEint, klow=klow, rho_sums=rho_sums,
                    Tl=Tl, outf=None)
    rv = np.zeros(len(rho))
    kv = np.zeros(len(rho))
    for i in range(len(rho)):
        kvl = [k[i] for k in kl]
        if min(kvl) < klow:
            rv[i] = rcoordl[kvl.index(min(kvl))]
            continue
        rmin, kmin = findmin(rcoordl, kvl)
        rv[i] = rmin
        kv[i] = kmin
    
    if ofp is not None:
        # thermal rates
        if Tl is None: Tl = [300., 500., 750., 1000., 1500., 2000.]
        ktstl, kvtst, rvtst = cvtrates(Tl, rovibm, rovibcl, E0l, deltaH0l, rcoordl,
                                       convK=convK, convJ=(convJ or centcorr))
        knuml = numrates(Tl, rho, [kv], np.arange(nbin, dtype=float)*dE)
        ofp.write("# HPL rate constants [s-1]:\n")
        if full_output:
            ofp.write("#    T[K]  k(num)     k(cvt)     r(cvt)  %s\n"
                      % (" ".join("k@r=%-6g" % x for x in rcoordl)))
            for iT in range(len(Tl)):
                ofp.write("#%8.1f %10.4e %10.4e %-8s %s\n"
                          % (Tl[iT], knuml[0][iT], kvtst[iT], "%-6g" % rvtst[iT],
                             " ".join("%10.4e" % ktstl[iv][iT] for iv in range(len(ktstl)))))
        else:
            ofp.write("#    T[K]  k(num)     k(cvt)     r(cvt)\n")
            for iT in range(len(Tl)):
                ofp.write("#%8.1f %10.4e %10.4e %g\n"
                          % (Tl[iT], knuml[0][iT], kvtst[iT], rvtst[iT]))
            
        # rho(E) and k(E)
        ofp.write("# rho(E) and k(E):\n")
        if full_output:
            ofp.write("# E[cm-1] rho[cm]     k[s-1]      r       %s\n"
                      % (" ".join("  k@r=%-6g" % x for x in rcoordl)))
            for i in range(nbin):
                ofp.write(" %6d %12.6e %12.6e %-8s %s\n"
                          % (i*dE, rho[i], kv[i], "%-6g" % rv[i],
                             " ".join("%12.6e" % kl[iv][i] for iv in range(len(kl)))))
        else:
            ofp.write("# E[cm-1] rho[cm]     k[s-1]      r\n")
            for i in range(nbin):
                ofp.write(" %6d %12.6e %12.6e %g\n" % (i*dE, rho[i], kv[i], rv[i]))
        
    return rho, kl, kv, rv


def vrrkmEJ(maxE, dE, maxJ, rovibm, rovibcl, E0l, deltaH0l, rotB2Dprodl, rcoordl,
           dEint=1, Tl=None, outf=sys.stdout, full_output=False, verbose=False):
    """ 2D variational RRKM calculation """
    t0 = time.time()
    if outf is None: ofp = None
    elif hasattr(outf, "write"): ofp = outf
    else: ofp = open(outf, "w")
    nbin = int(1. + maxE/dE)
    nchan = 1
    
    if ofp is not None:
        ofp.write("# reactant:\n")
        rovibm.write_to(ofp, prefix="#   ")
        ofp.write("# variational RRKM(E,J) calculation (%d points):\n" % len(rovibcl))
        ofp.write("#   nbin = %d\n" % (nbin))
        ofp.write("#   dE = %d\n" % (dE))
        ofp.write("#   maxJ = %d\n" % (maxJ))
        ofp.write("#   nchan = %d\n" % (nchan))
        ofp.write("#   B2D = %g\n" % (rovibm.rotB2D))
        if full_output:
            ofp.write("# E[cm-1] rho[cm]     k[s-1]      r       %s\n"
                      % (" ".join("  k@r=%-6g" % x for x in rcoordl)))
        else:
            ofp.write("# E[cm-1] rho[cm]     k[s-1]      r\n")
    
    rho = rovibm.dens(nbin, dE, convK=True, convJ=False, dEint=dEint)
    sumsl = []
    nsums = nbin + max(int(nbin/100), 10)
    for rovibc in rovibcl:
        sumsl.append(rovibc.sums(nsums, dE, convK=True, convJ=False, dEint=dEint))
    rho_sums = (rho, sumsl)
    
    J = 0
    totcount = 0
    offsetl, rhol, kvl, rvl = [], [], [], []
    while J <= maxJ:
        Erot = rovibm.rotB2D * J * (J + 1.)
        offset = int(Erot/dE)
        if offset >= nbin: break
        
        E0Jl, deltaH0Jl = [], []
        for iv in range(len(rovibcl)):
            Erotts = rovibcl[iv].rotB2D * J * (J + 1.)
            if rotB2Dprodl[iv] is None: Erotprod = 0.
            else: Erotprod = rotB2Dprodl[iv] * J * (J + 1.)
            E0J = E0l[iv] + (Erotts - Erot)
            if E0J < 0.: E0J = 0.
            E0Jl.append(E0J)
            deltaH0Jl.append(deltaH0l[iv] + (Erotprod - Erot))
        
        rho, kl, kv, rv = vrrkmE(maxE, dE, rovibm, rovibcl, E0Jl, deltaH0Jl, rcoordl,
                                 rho_sums=rho_sums, outf=None)
        rho = rho * (2.*J + 1.)  # 2J+1 factor is applied to rho
        offsetl.append(offset)
        rhol.append(rho)
        kvl.append(kv)
        rvl.append(rv)
        
        if ofp is not None:
            ofp.write("# J = %d\n" % J)
            if full_output:
                for ig in range(nbin-offset):
                    ofp.write(" %6d %12.6e %12.6e %-8s %s\n" 
                              % (ig*dE, rho[ig], kv[ig], "%-6g" % rv[ig],
                                 " ".join("%12.6e" % kl[iv][ig] for iv in range(len(kl)))))
            else:
                for ig in range(nbin-offset):
                    ofp.write(" %6d %12.6e %12.6e %g\n" % (ig*dE, rho[ig], kv[ig], rv[ig]))
            ofp.write("\n\n")
            if hasattr(ofp, "flush"): ofp.flush()
        totcount += nbin-offset
        if verbose and (J % 20) == 0:
            print("# vrrkmEJ: J = %d, %.1f s elapsed" % (J, time.time()-t0))
            sys.stdout.flush()
        J += 1
    
    if ofp is not None:
        ofp.write("# total count = %d\n" % totcount)

        # thermal rates
        if Tl is None: Tl = [300., 500., 750., 1000., 1500., 2000.]
        ktstl, kvtst, rvtst = cvtrates(Tl, rovibm, rovibcl, E0l, deltaH0l, rcoordl)
        El = []
        for i in range(len(offsetl)):
            El.append((np.arange(nbin, dtype=float)+offsetl[i])*dE)
        Ea = np.ravel(El)
        rhoa = np.ravel(rhol)
        kva = np.ravel(kvl)
        knuml = numrates(Tl, rhoa, [kva], Ea)
        ofp.write("# HPL rate constants [s-1]:\n")
        if full_output:
            ofp.write("#    T[K]  k(num)     k(cvt)     r(cvt)  %s\n"
                      % (" ".join("k@r=%-6g" % x for x in rcoordl)))
            for iT in range(len(Tl)):
                ofp.write("#%8.1f %10.4e %10.4e %-8s %s\n"
                          % (Tl[iT], knuml[0][iT], kvtst[iT], "%-6g" % rvtst[iT],
                             " ".join("%10.4e" % ktstl[iv][iT] for iv in range(len(ktstl)))))
        else:
            ofp.write("#    T[K]  k(num)     k(cvt)    r(cvt)\n")
            for iT in range(len(Tl)):
                ofp.write("#%8.1f %10.4e %10.4e %g\n" % (Tl[iT], knuml[0][iT], kvtst[iT], rvtst[iT]))
        if hasattr(ofp, "flush"): ofp.flush()
    if verbose:
        print("# vrrkmEJ: done, %.1f s elapsed" % (time.time()-t0))
        sys.stdout.flush()

    return offsetl, rhol, kvl, rvl


def ilt(maxE, dE, rovibm, ilt_A, ilt_E, convK=True, convJ=True,
        dEint=1, outf=sys.stdout):
    """ ILT: k(T) = A exp(-Ea/RT) => k(E) = A rho(E-Ea) / rho(E) """
    if outf is None: ofp = None
    elif hasattr(outf, "write"): ofp = outf
    else: ofp = open(outf, "w")
    nbin = int(1. + maxE/dE)
    nchan = 1
    
    if ofp is not None:
        ofp.write("# reactant:\n")
        rovibm.write_to(ofp, prefix="#   ")
        ofp.write("# ILT calculation:\n")
        ofp.write("#   convK, convJ = %s, %s\n" % (convK, convJ))
        ofp.write("#   A = %g\n" % (ilt_A))
        ofp.write("#   E = %g\n" % (ilt_E))
        ofp.write("#   nchan = %d\n" % (nchan))
    
    rho = rovibm.dens(nbin, dE, convK=convK, convJ=convJ, dEint=dEint)
    
    iE = int(ilt_E / dE)
    k = np.zeros(nbin)
    for i in range(nbin):
        if rho[i] <= 0.: continue
        if i >= iE: k[i] = ilt_A * rho[i-iE] / rho[i]
        else: k[i] = 0.
    if ofp is not None:
        kl = [k]
        ofp.write("# E[cm-1] rho[cm] %s\n" %
                  (" ".join("    k%02d[s-1]" % (ich+1) for ich in range(nchan))))
        for i in range(nbin):
            ofp.write(" %6d %12.6e %s\n" %
                      (i*dE, rho[i], " ".join("%12.6e" % k[i] for k in kl)))
    return rho, k



def rrkmth(infn=None, outfn=None, masfn=None):
    """ Compatible with the RRKMTH (rev.0.15) program of SSUMES software """
    if infn is None: fp = sys.stdin
    else: fp = open(infn)
    title = fp.readline().rstrip()
    iterx = iter(fp.read().strip().split())
    def nextf(): return float(next(iterx))
    def nexti(): return int(next(iterx))
    
    nn, inc = nexti(), nexti()
    npress, ntemp = nexti(), nexti()
    nchan, jav = nexti(), nexti()
    if jav != 0: raise ValueError("JAV should be 0.")
    rch = range(nchan)
    rovibm = RoVib()
    rovibcl = [RoVib() for ich in rch]
    
    jf, nf = [nexti() for ich in rch], nexti()
    eok = [nextf() for ich in rch]
    src, srm = [nextf() for ich in rch], nextf()
    for ich in rch: rovibcl[ich].nsym = src[ich]
    rovibm.nsym = srm
    
    bcmplx, bmolec = [nextf() for ich in rch], nextf()
    iexmd, iexrtd = 3, [3 for ich in rch]
    if bmolec < 0.: iexmd = 2
    for ich in rch:
        if bcmplx[ich] < 0: iexrtd[ich] = 2
    for ich in rch: rovibcl[ich].rotB2D = abs(bcmplx[ich])
    rovibm.rotB2D = abs(bmolec)
    
    n, nintr = [nexti() for ich in rch], nexti()
    for ich in range(nchan+1):
        if ich == nchan: nthis, rvthis = nintr, rovibm
        else: nthis, rvthis = n[ich], rovibcl[ich]
        introt = []
        for j in range(nthis):
            b, sig, irtd = nextf(), nextf(), nexti()
            hohnd, v0 = nextf(), nextf()
            if (irtd != 1) and (hohnd >= 0.):
                raise ValueError("hind rot must be 1d.")
            introt.append((b, sig, hohnd, v0))
            if hohnd <= 0. and irtd > 1: # >1d free rotor
                for ir in range(irtd-1): introt.append((b, 1, hohnd, v0))
            rvthis.introt = introt
    
    sgma, wt1, wt2, eps = nextf(), nextf(), nextf(), nextf()
    for ich in range(nchan+1):
        if ich == nchan: nthis, rvthis = nf, rovibm
        else: nthis, rvthis = jf[ich], rovibcl[ich]
        freq = []
        for i in range(nthis):
            nc, jc = nexti(), nexti()
            for j in range(jc): freq.append(nc)
        rvthis.freq = freq
    
    press = [nextf() for i in range(npress)]
    temp = [nextf() for i in range(ntemp)]
    itptun = [0 for ich in rch]
    deltah = [None for ich in rch]
    rdmirc = [None for ich in rch]
    try:
        itptun = [nexti() for ich in rch]
        for ich in rch:
            if itptun != 0:
                deltah[ich], rovibcl[ich].freqimg = nextf(), nextf()
                rdmirc[ich] = nextf()
    except StopIteration: pass
    try:
        ncnfc, ncnfm = [nexti() for ich in rch], nexti()
        for ich in range(nchan+1):
            if ich == nchan: nthis, rvthis = ncnfm, rovibm
            else: nthis, rvthis = ncnfc[ich], rovibcl[ich]
            rvthis.states = [(nexti(), nextf()) for i in range(nthis)]
    except StopIteration: pass
    if infn is not None: fp.close()
    
    wka = 349.7751 # kcal/mol => cm^-1
    E0l = [x*wka for x in eok]
    deltaH0l = [None for ich in rch]
    for ich in rch:
        if deltah[ich] is not None: deltaH0l[ich] = deltah[ich] * wka
    dE = inc
    Ethresh = min(E0l)
    for x in deltaH0l:
        if x is not None and 0. <= x < Ethresh: Ethresh = x
    nreact = int(Ethresh/float(dE))
    nbin = nreact + nn
    maxE = (nbin - 1) * dE
    
    # output
    if outfn is None: ofp = sys.stdout
    else: ofp = open(outfn, "w")
    
    dEint = 10 # compatible with rrkmth
    # rrkm calculation (centrifugal correction, scaled rho)
    rho, kl = rrkmE(maxE, dE, rovibm, rovibcl, E0l, deltaH0l, convJ=False,
                    centcorr=iexmd, dEint=dEint, Tl=temp, outf=ofp)
    
    # output master file (only for ssumes)
    if masfn is None: mfp = sys.stdout
    else: mfp = open(masfn, "w")
    def write7(w7fp, arr):
        div, mod = divmod(len(arr),7)
        for i in range(div+1):
            if i == div and mod == 0: break
            arr1 = arr[i*7 : min((i+1)*7, len(arr))]
            w7fp.write("  %s\n" % (" ".join("%10.4E" % x for x in arr1)))
    mfp.write(" %s\n %5d        %5d %4d\n" % (title, 400, nchan, inc))
    write7(mfp, (1e-6,1e-6,1e-6))
    mfp.write(" %12.3f\n -2\n  400.\n  0 1\n" % ((nreact+0.5)*inc/wka))
    mfp.write(" %4d\n" % (npress))
    write7(mfp, press)
    mfp.write("  0.0\n")
    mfp.write(" %6d\n" % (jav))
    mfp.write(" %4d\n %s\n" % (ntemp, "  ".join("%8.2f" % x for x in temp)))
    mfp.write("  %11.3f %11.3f %11.3f %11.3f\n" % (sgma, wt1, wt2, eps))
    mfp.write(" %2d %1d\n" % (1, max(nchan, 2)))
    mfp.write(" %6d\n" % (nbin))
    write7(mfp, rho)
    for itemp in range(ntemp):
        mfp.write(" %9.3f  %8.3f\n" % (temp[itemp], 1.)) # CORRAT is dummy
    mfp.write(" %6d\n" % (nbin-nreact))
    for ich in rch: write7(mfp, kl[ich][nreact:])
    if nchan == 1: write7(mfp, np.zeros(nbin-nreact))
    if masfn is not None: mfp.close()
    
    if outfn is not None: ofp.close()
    return



