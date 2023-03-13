#! /usr/bin/env python3

"""
TST rate constants
"""

import os, sys, time
import numpy as np

from . import constants
from .utils import findmin
from .utils import name2weight


def tunnel_eck_T(T, E0, deltaH0, freqimg, fE=50., fT=20.):
    """ One-dimensional correction for tunneling through Eckart potential
    [Refs: C. Eckart, Phys. Rev. 35 (1930) 1303;
           B. C. Garrett and D. G. Truhlar, J. Phys. Chem. 83 (1979) 2921.]
    arguments:
      T: temperature(s)
      E0, deltaH0, freqimg: potential parameters
      fE: integration grain size = kB*min(T) / fE
      fT: integration upper limit = E0 + kB*max(T)*fT
    """
    T = np.atleast_1d(T)
    Bv = (np.sqrt(E0) + np.sqrt(E0-deltaH0))**2
    alpha_hb = freqimg * np.sqrt(Bv / (2.*E0*(E0-deltaH0)))
    tmpabd = 2. * np.pi * np.sqrt(2.) / alpha_hb
    d = tmpabd * np.sqrt(Bv - 0.5 * (alpha_hb/2.)**2)

    dEint = min(T) / constants.cm2k / fE
    Emin = max(0., deltaH0)
    Emax = E0 + max(T) / constants.cm2k * fT
    Ec = np.arange(Emin, Emax, dEint, dtype=float)
    a = tmpabd * np.sqrt(Ec)
    b = tmpabd * np.sqrt(Ec-deltaH0)
    ex2a, ex2b = np.exp(-2.*a), np.exp(-2.*b)
    ex2ab = ex2a*ex2b
    ptran = (1. + ex2ab - ex2a - ex2b) \
            / (1. + ex2ab + np.exp(d-(a+b)) + np.exp(-d-(a+b)))
    
    corr = np.ones(len(T))
    for i in range(len(T)):
        intg = (ptran * np.exp(-Ec * constants.cm2k / T[i])).sum() * dEint
        corr[i] = intg / (T[i] / constants.cm2k)  / np.exp(-E0 * constants.cm2k / T[i])
    return corr



def tstrates(T, rovibm, rovibcl, E0l, deltaH0l, convK=True, convJ=True):
    """ thermal rate constants for unimolecular reactant from transition state theory """
    T = np.atleast_1d(T)
    nchan = len(rovibcl)
    qm = rovibm.part(T, convK=convK, convJ=convJ)
    kl = []
    for ich in range(nchan):
        E0, deltaH0, rovibc = E0l[ich], deltaH0l[ich], rovibcl[ich]
        qc = rovibc.part(T, convK=convK, convJ=convJ)
        k = (constants.kb * T / constants.h) * (qc / qm) * \
            np.exp(- E0 * constants.cm2k / T)
        if ((rovibc.freqimg is not None) and (rovibc.freqimg > 0.) and (deltaH0 is not None)
            and (E0 > deltaH0) and (E0 > 0.)):
            k *= tunnel_eck_T(T, E0, deltaH0, rovibc.freqimg)
        kl.append(k)
    return kl


def cvtrates(T, rovibm, rovibcl, E0l, deltaH0l, rcoordl, convK=True, convJ=True):
    """ thermal rate constants for unimolecular reactant from canonical variational transition state theory """
    T = np.atleast_1d(T)
    ktstl = tstrates(T, rovibm, rovibcl, E0l, deltaH0l, convK=convK, convJ=convJ)
    kv = np.zeros(len(T))
    rv = np.zeros(len(T))
    for iT in range(len(T)):
        ktstv = [ktst[iT] for ktst in ktstl]
        rmin, kmin = findmin(rcoordl, ktstv)
        kv[iT] = kmin
        rv[iT] = rmin
    return ktstl, kv, rv


def bimol_tstrate(T, rovib_reac1, rovib_reac2, rovib_ts, gelec_reac1, gelec_reac2, gelec_ts,
                  wt_or_name_1, wt_or_name_2, E0, deltaH0):
    """ thermal rate constant for a bimolecular reaction from transition state theory"""
    T = np.atleast_1d(T)
    if isinstance(wt_or_name_1, str): wt_or_name_1 = name2weight(wt_or_name_1)
    if isinstance(wt_or_name_2, str): wt_or_name_2 = name2weight(wt_or_name_2)
    
    redmass = wt_or_name_1 * wt_or_name_2 / (wt_or_name_1 + wt_or_name_2)
    conv = (2. * np.pi * constants.amu * constants.kb / constants.h / constants.h)**1.5 * 1e-6
    qtrans_ratio = 1. / (conv * redmass**1.5 * T**1.5)
    qelec_ratio = gelec_ts / (gelec_reac1 * gelec_reac2)
    qrovib_ratio = rovib_ts.part(T) / (rovib_reac1.part(T) * rovib_reac2.part(T))
    
    k = (constants.kb * T / constants.h) * qtrans_ratio * qelec_ratio * qrovib_ratio * np.exp(- E0 * constants.cm2k / T)
    if ((rovib_ts.freqimg is not None) and (deltaH0 is not None)
        and (E0 > deltaH0) and (E0 > 0.)):
        k *= tunnel_eck_T(T, E0, deltaH0, rovib_ts.freqimg)
    return k

