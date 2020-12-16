#! /usr/bin/env python3

"""
collision frequency
"""

import numpy as np
from scipy import interpolate

from . import constants

def name2weight(name):
    """ returns moleculer weights [amu] of given molecular formula """
    name = name.strip()
    n = len(name)
    l = []
    s = ""
    state = 0 # 0  1 
    for i in range(n):
        if len(s) == 0:
            s += name[i]
            continue
        if name[i].isdigit():
            if s[-1].isdigit(): s += name[i]
            else: l.append(s); s = name[i]
        else:
            if not s[-1].isdigit(): s += name[i]
            else: l.append(s); s = name[i]
    if len(s) > 0: l.append(s)
    
    l2 = []
    for s in l:
        if s.isdigit():
            for i in range(int(s)-1):
                l2.append(l2[-1])
            continue
        while len(s) > 0:
            if len(s) == 1:
                l2.append(s)
                s = ""
                break
            if s[:2] in constants.atom_name:
                l2.append(s[:2])
                s = s[2:]
            else:
                l2.append(s[0])
                s = s[1:]
    return sum(constants.amass[x] for x in l2)
    
    
def collfreq_hs(wt_or_name_1, wt_or_name_2, sigmaA, T):
    """ hard-sphere collision frequency [cm^3 molecule^-1 s^-1]
    arguments:
      wt_or_name_1, wt_or_name_2: molecular weights [amu] or formula
      sigmaA: collision diameter [A]
      T: temperature [K]
    """
    if isinstance(wt_or_name_1, str): wt_or_name_1 = name2weight(wt_or_name_1)
    if isinstance(wt_or_name_2, str): wt_or_name_2 = name2weight(wt_or_name_2)
    redmass = wt_or_name_1 * wt_or_name_2 / (wt_or_name_1 + wt_or_name_2)
    cross_section = np.pi * sigmaA * sigmaA * 1.e-16
    mean_vel = np.sqrt(8*constants.kb*T / (np.pi*redmass*constants.amu))*100
    return cross_section * mean_vel


def collfreq_capture(wt_or_name_1, wt_or_name_2, sigmaA, epsK, T):
    """ dispersion capture rate constant [cm^3 molecule^-1 s^-1]
    [Ref: e.g., A. Matsugi, J. Phys. Chem. A 122 (2018) 1972]
    arguments:
      wt_or_name_1, wt_or_name_2: molecular weights [amu] or formula
      sigmaA: L-J diameter parameter [A]
      epsK: L-J well depth parameter [K]
      T: temperature [K]
    """
    g23 = 1.35411793943  # gamma(2./3.)
    return 2.*g23 * (T/epsK)**(-1./3.) * \
           collfreq_hs(wt_or_name_1, wt_or_name_2, sigmaA, T)


def collfreq_lj_approx(wt_or_name_1, wt_or_name_2, sigmaA, epsK, T):
    """ Lennard-Jones collision frequency [cm^3 molecule^-1 s^-1]
    Approximate expression of [Ref: J. Troe, J. Chem. Phys. 66 (1977) 4758].
    arguments:
      wt_or_name_1, wt_or_name_2: molecular weights [amu] or formula
      sigmaA: L-J diameter parameter [A]
      epsK: L-J well depth parameter [K]
      T: temperature [K]
    """
    omega22 = 1. / (0.636 + 0.567 * np.log10(T / epsK))
    return omega22 * collfreq_hs(wt_or_name_1, wt_or_name_2, sigmaA, T)


def collfreq_lj(wt_or_name_1, wt_or_name_2, sigmaA, epsK, T):
    """ Lennard-Jones collision frequency [cm^3 molecule^-1 s^-1]
    Collision integral interpolated using the tablulated values.
    [Ref: M. Klein and F.J. Smith, J Res Natl Bur Stand A Phys Chem.
          1968 Jul-Aug; 72A(4): 359; doi: 10.6028/jres.072A.033]
    arguments:
      wt_or_name_1, wt_or_name_2: molecular weights [amu] or formula
      sigmaA: L-J diameter parameter [A]
      epsK: L-J well depth parameter [K]
      T: temperature [K]
    """
    global _omega22_table_interp
    if _omega22_table_interp is None:
        _omega22_table_interp = interpolate.interp1d(_omega22_table_Ts, _omega22_table)
    omega22 = _omega22_table_interp(T / epsK)
    return omega22 * collfreq_hs(wt_or_name_1, wt_or_name_2, sigmaA, T)

_omega22_table_interp = None

# collision integral table
#  [Ref: M. Klein and F.J. Smith, J Res Natl Bur Stand A Phys Chem.
#        1968 Jul-Aug; 72A(4): 359; doi: 10.6028/jres.072A.033]
_omega22_table_Ts = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                     0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                     1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                     2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5,
                     5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
                     11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                     22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0,
                     45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0,
                     95.0, 100.0, 125.0, 150.0, 175.0, 200.0]
_omega22_table = [4.1039, 3.5879, 3.2686, 3.0336, 2.8456, 2.6777, 2.5327,
                  2.4003, 2.2867, 2.1774, 2.0828, 1.9994, 1.9230, 1.8521,
                  1.7885, 1.7327, 1.6821, 1.6360, 1.5933, 1.5175, 1.4542,
                  1.4010, 1.3550, 1.3151, 1.2801, 1.2491, 1.2215, 1.1971,
                  1.1753, 1.1379, 1.1069, 1.0807, 1.0581, 1.0386, 1.0215,
                  1.0064, 0.9929, 0.9807, 0.9697, 0.9461, 0.9266, 0.9102,
                  0.8961, 0.8837, 0.8726, 0.8627, 0.8537, 0.8454, 0.8378,
                  0.8308, 0.8242, 0.8123, 0.8017, 0.7921, 0.7834, 0.7754,
                  0.7681, 0.7613, 0.7549, 0.7489, 0.7433, 0.7330, 0.7238,
                  0.7153, 0.7076, 0.7005, 0.6939, 0.6878, 0.6820, 0.6766,
                  0.6715, 0.6599, 0.6496, 0.6405, 0.6322, 0.6246, 0.6177,
                  0.6113, 0.6053, 0.5998, 0.5946, 0.5897, 0.5851, 0.5654,
                  0.5497, 0.5366, 0.5256]

