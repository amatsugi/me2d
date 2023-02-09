#! /usr/bin/env python3

"""
thermodynamic functions
"""

import os, sys, time
import numpy as np
from scipy.special import iv as bessel_iv
from scipy.optimize import leastsq

from . import constants
from .utils import name2weight
from .utils import name2atoms


def therm(T, rovib, gelec, wt_or_name):
    """ Calculate DF=H(T)-H0, S(T), and Cp(T),
    gelec: electronic degeneracy
    wt_or_name: molecular weight [amu] or formula
    """
    if isinstance(wt_or_name, str): mw = name2weight(wt_or_name)
    else: mw = wt_or_name
    freerot, hindrot = rovib.findrot(convK=True, convJ=True)
    freq  = [x*rovib.fscale for x in rovib.freq]
    nsym = rovib.nsym
    states = rovib.states
    
    T = np.atleast_1d(T)
    DH = np.zeros(len(T))
    S = np.zeros(len(T))
    Cp = np.zeros(len(T))
    kT = T / constants.cm2k
    def bol(E):
        return np.exp(- E / kT)
    
    for x in freerot:
        B, sig, dim = x
        DH += dim/2.
        Cp += dim/2.
        q = (kT/B)**(0.5*dim) / float(sig)
        if dim % 2 == 1: q *= np.sqrt(np.pi)
        S += dim/2. + np.log(q)
    
    for x in hindrot:
        B, sig, hofreq, V0 = x
        bolf = bol(hofreq)
        # H.O
        H_ho = (hofreq/kT) * bolf / (1. - bolf)
        S_ho = H_ho - np.log(1. - bolf)
        Cp_ho = H_ho**2 / bolf
        # classical H.O.
        H_hoc = 1.
        S_hoc = 1. - np.log(hofreq/kT)
        Cp_hoc = 1.
        # classical H.R.
        V0half = V0/2./ kT
        b0 = bessel_iv(0, V0half)
        b1 = bessel_iv(1, V0half)
        b2 = bessel_iv(2, V0half)
        H_hrc = 0.5 + V0half * (1. - b1 / b0)
        q_fr = (np.pi * kT / B)**0.5 / float(sig)
        q_hrc = q_fr * bol(V0/2.) * b0
        S_hrc = H_hrc + np.log(q_hrc)
        Cp_hrc = 0.5 + V0half*b1/b0 - (V0half*b1/b0)**2 + V0half**2*b2/b0
        # P.G. approx
        DH += H_hrc + H_ho - H_hoc
        S += S_hrc + S_ho - S_hoc
        Cp += Cp_hrc + Cp_ho - Cp_hoc
    
    if len(states) > 0:
        sum0, sum1, sum2 = 0., 0., 0.
        for x in states:
            degen, level = x
            sum0 += degen * bol(level)
            sum1 += degen * (level/kT) * bol(level)
            sum2 += degen * (level/kT)**2 * bol(level)
        DH += sum1 / sum0
        S += np.log(sum0) + sum1 / sum0
        Cp += sum2 / sum0 - (sum1 / sum0)**2
    
    if gelec is not None:
        S += np.log(gelec)
    
    for hofreq in freq:
        bolf = bol(hofreq)
        H_ho = (hofreq/kT) * bolf / (1. - bolf)
        DH += H_ho
        S += H_ho - np.log(1. - bolf)
        Cp += H_ho**2 / bolf

    # trans (1bar)
    DH += 2.5  # 3/2 + (H - U)
    Cp += 2.5  # 3/2 + (Cp - Cv)
    S += 1.5 * np.log(2. * np.pi * constants.amu*mw / constants.h**2) \
         + 2.5 * np.log(constants.kb*T) - np.log(1e5) + 2.5
    
    # symmetry number
    S -= np.log(nsym)

    # units (kJ/mol or J/K/mol
    DH *= constants.r * T / 1000.
    S *= constants.r
    Cp *= constants.r
    
    return DH, S, Cp


def therm_nasa7fit(T, Cp, H298, S298, midT=1000., opt_midT=True, woffT=50., verbose=False):
    """ Calculate NASA-7 polynomial coefficients. """
    T = np.atleast_1d(T)
    Cp = np.atleast_1d(Cp)

    midCp_guess = Cp[(np.abs(T - midT)).argmin()]
    if opt_midT:
        midT_guess = midT
        guess = [midT_guess, midCp_guess]
    else:
        midT_fix = midT
        guess = [midCp_guess]
    
    def fc(p, rescoefs=False):
        if opt_midT: midT, midCp = p
        else: midCp = p; midT = midT_fix
        w = np.exp(-abs(T-midT)/woffT)
        wlow = np.where(T <= midT, 1., w)
        whigh = np.where(T >= midT, 1., w)

        # linear least sq. fit using svd to derive coefficients with given Cp and/at midT
        obj = (Cp - midCp) / constants.r
        a = np.asarray([T**i - midT**i for i in range(1,5)])
        cutoff = len(T) * np.finfo(T.dtype).eps
        
        u, s, v = np.linalg.svd((a*wlow).T, full_matrices=0)
        pinv_s = np.array([1./ss if ss > cutoff else 0. for ss in s])
        plow = np.dot(np.dot(u.T, obj*wlow)*pinv_s, v)
        clow1 = midCp / constants.r - sum(plow[i-1]*midT**i for i in range(1,5))
        Cplow = (clow1 + sum(plow[i-1]*T**i for i in range(1,5))) * constants.r
        
        u, s, v = np.linalg.svd((a*whigh).T, full_matrices=0)
        pinv_s = np.array([1./ss if ss > cutoff else 0. for ss in s])
        phigh = np.dot(np.dot(u.T, obj*whigh)*pinv_s, v)
        chigh1 = midCp / constants.r - sum(phigh[i-1]*midT**i for i in range(1,5))
        Cphigh = (chigh1 + sum(phigh[i-1]*T**i for i in range(1,5))) * constants.r

        fitCp = np.where(T <= midT, Cplow, Cphigh)
        if rescoefs: return fitCp, (clow1, *plow), (chigh1, *phigh)
        else: return Cp -fitCp

    # linear least sq. fit to determine Cp and/at midT
    params, cov, infodict, mesg, ier = leastsq(fc, guess, full_output=True, epsfcn=1e-8)
    if verbose: 
        print("Number of function calls: %d" % (infodict['nfev'],))
        print("ier = %d: %s" % (ier, mesg))
    if ier not in [1,2,3,4]:
        raise ValueError("leastsq solution not found (ier = %d)." % (ier))
    
    if opt_midT: midT, midCp = params
    else: midCp = params; midT = midT_fix
    fitCp, clow, chigh = fc(params, rescoefs=True)
    if verbose:
        dCpdT_low = sum(i*clow[i]*midT**(i-1.) for i in range(1,5)) * constants.r
        dCpdT_high = sum(i*chigh[i]*midT**(i-1.) for i in range(1,5)) * constants.r
        print("dCp/dT at midT(=%.2fK) = %.4f, %.4f [J/K2/mol] (%.2f%% diff.)" 
              % (midT, dCpdT_low, dCpdT_high, 100*abs(dCpdT_low-dCpdT_high)/dCpdT_low))
    
    T298 = 298.15
    a6low = H298*1000./constants.r - sum(clow[i]/(i+1.)*T298**(i+1) for i in range(5))
    a7low = S298/constants.r - (clow[0]*np.log(T298) + sum(clow[i]/i*T298**i for i in range(1,5)))
    a6high = a6low + sum((clow[i]-chigh[i])/(i+1.)*midT**(i+1) for i in range(5))
    a7high = a7low + (clow[0]-chigh[0])*np.log(midT) + sum((clow[i]-chigh[i])/i*midT**i for i in range(1,5))
    
    coeffs = [*chigh, a6high, a7high, *clow, a6low, a7low]
    return min(T), max(T), midT, coeffs


def therm_nasa7stat(T, H, S, Cp, lowT, highT, midT, coeffs, verbose=False):
    """ Calculate stats (MSD, MAD, RMSD, MAXPOSD, MAXNEGD) of NASA-7 polynomial coefficients """
    T = np.atleast_1d(T)
    Hdevs, Sdevs, Cpdevs = None, None, None

    if verbose: print("               MSD     MAD     RMSD   MAXPOSD MAXNEGD")
    def devstat(dev, name=""):
        n = float(len(dev))
        msd = sum(dev)  / n
        mad = sum(np.abs(dev))  / n
        rmsd = np.sqrt(sum(np.square(dev)) / n)
        maxposd = max(0, max(dev))
        maxnegd = min(0, min(dev))
        if verbose:
            print("%-12s %7.3f %7.3f %7.3f %7.3f %7.3f" % (name, msd, mad, rmsd, maxposd, maxnegd))
        return msd, mad, rmsd, maxposd, maxnegd
    
    if H is not None:
        H = np.atleast_1d(H)
        Hhigh = constants.r / 1000. * (sum(coeffs[i]/(i+1.)*T**(i+1) for i in range(5)) + coeffs[5])
        Hlow = constants.r / 1000. * (sum(coeffs[i+7]/(i+1.)*T**(i+1) for i in range(5)) + coeffs[5+7])
        Hdevs = devstat(np.where(T<=midT, Hlow, Hhigh) - H, "H[kJ/mol]")

    if S is not None:
        S = np.atleast_1d(S)
        Shigh = constants.r * (coeffs[0]*np.log(T) + sum(coeffs[i]/i*T**i for i in range(1,5)) + coeffs[6])
        Slow = constants.r * (coeffs[7]*np.log(T) + sum(coeffs[i+7]/i*T**i for i in range(1,5)) + coeffs[6+7])
        Sdevs = devstat(np.where(T<=midT, Slow, Shigh) - S, "S[J/K/mol]")

    if Cp is not None:
        Cp = np.atleast_1d(Cp)
        Cphigh = constants.r * sum(coeffs[i]*T**i for i in range(5))
        Cplow = constants.r * sum(coeffs[i+7]*T**i for i in range(5))
        Cpdevs = devstat(np.where(T<=midT, Cplow, Cphigh) - Cp, "Cp[J/K/mol]")

    return Hdevs, Sdevs, Cpdevs


def therm_nasa7str(name, comment, atoms_or_formula, phase, lowT, highT, midT, coeffs):
    """ Format NASA-7 polynomial coefficients """
    if isinstance(atoms_or_formula, str): atoms = name2atoms(atoms_or_formula)
    else: atoms = atoms_or_formula
    elements = []
    for x in atoms:
        if x not in elements: elements.append(x)
    if len(elements) > 5: raise ValueError("ERROR: No. of elements exceeds 5: %s" % (str(elements)))
    
    def atom2str(idx):
        if idx > len(elements) - 1: return "     "
        else: return "%-2s%3d" % (elements[idx], atoms.count(elements[idx]))
    def formcoeff(x):
        s = "% 16.9E" % x
        if s[-3] == "0": s = s[:-3] + s[-2:]
        else: s = s[:-5] + s[-4:]
        return s
    
    if len(name) <= 10:
        namecomment = "%-10s%14s" % (name, comment)
    else:
        cfmt = "%%%ds" % (24 - len(name) - 1)
        namecomment = name + " " + (cfmt % comment)
    if len(namecomment) > 24: namecomment = namecomment[:24]
    
    s = "%s%s%s%s%s%s% 10.3f% 10.3f% 8.2f%s 1\n" \
        % (namecomment, atom2str(0), atom2str(1), atom2str(2),
           atom2str(3), phase, lowT, highT, midT, atom2str(4))
    s += "".join(map(formcoeff, coeffs[0:5])) + "    2\n"
    s += "".join(map(formcoeff, coeffs[5:10])) + "    3\n"
    s += "".join(map(formcoeff, coeffs[10:])) + "                   4\n"
    return s

