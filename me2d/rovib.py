#! /usr/bin/env python3

"""
Rovibrational properties
"""

import os, sys, time
from io import StringIO
import numpy as np
from scipy.special import gamma, ellipk, ellipe
from scipy.special import iv as bessel_iv

from . import constants


class RoVib(object):
    """ Rovibrational properties of molecule """
    
    def __init__(self, nsym=1, rotA=None, rotB2D=None, freq=None, fscale=1.,
                 introt=None, states=None, freqimg=None):
        self.nsym = nsym     # symmetry number
        self.rotA = rotA     # rot. const. for active 1D rotor [cm^-1]
        self.rotB2D = rotB2D # rot. const. for 2D J-rotor [cm^-1]

        # list of frequencies [cm^-1] and scaling factor
        if freq is None: self.freq = []
        else: self.freq = freq
        self.fscale = fscale

        # internal rotor [(B, sig, hofreq, V0), ...] (B, hofreq, V0 in cm^-1)
        # set hofreq=0 for free rotor
        # set V0<0 to estimate V0 from hofreq
        if introt is None: self.introt = []
        else: self.introt = introt

        # states [(degen, energy), ....]  (energy in cm^-1)
        # None for single state
        if states is None: self.states = [(1, 0.)]
        else: self.states = states
        
        self.freqimg = freqimg  # imaginary freq for TS [cm^-1]
    
    def findrot(self, convK=True, convJ=True):
        freerot = []  # [(B, sigma, dimen), ...]
        hindrot = []  # [(B, sigma, hofreq, V0), ...]
        if convK and (self.rotA is not None):
            freerot.append((self.rotA, 1, 1))
        if convJ and (self.rotB2D is not None):
            freerot.append((self.rotB2D, 1, 2))
        for x in self.introt:
            B, sigma, hofreq, V0 = x
            if hofreq <= 0.:
                freerot.append((B, sigma, 1))
            else:
                if V0 < 0.: # use estimated value if V0<0
                    V0 = float(hofreq * hofreq) / (sigma * sigma * B)
                hindrot.append((B, sigma, hofreq, V0))
        return freerot, hindrot
    
    def dens(self, nbin, dE, convK=True, convJ=True, dEint=1):
        freerot, hindrot = self.findrot(convK, convJ)
        freq  = [x*self.fscale for x in self.freq]
        return mod_bs(nbin, dE, self.nsym, freq, freerot, hindrot,
                      self.states, issum=False, dEint=dEint)
    
    def sums(self, nbin, dE, convK=True, convJ=True, dEint=1):
        freerot, hindrot = self.findrot(convK, convJ)
        freq  = [x*self.fscale for x in self.freq]
        return mod_bs(nbin, dE, self.nsym, freq, freerot, hindrot,
                      self.states, issum=True, dEint=dEint)

    def part(self, T, convK=True, convJ=True):
        freerot, hindrot = self.findrot(convK, convJ)
        freq  = [x*self.fscale for x in self.freq]
        return partfunc_rv(T, self.nsym, freq, freerot, hindrot,
                           self.states)

    def __str__(self):
        outf = StringIO()
        self.write_to(outf=outf)
        return outf.getvalue()

    def write_to(self, outf=None, prefix=None, E0=None, deltaH0=None, rotB2Dprod=None):
        if outf is None: fp = sys.stdout
        elif hasattr(outf, "write"): fp = outf
        else: fp = open(outf, "w")
        
        if prefix is None: prefix = ""
        if self.nsym is not None:
            fp.write("%s nsym = %g\n" % (prefix, self.nsym))
        if self.rotA is not None:
            fp.write("%s rotA = %g\n" % (prefix, self.rotA))
        if self.rotB2D is not None:
            fp.write("%s rotB2D = %g\n" % (prefix, self.rotB2D))
        
        if self.freq is not None and len(self.freq) > 0:
            n1 = 8
            div, mod = divmod(len(self.freq), n1)
            spref = "%s freq = [" % (prefix)
            for i in range(div+1):
                tmp = self.freq[i*n1 :min((i+1)*n1, len(self.freq))]
                tmp = ", ".join("%g" % x for x in tmp)
                if (i == div-1 and mod == 0) or (i == div):
                    fp.write("%s%s]\n" % (spref, tmp))
                    break
                fp.write("%s%s,\n" % (spref, tmp))
                spref = "%s%s" % (prefix, " "*(len(spref)-len(prefix)))

        if self.fscale is not None:
            fp.write("%s fscale = %g\n" % (prefix, self.fscale))
                
        if self.introt is not None and len(self.introt) > 0:
            spref = "%s introt = [" % (prefix)
            for i in range(len(self.introt)):
                x = self.introt[i]
                addstr = ""
                if x[2] <= 0.:
                    addstr = "  # free rotor"
                if x[2] > 0. and x[3] < 0.:
                    V0 = float(x[2]*x[2]) / (x[1]*x[1]*x[0])
                    addstr = "  # V0 = %g" % (V0)
                if i == len(self.introt) - 1:
                    fp.write("%s(%g, %g, %g, %g)]%s\n" % (spref, *x, addstr))
                    break
                fp.write("%s(%g, %g, %g, %g),%s\n" % (spref, *x, addstr))
                spref = "%s%s" % (prefix, " "*(len(spref)-len(prefix)))
                
        if self.states is not None and len(self.states) > 0:
            n1 = 6
            div, mod = divmod(len(self.states), n1)
            spref = "%s states = [" % (prefix)
            for i in range(div+1):
                tmp = self.states[i*n1 : min((i+1)*n1, len(self.states))]
                tmp = ", ".join("(%g, %g)" % (x[0], x[1]) for x in tmp)
                if (i == div-1 and mod == 0) or (i == div):
                    fp.write("%s%s]\n" % (spref, tmp))
                    break
                fp.write("%s%s,\n" % (spref, tmp))
                spref = "%s%s" % (prefix, " "*(len(spref)-len(prefix)))
            
        if self.freqimg is not None:
            fp.write("%s freqimg = %g\n" % (prefix, self.freqimg))
        if E0 is not None:
            fp.write("%s E0 = %g\n" % (prefix, E0))
        if deltaH0 is not None:
            fp.write("%s deltaH0 = %g\n" % (prefix, deltaH0))
        if rotB2Dprod is not None:
            fp.write("%s rotB2Dprod = %g\n" % (prefix, rotB2Dprod))


def mod_bs(nbin, dE, nsym, freq, freerot, hindrot, states,
           issum=False, dEint=1):
    """ Modified Beyer-Swinehart algorithm for density or sum of states.
    [Ref: R.G. Gilbert and S.C. Smith, Theory of Unimolecular and
    Recombination Reactions, Blackwell: Oxford, UK (1990), page 157.]
    - classical free rotor and quasi-quantum sinusoidally hindered rotor
    - arbitrarily given states (electronic states / isomers / conformers)
    - Beyer-Swinehart direct count
    arguments:
      nbin:  number of output energy bins
      dE:    bin size for the output array
      nsym:  symmetru number
      freq:  list of frequencies
      freerot: [(B, sigma, dimen), ...]
      hindrot:[(B, sigma, hofreq, V0), ...]
      states: [(degen, energy), ...]
      issum: True to calculate sum (otherwise density) of states
      dEint: internal energy grain size
    """
    
    n = int(float((nbin+1)*dE)/dEint + 1)  # No. of internal grains
    
    # initialize work array with free-rotor rho*dEint
    prod = 1.0
    hdimsum = 0.
    if len(freerot) > 0:
        for x in freerot:
            B, sig, dim = x
            hdim = 0.5 * dim
            prod *= gamma(hdim) / float(sig) / B**hdim
            hdimsum += hdim
        work = (np.arange(1, n+1, dtype=float)*dEint)**(hdimsum-1.) \
               * prod / gamma(hdimsum) * dEint
    else:
        # no free rotor
        work = np.zeros(n)
        work[0] = 1.0
    
    # convolve hindered rotor density of states
    for x in hindrot:
        B, sig, hofreq, V0 = x
        rhohr = hindrot_ds(n, dEint, hofreq, B, sig, V0, issum=False)
        for i in reversed(range(n)): # convolve into work
            work[i] = (work[i::-1] * rhohr[:i+1]).sum()
    
    # convolve states (electronic states / isomers / conformers)
    if len(states) == 0: pass
    elif len(states) == 1 and int(float(states[0][1])/dEint+0.5) == 0:
        work *= states[0][0]
    else:
        rhost = np.zeros(n)
        for x in states:
            degen, level = x
            i = int(float(level)/dEint + 0.5)
            if i > n - 1: continue
            rhost[i] += float(degen)
        for i in reversed(range(n)): # convolve into work
            work[i] = (work[i::-1] * rhost[:i+1]).sum()
    
    # numerical integration for sum of states
    if issum:
        for i in range(1,n): work[i] += work[i-1]
        work *= dEint
    
    # Beyer-Swinehart direct count
    for fr in freq:
        ifr = int(float(fr)/dEint)
        if float(fr)/dEint - ifr >= 0.5: ifr += 1 # round to grain
        for j in range(ifr, n):
            work[j] += work[j-ifr]
    
    # reduce to dE*nbin
    res = np.zeros(nbin)
    istep = int(round(float(dE)/dEint))
    for i in range(nbin):
        res[i] = work[i*istep:(i+1)*istep].sum()
    res /= float(istep) * dEint * float(nsym)
    
    return res


def hindrot_ds(n, dE, hofreq, B, sigma, V0, issum=False):
    """ Density or sum of states for a sinusoidally hindered rotor
    Quantum effects are corrected by the Pitzer-Gwinn approximation.
    [Ref: V.D. Knyazev, J. Phys. Chem. A 102 (1998) 3916.]
    arguments:
      n: number of energy grains
      dE: energy grain size
      hofreq, B, sigma, V0: hindered rotor parameters
      issum: True to calculate sum (otherwise density) of states
    """
    w = np.zeros(n)
    Eq = 0.5*hofreq + dE * np.arange(1, n+1, dtype=float)
    ic = min(n, max(0, int((V0 - 0.5*hofreq) / dE))) # index for V0
    
    # E >= V0
    w[ic:] = 4./(np.pi*sigma) * np.sqrt(Eq[ic:]/B) * ellipe(V0/Eq[ic:])
    
    # E < V0
    Er = hofreq * (np.arange(int(V0/hofreq+0.5), dtype=float) + 0.5)
    for i in range(ic):
        w[i] = ellipk((Eq[i] - Er[:int(Eq[i]/hofreq+0.5)]) / V0).sum()
    w[:ic] *= hofreq * 2. / (np.pi*sigma*np.sqrt(B*V0))
        
    if not issum: w[1:] -= w[:n-1] # differentiation to give density
    return w


def partfunc_rv(T, nsym, freq, freerot, hindrot, states):
    """ Rovibrational partition function """
    T = np.atleast_1d(T)
    q = np.ones(len(T))
    kT = T / constants.cm2k
    def bol(E):
        return np.exp(- E / kT)
    
    for x in freerot:
        B, sig, dim = x
        q *= (kT/B)**(0.5*dim) / float(sig)
        if dim % 2 == 1: q *= np.sqrt(np.pi)
    
    for x in hindrot:
        B, sig, hofreq, V0 = x
        qho = 1. / (1. - bol(hofreq))
        qhocl = kT / hofreq
        qfr = (np.pi * kT / B)**0.5 / float(sig)
        qhrcl = qfr * bol(V0/2.) * bessel_iv(0, V0/2. / kT)
        q *= qhrcl * (qho / qhocl)

    if len(states) > 0:
        qstates = 0.
        for x in states:
            degen, level = x
            qstates += degen * bol(level)
        q *= qstates
    
    for fr in freq:
        q *= 1. / (1. - bol(fr))
    
    q /= float(nsym)
    return q


def equilibrium_consts(T, rovib_reac, rovib_prod, deltaH0):
    """ Equilibrium constants for unimolecular reactant and product """
    T = np.atleast_1d(T)
    qrovib_reac = rovib_reac.part(T)
    qrovib_prod = rovib_prod.part(T)
    keq = qrovib_prod / qrovib_reac * np.exp( - deltaH0 * constants.cm2k / T)
    return keq

def equilibrium_consts_12(T, rovib_reac, rovib_prod1, rovib_prod2,
                          redmass, gelec_ratio, deltaH0):
    """ Equilibrium constants for unimolecular reactant and bimolecular products, A => B + C.
    redmass: reduced mass (m_B * m_C / m_A)
    gelec_ratio: ratio of electronic degeneracy (gB * gC / gA)
    """
    T = np.atleast_1d(T)
    qrovib_reac = rovib_reac.part(T)
    qrovib_prod1 = rovib_prod1.part(T)
    qrovib_prod2 = rovib_prod2.part(T)
    conv = (2. * np.pi * constants.amu * constants.kb / constants.h / constants.h)**1.5 * 1e-6
    qtrans_ratio = conv * redmass**1.5 * T**1.5
    keq = qtrans_ratio * gelec_ratio * qrovib_prod1 * qrovib_prod2 / qrovib_reac \
          * np.exp( - deltaH0 * constants.cm2k / T)
    return keq


def read_rovib_gaussian(fn):
    """ Read rovibrational properties from gaussian output file """
    fp = open(fn)
    rovib = RoVib()
    
    freql = []
    freqimagl = []
    nsym = None
    rotl = []
    for l in fp:
        if l.startswith(" Frequencies --"):
            ls = l.split()
            if len(ls) > 2:
                for x in ls[2:]:
                    if float(x) > 0.0: freql.append(float(x))
                    else: freqimagl.append(float(x))
        if l.startswith(" Rotational symmetry number"):
            nsym = float(l.split()[3])
        if l.startswith(" Rotational constants (GHZ):"):
            ls = l.split()
            rotl = [float(x)*1e9/(constants.c0*100.) for x in ls[3:]]

    rovib.freq = freql
    if len(freqimagl) == 1:
        rovib.freqimg = abs(freqimagl[0])
    elif len(freqimagl) > 1:
        raise ValueError("multiple imag. freqs")

    if nsym is not None: rovib.nsym = nsym
    
    if len(rotl) == 3:
        rovib.rotA = rotl[0]
        rovib.rotB2D = np.sqrt(rotl[1]*rotl[2])
    else:
        rovib.rotB2D = rotl[0]
    return rovib


def read_rovib_gpo(fn):
    """ Read rovibrational properties from a GPOP (.gpo) file """
    fp = open(fn)
    rovib = RoVib()
    freq = None
    introt_idxVibl = []
    isTS = False
    for l in fp:
        ls = l.strip().split()
        if len(ls) == 0: continue
        if ls[0].startswith("#"): continue
        elif ls[0] == "scaleFactVib": rovib.fscale = float(ls[1])
        elif ls[0] == "rotSymNbr": rovib.nsym = float(ls[1])
        elif ls[0] == "numIsomers": rovib.states = [(int(ls[1]), 0.)]
        elif ls[0] == "rotcGHz":
            rotl = [float(x)*1e9/(constants.c0*100.) for x in ls[1:]]
            if len(rotl) == 3:
                rovib.rotA = rotl[0]
                rovib.rotB2D = np.sqrt(rotl[1]*rotl[2])
            else:
                rovib.rotB2D = rotl[0]
        elif ls[0] == "frequencies{":
            nf = int(ls[1])
            freq = []
            while True:
                ls = next(fp).strip().split()
                if ls[0] == "}frequencies": break
                freq.append(float(ls[1]))
                nnrm = int(next(fp).strip().split()[1])
                while True:
                    ls = next(fp).strip().split()
                    if ls[0] == "}nrm": break
            if len(freq) != nf: raise ValueError("len(freq) != nf")
        elif ls[0] == "isTS" and ls[1].lower() == "true": isTS = True
        elif ls[0] == "isomStates{": 
            ni = int(ls[1])
            isomstates = []
            while True:
                ls = next(fp).strip().split()
                if ls[0] == "}isomStates": break
                isomstates.append((int(ls[0]), float(ls[1])))
            if len(isomstates) != ni: raise ValueError("len(isomstates) != ni")
            rovib.states = isomstates
        elif ls[0] == "intRotors{":
            nr = int(ls[1])
            introt = []
            while True:
                ls = next(fp).strip().split()
                if ls[0] == "}intRotors": break
                idxVib = int(ls[0]) - 1
                hofreq = freq[idxVib]
                introt_idxVibl.append(idxVib)
                sig = int(ls[1])
                redB = float(ls[2])
                V0 = float(ls[5])
                if ls[3].lower() == "true": hofreq, V0 = 0., 0.
                elif ls[4].lower() == "false": V0 = -1.
                introt.append((redB, sig, hofreq, V0))
            rovib.introt = introt
            if len(introt) != nr: raise ValueError("len(introt) != nr")
    freq2 = []
    for i in range(len(freq)):
        if i not in introt_idxVibl: freq2.append(freq[i])
    freq = freq2[:]
    if isTS:
        freqimg = freq.pop(0)
        if freqimg >= 0.: raise ValueError("freqimg >= 0.")
        rovib.freqimg = abs(freqimg)
    for fr in freq:
        if fr <= 0.: raise ValueError("negative or zero freq.")
    rovib.freq = freq
    return rovib


def read_rovib(fn):
    if not os.path.exists(fn):
        raise ValueError("# %s: file not exist." % (fn))
    if fn.endswith(".gpo"):
        return read_rovib_gpo(fn)
    elif fn.endswith(".log"):
        return read_rovib_gaussian(fn)
    raise ValueError("# %s: file not supported." % (fn))


