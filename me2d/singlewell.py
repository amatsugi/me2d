#! /usr/bin/env python3

"""
single well master equations
"""

import sys, os, time
import numpy as np

from . import __version__
from . import constants
from .solverlib import load_library
from .solverlib import set_num_threads
from .solverlib import restore_num_threads
from .solverlib import get_num_threads

from .readfile import read1d
from .readfile import read2d


class MEBase(object):
    """ Base class for master equations """
    
    def __init__(self):
        self.nsiz = None
        self.nchan = None
        self.rhoa = None
        self.ka = None
        self.kchl = None
        self.Ea = None
        self.is1d = False
        self.lib, self.libfn = load_library()
    
    
    def set_k(self, dens):
        """ virtual method reserved for strong-collider-in-J model """
        pass
    
    def get_channel_strings(self, name="k", unit="s-1"):
        return ["%s%d[%s]" % (name, ich+1, unit) for ich in range(self.nchan)]
    
    def k_chemical_activation(self, khpl, kcal, kdl):
        """ phenomenological rate constants for chemical activation reactions
        calculated from steady-state solutions
        khpl: HPL bimolecular rate constant
        kcal: apparent decomposition rate constants in the chemical activation steady state
        kdl: decomposition rate constants
        returns: kr(->recombination), kbl(->bimolecular products)
        """
        ndisch = len(kdl)
        sumkcal = sum(kcal)
        kr = khpl * sum(kdl) / sumkcal
        kbl = [khpl * (kcal[i] - kdl[i]) / sumkcal for i in range(ndisch)]
        return kr, kbl
    
    def hpl(self, T):
        self.set_k(None)
        ga = self.rhoa * np.exp(- self.Ea * constants.cm2k / T)
        ga /= ga.sum()
        # ktot, [kch], ga (thermal distrib.)
        return (self.ka * ga).sum(), [(kch * ga).sum() for kch in self.kchl], ga

    
    def solve(self, T, p, gguess=None, solver="", bandpcrit=1e-9, neig=1,
              chemact_ch=None, verbose=False, nthreads=None, maxmemGB=None, chkfn=None):
        """ solve ME by calling solve1d or solve2d function of the library
        T: temperature in K
        p: pressure in bar
        gguess: initial guess for iterative solver
        solver: see me2d.show_solvers()
        bandpcrit: truncation threshold for banded matrix (None to use dense matrix)
        neig: number of eigenpairs to be computed
        chemact_ch: recombination channel (for chemical activation with solver=LinEq; gguess has to be None)
        verbose: verbose flag (True/False or integer)
        nthreads: number of threads to be used in the computation
        maxmemGB: max memory size used by the solver in GB
        chkfn: file name for storing matrix data (optional; used with the solver options save/load)
        """
        logfp = sys.stdout
        
        if bandpcrit is None: bandpcrit = -1.
        if nthreads is not None: 
            max_threads_orig = get_num_threads()
            set_num_threads(nthreads)
        if maxmemGB is not None: self.lib.set_me_maxmem_GB(maxmemGB)
        if chkfn is None: chkfn = ""
        
        if verbose:
            logfp.write("%s ver.%s: %s.solve started at %s\n"
                        % (__package__, __version__, self.__class__.__name__,
                           time.strftime("%a, %d %b %Y %H:%M:%S")))
            logfp.write("Library: %s\n" % (self.libfn))
            max_threads = get_num_threads()
            if len(max_threads) > 0:
                for x in max_threads: logfp.write("%s max threads = %s\n" % (x[0], x[1]))
            if maxmemGB is not None:
                logfp.write("Max memory size = %s GB\n" % (self.lib.get_me_maxmem_GB()))
            
            solverstr = solver
            if bandpcrit >= 0: solverstr += " (banded, pcrit=%.1e)" % (bandpcrit)
            if chemact_ch is not None: solverstr += " (chemact_ch = %s)" % (chemact_ch)
            logfp.write("T=%.0f, p=%.2e, solver=%s\n" % (T, p, solverstr))
            logfp.flush()
            
        # p is given in bar
        dens = p * 0.1 / constants.kb / T  # molecule/cm3
        ZM = self.Z * dens # s^-1
        kbT = T / constants.cm2k # cm^-1
        
        self.set_k(dens)
        
        if verbose is True: verbose = 1
        elif verbose is False: verbose = 0
        
        if gguess is not None: vec = np.array(gguess)
        elif chemact_ch is not None:
            # chemical activation flux
            vec = self.rhoa * np.exp(- self.Ea * constants.cm2k / T)
            vec *= self.kchl[chemact_ch-1]
        else: vec = self.rhoa * np.exp(- self.Ea * constants.cm2k / T) # thermal distrib
        
        vals = np.zeros(neig)
        if len(vec) < neig*self.nsiz: vec = np.append(vec, np.zeros(neig*self.nsiz - len(vec)))
        
        if self.is1d:
            ptype = -1 # downward prob. given
            res = self.lib.solve1d(self.nsiz, neig, vals, vec,
                                   self.Ea, self.rhoa, self.ka,
                                   self.y_e, self.ainv_ea, ptype,
                                   bandpcrit, ZM, kbT, solver.encode(), chkfn.encode(), verbose)
        else:
            ptype = 0 # symmetrized prob. given
            res = self.lib.solve2d(self.nsiz, neig, vals, vec,
                                   self.Ea, self.ea, self.Ja, self.rhoa, self.ka,
                                   self.y_e, self.y_J, self.ainv_ea, self.ainv_Ja, ptype,
                                   bandpcrit, ZM, kbT, solver.encode(), chkfn.encode(), verbose)
            
        if nthreads is not None:
            restore_num_threads(max_threads_orig)
        
        if res < 0: raise ValueError("ERROR in solver: res = %g" % res)
        ksol = -vals[0]
        ga = vec[:self.nsiz]
        
        k = (self.ka * ga).sum()
        kchl = [(kch * ga).sum() for kch in self.kchl]
        kdiff = abs((k - ksol) / k)
        if (kdiff > 0.01) and (not solver.startswith("LinEq")):
            logfp.write("WARNING: |k-ksol|/|k| = %.2e (k=%.6e, ksol=%.6e)\n" % (kdiff, k, ksol))
            logfp.flush()
        return k, kchl, ga, vals, vec



class ME1D(MEBase):
    """ 1D master equation """
    
    @classmethod
    def read_from(cls, fn, maxE=None, offset0=None):
        """ read E-resolved rho and k[ch] from rrkmE file, and return an instance of the class """
        dE, rho, kl = read1d(fn, maxE=maxE)
        return cls(dE, rho, kl, offset0=offset0)
    
    
    def __init__(self, dE, rhoa, kl, offset0=None):
        """ arguments: E-resolved rho and k
        dE: scalar
        rhol: rho array
        kl: list(ch) of k arrays
        offset0: offset for E (reserved for multiwell ME)
        """
        super().__init__()
        
        self.nchan = len(kl)
        self.nsiz = len(rhoa)
        self.is1d = True
        self.dE = dE
        self.rhoa = np.array(rhoa)
        self.kl = kl
        self.offset0 = offset0
        if self.offset0 is None: self.offset0 = 0
        
        self.ka = np.zeros(self.nsiz)
        self.kchl = []
        for i in range(self.nchan):
            ka = np.array(self.kl[i])
            self.kchl.append(ka)
            self.ka += ka
            
        self.Ea = (np.arange(self.nsiz, dtype=np.float64) + self.offset0) * self.dE
    
    
    def set_params(self, Z, y_e, a_e):
        """ set energy transfer parameters
        P(E,E') = C(E') * exp[- (|E-E'|/a_e(E'))^y_e ] for downward transition
        Z: scalar
        y: scalar
        a: scalar(constant value) or callable(function of (E0))
        """
        self.Z = Z
        self.y_e = y_e
        # store a^-1
        if callable(a_e): self.ainv_ea = 1. / a_e(self.Ea)
        else: self.ainv_ea = 1. / a_e
        if np.isscalar(self.ainv_ea):
            self.ainv_ea = np.full(self.nsiz, self.ainv_ea)



class ME1DuJ(ME1D):
    """ Microcanonical strong-collider-in-J model
    [Ref: J.A. Miller et al., J. Phys. Chem. A 106 (2002) 4904.]
    """
    
    @classmethod
    def read_from(cls, fn, dJ, maxE=None, maxJ=None):
        """ read E,J-resolved rho and k[ch] from rrkmEJ file, and return an instance of the class """
        dE, B2D, Jl, offsetl, rhol, kll = read2d(fn, dJ=dJ, maxE=maxE, maxJ=maxJ)
        return cls(dE, Jl, offsetl, rhol, kll)
    
    
    def __init__(self, dE, Jl, offsetl, rhol, kll):
        """ arguments: E,J-resolved rho and k
        dE: scalar
        Jl: list of J; J = Jl[iJ]
        offsetl: list(J) of offset; offset = offsetl[iJ]
        rhol: list(J) of rho arrays, rho = rhol[iJ][:]  (2J+1 degeneracy is included in rho)
        kll: list(ch,J) of k arrays, k = kll[ich][iJ][:]
        """
        super(ME1D, self).__init__()
        
        self.nsiz = len(rhol[0])
        self.nchan = len(kll)
        self.is1d = True
        self.dE = dE
        self.Jl = Jl  # [iJ]
        self.offsetl = offsetl  # [iJ]
        self.rhol = rhol # [iJ][iE]
        self.kll = kll # [ichan][iJ][iE]
        self.Ea = np.arange(self.nsiz, dtype=np.float64) * self.dE
        
        self.kal = [] # J, E
        self.rhoa = np.zeros(self.nsiz) # E
        for iJ in range(len(self.rhol)):
            ka = np.zeros(len(self.rhol[iJ]))
            for ich in range(self.nchan): ka += self.kll[ich][iJ]
            self.kal.append(ka)
            self.rhoa[self.offsetl[iJ]:] += self.rhol[iJ]
        self.ka = None # E
        self.kchl = None # ch, E
        
        
    def set_k(self, dens):
        """ set J-averaged k of microcanonical strong-collider-in-J model
        dens: [M];  None for thermal average
        """
        self.ka = np.zeros(self.nsiz)
        self.kchl = [np.zeros(self.nsiz) for ich in range(self.nchan)]
        yl = np.zeros(self.nsiz)
        
        for iJ in range(len(self.Jl)):
            rho = self.rhol[iJ]
            kchl = [self.kll[ich][iJ] for ich in range(self.nchan)]
            igsh = self.offsetl[iJ]
            y = np.copy(rho)
            if dens is not None: y /= self.Z*dens + self.kal[iJ]
            yl[igsh:] += y
            for ich in range(self.nchan):
                self.kchl[ich][igsh:] += kchl[ich] * y
        
        for ich in range(self.nchan):
            self.kchl[ich] /= yl
            self.ka += self.kchl[ich]



class ME2D(MEBase):
    """ 2D master equation """
    
    @classmethod
    def read_from(cls, fn, dJ, maxE=None, maxJ=None, offset0=None):
        """ read E,J-resolved rho and k[ch] from rrkmEJ file, and return an instance of the class """
        dE, B2D, Jl, offsetl, rhol, kll = read2d(fn, dJ=dJ, maxE=maxE, maxJ=maxJ)
        return cls(dE, B2D, Jl, offsetl, rhol, kll, offset0=offset0)
    
    
    def __init__(self, dE, B2D, Jl, offsetl, rhol, kll, offset0=None):
        """ arguments: E,J-resolved rho and k
        dE: scalar, energy grain size
        B2D: scalar, B2D of reactant
        Jl: list of J; J = Jl[iJ]
        offsetl: list(J) of offset; offset = offsetl[iJ]
        rhol: list(J) of rho arrays, rho = rhol[iJ][:]  (2J+1 degeneracy is included in rho)
        kll: list(ch,J) of k arrays, k = kll[ich][iJ][:]
        offset0: offset for E (reserved for multiwell ME)
        """
        super().__init__()
        
        self.nchan = len(kll)
        self.is1d = False
        self.dE = dE
        self.Jl = Jl  # [iJ]
        self.offsetl = offsetl  # [iJ]
        self.rhol = rhol # [iJ][iE]
        self.kll = kll # [ichan][iJ][iE]
        self.offset0 = offset0
        if self.offset0 is None: self.offset0 = 0
        
        self.sizE = len(self.rhol[0]) # size for J=0
        self.sizJl = [0 for iE in range(self.sizE)]  # [iE], size list for E=0,dE,2dE,...
        for iJ in range(len(self.Jl)):
            offset = self.offsetl[iJ]
            if self.sizE <= offset: break
            for iE in range(self.sizE - offset):
                self.sizJl[iE+offset] += 1
        self.nsiz = sum(self.sizJl)
        
        # a: flattened 1d arrays [(E=0)J=Jmax,...,J=0 (E=dE)... (E=2dE)... ...]
        # outer: ascending in E, inner: descending in J
        self.Ea = np.zeros(self.nsiz)
        self.Ja = np.zeros(self.nsiz)
        self.rhoa = np.zeros(self.nsiz)
        self.ka = np.zeros(self.nsiz)
        self.kchl = [np.zeros(self.nsiz) for ich in range(self.nchan)]
        iEstart = 0
        for iE in range(self.sizE):
            for iJ in range(self.sizJl[iE]):
                iJr = self.sizJl[iE] -1 - iJ # descending in J
                offset = self.offsetl[iJ]
                self.Ea[iEstart+iJr] = self.dE*(self.offset0 + iE)
                self.Ja[iEstart+iJr] = self.Jl[iJ]
                self.rhoa[iEstart+iJr] = self.rhol[iJ][iE-offset]
                for ich in range(self.nchan):
                    self.kchl[ich][iEstart+iJr] = self.kll[ich][iJ][iE-offset]
            iEstart += self.sizJl[iE]
        for ich in range(self.nchan):
            self.ka += self.kchl[ich]

        self.ea = self.Ea - B2D * self.Ja * (self.Ja + 1)
        emin = self.dE * self.offset0
        for i in range(self.nsiz):
            if self.ea[i] < emin: self.ea[i] = emin

        return
    
    def set_params(self, Z, y_e, y_J, a_e, a_J):
        """ set collisional transition parameters
        P(E,J;E',J') = C(E',J') * [S(E,J)/S(E',J')] * F_s(E,J;E',J')
        F_s = exp[-(|eps-eps'|/a_e(E',J'))^y_e] * exp[-(|J-J'|/a_J(E',J'))^y_J]
        See [Ref: A. Matsugi, J. Phys. Chem. A 124 (2020) 6645]
        Z: scalar
        y: scalar
        a: scalar(constant value) or callable(function of (E0,J0))
        """
        self.Z = Z
        self.y_e = y_e
        self.y_J = y_J
        # store a^-1
        if callable(a_e): self.ainv_ea = 1. / a_e(self.Ea, self.Ja)
        else: self.ainv_ea = 1. / a_e
        if np.isscalar(self.ainv_ea):
            self.ainv_ea = np.full(self.nsiz, self.ainv_ea)
            
        if callable(a_J): self.ainv_Ja = 1. / a_J(self.Ea, self.Ja)
        else: self.ainv_Ja = 1. / a_J
        if np.isscalar(self.ainv_Ja):
            self.ainv_Ja = np.full(self.nsiz, self.ainv_Ja)
    
    
