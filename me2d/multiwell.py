#! /usr/bin/env python3

"""
multiple-well master equations
"""

import sys, os, time
import ctypes
import numpy as np
import scipy.linalg

from . import __version__
from . import constants
from .solverlib import load_library
from .solverlib import set_num_threads
from .solverlib import restore_num_threads
from .solverlib import get_num_threads

from .readfile import read1d
from .readfile import read2d
from .singlewell import ME1D
from .singlewell import ME2D


def prepare_multi1d(well_list, maxE=None):
    """ prepare well objects with aligned top-energy grains for multiple-well 1D ME
    well_list: list of tuple (well_name, rrkm_file, relative_energy)
    """
    names = []
    tmpl = []
    dE = None
    topEbin = None
    for well in well_list:
        name, fn, relE = well
        if name in names: raise ValueError("Duplicated well name: %s" % name)
        names.append(name)
        dE_this, rho, kl = read1d(fn)
        if dE is None: dE = dE_this
        elif dE != dE_this: raise ValueError("inconsistent dE: %s, %s" % (dE, dE_this))
        offset0 = int(relE/dE + 0.5)
        if (topEbin is None) or (topEbin > offset0 + len(rho)): topEbin = offset0 + len(rho)
        tmpl.append([rho, kl, offset0])
    
    if (maxE is not None) and ((topEbin-1) * dE > maxE): topEbin = int(maxE / dE) + 1
    
    wells = []
    for x in tmpl:
        rho, kl, offset0 = x
        nsiz_this = topEbin - offset0
        rho = rho[:nsiz_this]
        for ich in range(len(kl)): kl[ich] = kl[ich][:nsiz_this]
        wells.append(ME1D(dE, rho, kl, offset0))
    return names, wells


def prepare_multi2d(well_list, dJ, maxE=None, maxJ=None):
    """ prepare well objects with aligned top-energy grains for multiple-well 2D ME
    well_list: list of tuple (well_name, rrkmEJ_file, relative_energy)
    """
    names = []
    tmpl = []
    dE = None
    topEbin = None
    for well in well_list:
        name, fn, relE = well
        if name in names: raise ValueError("Duplicated well name: %s" % name)
        names.append(name)
        dE_this, B2D, Jl, offsetl, rhol, kll = read2d(fn, dJ, maxJ=maxJ)
        if dE is None: dE = dE_this
        elif dE != dE_this: raise ValueError("inconsistent dE: %s, %s" % (dE, dE_this))
        offset0 = int(relE/dE + 0.5)
        sizE = len(rhol[0])
        if (topEbin is None) or (topEbin > offset0 + sizE): topEbin = offset0 + sizE
        tmpl.append([B2D, Jl, offsetl, rhol, kll, offset0])
    
    if (maxE is not None) and ((topEbin-1) * dE > maxE): topEbin = int(maxE / dE) + 1
    
    wells = []
    for x in tmpl:
        B2D, Jl, offsetl, rhol, kll, offset0 = x
        for iJ in range(len(Jl)):
            sizE_this = topEbin - offset0 - offsetl[iJ]
            rhol[iJ] = rhol[iJ][:sizE_this]
            for ich in range(len(kll)): kll[ich][iJ] = kll[ich][iJ][:sizE_this]
        wells.append(ME2D(dE, B2D, Jl, offsetl, rhol, kll, offset0))
    return names, wells
    


class MEBaseMW(object):
    """ Base class for multiple-well master equations """
    
    def __init__(self):
        self.nwell = None
        self.names = None
        self.wells = None
        self.connections = None
        self.nsiz = None
        self.nsizl = None
        self.posl = None
        self.rhoa = None
        self.ka = None
        self.Ea = None
        self.channels = None
        self.kisom_sym = None
        self.kisom_i = None
        self.kisom_j = None
        self.is1d = False
        self.lib, self.libfn = load_library()

    def __getitem__(self, name):
        """ returns well object """
        if name in self.names:
            return self.wells[self.names.index(name)]
        else:
            raise IndexError("Unknown well name: %s" % (name))
        
    def set_kisom(self, iwell1, ich1, iwell2, ich2):
        """ virtual method """
        raise TypeError("method set_kisom not provided")
        
    
    def set_channels(self):
        # self.channels[iwell_from][ich] = iwell_to; None for dissociation
        self.channels = [[None for ich in range(well.nchan)] for well in self.wells]
        self.kisom_sym = np.array([]) # array of symmetrized k for isomerization
        self.kisom_i = np.array([], dtype=np.int64) # array of index i
        self.kisom_j = np.array([], dtype=np.int64) # array of index j (j>i)
        for icon in range(len(self.connections)):
            (name1, ch1), (name2, ch2) = self.connections[icon]
            ich1 = ch1 - 1
            ich2 = ch2 - 1
            if name1 not in self.names: raise ValueError("Well not found: %s" % name1)
            if name2 not in self.names: raise ValueError("Well not found: %s" % name2)
            if name1 == name2: raise ValueError("Invalid connection: %s" % (name1, name2))
            iwell1 = self.names.index(name1)
            iwell2 = self.names.index(name2)
            if ich1 < 0 or ich1 >= self.wells[iwell1].nchan:
                raise ValueError("Invalid channel %s (well %s)" % (ich1+1, name1))
            if ich2 < 0 or ich2 >= self.wells[iwell2].nchan:
                raise ValueError("Invalid channel %s (well %s)" % (ich2+1, name2))
            if self.channels[iwell1][ich1] is not None:
                raise ValueError("Duplicated channel %s (well %s)" % (ich1+1, name1))
            if self.channels[iwell2][ich2] is not None:
                raise ValueError("Duplicated channel %s (well %s)" % (ich2+1, name2))
            self.channels[iwell1][ich1] = iwell2
            self.channels[iwell2][ich2] = iwell1
            self.set_kisom(iwell1, ich1, iwell2, ich2)
        if self.nwell > 1:
            isolated = True
            for iwell in range(self.nwell):
                if all(x is None for x in self.channels[iwell]):
                    raise ValueError("isolated well: %s" % (self.names[iwell]))

    def get_channel_strings(self):
        kstrs = []
        for iwell in range(self.nwell):
            well = self.wells[iwell]
            for ich in range(well.nchan):
                if self.channels[iwell][ich] is not None: continue
                s = "%s-k%d" % (self.names[iwell], ich+1)
                if self.channels[iwell][ich] is None: s += "(dis)"
                else: s += "(to-%s)" % (self.names[self.channels[iwell][ich]])
                kstrs.append(s)
        xstrs = []
        for iwell in range(self.nwell):
            s = "x(%s)" % (self.names[iwell])
            xstrs.append(s)
        return kstrs, xstrs

    def get_channel_strings_phnm(self):
        dischl = []
        for iwell in range(self.nwell):
            well = self.wells[iwell]
            for ich in range(well.nchan):
                if self.channels[iwell][ich] is None:
                    dischl.append("%s-ch%d" % (self.names[iwell], ich+1))
        kdstrl = [[] for x in dischl]
        for i in range(len(dischl)):
            for j in range(self.nwell):
                kdstrl[i].append("%s->%s" % (self.names[j], dischl[i]))
        kwstrl = [[] for i in range(self.nwell)]
        for i in range(self.nwell):
            jc = 0
            for j in range(self.nwell):
                if i == j: continue
                kwstrl[i].append("%s->%s" % (self.names[j], self.names[i]))
        return kdstrl, kwstrl

    def kphnm_from_ss(self, kll, popll):
        """ phenomenological rate constants from steady-state solutions
        kll: channel-specific overall decomposition rate constants in
             the steady-state decomposition of wells
        popll: steady-state populations during the steady-state
               decomposition of wells
        returns lists of the rate constants for dissociation and isomerization,
        kdl and kwl, corresponding to kdstrl and kwstrl, respectively, of the
        get_channel_strings_phnm() method
        """
        ndisch = len(kll[0])
        G = np.zeros((self.nwell, self.nwell))
        for i in range(self.nwell):
            for j in range(self.nwell): G[i,j] = popll[i][j]
        kdl = [np.zeros(self.nwell) for ich in range(ndisch)]
        G_LU = scipy.linalg.lu_factor(G)
        for i in range(ndisch):
            kssl = np.zeros(self.nwell)
            for j in range(self.nwell): kssl[j] = kll[j][i]
            kdl[i] = scipy.linalg.lu_solve(G_LU, kssl)
        
        nm1 = self.nwell-1
        GW = np.zeros((self.nwell*nm1, self.nwell*nm1))
        dvec = np.zeros(self.nwell*nm1)
        for i in range(self.nwell):
            jc = 0
            for j in range(self.nwell):
                if i == j: continue
                dvec[i*nm1+jc] = G[j][i] * sum(kdl[ich][i] for ich in range(ndisch))
                kc = 0
                for k in range(self.nwell):
                    if k == i: continue
                    GW[i*nm1+jc][i*nm1+kc] = G[j][k]
                    ic = i
                    if i > k: ic -= 1
                    GW[i*nm1+jc][k*nm1+ic] = -G[j][i]
                    kc += 1
                jc += 1
        GW_LU = scipy.linalg.lu_factor(GW)
        kw_all = scipy.linalg.lu_solve(GW_LU, dvec)
        kwl = [np.copy(kw_all[iwell*nm1:(iwell+1)*nm1]) for iwell in range(self.nwell)]
        return kdl, kwl

    def get_channel_strings_phnm_ca(self, chemact_well_ch):
        chemact_well = chemact_well_ch[0]
        chemact_ch = chemact_well_ch[1]
        if chemact_well in self.names: chemact_well = self.names.index(chemact_well)
        chemact = "%s-ch%d" % (self.names[chemact_well], chemact_ch)
        dischl = []
        for iwell in range(self.nwell):
            well = self.wells[iwell]
            for ich in range(well.nchan):
                if self.channels[iwell][ich] is None:
                    dischl.append("%s-ch%d" % (self.names[iwell], ich+1))
        krstrl = ["%s->%s" % (chemact, self.names[i]) for i in range(self.nwell)]
        kbstrl = []
        for x in dischl:
            if x == chemact: kbstrl.append("%s(no-rxn)" % (chemact))
            else: kbstrl.append("%s->%s" % (chemact, x))
        
        return krstrl, kbstrl
    
    def kphnm_from_cass(self, khpl, kl, popl, kdl, kwl):
        """ phenomenological rate constants from chemical activation steady-state solution
        khpl: HPL bimolecular rate constant
        kl: channel-specific apparent decomposition rate constants in
            the chemical activation steady state
        popl: steady-state populations during the chemical activation steady state
        kdl, kwl: outputs of kphnm_from_ss()
        returns lists of the rate constants for reactant-to-well and reactant-to-fragments,
        krl and kbl, corresponding to krstrl and kbstrl, respectively, of the
        get_channel_strings_phnm_ca() method
        """
        ndisch = len(kl)
        krl = [0. for i in range(self.nwell)]
        kbl = [0. for i in range(ndisch)]
        sumkl = sum(kl)
        for j in range(self.nwell):
            for i in range(self.nwell):
                if i == j: continue
                jc = j
                ic = i
                if j > i: jc -= 1
                if i > j: ic -= 1
                krl[j] += popl[j]*kwl[i][jc] - popl[i]*kwl[j][ic]
            for l in range(ndisch): krl[j] += popl[j] * kdl[l][j]
            krl[j] *= khpl / sumkl
    
        for l in range(ndisch):
            kbl[l] = kl[l]
            for j in range(self.nwell): kbl[l] -= popl[j] * kdl[l][j]
            kbl[l] *= khpl / sumkl
        return krl, kbl
    
    def hpl(self, T):
        ga = self.rhoa * np.exp(- self.Ea * constants.cm2k / T)
        ga /= ga.sum()
        kdis = 0.
        kl = []
        popl = []
        for iwell in range(self.nwell):
            well = self.wells[iwell]
            ga_this = ga[self.posl[iwell]:self.posl[iwell]+self.nsizl[iwell]]
            popl.append(ga_this.sum())
            for ich in range(well.nchan):
                if self.channels[iwell][ich] is not None: continue
                k = (well.kchl[ich] * ga_this).sum()
                kl.append(k)
                if self.channels[iwell][ich] is None: # dissoc
                    kdis += k
        return kdis, kl, ga, popl


    def solve(self, T, p, gguess=None, solver="", bandpcrit=1e-9, neig=1,
              reactant=None, chemact_well_ch=None,
              verbose=False, nthreads=None, maxmemGB=None):
        """ solve ME by calling solve1d or solve2d function of the library
        T: temperature in K
        p: pressure in bar
        gguess: initial guess for iterative solver
        solver: see me2d.show_solvers()
        bandpcrit: truncation threshold for banded matrix (None to use dense matrix)
        neig: number of eigenpairs to be computed
        reactant: name of the reactant well (only for InvIter solver for strady-state decomposition)
        chemact_well_ch: recombination (well-name, channel) (for chemical activation with solver=LinEq; gguess has to be None)
        verbose: verbose flag (True/False or integer)
        nthreads: number of threads to be used in the computation
        maxmemGB: max memory size used by the solver in GB
        """
        logfp = sys.stdout
        
        if bandpcrit is None: bandpcrit = -1.
        if reactant is None: reactant = -1
        elif reactant in self.names: reactant = self.names.index(reactant)
        if chemact_well_ch is not None:
            chemact_well = chemact_well_ch[0]
            chemact_ch = chemact_well_ch[1]
            if chemact_well in self.names: chemact_well = self.names.index(chemact_well)
        
        if nthreads is not None: 
            max_threads_orig = get_num_threads()
            set_num_threads(nthreads)
        if maxmemGB is not None: self.lib.set_me_maxmem_GB(maxmemGB)
        
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
            if reactant >= 0: solverstr += " (reactant = %s)" % (self.names[reactant])
            if chemact_well_ch is not None:
                solverstr += " (chemact: well = %s, ch = %s)" % (self.names[chemact_well], chemact_ch)
            logfp.write("T=%.0f, p=%.2e, solver=%s\n" % (T, p, solverstr))
            logfp.flush()
        
        nsiz = np.array(self.nsizl, dtype=np.int64)
        bandpcrit = np.full(self.nsiz, bandpcrit)
        
        # p is given in bar
        dens = p * 0.1 / constants.kb / T  # molecule/cm3
        ZM = np.array([well.Z * dens for well in self.wells]) # s^-1
        kbT = T / constants.cm2k # cm^-1
        
        if verbose is True: verbose = 1
        elif verbose is False: verbose = 0
        
        if gguess is not None: vec = np.array(gguess)
        elif chemact_well_ch is not None:
            # chemical activation flux
            if self.channels[chemact_well][chemact_ch-1] is not None:
                logfp.write("WARNING: THIS IS AN ISOMERIZATION CHANNEL!\n")
                logfp.flush()
            ga = self.rhoa * np.exp(- self.Ea * constants.cm2k / T)
            ga_sel = ga[self.posl[chemact_well]:self.posl[chemact_well]+self.nsizl[chemact_well]]
            flx = ga_sel * self.wells[chemact_well].kchl[chemact_ch-1]
            vec = np.zeros(len(ga))
            vec[self.posl[chemact_well]:self.posl[chemact_well]+self.nsizl[chemact_well]] = flx
        else: vec = self.rhoa * np.exp(- self.Ea * constants.cm2k / T) # thermal distrib
        
        vals = np.zeros(neig)
        if len(vec) < neig*self.nsiz: vec = np.append(vec, np.zeros(neig*self.nsiz - len(vec)))
        
        if self.is1d:
            y_e = np.array([well.y_e for well in self.wells])
            ainv_ea = np.concatenate([well.ainv_ea for well in self.wells])
            ptype = np.full(self.nsiz, -1, dtype=np.int64) # downward prob. given
            res = self.lib.solve1d_mw(self.nwell, nsiz, neig, vals, vec,
                                      self.Ea, self.rhoa, self.ka,
                                      y_e, ainv_ea, ptype,
                                      len(self.kisom_sym), self.kisom_sym,
                                      self.kisom_i, self.kisom_j,
                                      bandpcrit, ZM, kbT, solver.encode(), reactant, verbose)
        else:
            y_e = np.array([well.y_e for well in self.wells])
            y_J = np.array([well.y_J for well in self.wells])
            ainv_ea = np.concatenate([well.ainv_ea for well in self.wells])
            ainv_Ja = np.concatenate([well.ainv_Ja for well in self.wells])
            ptype = np.full(self.nsiz, 0, dtype=np.int64) # symmetrized prob. given
            res = self.lib.solve2d_mw(self.nwell, nsiz, neig, vals, vec,
                                      self.Ea, self.ea, self.Ja, self.rhoa, self.ka,
                                      y_e, y_J, ainv_ea, ainv_Ja, ptype,
                                      len(self.kisom_sym), self.kisom_sym,
                                      self.kisom_i, self.kisom_j,
                                      bandpcrit, ZM, kbT, solver.encode(), reactant, verbose)
        
        if nthreads is not None:
            restore_num_threads(max_threads_orig)
        
        if res < 0: raise ValueError("ERROR in solver: res = %g" % res)
        ksol = -vals[0]
        ga = vec[:self.nsiz]

        kdis = 0.
        kl = []
        popl = []
        for iwell in range(self.nwell):
            well = self.wells[iwell]
            ga_this = ga[self.posl[iwell]:self.posl[iwell]+self.nsizl[iwell]]
            popl.append(ga_this.sum())
            for ich in range(well.nchan):
                if self.channels[iwell][ich] is not None: continue
                k = (well.kchl[ich] * ga_this).sum()
                kl.append(k)
                if self.channels[iwell][ich] is None: # dissoc
                    kdis += k
        kdiff = abs((kdis - ksol) / kdis)
        if (kdiff > 0.01) and (not solver.startswith("LinEq")):
            logfp.write("WARNING: |kdis-ksol|/|kdis| = %.2e (kdis=%.6e, ksol=%.6e)\n"
                        % (kdiff, kdis, ksol))
            logfp.flush()
        return kdis, kl, ga, popl, vals, vec



class ME1DMW(MEBaseMW):
    """ multiple-well 1D master equation """
    
    @classmethod
    def read_from(cls, well_list, connections, maxE=None):
        """ read well_list and return an instance of the class
        arguments:
          well_list: list of tuple (well_name, rrkm_file, relative_energy)
          connections: list of isomerization channels ((well_name, ichan), (well_name, ichan))
                       (note: channel index starts from 1)
        """
        names, wells = prepare_multi1d(well_list, maxE=maxE)
        return cls(names, wells, connections)
    
    
    def __init__(self, names, wells, connections):
        """
        names: list of well names
        wells: list of ME1D objects
        connections: list of isomerization channels ((well_name, ichan), (well_name, ichan))
                     (note: channel index starts from 1)
        """
        super().__init__()
        
        self.is1d = True
        self.nwell = len(names)
        self.names = names
        self.wells = wells
        self.connections = connections
        
        self.dE = None
        self.topE = None
        self.nsizl = []
        self.posl = []
        for iwell in range(self.nwell):
            dE_this = self.wells[iwell].dE
            if self.dE is None: self.dE = dE_this
            elif self.dE != dE_this: 
                raise ValueError("inconsistent dE: %s, %s" % (self.dE, dE_this))
            topE_this = self.wells[iwell].Ea[-1]
            if self.topE is None: self.topE = topE_this
            elif self.topE != topE_this: 
                raise ValueError("inconsistent topE: %s, %s" % (self.topE, topE_this))
            
            self.nsizl.append(self.wells[iwell].nsiz)
            if iwell == 0: self.posl.append(0)
            else: self.posl.append(self.posl[iwell-1] + self.wells[iwell-1].nsiz)
        
        self.nsiz = sum(self.nsizl)
        self.Ea = np.concatenate([well.Ea for well in self.wells])
        self.rhoa = np.concatenate([well.rhoa for well in self.wells])
        self.ka = np.concatenate([well.ka for well in self.wells])
        
        self.set_channels()
        
    
    def set_kisom(self, iwell1, ich1, iwell2, ich2):
        well1 = self.wells[iwell1]
        well2 = self.wells[iwell2]
        if (well1.offset0 > well2.offset0):
            start1 = 0
            start2 = well1.offset0 - well2.offset0
            nk = well1.nsiz
        else:
            start1 = well2.offset0 - well1.offset0
            start2 = 0
            nk = well2.nsiz

        rho1 = well1.rhoa[start1:]
        rho2 = well2.rhoa[start2:]
        k1 = well1.kchl[ich1][start1:]
        k2 = well2.kchl[ich2][start2:]
        nh1 = rho1 * k1
        nh2 = rho2 * k2
        nh = np.sqrt(nh1) * np.sqrt(nh2)
        
        # check symmetry
        rdiffmax = 0.
        for ik in range(nk):
            # skip check for small k
            if (k1[ik] < 1e-6) and (k2[ik] < 1e-6): continue
            
            # +/- 1 grain tolerance (for numerical discretization error)
            diff = None
            if ik == 0:
                if nh1[ik] > nh2[ik+1]: diff = nh1[ik] - nh2[ik+1]
            elif ik == nk-1:
                if nh1[ik] < nh2[ik-1]: diff = nh2[ik-1] - nh1[ik]
            else:
                if nh1[ik] < nh2[ik-1]: diff = nh2[ik-1] - nh1[ik]
                elif nh1[ik] > nh2[ik+1]: diff = nh1[ik] - nh2[ik+1]
            if diff is not None:
                rdiff = abs(diff) / nh[ik]
                if rdiff > rdiffmax: rdiffmax = rdiff

        if rdiffmax > 0.3:
            raise ValueError("asymmetry detected: %s %% between %s and %s"
                             % (rdiffmax*100, self.names[iwell1], self.names[iwell2]))

        ksym = nh / (np.sqrt(rho1) * np.sqrt(rho2)) # store symmetrized k (= k_i * sqrt(rho_i/rho_j))
        self.kisom_sym = np.append(self.kisom_sym, ksym)

        pos1 = self.posl[iwell1] + start1
        pos2 = self.posl[iwell2] + start2
        pos1a = pos1 + np.arange(nk, dtype=np.int64)
        pos2a = pos2 + np.arange(nk, dtype=np.int64)
        if pos1 < pos2:
            self.kisom_i = np.append(self.kisom_i, pos1a) # array of index i
            self.kisom_j = np.append(self.kisom_j, pos2a) # array of index j (j>i)
        else:
            self.kisom_i = np.append(self.kisom_i, pos2a) # array of index i
            self.kisom_j = np.append(self.kisom_j, pos1a) # array of index j (j>i)
        return


class ME2DMW(MEBaseMW):
    """ multiple-well 2D master equation """
    
    @classmethod
    def read_from(cls, well_list, connections, dJ, maxE=None, maxJ=None):
        """ read well_list and return an instance of the class
        arguments:
          well_list: list of tuple (well_name, rrkmEJ_file, relative_energy)
          connections: list of isomerization channels ((well_name, ichan), (well_name, ichan))
                       (note: channel index starts from 1)
        """
        names, wells = prepare_multi2d(well_list, dJ, maxE=maxE, maxJ=maxJ)
        return cls(names, wells, connections)
    
    
    def __init__(self, names, wells, connections):
        """
        names: list of well names
        wells: list of ME2D objects
        connections: list of isomerization channels ((well_name, ichan), (well_name, ichan))
                     (note: channel index starts from 1)
        """
        super().__init__()
        
        self.is1d = False
        self.nwell = len(names)
        self.names = names
        self.wells = wells
        self.connections = connections
        
        self.dE = None
        self.topE = None
        self.Jl = None
        self.nsizl = []
        self.posl = []
        for iwell in range(self.nwell):
            dE_this = self.wells[iwell].dE
            if self.dE is None: self.dE = dE_this
            elif self.dE != dE_this: 
                raise ValueError("inconsistent dE: %s, %s" % (self.dE, dE_this))
            topE_this = max(self.wells[iwell].Ea)
            if self.topE is None: self.topE = topE_this
            elif self.topE != topE_this: 
                raise ValueError("inconsistent topE: %s, %s" % (self.topE, topE_this))

            Jl_this = self.wells[iwell].Jl
            if self.Jl is None: self.Jl = Jl_this[:]
            else:
                Jlen = min(len(self.Jl), len(Jl_this))
                for iJ in range(Jlen):
                    if self.Jl[iJ] != Jl_this[iJ]:
                        raise ValueError("inconsistent Jl: %s, %s" % (self.Jl, Jl_this))
                if len(Jl_this) > len(self.Jl): self.Jl = Jl_this[:]
            
            self.nsizl.append(self.wells[iwell].nsiz)
            if iwell == 0: self.posl.append(0)
            else: self.posl.append(self.posl[iwell-1] + self.wells[iwell-1].nsiz)
        
        self.nsiz = sum(self.nsizl)
        self.Ea = np.concatenate([well.Ea for well in self.wells])
        self.ea = np.concatenate([well.ea for well in self.wells])
        self.Ja = np.concatenate([well.Ja for well in self.wells])
        self.rhoa = np.concatenate([well.rhoa for well in self.wells])
        self.ka = np.concatenate([well.ka for well in self.wells])
        
        self.set_channels()
    
    
    def set_kisom(self, iwell1, ich1, iwell2, ich2):
        well1 = self.wells[iwell1]
        well2 = self.wells[iwell2]
        Jlen = min(len(well1.Jl), len(well2.Jl))
        for iJ in range(Jlen):
            J = well1.Jl[iJ]
            off1 = well1.offset0 + well1.offsetl[iJ]
            off2 = well2.offset0 + well2.offsetl[iJ]
            if (off1 > off2):
                start1 = 0
                start2 = off1 - off2
                nk = len(well1.rhol[iJ])
            else:
                start1 = off2 - off1
                start2 = 0
                nk = len(well2.rhol[iJ])
            
            rho1 = well1.rhol[iJ][start1:]
            rho2 = well2.rhol[iJ][start2:]
            k1 = well1.kll[ich1][iJ][start1:]
            k2 = well2.kll[ich2][iJ][start2:]
            nh1 = rho1 * k1
            nh2 = rho2 * k2
            nh = np.sqrt(nh1) * np.sqrt(nh2)
        
            # check symmetry
            rdiffmax = 0.
            for ik in range(nk):
                # skip check for small k
                if (k1[ik] < 1e-6) and (k2[ik] < 1e-6): continue
                
                # +/- 1 grain tolerance (for numerical discretization error)
                diff = None
                if ik == 0:
                    if nh1[ik] > nh2[ik+1]: diff = nh1[ik] - nh2[ik+1]
                elif ik == nk-1:
                    if nh1[ik] < nh2[ik-1]: diff = nh2[ik-1] - nh1[ik]
                else:
                    if nh1[ik] < nh2[ik-1]: diff = nh2[ik-1] - nh1[ik]
                    elif nh1[ik] > nh2[ik+1]: diff = nh1[ik] - nh2[ik+1]
                if diff is not None:
                    rdiff = abs(diff) / nh[ik]
                    if rdiff > rdiffmax: rdiffmax = rdiff
            
            if rdiffmax > 0.3:
                raise ValueError("asymmetry detected: %s %% between %s and %s"
                                 % (rdiffmax*100, self.names[iwell1], self.names[iwell2]))
            
            ksym = nh / (np.sqrt(rho1) * np.sqrt(rho2)) # store symmetrized k (= k_i * sqrt(rho_i/rho_j))
            self.kisom_sym = np.append(self.kisom_sym, ksym)
            
            pos1a = np.zeros(nk, dtype=np.int64)
            iEstart = 0
            for iE in range(well1.offsetl[iJ] + start1):
                iEstart += well1.sizJl[iE]
            for ik in range(nk):
                iE = well1.offsetl[iJ] + start1 + ik
                iJr = well1.sizJl[iE] - 1 - iJ
                pos1a[ik] = self.posl[iwell1] + iEstart + iJr
                iEstart += well1.sizJl[iE]
            
            pos2a = np.zeros(nk, dtype=np.int64)
            iEstart = 0
            for iE in range(well2.offsetl[iJ] + start2):
                iEstart += well2.sizJl[iE]
            for ik in range(nk):
                iE = well2.offsetl[iJ] + start2 + ik
                iJr = well2.sizJl[iE] - 1 - iJ
                pos2a[ik] = self.posl[iwell2] + iEstart + iJr
                iEstart += well2.sizJl[iE]

            for ik in range(nk):
                if well1.Ea[pos1a[ik]-self.posl[iwell1]] != well2.Ea[pos2a[ik]-self.posl[iwell2]]:
                    raise ValueError("Error in set_kisom: inconsistent Ea")
                if well1.Ja[pos1a[ik]-self.posl[iwell1]] != well2.Ja[pos2a[ik]-self.posl[iwell2]]:
                    raise ValueError("Error in set_kisom: inconsistent Ja")
                    
            if self.posl[iwell1] < self.posl[iwell2]:
                self.kisom_i = np.append(self.kisom_i, pos1a) # array of index i
                self.kisom_j = np.append(self.kisom_j, pos2a) # array of index j (j>i)
            else:
                self.kisom_i = np.append(self.kisom_i, pos2a) # array of index i
                self.kisom_j = np.append(self.kisom_j, pos1a) # array of index j (j>i)
        return



