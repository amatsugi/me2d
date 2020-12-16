#! /usr/bin/env python3

"""
calculate rotational constants
"""

import os
import numpy as np

from . import constants


def get_com(coord, weights):
    com = np.zeros(3)
    for i in range(len(weights)):
        com += weights[i] * coord[i]
    return com / sum(weights)

def centerized(coord, weights):
    crd = []
    com = get_com(coord, weights)
    for i in range(len(weights)):
        crd.append(coord[i] - com)
    return crd

def get_prcaxes(weights, coord):
    itensor = np.zeros((3,3))
    for i, x in enumerate(coord):
        itensor[0][0] += weights[i] * (x[1]*x[1] + x[2]*x[2])
        itensor[1][1] += weights[i] * (x[2]*x[2] + x[0]*x[0])
        itensor[2][2] += weights[i] * (x[0]*x[0] + x[1]*x[1])
        itensor[0][1] -= weights[i] * x[0]*x[1]
        itensor[1][2] -= weights[i] * x[1]*x[2]
        itensor[0][2] -= weights[i] * x[2]*x[0]
    itensor[1][0] = itensor[0][1]
    itensor[2][0] = itensor[0][2]
    itensor[2][1] = itensor[1][2]
    imom, diagonalmat = np.linalg.eigh(itensor)
    prcaxes = diagonalmat.T
    
    imom_l = []
    prcaxes_l = []
    tmp = [(imom.real[i], prcaxes[i]) for i in range(3)]
    for x in sorted(tmp):
        imom_l.append(x[0])
        prcaxes_l.append(x[1])
    if np.linalg.det(prcaxes) < 0.: prcaxes[2] = - prcaxes[2]
    return imom_l, prcaxes_l


def external_rotc(atoms, coordA):
    """ rotational constants
    arguments:
      atoms: list of atoms
      coordA: coordinates [len(atoms)][3] in angstroms
    """
    coord = [np.array(x)/(constants.bohr/constants.ang) for x in coordA]  # bohr
    weights = [constants.amass[x]/(constants.au/constants.amu) for x in atoms]
    coord = centerized(coord, weights)
    imom, prcaxes = get_prcaxes(weights, coord)
    rotc = [1. / (2.*x) / constants.cm2eh for x in imom] # cm^-1
    return rotc


def internal_rotc(atoms, coordA, ipiv, itop, checktop=True):
    """ Reduced rotational constant internal rotation
    [Ref: K.S. Pitzer, J. Chem. Phys. 14 (1946) 239.]
    arguments:
      atoms: list of atoms
      coordA: coordinates [len(atoms)][3] in angstroms
      ipiv: index (starting from 1) of pivot atom
      itop: if int: index of pivot atom of the top moiety
            otherwise: tuple of indexes (starting from 1) of the top moiety
    """
    if isinstance(itop, int): itop = find_top(atoms, coordA, ipiv, itop)
    if checktop: check_top(atoms, coordA, ipiv, itop)
    coord = [np.array(x)/(constants.bohr/constants.ang) for x in coordA]  # bohr
    weights = [constants.amass[x]/(constants.au/constants.amu) for x in atoms]
    coord = centerized(coord, weights)
    imom, prcaxes = get_prcaxes(weights, coord)
    
    piv = coord[ipiv-1]
    coord_top = [coord[i-1] for i in itop]
    wt_top = [weights[i-1] for i in itop]
    
    za = piv - coord_top[0]
    za /= np.linalg.norm(za)
    topcom = get_com(coord_top[1:], wt_top[1:])
    pt = np.dot(za, topcom - piv)
    origin = piv + pt * za
    
    xa = topcom - origin
    xlen = np.linalg.norm(xa)
    while xlen < 1e-6:
        dummycom = np.random.uniform(-2., 2., 3)
        pt = np.dot(za, dummycom - piv)
        dummyorigin = piv + pt * za
        xa = dummycom - dummyorigin
        xlen = np.linalg.norm(xa)
    xa /= xlen
    
    ya = np.cross(za, xa)
    ya /= np.linalg.norm(ya)
    maxes = [xa, ya, za]
    
    alpha = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            alpha[i][j] = np.dot(maxes[i], prcaxes[j])
    
    A, B, C, U = 0., 0., 0., 0.
    for i in range(1, len(coord_top)):
        w = wt_top[i]
        rel_coord = coord_top[i] - origin
        x = np.dot(maxes[0], rel_coord)
        y = np.dot(maxes[1], rel_coord)
        z = np.dot(maxes[2], rel_coord)
        A += w * (x*x + y*y)
        B += w * x * z
        C += w * y * z
        U += w * x
    
    r = [np.dot(origin, prcaxes[0]),
         np.dot(origin, prcaxes[1]),
         np.dot(origin, prcaxes[2])]
    beta = np.zeros(3)
    for i in range(3):
        im1 = i - 1
        if im1 == -1: im1 = 2
        ip1 = i + 1
        if ip1 == 3: ip1 = 0
        beta[i] = alpha[2][i] * A - alpha[0][i] * B - alpha[1][i] * C \
                  + U * (alpha[1][im1] * r[ip1] - alpha[1][ip1] * r[im1])
    
    lam = 0.
    M = sum(weights)
    for i in range(3):
        lam += (alpha[1][i] * U)**2 / M + (beta[i])**2 / imom[i]
    Im = A - lam
    B = 1. / (2.*Im) / constants.cm2eh # cm^-1
    return B


max_bond_X = 1.8  # for X-Y (X != H and Y != H)
max_bond_H1 = 1.4 # for X-H (X != H)
max_bond_H2 = 1.0 # for H-H
max_bond_piv = 3.0  # between pivot atoms

def max_bond(atom1, atom2):
    count = 0
    if atom1 in ("H", "h"): count += 1
    if atom2 in ("H", "h"): count += 1
    if count == 0: return max_bond_X
    elif count == 1: return max_bond_H1
    elif count == 2: return max_bond_H2


def find_top(atoms, coordA, ipiv1, ipiv2):
    """ find top moiety for ipiv2 """
    coordA = [np.array(x) for x in coordA]
    n = len(coordA)
    bmat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            bmat[i][j] = np.linalg.norm(coordA[i] - coordA[j])
    
    rest = []
    for i in range(len(atoms)):
        if (i != ipiv1-1) and (i != ipiv2-1): rest.append(i+1)
    itop = [ipiv2]
    num_prev = len(itop)
    while True:
        for ir in rest:
            for it in itop:
                if bmat[ir-1][it-1] < max_bond(atoms[ir-1], atoms[it-1]):
                    itop.append(ir)
                    rest.remove(ir)
                    break
        if len(itop) == num_prev: break
        num_prev = len(itop)
    return tuple(itop)

def check_top(atoms, coordA, ipiv, itop):
    """ check ipiv and itop """
    for i in range(len(itop)):
        if itop[i] == ipiv:
            raise ValueError("ERROR: itop and ipiv duplicated")
        for j in range(i+1, len(itop)):
            if itop[i] == itop[j]:
                raise ValueError("ERROR: itop duplicated")
    
    coordA = [np.array(x) for x in coordA]
    n = len(coordA)
    bmat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            bmat[i][j] = np.linalg.norm(coordA[i] - coordA[j])
    
    if bmat[ipiv-1][itop[0]-1] > max_bond_piv:
        print("WARNING: check_top: large separation between pivot atoms, %s A"
              % (bmat[ipiv-1][itop[0]-1]))
    
    for i in range(len(itop)):
        isolated = True
        for j in range(len(itop)):
            if i == j: continue
            if bmat[itop[i]-1][itop[j]-1] <= max_bond(atoms[itop[i]-1], atoms[itop[j]-1]):
                isolated = False
                break
        if isolated:
            print("WARNING: check_top: atom #%d separated from top" % (itop[i]))
    
    for i in range(n):
        if i + 1 == ipiv: continue
        if i + 1 in itop: continue
        for j in range(len(itop)):
            if bmat[i][itop[j]-1] <= max_bond(atoms[i], atoms[itop[j]-1]):
                print(("WARNING: check_top: atom #%d placed close "
                       "to atom #%d in top") % (i+1, itop[j]))
    
    return


def read_geom_xyzfile(xyzfn):
    fp = open(xyzfn)
    natoms = int(next(fp).strip())
    title = next(fp).rstrip()
    atoms, coordA = [], []
    for i in range(natoms):
        ls = next(fp).split()
        atm = ls[0]
        if atm in constants.anum2atom:
            atm = constants.anum2atom[atm]
        atoms.append(atm)
        coordA.append(np.array([float(x) for x in ls[1:4]]))
    fp.close()
    return atoms, coordA


def read_geom_xyzstr(s):
    s = s.strip().splitlines()
    natoms = int(s[0].strip())
    title = s[1].rstrip()
    atoms, coordA = [], []
    for i in range(natoms):
        ls = s[i+2].split()
        atm = ls[0]
        if atm in constants.anum2atom:
            atm = constants.anum2atom[atm]
        atoms.append(atm)
        coordA.append(np.array([float(x) for x in ls[1:4]]))
    return atoms, coordA


def read_geom_gaussian(fn):
    fp = open(fn)
    atoms = None
    coordA = None
    for l in fp:
        if (l.find("Z-Matrix orientation:") > -1
            or l.find("Input orientation:") > -1
            or l.find("Standard orientation:") > -1):
            for i in range(4): next(fp)
            atoms = []
            coordA = []
            while True:
                ls = next(fp).split()
                if len(ls) < 6: break
                anum = int(ls[1])
                atoms.append(constants.anum2atom[anum])
                coordA.append([float(ls[j+3]) for j in range(3)])
    return atoms, coordA
    

def read_geom(fn_or_str):
    if os.path.exists(fn_or_str):
        if fn_or_str.endswith(".xyz"):
            return read_geom_xyzfile(fn_or_str)
        elif fn_or_str.endswith(".log"):
            return read_geom_gaussian(fn_or_str)
    else:
        l = fn_or_str.strip().splitlines()
        if len(l) > 2 and l[0].strip().isdigit():
            return read_geom_xyzstr(fn_or_str)
    raise ValueError("# %s: cannot read geometry" % (fn_or_str))

