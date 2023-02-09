#! /usr/bin/env python3

"""
utils
"""

import numpy as np

from . import constants


def name2weight(name):
    """ moleculer weights [amu] of given molecular formula """
    return sum(constants.amass[x] for x in name2atoms(name))

def name2atoms(name):
    """ list of atoms from given molecular formula """
    name = name.strip()
    n = len(name)
    l = []
    s = ""
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
    return l2


def findmin(x, y):
    minind = list(y).index(min(y))
    if minind == 0:
        x1, x2, x3 = x[minind], x[minind+1], x[minind+2]
        y1, y2, y3 = y[minind], y[minind+1], y[minind+2]
    elif minind == len(x)-1:
        x1, x2, x3 = x[minind-2], x[minind-1], x[minind]
        y1, y2, y3 = y[minind-2], y[minind-1], y[minind]
    else:
        x1, x2, x3 = x[minind-1], x[minind], x[minind+1]
        y1, y2, y3 = y[minind-1], y[minind], y[minind+1]
    a = ((y1-y2)*(x1-x3) - (y1-y3)*(x1-x2)) / ((x1-x2)*(x1-x3)*(x2-x3))
    b = (y1-y2) / (x1-x2) - a * (x1+x2)
    c = y1 - a*x1*x1 - b*x1
    xmin = -b / (2. * a)
    ymin = a * xmin * xmin + b * xmin + c
    if (ymin > y[minind]) or (xmin < min(x)) or (xmin > max(x)):
        xmin = x[minind]
        ymin = y[minind]
    return xmin, ymin



def fit_arrhenius(T, k, mod=True, verbose=False):
    """ (modified) Arrhenius fit of rate constant """
    T = np.atleast_1d(T)
    k = np.atleast_1d(k)
    # mod: log(k) = log(A) + n*log(T) + (E/R)*(-1/T)
    # not mod: log(k) = log(A) + (E/R)*(-1/T)
    
    if mod: a = np.asarray([np.ones(len(T)), np.log(T), -1./T])
    else: a = np.asarray([np.ones(len(T)), -1./T])
    
    cutoff = len(T) * np.finfo(T.dtype).eps
        
    u, s, v = np.linalg.svd(a.T, full_matrices=0)
    pinv_s = np.array([1./ss if ss > cutoff else 0. for ss in s])
    res = np.dot(np.dot(u.T, np.log(k))*pinv_s, v)
    
    if mod: 
        logA, n, ER = res
        A = np.exp(logA)
        params = (A, n, ER)
        kfit = A * T**n * np.exp(-ER/T)
    else: 
        logA, ER = res
        A = np.exp(logA)
        params = (A, ER)
        kfit = A * np.exp(-ER/T)
    
    # relative deviation stats (in %)
    rdev = (kfit - k) / k
    num = float(len(rdev))
    mad = 100*sum(np.abs(rdev))  / num
    rmsd = 100*np.sqrt(sum(np.square(rdev)) / num)
    maxposd = 100*max(0, max(rdev))
    maxnegd = 100*min(0, min(rdev))
    stats = (mad, rmsd, maxposd, maxnegd)
    if verbose:
        if mod: print("k = %.3e * T**%.4f * exp(- %.1f/T)" % (A, n, ER))
        else:  print("k = %.3e * exp(- %.1f/T)" % (A, ER))
        print("MAD = %.1f%%, RMSD = %.1f%%, MAXPOSD = %.1f%%, MAXNEGD = %.1f%%" % stats)
    
    return params, kfit, stats
