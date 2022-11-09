#! /usr/bin/env python3

"""
utils
"""

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


