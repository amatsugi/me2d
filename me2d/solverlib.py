#! /usr/bin/env python3

"""
solver library (through ctypes)
"""

import sys, os
import ctypes
import ctypes.util
import numpy as np

_lib = None
_libfn = None

_threads_funcs = [("OpenMP", "omp_set_num_threads", "omp_get_max_threads"),
                 ("OpenBLAS", "openblas_set_num_threads", "openblas_get_num_threads"),
                 ("MKL", "MKL_Set_Num_Threads", "MKL_Get_Max_Threads")]


def load_library(libname="libme2d"):
    global _lib
    global _libfn
    if _lib is not None: return _lib, _libfn
    if sys.platform.startswith("win"): ext = ".dll"
    elif sys.platform.startswith("darwin"): ext = ".dylib"
    else: ext = ".so"
    lib = None
    libfn = None
    fdir = os.path.dirname(__file__)
    libfnl = [os.path.join(fdir, libname+ext),
              os.path.join(fdir, "lib", libname+ext),
              os.path.join(fdir, "bin", libname+ext),
              os.path.join(".", libname+ext),
              os.path.join("./lib", libname+ext),
              os.path.join("./bin", libname+ext)]
    for x in libfnl:
        if os.path.exists(x):
            libfn = x
            break
    if libfn is None:
        x = ctypes.util.find_library(libname)
        if x is not None: libfn = x
    
    if libfn is None:
        raise ValueError("Library not found: %s" % (libname))
    lib = ctypes.cdll.LoadLibrary(libfn)
    CDBL = ctypes.c_double
    CI64 = ctypes.c_int64
    CCHP = ctypes.c_char_p
    DARR = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C')
    IARR = np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags='C')
    lib.set_me_maxmem_GB.restype = None
    lib.set_me_maxmem_GB.argtypes = [CDBL]
    lib.get_me_maxmem_GB.restype = CDBL
    lib.get_me_maxmem_GB.argtypes = None
    lib.show_solvers.restype = None
    lib.show_solvers.argtypes = None
    lib.solve1d.restype = CI64
    lib.solve1d.argtypes = [CI64, CI64, DARR, DARR,
                            DARR, DARR, DARR,
                            CDBL, DARR, CI64,
                            CDBL, CDBL, CDBL, CCHP, CCHP, CI64]
    lib.solve2d.restype = CI64
    lib.solve2d.argtypes = [CI64, CI64, DARR, DARR,
                            DARR, DARR, DARR, DARR, DARR,
                            CDBL, CDBL, DARR, DARR, CI64,
                            CDBL, CDBL, CDBL, CCHP, CCHP, CI64]
    
    lib.solve1d_mw.restype = CI64
    lib.solve1d_mw.argtypes = [CI64, IARR, CI64, DARR, DARR,
                               DARR, DARR, DARR,
                               DARR, DARR, IARR,
                               CI64, DARR, IARR, IARR,
                               DARR, DARR, CDBL, CCHP, CI64, CCHP, CI64]
    lib.solve2d_mw.restype = CI64
    lib.solve2d_mw.argtypes = [CI64, IARR, CI64, DARR, DARR,
                               DARR, DARR, DARR, DARR, DARR,
                               DARR, DARR, DARR, DARR, IARR,
                               CI64, DARR, IARR, IARR,
                               DARR, DARR, CDBL, CCHP, CI64, CCHP, CI64]
    _lib = lib
    _libfn = libfn
    return _lib, _libfn



def show_solvers():
    lib, libfn = load_library()
    lib.show_solvers()
    return

def set_num_threads(n):
    lib, libfn = load_library()
    for x in _threads_funcs:
        func = getattr(lib, x[1], None)
        if func is not None:
            try: func(int(n))
            except: pass
    return
    
def restore_num_threads(max_threads):
    lib, libfn = load_library()
    for x in _threads_funcs:
        for y in max_threads:
            func = getattr(lib, x[1], None)
            if func is not None:
                try: func(int(y[1]))
                except: pass
    return

def get_num_threads():
    lib, libfn = load_library()
    max_threads = []
    for x in _threads_funcs:
        func = getattr(lib, x[2], None)
        if func is not None:
            try: max_threads.append((x[0], func()))
            except: pass
    return max_threads

