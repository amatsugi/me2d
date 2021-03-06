# ME2D
Python/C++ code for solving two-dimensional master equations of unimolecular reactions

## Requirements
  - Python 3 with Numpy and Scipy packages
  - C++ compiler
  - BLAS and LAPACK libraries (such as [OpenBLAS](https://github.com/xianyi/OpenBLAS))
  - (optional) ARPACK library ([ARPACK-NG](https://github.com/opencollab/arpack-ng))

## Installation
- Compile the solver library
  - change directory to `me2d/lib`.
  - create `makefile` referring to the sample files `makefile.*`.
  - compile with `make`.

- Add the installation directory to the environment variable `PYTHONPATH`
 (e.g., `export PYTHONPATH=/home/user/lib/me2d:$PYTHONPATH`)

## Examples and References
Examples of use can be found in the `examples` directory.
- [examples/chf3](examples/chf3): Single-channel thermal decomposition of trifluoromethane [1]
- [examples/nc3h7](examples/nc3h7): Two-channel thermal decomposition of n-propyl [2]
- [examples/c5h11](examples/c5h11): Thermal decomposition in the multiple-well (1-,2-pentyl) system [2]

References  
 [1] [J. Phys. Chem. A 124 (2020) 6645](https://doi.org/10.1021/acs.jpca.0c05906)  
 [2] [J. Phys. Chem. A 125 (2021) 2532](https://doi.org/10.1021/acs.jpca.1c00666)  
