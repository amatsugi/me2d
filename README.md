# ME2D
Python/C++ code for solving two-dimensional master equations of unimolecular reactions

Reference
  - [J. Phys. Chem. A 124(2020) 6645](https://doi.org/10.1021/acs.jpca.0c05906)

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

## Usage
Examples of use can be found in the `examples` directory.

