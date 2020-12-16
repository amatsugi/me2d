#include "extlib.h"
#include <iostream>
using std::cout;
using std::endl;

#if NO_ARPACK
// define dummy routines
void ARPACK(dsaupd)(arpack_int &ido, char &bmat, arpack_int &n, char *which,
                    arpack_int &nev, double &tol, double *resid, arpack_int &ncv,
                    double *v, arpack_int &ldv, arpack_int *iparam, arpack_int *ipntr,
                    double *workd, double *workl, arpack_int &lworkl, arpack_int &info) {
  cout << "ARPACK not linked." << endl;
  ido = 99;
  info = -99999;
}
void ARPACK(dseupd)(arpack_int &rvec, char &howmny, arpack_int *select, double *d, double *z,
                    arpack_int &ldz, double &sigma, char &bmat, arpack_int &n, char *which,
                    arpack_int &nev, double &tol, double *resid, arpack_int &ncv,
                    double *v, arpack_int &ldv, arpack_int *iparam, arpack_int *ipntr,
                    double *workd, double *workl, arpack_int &lworkl, arpack_int &info) {
  cout << "ARPACK not linked." << endl;
  info = -99999;
}
#endif


