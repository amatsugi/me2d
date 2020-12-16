#ifndef _EXTLIB_H_
#define _EXTLIB_H_

#include "defs.h"

extern "C" {

void BLAS(dsymv)(char &uplo, blas_int &n, double &alpha, double *a, blas_int &lda,
                 double *x, blas_int &incx, double &beta, double *y, blas_int &incy);
void BLAS(dsbmv)(char &uplo, blas_int &n, blas_int &k, double &alpha, double *a, blas_int &lda,
                 double *x, blas_int &incx, double &beta, double *y, blas_int &incy);

void LAPACK(dpotrf)(char &uplo, lapack_int &n, double *a, lapack_int &lda,
                    lapack_int &info);
void LAPACK(dpotrs)(char &uplo, lapack_int &n, lapack_int &nrhs, double *a,
                    lapack_int &lda, double *b, lapack_int &ldb, lapack_int &info);
void LAPACK(dgetrf)(lapack_int &m, lapack_int &n, double *a, lapack_int &lda,
                    lapack_int *ipiv, lapack_int &info);
void LAPACK(dgetrs)(char &trans, lapack_int &n, lapack_int &nrhs, double *a,
                    lapack_int &lda, lapack_int *ipiv, double *b, lapack_int &ldb,
                    lapack_int &info);
void LAPACK(dsytrf)(char &uplo, lapack_int &n, double *a, lapack_int &lda,
                    lapack_int *ipiv, double *work, lapack_int &lwork,
                    lapack_int &info);
void LAPACK(dsytrs)(char &uplo, lapack_int &n, lapack_int &nrhs, double *a,
                    lapack_int &lda, lapack_int *ipiv, double *b, lapack_int &ldb,
                    lapack_int &info);

void LAPACK(dpbtrf)(char &uplo, lapack_int &n, lapack_int &kd, double *ab,
                    lapack_int &ldab, lapack_int &info);
void LAPACK(dpbtrs)(char &uplo, lapack_int &n, lapack_int &kd, lapack_int &nrhs,
                    double *ab, lapack_int &ldab, double *b, lapack_int &ldb,
                    lapack_int &info);

void LAPACK(dpttrf)(lapack_int &n, double *d, double *e, lapack_int&info);
void LAPACK(dpttrs)(lapack_int &n, lapack_int &nrhs, double *d, double *e,
                    double *b, lapack_int &ldb, lapack_int&info);

void LAPACK(dsyevr)(char &jobz, char &range, char &uplo, lapack_int &n, double *a,
                    lapack_int &lda, double &vl, double &vu, lapack_int &il,
                    lapack_int &iu, double &abstol, lapack_int &m, double *w,
                    double *z, lapack_int &ldz, lapack_int *isuppz, double *work,
                    lapack_int &lwork, lapack_int *iwork, lapack_int &liwork,
                    lapack_int &info);

void ARPACK(dsaupd)(arpack_int &ido, char &bmat, arpack_int &n, char *which,
                    arpack_int &nev, double &tol, double *resid, arpack_int &ncv,
                    double *v, arpack_int &ldv, arpack_int *iparam, arpack_int *ipntr,
                    double *workd, double *workl, arpack_int &lworkl, arpack_int &info);

void ARPACK(dseupd)(arpack_int &rvec, char &howmny, arpack_int *select, double *d, double *z,
                    arpack_int &ldz, double &sigma, char &bmat, arpack_int &n, char *which,
                    arpack_int &nev, double &tol, double *resid, arpack_int &ncv,
                    double *v, arpack_int &ldv, arpack_int *iparam, arpack_int *ipntr,
                    double *workd, double *workl, arpack_int &lworkl, arpack_int &info);

}


#endif

