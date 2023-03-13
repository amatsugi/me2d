#include "mesingle.h"
#include "extlib.h"
#include "funcs.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using std::cout;
using std::endl;


int MESingle::init_dens(int64_t nsiz) {
  int res;
  this->nsiz = nsiz;
  me_type = Dens;
  
  if (verbose > 0) {
    cout << "init_dens: matrix size = (" << nsiz << "," << nsiz << "), "
      << fmt_size(sizeof(double)*nsiz*nsiz) << endl;
  }
  delete_array(U, U_size);
  delete_array(ka, ka_size);
  delete_array(sa, sa_size);
  if ((res = allocate_array(U, nsiz*nsiz)) < 0) { return res; }
  U_size = nsiz*nsiz;
  if ((res = allocate_array(ka, nsiz)) < 0) { return res; }
  ka_size = nsiz;
  if ((res = allocate_array(sa, nsiz)) < 0) { return res; }
  sa_size = nsiz;
  return 0;
}

int MESingle::init_band(int64_t nsiz, double *ea, double *Ja, double y_e, double y_J,
                  double *ainv_ea, double *ainv_Ja, double bandpcrit) {
  int64_t i, j, checksiz;
  int res, nthread=1, ithread=0;
  double prob, sumfs;
  this->nsiz = nsiz;
  me_type = Band;
#ifdef _OPENMP
  nthread = omp_get_max_threads();
#endif
  
  if (verbose > 0) { cout << "init_band: calculate bandwidth..." << endl; }
  // temporarily use U to store thread-dependent prob.
  if ((res = allocate_array(U, nsiz*nthread)) < 0) { return res; }
  U_size = nsiz*nthread;
  
  bwidth = 1;
#pragma omp parallel for private(ithread,checksiz,sumfs,j) schedule(dynamic,16)
  for (i=0; i<nsiz; i++) {
#pragma omp critical
    checksiz = i + bwidth;
    if (checksiz >= nsiz - 1) { continue; }
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#endif
    sumfs = 0.;
    for (j=0; j<nsiz; j++) {
      if (Ja == nullptr) { prob = probfunc1d(ea[i] - ea[j], y_e, ainv_ea[j]); }
      else { prob = probfunc2d(ea[i] - ea[j], Ja[i] - Ja[j], y_e, y_J, ainv_ea[j], ainv_Ja[j]); }
      U[nsiz*ithread + j] = prob;
      sumfs += prob;
    }
    sumfs *= bandpcrit;
#pragma omp critical
    for (j=nsiz-1; j>=i+bwidth; j--) {
      if (U[nsiz*ithread + j] > sumfs) { bwidth = j - i + 1; break; }
    }
  }
  if (verbose > 0) {
    cout << "init_band: bandwidth = " << bwidth << endl;
    cout << "init_band: matrix size = (" << nsiz << "," << bwidth+1 << "), "
      << fmt_size(sizeof(double)*nsiz*(bwidth+1)) << endl;
  }
  delete_array(U, U_size);
  delete_array(sa, sa_size);
  delete_array(ka, ka_size);
  if ((res = allocate_array(U, nsiz*(bwidth+1))) < 0) { return res; }
  U_size = nsiz*(bwidth+1);
  if ((res = allocate_array(ka, nsiz)) < 0) { return res; }
  ka_size = nsiz;
  if ((res = allocate_array(sa, nsiz)) < 0) { return res; }
  sa_size = nsiz;
  return 0;
}


void MESingle::set_ka(double *ka_) {
  for (int64_t i=0; i<nsiz; i++) { ka[i] = ka_[i]; }
}

void MESingle::set_sa(double *rhoa, double *Ea, double kbT) {
  for (int64_t i=0; i<nsiz; i++) { sa[i] = std::sqrt(rhoa[i] * std::exp( - Ea[i] / kbT )); }
}

void MESingle::set_prob(double *ea, double *Ja, double y_e, double y_J,
                  double *ainv_ea, double *ainv_Ja, int64_t ptype) {
  int64_t i, j, ihigh, bp1=bwidth+1;
  double prob;
  
#pragma omp parallel for private(ihigh,prob,j) schedule(dynamic,16)
  for (i=0; i<nsiz; i++) {
    if (me_type == Band) { ihigh = std::min(nsiz, i+bp1); }
    else { ihigh = nsiz; }
    for (j=i; j<ihigh; j++) {
      if (Ja == nullptr) { prob = probfunc1d(ea[i] - ea[j], y_e, ainv_ea[j]); }
      else { prob = probfunc2d(ea[i] - ea[j], Ja[i] - Ja[j], y_e, y_J, ainv_ea[j], ainv_Ja[j]); }
      
      // ptype: 0 = symmetric, -1 = down-prob given, 1 = up-prob given.
      if (me_type == Band) {
        if (ptype == 0) { U[bwidth+i-j+j*bp1] = prob * sa[i] / sa[j]; }
        else if (ptype < 0) { U[bwidth+i-j+j*bp1] = prob; }
        else { U[bwidth+i-j+j*bp1] = prob * (sa[i]*sa[i]) / (sa[j]*sa[j]); }
      }
      else if (me_type == Dens) {
        if (ptype == 0) {
          U[i*nsiz+j] = prob * sa[i] / sa[j];
          if (j > i) { U[i+j*nsiz] = prob * sa[j] / sa[i]; }
        }
        else if (ptype < 0) {
          U[i*nsiz+j] = prob;
          if (j > i) { U[i+j*nsiz] = prob * (sa[j]*sa[j]) / (sa[i]*sa[i]); }
        }
        else {
          U[i*nsiz+j] = prob * (sa[i]*sa[i]) / (sa[j]*sa[j]);
          if (j > i) { U[i+j*nsiz] = prob; }
        }
      }
    }
  }
  return;
}

int MESingle::normalize() {
  int res;
  int64_t i, j, ip1, ilow, ihigh, bp1=bwidth+1;
  double downUsum, upPsum;
  double *norm = nullptr;
  if ((res = allocate_array(norm, nsiz)) < 0) { return res; }
  
  // set normalization coeffs.
#pragma omp parallel for private(ilow,downUsum,j) schedule(dynamic,16)
  for (i=0; i<nsiz; i++) {
    downUsum = 0.;
    if (me_type == Band) {
      if (i > bwidth) { ilow = i - bwidth; }
      else { ilow = 0; }
      for (j=ilow; j<i+1; j++) { downUsum += U[bwidth+j-i+i*bp1]; }
    }
    else if (me_type == Dens) { for (j=0; j<i+1; j++) { downUsum += U[j*nsiz + i]; } }
    norm[i] = 1. / downUsum;
  }
  // top to down scheme, not parallelized
  for (ip1=nsiz; ip1>0; ip1--) {
    i = ip1 - 1;
    upPsum =  0.;
    if (me_type == Band) {
      ihigh = std::min(nsiz, i+bp1);
      for (j=ihigh-1; j>=i+1; j--) { upPsum += U[bwidth+i-j+j*bp1] * norm[j]; }
      norm[i] *= (sa[i]*sa[i] - upPsum);
    }
    else if (me_type == Dens) {
      for (j=nsiz-1; j>=i+1; j--) { upPsum += U[j*nsiz + i] * norm[j]; }
      norm[i] *= (1. - upPsum);
    }
  }
  if (me_type == Band) { for (i=0; i<nsiz; i++) { norm[i] /= (sa[i]*sa[i]); } }
  // normalize
#pragma omp parallel for private(ihigh,j) schedule(dynamic,16)
  for (i=0; i<nsiz; i++) {
    if (me_type == Band) {
      ihigh = std::min(nsiz, i+bp1);
      for (j=i; j<ihigh; j++) { U[bwidth+i-j+j*bp1] *= norm[j]; }
    }
    else if (me_type == Dens) {
      for (j=i; j<nsiz; j++) { U[i*nsiz+j] *= norm[j]; }
      for (j=i+1; j<nsiz; j++) { U[j*nsiz+i] *= norm[j]; }
    }
  }
  delete_array(norm, nsiz);
  return check_sum();
}

int MESingle::check_sum() {
  int res;
  int64_t i, j, ilow, ihigh, bp1=bwidth+1;
  double Usum, maxdev=0.;
  double *sums = nullptr;
  if ((res = allocate_array(sums, nsiz)) < 0) { return res; }
  
#pragma omp parallel for private(Usum,ilow,ihigh,j) schedule(dynamic,16)
  for (i=0; i<nsiz; i++) {
    Usum = 0;
    if (me_type == Band) {
      if (i > bwidth) { ilow = i - bwidth; }
      else { ilow = 0; }
      ihigh = std::min(nsiz, i+bp1);
      for (j=ilow; j<=i; j++) { Usum += U[bwidth+j-i+i*bp1]; }
      for (j=i+1; j<ihigh; j++) { Usum += U[bwidth+i-j+j*bp1] * sa[j]*sa[j] / (sa[i]*sa[i]); }
    }
    else if (me_type == Dens) {
      for (j=0; j<nsiz; j++) { Usum += U[i + j*nsiz]; }
    }
    sums[i] = Usum;
  }
  for (i=0; i<nsiz; i++) {
    if (fabs(sums[i]-1.) > fabs(maxdev)) { maxdev = sums[i]-1.; }
  }
  if (maxdev > TOL_NormSum) {
    cout << "ERROR: check_sum: normalized sum != 1 (max. deviation = " << maxdev << ")" << endl;
    delete_array(sums, nsiz);
    return ErrChkSum;
  } else if (maxdev > TOL_NormSum_Warn) {
    cout << "WARNING: check_sum: max. deviation in normalized sums = " << maxdev << endl;
  }
  delete_array(sums, nsiz);
  return 0;
}


int MESingle::symmetrize() {
  int64_t i, j, ihigh, bp1=bwidth+1;
  if (me_type == Band) {
#pragma omp parallel for private(ihigh,j) schedule(dynamic,16)
    for (i=0; i<nsiz; i++) {
    ihigh = std::min(nsiz, i+bp1);
      for (j=i; j<ihigh; j++) { U[bwidth+i-j+j*bp1] *= sa[j] / sa[i]; }
    }
  }
  else if (me_type == Dens) {
#pragma omp parallel for private(j)
    for (i=0; i<nsiz; i++) {
      for (j=0; j<nsiz; j++) { U[i*nsiz + j] *= sa[j] / sa[i]; }
    }
  }
  return check_symmetry();
}

int MESingle::check_symmetry() {
  int64_t i, j, iasym = 0;
  double Uij, Uji;
  if (me_type == Band) { return 0; } // symmetric storage
  else if (me_type == Dens) {
#pragma omp parallel for private(Uij,Uji,j) schedule(dynamic,16)
    for (i=0; i<nsiz; i++) {
      for (j=i+1; j<nsiz; j++) {
        Uij = U[i*nsiz + j];
        Uji = U[i + j*nsiz];
        if ((fabs(Uij) + fabs(Uji) > TOL_SymBase) && (fabs(Uij/Uji - 1.) > TOL_Sym)) {
#pragma omp atomic
          iasym++;
        }
      }
    }
  }
  if (iasym != 0) {
    cout << "ERROR: check_symmetry: asymmetry detected." << endl;
    return ErrChkSym;
  }
  return 0;
}


void MESingle::apply_collision_rate(double ZM) {
  int64_t i, j, ihigh, bp1=bwidth+1;
  if (me_type == Band) {
#pragma omp parallel for private(ihigh,j) schedule(dynamic,16)
    for (i=0; i<nsiz; i++) {
      ihigh = std::min(nsiz, i+bp1);
      for (j=i; j<ihigh; j++) { U[bwidth+i-j+j*bp1] *= ZM; }
      U[bwidth+i*bp1] -= ZM;
    }
  }
  else if (me_type == Dens) {
#pragma omp parallel for private(j)
    for (i=0; i<nsiz; i++) {
      for (j=0; j<nsiz; j++) { U[i*nsiz + j] *= ZM; }
      U[i*nsiz + i] -= ZM;
    }
  }
  return;
}


void MESingle::apply_rates() {
  int64_t i, bp1=bwidth+1;
  if (me_type == Band) { for (i=0; i<nsiz; i++) { U[bwidth+i*bp1] -= ka[i]; } }
  else if (me_type == Dens) { for (i=0; i<nsiz; i++) { U[i*nsiz + i] -= ka[i]; } }
  return;
}

void MESingle::invert_sign() {
  int64_t i, j, ihigh, bp1=bwidth+1;
  if (me_type == Band) {
#pragma omp parallel for private(ihigh,j) schedule(dynamic,16)
    for (i=0; i<nsiz; i++) {
      ka[i] = - ka[i];
      ihigh = std::min(nsiz, i+bp1);
      for (j=i; j<ihigh; j++) { U[bwidth+i-j+j*bp1] = -U[bwidth+i-j+j*bp1]; }
    }
  }
  else if (me_type == Dens) {
#pragma omp parallel for private(j)
    for (i=0; i<nsiz; i++) {
      ka[i] = - ka[i];
      for (j=0; j<nsiz; j++) { U[i*nsiz + j] = -U[i*nsiz + j]; }
    }
  }
  return;
}


void MESingle::multiply_U(double *x, double *y) {
  int64_t i;
  blas_int N=(blas_int)nsiz, K=(blas_int)bwidth, lda=K+1, inc=1;
  double alpha = 1., beta = 0.;
  char uplo='U';
  if (me_type == Band)
    { BLAS(dsbmv)(uplo, N, K, alpha, U, lda, x, inc, beta, y, inc); }
  else if (me_type == Dens)
    { BLAS(dsymv)(uplo, N, alpha, U, N, x, inc, beta, y, inc); }
  for (i=0; i<nsiz; i++) { y[i] -= ka[i] * x[i]; }
  return;
}


void MESingle::multiply_U_semiquad(quad_float *x, quad_float *y) {
  int64_t i, j, ilow, bp1=bwidth+1;;
  quad_float tmp1, tmp2;
  if (me_type == Band) {
    // not parallelized
    for (i=0; i<nsiz; i++) { y[i] = 0.; }
    for (i=0; i<nsiz; i++) {
      if (i > bwidth) { ilow = i - bwidth; }
      else { ilow = 0; }
      tmp1 = 0.;
      for (j=ilow; j<i; j++) {
        tmp2 = U[bwidth+j-i+i*bp1];
        y[j] += tmp2 * x[i];
        tmp1 += tmp2 * x[j];
      }
      y[i] += U[bwidth+i*bp1] * x[i] + tmp1;
    }
  }
  else if (me_type == Dens) {
#pragma omp parallel for private(j)
    for (i=0; i<nsiz; i++) {
      y[i] = 0.;
      for (j=0; j<nsiz; j++) { y[i] += U[i*nsiz+j] * x[j]; }
    }
  }
  for (i=0; i<nsiz; i++) { y[i] -= ka[i] * x[i]; }
  return;
}


void MESingle::diagonal_U(double *d) {
  int64_t i, bp1=bwidth+1;
  if (me_type == Band) { for (i=0; i<nsiz; i++) { d[i] = U[bwidth+i*bp1] - ka[i]; } }
  else if (me_type == Dens) { for (i=0; i<nsiz; i++) { d[i] = U[i*nsiz + i] - ka[i]; } }
  return;
}

void MESingle::reduced_banded_U(double *bm, int64_t red_bw) {
  int64_t i, j, ihigh, bp1=bwidth+1, red_bp1=red_bw+1;
  if (me_type == Band) {
#pragma omp parallel for private(ihigh,j) schedule(dynamic,16)
    for (i=0; i<nsiz; i++) {
      ihigh = std::min(nsiz, i+red_bp1);
      for (j=i; j<ihigh; j++) {
        if ((j-i) < bp1) { bm[red_bw+i-j+j*red_bp1] = U[bwidth+i-j+j*bp1]; }
        else { bm[red_bw+i-j+j*red_bp1] = 0; }
      }
      bm[red_bw+i*red_bp1] -= ka[i];
    }
  }
  else if (me_type == Dens) {
#pragma omp parallel for private(ihigh,j) schedule(dynamic,16)
    for (i=0; i<nsiz; i++) {
      ihigh = std::min(nsiz, i+red_bp1);
      for (j=i; j<ihigh; j++) {
        if (j < nsiz) { bm[red_bw+i-j+j*red_bp1] = U[i*nsiz+j]; }
        else { bm[red_bw+i-j+j*red_bp1] = 0; }
      }
      bm[red_bw+i*red_bp1] -= ka[i];
    }
  }
  return;
}

int64_t MESingle::get_reduced_bwidth(double relb) {
  int64_t red_bwidth=0;
  if (relb < 0.) { red_bwidth = (int64_t)fabs(relb); }
  else {
    if (me_type == Band) { red_bwidth = (int64_t)((double)bwidth * relb); }
    else if (me_type == Dens) { red_bwidth = (int64_t)((double)nsiz * relb); }
  }
  if ((me_type == Band) && (red_bwidth > bwidth)) { red_bwidth = bwidth; }
  else if ((me_type == Dens) && (red_bwidth > nsiz)) { red_bwidth = nsiz; }
  return red_bwidth;
}



