#include "memulti.h"
#include "funcs.h"
#include <iostream>
#include <cmath>
#include <algorithm>
using std::cout;
using std::endl;


int MEMulti::init(int64_t nwell) {
  int res;
  int64_t iwell;
  this->nwell = nwell;
  me_type = MW;
  
  delete_array(wells, wells_size);
  if ((res = allocate_array(wells, nwell)) < 0) { return res; }
  wells_size = nwell;
  
  for (iwell=0; iwell<nwell; iwell++) { wells[iwell].verbose = verbose; }

  if (verbose > 0) { cout << "init: nwell = " << nwell << endl; }
  return 0;
}


int MEMulti::post_set_wells() {
  int res;
  int64_t iwell, i;
  delete_array(poswell, poswell_size);
  if ((res = allocate_array(poswell, nwell)) < 0) { return res; }
  poswell_size = nwell;
  
  nsiz = 0;
  for(iwell=0; iwell<nwell; iwell++) {
    nsiz += wells[iwell].nsiz;
    if (iwell == 0) { poswell[iwell] = 0; }
    else { poswell[iwell] = poswell[iwell-1] + wells[iwell-1].nsiz; }
  }
  if (verbose > 0) { cout << "post_set_wells: total nsiz = " << nsiz << endl; }
  
  delete_array(ka, ka_size);
  delete_array(sa, sa_size);
  if ((res = allocate_array(ka, nsiz)) < 0) { return res; }
  ka_size = nsiz;
  if ((res = allocate_array(sa, nsiz)) < 0) { return res; }
  sa_size = nsiz;
  
#pragma omp parallel for private(i) schedule(dynamic,1)
  for(iwell=0; iwell<nwell; iwell++) {
    for(i=0; i<wells[iwell].nsiz; i++) {
      ka[poswell[iwell]+i] = wells[iwell].ka[i];
      sa[poswell[iwell]+i] = wells[iwell].sa[i];
    }
  }
  
  return 0;
}

int MEMulti::set_kisom(int64_t nkisom_, double *kisom_,
                    int64_t *kisom_i_, int64_t *kisom_j_) {
  int res;
  int64_t i;
  nkisom = nkisom_;
  
  delete_array(kisom, kisom_size);
  delete_array(kisom_i, kisom_i_size);
  delete_array(kisom_j, kisom_j_size);
  if ((res = allocate_array(kisom, nkisom)) < 0) { return res; }
  kisom_size = nkisom;
  if ((res = allocate_array(kisom_i, nkisom)) < 0) { return res; }
  kisom_i_size = nkisom;
  if ((res = allocate_array(kisom_j, nkisom)) < 0) { return res; }
  kisom_j_size = nkisom;

#pragma omp parallel for
  for (i=0; i<nkisom; i++) {
    kisom[i] = kisom_[i];
    kisom_i[i] = kisom_i_[i];
    kisom_j[i] = kisom_j_[i];
  }
  if (verbose > 0) { cout << "set_kisom: nkisom = " << nkisom << endl; }
  return 0;
}


int MEMulti::set_full_matrix() {
  int res;
  int64_t iwell, i, j, posw, nsizw, ihigh, bwidthw, bp1w;
  
  if (verbose > 0) {
    cout << "set_full_matrix: full matrix size = (" << nsiz << "," << nsiz << "), "
      << fmt_size(sizeof(double)*nsiz*nsiz) << endl;
  }
  delete_array(U, U_size);
  if ((res = allocate_array(U, nsiz*nsiz)) < 0) { return res; }
  U_size = nsiz*nsiz;

#pragma omp parallel for
  for (i=0; i<U_size; i++) { U[i] = 0.; }
  
  for (iwell=0; iwell<nwell; iwell++) {
    nsizw = wells[iwell].nsiz;
    posw = poswell[iwell];
    if (wells[iwell].me_type == wells[iwell].Band) {
      bwidthw = wells[iwell].bwidth;
      bp1w = bwidthw + 1;
#pragma omp parallel for private(ihigh,j) schedule(dynamic,16)
      for(i=0; i<nsizw; i++) {
        ihigh = std::min(nsizw, i+bp1w);
        U[(posw+i)*nsiz + (posw+i)] = wells[iwell].U[bwidthw+i*bp1w];
        for (j=i+1; j<ihigh; j++) {
          U[(posw+i)*nsiz + (posw+j)] = wells[iwell].U[bwidthw+i-j+j*bp1w];
          U[(posw+j)*nsiz + (posw+i)] = wells[iwell].U[bwidthw+i-j+j*bp1w];
        }
      }
    }
    else if (wells[iwell].me_type == wells[iwell].Dens) {
#pragma omp parallel for private(j)
      for(i=0; i<nsizw; i++) {
        for(j=0; j<nsizw; j++) {
          U[(posw+i)*nsiz + (posw+j)] = wells[iwell].U[i*nsizw + j];
        }
      }
    }
  }
  return 0;
}


int MEMulti::set_reactant(int64_t reactant) {
  int64_t iwell, i;
  int res;
  if ((reactant < 0) || (reactant >= nwell)) { return 0; }
  if (solver != InvIter) {
    cout << "WARNING: set_reactant: reactant designation (only for InvIter) ignored." << endl;
  }
  delete_array(source, source_size);
  if ((res = allocate_array(source, nsiz)) < 0) { return res; }
  source_size = nsiz;
  
#pragma omp parallel for private(i) schedule(dynamic,1)
  for(iwell=0; iwell<nwell; iwell++) {
    if (reactant == iwell) { for(i=0; i<wells[iwell].nsiz; i++) { source[poswell[iwell]+i] = 1.; } }
    else { for(i=0; i<wells[iwell].nsiz; i++) { source[poswell[iwell]+i] = 0.; } }
  }
  return 0;
}


void MEMulti::apply_rates() {
  int64_t iwell, i;
  for (iwell=0; iwell<nwell; iwell++) { wells[iwell].apply_rates(); }
  
  if (U != nullptr) {
#pragma omp parallel for
    for (i=0; i<nsiz; i++) { U[i*nsiz + i] -= ka[i]; }
#pragma omp parallel for
    for (i=0; i<nkisom; i++) {
      U[kisom_i[i]*nsiz+kisom_j[i]] = kisom[i];
      U[kisom_j[i]*nsiz+kisom_i[i]] = kisom[i];
    }
  }
  return;
}


void MEMulti::invert_sign() {
  int64_t iwell, i;
  for (iwell=0; iwell<nwell; iwell++) { wells[iwell].invert_sign(); }
  for (i=0; i<nkisom; i++) { kisom[i] = - kisom[i]; }
  for (i=0; i<nsiz; i++) { ka[i] = - ka[i]; }

  if (U != nullptr) {
#pragma omp parallel for
    for (i=0; i<U_size; i++) { U[i] = - U[i]; }
  }
}


void MEMulti::multiply_U(double *x, double *y) {
  int64_t iwell, i;
  
#pragma omp parallel for schedule(dynamic,1)
  for (iwell=0; iwell<nwell; iwell++) {
    wells[iwell].multiply_U(x+poswell[iwell], y+poswell[iwell]);
  }
  for (i=0; i<nkisom; i++) {
    y[kisom_i[i]] += kisom[i] * x[kisom_j[i]];
    y[kisom_j[i]] += kisom[i] * x[kisom_i[i]];
  }
}

void MEMulti::multiply_U_semiquad(quad_float *x, quad_float *y) {
  int64_t iwell, i;
#pragma omp parallel for schedule(dynamic,1)
  for (iwell=0; iwell<nwell; iwell++) {
    wells[iwell].multiply_U_semiquad(x+poswell[iwell], y+poswell[iwell]);
  }
  for (i=0; i<nkisom; i++) {
    y[kisom_i[i]] += kisom[i] * x[kisom_j[i]];
    y[kisom_j[i]] += kisom[i] * x[kisom_i[i]];
  }
}


void MEMulti::diagonal_U(double *d) {
  int64_t iwell;
#pragma omp parallel for schedule(dynamic,1)
  for (iwell=0; iwell<nwell; iwell++) {
    wells[iwell].diagonal_U(d+poswell[iwell]);
  }
}

void MEMulti::reduced_banded_U(double *bm, int64_t red_bw) {
  int64_t iwell, i, j, posbw, red_bp1=red_bw+1;
  
#pragma omp parallel for private(j)
  for (i=0; i<nsiz; i++) {
    for (j=0; j<red_bp1; j++) { bm[i*red_bp1 + j] = 0.; }
  }
  for (iwell=0; iwell<nwell; iwell++) {
    posbw = poswell[iwell] * red_bp1;
    wells[iwell].reduced_banded_U(bm+posbw, red_bw);
  }
}

int64_t MEMulti::get_reduced_bwidth(double relb) {
  int64_t iwell, maxsiz=0, red_bwidth=0;
  double average=0.;
  
  if (relb < 0.) { red_bwidth = (int64_t)fabs(relb); }
  else {
    for (iwell=0; iwell<nwell; iwell++) {
      if (wells[iwell].me_type == wells[iwell].Band) { average += (double)wells[iwell].bwidth; }
      else if (wells[iwell].me_type == wells[iwell].Dens) { average += (double)wells[iwell].nsiz; }
    }
    if (nwell > 0) { average = average / (double)nwell; }
    red_bwidth = (int64_t)(average * relb);
  }
  
  for (iwell=0; iwell<nwell; iwell++) {
    if ((wells[iwell].me_type == wells[iwell].Band) && (wells[iwell].bwidth > maxsiz))
      { maxsiz = wells[iwell].bwidth; }
    else if ((wells[iwell].me_type == wells[iwell].Dens) && (wells[iwell].nsiz > maxsiz))
      { maxsiz = wells[iwell].nsiz; }
  }
  if (red_bwidth > maxsiz) { red_bwidth = maxsiz; }
  return red_bwidth;
}

