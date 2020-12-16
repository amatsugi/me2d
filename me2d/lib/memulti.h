#ifndef _MEMULTI_H_
#define _MEMULTI_H_

#include "defs.h"
#include "mesolver.h"
#include "mesingle.h"

/*
 Multiple-Well ME class
  block diagonal: nwell * single-well-ME
  off-block-diagonal: U[kisom_i][kisom_j] = kisom
 */

class MEMulti: public MESolver {
public:
  int64_t nwell;
  MESingle *wells; int64_t wells_size;
  int64_t nkisom;
  double *kisom; int64_t kisom_size;
  int64_t *kisom_i; int64_t kisom_i_size;
  int64_t *kisom_j; int64_t kisom_j_size;
  int64_t *poswell; int64_t poswell_size;
  
public:
  MEMulti()
       : nwell(0), wells(nullptr), wells_size(0),
         nkisom(0), kisom(nullptr), kisom_size(0),
         kisom_i(nullptr), kisom_i_size(0),
         kisom_j(nullptr), kisom_j_size(0),
         poswell(nullptr), poswell_size(0) {}
  ~MEMulti() {
    delete_array(wells, wells_size);
    delete_array(kisom, kisom_size);
    delete_array(kisom_i, kisom_i_size);
    delete_array(kisom_j, kisom_j_size);
    delete_array(poswell, poswell_size);
  }
  
  int init(int64_t nwell);
  
  // well properties have to be set by calling member functions of wells[iwell]
  
  // post_set_wells: call after all wells are stored; set nsiz, poswell, ka, sa
  int post_set_wells();
  // set_kisom: call after post_set_wells
  int set_kisom(int64_t nkisom_, double *kisom_, int64_t *kisom_i_, int64_t *kisom_j_);
  // set_full_matrix: (if me_type == MWDens) call after set_kisom
  int set_full_matrix();
  // set_reactant: for inviter;  call after post_set_wells
  int set_reactant(int64_t reactant);
  
  void apply_rates() override;
  void invert_sign() override;
  void multiply_U(double *x, double *y) override;
  void multiply_U_semiquad(quad_float *x, quad_float *y) override;
  
  void diagonal_U(double *d) override;
  void reduced_banded_U(double *bm, int64_t red_bw) override;
  int64_t get_reduced_bwidth(double relb) override;

};


#endif

