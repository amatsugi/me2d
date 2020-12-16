#ifndef _MESINGLE_H_
#define _MESINGLE_H_

#include "defs.h"
#include "mesolver.h"
#include <cmath>

/*
 Single-Well ME class
 dense U: full nsiz*nsiz matrix
 banded U: nsiz*bp1 (bp1=bwidth+1), lapack symmetric (upper) band storage
 access to Uorig[i*nsiz + j] (i<=j):  U[bwidth+i-j+j*bp1]  (U[bwidth+i*bp1] for i=j)
 access to Uorig[i*nsiz + j] (i>j) :  U[bwidth+j-i+i*bp1] * sa[i]*sa[i]/sa[j]/sa[j]
 */

class MESingle: public MESolver {
public:
  enum {
    ErrChkSum = -121, ErrChkSym = 122
  };
  double TOL_NormSum, TOL_SymBase, TOL_Sym;
public:
  MESingle() : TOL_NormSum(1e-6), TOL_SymBase(1e-128), TOL_Sym(1e-6) {}
  ~MESingle() {}
  
  int init_dens(int64_t nsiz);
  int init_band(int64_t nsiz, double *ea, double *Ja, double y_e, double y_J,
                double *ainv_ea, double *ainv_Ja, double bandpcrit);
  
  void set_ka(double *ka_);
  void set_sa(double *rhoa, double *Ea, double kbT);
  void set_prob(double *ea, double *Ja, double y_e, double y_J,
                double *ainv_ea, double *ainv_Ja, int64_t ptype);
  int normalize();
  int symmetrize();
  void apply_collision_rate(double ZM);
  
  void apply_rates() override;
  void invert_sign() override;
  void multiply_U(double *x, double *y) override;
  void multiply_U_semiquad(quad_float *x, quad_float *y) override;
  
  void diagonal_U(double *d) override;
  void reduced_banded_U(double *bm, int64_t red_bw) override;
  int64_t get_reduced_bwidth(double relb) override;

private:
  int check_sum();
  int check_symmetry();
  
  // probability functions defined here
  inline double probfunc1d(double de, double y_e, double ainv_e) {
    double exponent = 0.;
    if (y_e == 1.) { exponent -= std::fabs(de) * ainv_e; }
    else if (y_e == 0.5) { exponent -= std::sqrt(fabs(de) * ainv_e); }
    else { exponent -= std::pow(std::fabs(de) * ainv_e, y_e); }
    return std::exp(exponent);
  }
  
  inline double probfunc2d(double de, double dJ, double y_e, double y_J,
                           double ainv_e, double ainv_J) {
    double exponent = 0.;
    if (y_e == 1.) { exponent -= std::fabs(de) * ainv_e; }
    else if (y_e == 0.5) { exponent -= std::sqrt(fabs(de) * ainv_e); }
    else { exponent -= std::pow(std::fabs(de) * ainv_e, y_e); }
    if (y_J == 1.) { exponent -= std::fabs(dJ) * ainv_J; }
    else if (y_J == 0.5) { exponent -= std::sqrt(std::fabs(dJ) * ainv_J); }
    else { exponent -= std::pow(std::fabs(dJ) * ainv_J, y_J); }
    return std::exp(exponent);
  }
};


#endif

