#ifndef _MESOLVER_H_
#define _MESOLVER_H_

#include "defs.h"
#include <iostream>
#include <string>
#include <cmath>
#include <new>


class MESolver {
public:
  static int64_t max_array_size;
  static int64_t total_array_size;
  static int64_t total_array_size_used;
  enum me_type_t { Dens = 1, Band = 2, MW = 11, MWDens = Dens } me_type;
  enum solver_t { InvIter = 1, EigIter = 2, Eigen = 3, LinEq = 4 } solver;
  enum linear_solver_t { LsCHO = 1, LsLU = 2, LsLDLT = 3, LsCG = 11 } linear_solver;
  enum eigen_solver_t { EsDSYEVR = 1 } eigen_solver;
  enum precond_t { PcDiag = 1, PcBand = 2 } preconditioner;
  enum {
    ErrAlloc = -1, ErrMaxMem = -2,
    ErrSolver = -11, ErrMaxIter = -21, ErrLapack = -31,
    ErrCGMaxIter = -41, ErrArpack = -51,
  };
  int64_t InvIter_MaxIter;
  double InvIter_RTol;
  arpack_int Arpack_MaxIter;
  double Arpack_Tol;
  int64_t CG_MaxIter;
  double CG_RTol, CG_RTol_Semiquad;
  int64_t CG_guess_stored;
  int64_t nsiz, bwidth;
  double Pc_RelB;
  int64_t Pc_bwidth;
  int64_t verbose;
  double *U; int64_t U_size; // ME matrix (except for MW)
  double *ka; int64_t ka_size; // k(E)
  double *sa; int64_t sa_size; // symmetrize vector
  double *source; int64_t source_size; // source vector for specifing reactant in MW InvIter
private: // work arrays for linear solvers
  lapack_int *ipiv; int64_t ipiv_size;  // for LU and LDLT
  double *dwork_cg; int64_t dwork_cg_size;  // for CG
  quad_float *qwork_cg; int64_t qwork_cg_size;  // for CG semiquad
  double *dwork_pc; int64_t dwork_pc_size;  // for precond
  
public:
  MESolver()
       : me_type(Dens), solver(InvIter), linear_solver(LsCHO), eigen_solver(EsDSYEVR),
         preconditioner(PcBand),
         InvIter_MaxIter(100), InvIter_RTol(1e-6),
         Arpack_MaxIter(100), Arpack_Tol(1e-6),
         CG_MaxIter(2000), CG_RTol(1e-12), CG_RTol_Semiquad(-1.), CG_guess_stored(0),
         nsiz(0), bwidth(0),
         Pc_RelB(0.2), Pc_bwidth(100),
         verbose(0),
         U(nullptr), U_size(0), ka(nullptr), ka_size(0),
         sa(nullptr), sa_size(0), source(nullptr), source_size(0),
         ipiv(nullptr), ipiv_size(0),
         dwork_cg(nullptr), dwork_cg_size(0),
         qwork_cg(nullptr), qwork_cg_size(0),
         dwork_pc(nullptr), dwork_pc_size(0) {}
  virtual ~MESolver(){
    delete_array(U, U_size);
    delete_array(ka, ka_size);
    delete_array(sa, sa_size);
    delete_array(source, source_size);
    delete_array(ipiv, ipiv_size);
    delete_array(dwork_cg, dwork_cg_size);
    delete_array(qwork_cg, qwork_cg_size);
    delete_array(dwork_pc, dwork_pc_size);
  }
  
  static void show_solvers();
  int set_solver(std::string solver_str);
  int set_solver_option(std::string option_str);
  int solve(int64_t neig, double *vals, double *z);

  // virtual functions
  virtual void apply_rates() = 0; // apply ka (and kisom) to U
  virtual void invert_sign() = 0; // U = -U
  virtual void multiply_U(double *x, double *y) = 0; // y = Ux
  virtual void multiply_U_semiquad(quad_float *x, quad_float *y) = 0; // y = Ux

  // virtual functions for preconditioner
  virtual void diagonal_U(double *d) = 0; // d = diag(U)
  virtual void reduced_banded_U(double *bm, int64_t red_bw) = 0;
  virtual int64_t get_reduced_bwidth(double relb) = 0;
  
private:
  int inverse_iteration(double &val, double *ga);
  int eigen_iteration(int64_t neig, double *vals, double *z);
  int linear_equation(double &val, double *ga);
  
  int linear_solver_init();
  int linear_solver_solve(double *b);
  void linear_solver_clear();

  int pc_init();
  int pc_solve(double *b);
  int cg_init();
  int cg_solve(double *b);
  int cg_solve_semiquad(double *b, int x_as_guess=false);
  
  int eigen_solve(int64_t neig, double *vals, double *z);
  int eigen_solve_dsyevr(int64_t neig, double *vals, double *z);

protected:
  template<typename T> int allocate_array(T *&x, int64_t size) {
    int64_t required;
    double required_GB, max_GB;
    if (x != nullptr) {
      std::cout << "ERROR: allocate_array: already allocated (!=nullptr)" << std::endl;
      return ErrAlloc;
    }
    required = (total_array_size + sizeof(T)*size);
    if ((max_array_size > 0) && (required > max_array_size)) {
      required_GB = (double)required / 1024. / 1024. / 1024.;
      max_GB = (double)max_array_size / 1024. / 1024. / 1024.;
      std::cout << "ERROR: allocate_array: size limit exceeded: max_array_size="
        << max_array_size << " (" << max_GB << " GB)" << ", required="
          << required << " (" << required_GB << " GB)" << std::endl;
      return ErrMaxMem;
    }
    try { x = new T[size]; }
    catch (std::bad_alloc &err) {
      std::cout << "ERROR: allocate_array: failure to allocate " << err.what() << std::endl;
      delete[] x; x = nullptr; return ErrAlloc;
    }
    total_array_size += sizeof(T)*size;
    if (total_array_size > total_array_size_used) { total_array_size_used = total_array_size; }
    if (verbose > 1)
      { std::cout << "allocate_array:  total size = " << total_array_size << std::endl; }
    return 0;
  }
  
  template<typename T> void delete_array(T *&x, int64_t size) {
    if ((x == nullptr) || (size == 0)) { return; }
    delete[] x;
    x = nullptr;
    total_array_size -= sizeof(T) * size;
    if (verbose > 1)
      { std::cout << "delete_array:  total size = " << total_array_size << std::endl; }
    return;
  }
  
  template<typename T> void symmetrize_vector(T *x) {
    for (int64_t i=0; i<nsiz; i++) { x[i] /= (T)(sa[i]); }
  }
  template<typename T> void desymmetrize_vector(T *x) {
    for (int64_t i=0; i<nsiz; i++) { x[i] *= (T)(sa[i]); }
  }
  
  template<typename T> T nrm2(T *x) {
    T zero=0., one=1., scale=0., sumsq=1., absx;
    for (int64_t i=0; i<nsiz; i++) {
      absx = std::fabs(x[i]);
      if (scale < absx) {
        sumsq = one + sumsq * (scale/absx) * (scale/absx);
        scale = absx;
      }
      else if (absx > zero) {
        sumsq += (absx/scale) * (absx/scale);
      }
    }
    return scale * std::sqrt(sumsq);
  }
  
  template<typename T> T dot(T *a, T *b) {
    T dot=0.;
    for (int64_t i=0; i<nsiz; i++) { dot += a[i]*b[i]; }
    return dot;
  }
  template<typename T> T normalize_vector(T *x) {
    int64_t i;
    T sum=0.;
    for (i=0; i<nsiz; i++) { sum += x[i]; }
    for (i=0; i<nsiz; i++) { x[i] /= sum; }
    return sum;
  }
  
};



#endif

