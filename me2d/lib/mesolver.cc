#include "mesolver.h"
#include "funcs.h"
#include "extlib.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
using std::cout;
using std::endl;

int64_t MESolver::max_array_size = 0;
int64_t MESolver::total_array_size = 0;
int64_t MESolver::total_array_size_used = 0;


void MESolver::show_solvers() {
  cout <<
    R"(solver = "solver[,linear_or_eigen_solver[,option=value,...]]"
  solver
    - InvIter: inverse iteration with linear_solver
    - EigIter: iterative eigen solution (ARPACK, invert-mode with linear_solver)
    - Eigen: direct eigen solution with eigen_solver
  linear_solver
    - CHO: Cholesky decomposition (LAPACK DPOTR[FS] or DPBTR[FS])
    - LU: LU decomposition (LAPACK DGETR[FS]; not supported for banded matrix)
    - LDLT: LDLT decomposition (LAPACK DSYTR[FS]; not supported for banded matrix)
    - CG: conjugate gradient
    (default: CHO for single-well ME; CG for multiple-well ME)
  eigen_solver
    - DSYEVR: LAPACK DSYEVR (default)
  options
    for InvIter solver:
      - maxit: maximum number of iterations (default:100)
      - rtol: relative tolerance for convergence (default:1e-6)
    for EigIter solver:
      - maxit: maximum number of iterations (default:100)
      - rtol: relative tolerance for convergence (default:1e-6)
    for CG linear_solver:
      - CG_maxit: maximum number of iterations (default:2000)
      - CG_rtol: relative tolerance for convergence (default:1e-12)
      - CG_rtolsq: relative tolerance for convergence in the semiquad scheme
                   (partially use quadruple-precision-float) (default:-1)
          CG_rtolsq < 0: entirely use double float for CG
          CG_rtol >= 0 and CG_rtolsq >= 0: refine double-float solution with
                                           the semiquad scheme if needed
          CG_rtol < 0 and CG_rtolsq >= 0: entirely use the semiquad scheme
      - Pc_relb: relative bandwidth used for preconditioner (default:0.2)
)" << endl;
}

int MESolver::set_solver(std::string solver_str) {
  int res;
  int64_t i, n;
  std::string word;
  std::vector<std::string> words = split_string(solver_str, ",/ ");
  std::ostringstream oss;
  n = words.size();
  
  if (me_type == Dens) { oss << "dense matrix"; }
  else if (me_type == Band) { oss << "banded matrix"; }
  else if (me_type == MW) { oss << "multiple-well"; }
  else if (me_type == MWDens) { oss << "multiple-well (use full matrix)"; }
  
  if (n == 0) { // default solver = InvIter
    solver = InvIter; oss << ", inverse iteration";
  }
  else { // words[0] specifies solver
    word = words[0];
    if (keyeq(word, "InvIter")) { solver = InvIter; oss << ", inverse iteration"; }
    else if (keyeq(word, "EigIter")) { solver = EigIter; oss << ", iterative eigen solution (ARPACK)"; }
    else if (keyeq(word, "Eigen")) { solver = Eigen; oss << ", eigen solution";}
    else {
      cout << "ERROR: set_solver: invalid solver \"" + word + "\"" << endl;
      return ErrSolver;
    }
  }
  
  if (n <= 1) { // use default linear/eigen solvers 
    if ((solver == InvIter) || (solver == EigIter)) {
      if (me_type == MW) { linear_solver = LsCG; oss << " with conjugate gradient"; }
      else { linear_solver = LsCHO; oss << " with Cholesky decomposition"; }
    }
    else if (solver == Eigen) { eigen_solver = EsDSYEVR; oss << " with DSYEVR"; }
  }
  else { // words[1] specifies linear/eigen solvers
    word = words[1];
    if ((solver == InvIter) || (solver == EigIter)) {
      if (keyeq(word, "CHO")) { linear_solver = LsCHO; oss << " with Cholesky decomposition"; }
      else if (keyeq(word, "LU")) { linear_solver = LsLU; oss << " with LU decomposition"; }
      else if (keyeq(word, "LDLT")) { linear_solver = LsLDLT; oss << " with LDLT decomposition"; }
      else if (keyeq(word, "CG")) { linear_solver = LsCG; oss << " with conjugate gradient"; }
      else {
        cout << "ERROR: set_solver: invalid linear solver \"" + word + "\"" << endl;
        return ErrSolver;
      }
    }
    else if (solver == Eigen) {
      if (keyeq(word, "DSYEVR")) { eigen_solver = EsDSYEVR; oss << " with DSYEVR"; }
      else {
        cout << "ERROR: set_solver: invalid eigen solver \"" + word + "\"" << endl;
        return ErrSolver; }
    }
  }
  
  // check if solver implemented
  if (me_type == Band) {
    if ((solver == InvIter) || (solver == EigIter)) {
      if (linear_solver == LsLU) {
        cout << "ERROR: set_solver: banded LU solver not implemented" << endl;
        return ErrSolver;
      }
      else if (linear_solver == LsLDLT) {
        cout << "ERROR: set_solver: banded LDLT solver not implemented" << endl;
        return ErrSolver;
      }
    }
    else if (solver == Eigen) {
      cout << "ERROR: set_solver: banded eigen solver not implemented" << endl;
      return ErrSolver;
    }
  }
  
  if (verbose > 0) { cout << "set_solver: " << oss.str() << endl; }
  
  // check if full matrix is needed for MW
  if (me_type == MW) {
    if (((solver == InvIter) || (solver == EigIter)) && (linear_solver == LsCG)) {
      if (verbose > 0) { cout << "set_solver: the solver does not require full matrix." << endl; }
    }
    else {
      me_type = MWDens;
      if (verbose > 0) { cout << "set_solver: the solver requires full matrix." << endl; }
    }
  }
  
  if (n > 2) { // set options
    for (i=2; i<n; i++) {
      if ((res = set_solver_option(words[i])) < 0) { return res; }
    }
  }
  return 0;
}

int MESolver::set_solver_option(std::string option_str) {
  std::string key, val;
  std::vector<std::string> option = split_string(option_str, "=");
  std::ostringstream oss;
  
  try {
    if (option.size() == 2) {
      key = option[0]; val = option[1];
      if (solver == InvIter) {
        if (keyeq(key, "maxit"))
          { InvIter_MaxIter = stoi(val); oss << "InvIter_MaxIter set to " << InvIter_MaxIter; }
        else if (keyeq(key, "rtol"))
          { InvIter_RTol = stod(val); oss << "InvIter_RTol set to " << InvIter_RTol; }
      }
      if (solver == EigIter) {
        if (keyeq(key, "maxit"))
          { Arpack_MaxIter = stoi(val); oss << "Arpack_MaxIter set to " << Arpack_MaxIter; }
        else if (keyeq(key, "rtol"))
          { Arpack_Tol = stod(val); oss << "Arpack_Tol set to " << Arpack_Tol; }
      }
      if (((solver == InvIter) || (solver == EigIter)) && (linear_solver == LsCG)) {
        if (keyeq(key, "CG_maxit"))
          { CG_MaxIter = stoi(val); oss << "CG_MaxIter set to " << CG_MaxIter; }
        else if (keyeq(key, "CG_rtol"))
          { CG_RTol = stod(val); oss << "CG_RTol set to " << CG_RTol; }
        else if (keyeq(key, "CG_rtolsq"))
          { CG_RTol_Semiquad = stod(val); oss << "CG_RTol_Semiquad set to " << CG_RTol_Semiquad; }
        else if (keyeq(key, "Pc_relb")) {
          Pc_RelB = stod(val);
          oss << "Pc_RelB set to " << Pc_RelB;
          if (Pc_RelB == 0.) {
            preconditioner = PcDiag;
            oss << " (use diagonal preconditioner)";
          }
          else { preconditioner = PcBand; }
        }
      }
    }
  }
  catch (std::invalid_argument &err) {
    cout << "ERROR: set_solver_option: option " << key << " has invalid argument ("
      << err.what() << ")" << endl;
    return ErrSolver;
  }
  catch (std::out_of_range &err) {
    cout << "ERROR: set_solver_option: option " << key << " is out of range ("
      << err.what() << ")"  << endl;
    return ErrSolver;
  }
  
  if (oss.str().empty()) {
    cout << "WARNING: set_solver_option: option ignored: " << option_str << endl;
  }
  else if (verbose > 0) {
    cout << "set_solver_option: option " << oss.str() << endl;
  }
  return 0;
}


int MESolver::solve(int64_t neig, double *vals, double *z) {
  int res=0;
  if ((solver == InvIter) && (neig != 1)) {
    cout << "ERROR: solve: neig must be 1 for inverse iteration solver" << endl;
    return ErrSolver;
  }
  if (solver == InvIter) { res = inverse_iteration(vals[0], z); }
  else if (solver == EigIter) { res = eigen_iteration(neig, vals, z); }
  else if (solver == Eigen) { res = eigen_solve(neig, vals, z); }
  else {
    cout << "ERROR: solve: solver " << solver << " not implemented." << endl;
    return ErrSolver;
  }
  if (res < 0) { return res; }
  return 0;
}



// ---- iterative solvers ----


int MESolver::inverse_iteration(double &val, double *b) {
  int64_t i, it=0;
  int res;
  double rtol, sumb, relnorm, relval, val0=0.;
  double *b0 = nullptr;
  
  if (InvIter_RTol >= std::numeric_limits<double>::epsilon()) { rtol = InvIter_RTol; }
  else { rtol = std::numeric_limits<double>::epsilon(); }
  if (verbose > 0) { cout << "inverse_iteration: rtol = " << rtol << endl; }
  
  
  symmetrize_vector(b); normalize_vector(b);  // symmetrize initial guess vector
  invert_sign(); // invert sign to be positive definite
  
  if ((res = linear_solver_init()) < 0) { return res; }
  if ((res = allocate_array(b0, nsiz)) < 0) { return res; }
  
  it = 0;
  while (true) {
    it++;
    if (it > InvIter_MaxIter) {
      cout << "ERROR: inverse_iteration: max iterations exceeded." << endl;
      delete_array(b0, nsiz); return ErrMaxIter;
    }
    
    if (it == 1) { val0 = 0.; }
    else { val0 = val; }
    for (i=0; i<nsiz; i++) { b0[i] = b[i]; }
    if (source != nullptr) { for (i=0; i<nsiz; i++) { b[i] = b[i] * source[i]; } }
    
    res = linear_solver_solve(b);
    if (res < 0) { delete_array(b0, nsiz); return res; }
    val = - 1. / normalize_vector(b);
    
    for (i=0; i<nsiz; i++) { b0[i] = b[i] - b0[i]; } // residual
    relnorm = nrm2(b0) / nrm2(b);
    relval = std::fabs((val - val0) / val);
    if (verbose > 0) {
      cout << "inverse_iteration: it = " << it << ", relative residuals = "
        << relval << ", " << relnorm << endl;
    }
    if ((relval < rtol) && (relnorm < rtol)) { break; }
  }
  delete_array(b0, nsiz);
  linear_solver_clear();
  
  desymmetrize_vector(b); normalize_vector(b);
  if (source != nullptr) { // normarize with population ratio
    sumb = 0.;
    for (i=0; i<nsiz; i++) { sumb += b[i] * source[i]; }
    val *= sumb;
  }
  invert_sign();
  if (verbose > 0)
    { cout << "inverse_iteration: solved after " << it <<" iterations, val = " << val << endl; }
  return 0;
}


int MESolver::eigen_iteration(int64_t neig, double *vals, double *z) {
  int64_t i, lwork_all;
  int res;
  char bmat='I', which[]="LM", howmny='A';
  arpack_int ido=0, N=(arpack_int)nsiz, nev=(arpack_int)neig, ncv, ldv=N;
  arpack_int iparam[11], ipntr[11], lworkl, info;
  arpack_int rvec=1, ldz=N;
  double tol=Arpack_Tol, sigma=0.;
  double *work=nullptr; // for resid, v, workd, workl (tot. (4+ncv)*N * lworlk)
  double *resid, *v, *workd, *workl; // N, ncv*N, 3*N, lworkl
  double *x, *y; // used in rev. comm.
  arpack_int *select=nullptr; // ncv

  ncv = 2*nev;
  lworkl=ncv*ncv+8*ncv;
  info = 1; // info=0 use random guess; info!=0 use provided guess (stored in resid)
  iparam[1-1] = 1;
  iparam[3-1] = Arpack_MaxIter;  // max num of Arnoldi update iterations / return acutal num
  iparam[4-1] = 1;
  iparam[5-1] = 0; // number of converged eigenvalues found.
  iparam[7-1] = 3; // shift-and-invert mode
  
  invert_sign(); // invert sign to be positive definite (and vals[0] be the least-negative)
  
  if ((res = linear_solver_init()) < 0) { return res; }
  
  lwork_all = (int64_t)((4+ncv)*N + lworkl);
  if (verbose > 0) {
    cout << "eigen_iteration: ARPACK work size = " << lwork_all << " ("
      << fmt_size(sizeof(double)*lwork_all) << ")" << endl;
  }
  if ((res = allocate_array(work, lwork_all)) < 0) { return res; }
  if ((res = allocate_array(select, (int64_t)ncv)) < 0) { return res; }
  
  resid = work;
  v = resid + nsiz;
  workd = v + ((int64_t)ncv)*nsiz;
  workl = workd + 3*nsiz;
  for (i=0; i<(int64_t)ncv; i++) { select[i] = 1; }

  if (info != 0) { // symmetrize initial vector and use it as guess
    symmetrize_vector(z); normalize_vector(z);
    for (i=0; i<nsiz; i++) { resid[i] = z[i]; }
  }
  
  if (verbose > 0) {
    cout << "eigen_iteration: DSAUPD reverse communication ("
      << "neig=" << neig << ", ncv=" << ncv << ", tol=" << tol << ")" << endl;
  }
  while (true) {
    ARPACK(dsaupd)(ido, bmat, N, which, nev, tol, resid, ncv, v, ldv,
                   iparam, ipntr, workd, workl, lworkl, info);
    if ((ido == -1) || (ido == 1) || (ido == 2)) {
      x = workd + (ipntr[0] - 1);
      y = workd + (ipntr[1] - 1);
      for (i=0; i<nsiz; i++) { y[i] = x[i]; }
      if ((ido == -1) || (ido == 1)) {
        res = linear_solver_solve(y);
        if (res < 0) {
          delete_array(work, lwork_all);
          delete_array(select, (int64_t)ncv);
          return res;
        }
      }
    }
    else { break; }
  }
  if ((info != 0) || (ido != 99) || (iparam[5-1] < nev)) {
    cout << "ERROR: eigen_iteration: DSAUPD info = " << info << ", ido = "
      << ido << ", converged nev = " << iparam[5-1] << endl;
    delete_array(work, lwork_all);
    delete_array(select, (int64_t)ncv);
    return ErrArpack;
  }
  if (verbose > 0) {
    cout << "eigen_iteration: end of reverse communication ("
      << "niter=" << iparam[3-1] << ", nop=" << iparam[9-1]
        << ", nreo=" << iparam[11-1] << ")" << endl;
    for (i=0; i<neig; i++) {
      // these are invertly stored before calling dseupd
      cout << "eigen_iteration: Ritz value and relative error bound = "
        << -1. / workl[ipntr[6-1]-1+neig-1-i] << ", "
          << workl[ipntr[7-1]-1+neig-1-i] / workl[ipntr[6-1]-1+neig-1-i] << endl;
    }
  }
  
  if (verbose > 0) { cout << "eigen_iteration: calling DSEUPD..." << endl; }
  ARPACK(dseupd)(rvec, howmny, select, vals, z, ldz, sigma,
                 bmat, N, which, nev, tol, resid, ncv, v, ldv,
                 iparam, ipntr, workd, workl, lworkl, info);
  delete_array(work, lwork_all);
  delete_array(select, (int64_t)ncv);
  if (info != 0) {
    cout << "ERROR: eigen_iteration: DSEUPD info = " << info << endl;
    return ErrArpack;
  }
  linear_solver_clear();
  
  for (i=0; i<neig; i++) { vals[i] = -vals[i]; }
  for (i=0; i<neig; i++) { desymmetrize_vector(z+i*nsiz); normalize_vector(z+i*nsiz); }
  invert_sign();
  if (verbose > 0)
    { cout << "eigen_iteration: solved, least-negative eig. = " << vals[0] << endl; }
  
  return 0;
}


// ---- linear solvers ----

int MESolver::linear_solver_init() {
  lapack_int N=(lapack_int)nsiz, kd=(lapack_int)bwidth, ldab=kd+1;
  lapack_int lwork=-1, info=0;
  int res=0;
  double optimal_lwork;
  double *work = nullptr;
  char uplo='U';
  std::string routine;
  
  if ((me_type == Band) && (linear_solver == LsCHO)) { routine = "DPBTRF"; }
  else if ((me_type == Dens) && (linear_solver == LsCHO)) { routine = "DPOTRF"; }
  else if ((me_type == Dens) && (linear_solver == LsLU)) { routine = "DGETRF"; }
  else if ((me_type == Dens) && (linear_solver == LsLDLT)) { routine = "DSYTRF"; }
  else if ((me_type == Band) && (linear_solver == LsCG)) { routine = "CG_INIT(BAND)"; }
  else if ((me_type == Dens) && (linear_solver == LsCG)) { routine = "CG_INIT"; }
  else if ((me_type == MW) && (linear_solver == LsCG)) { routine = "CG_INIT(MW)"; }
  else {
    cout << "ERROR: linear_solver_init: (me_type,linear_solver) = ("
      << me_type << "," << linear_solver << ") not implemented." << endl;
    return ErrSolver;
  }

  linear_solver_clear();
  
  // allocate ipiv for LU/LDLT
  if ((linear_solver == LsLU) || (linear_solver == LsLDLT)) {
    if ((res = allocate_array(ipiv, nsiz)) < 0) { return res; }
    ipiv_size = nsiz;
  }
  
  // apply rate constants to U for direct solvers
  if ((linear_solver == LsCHO) || (linear_solver == LsLU)
      || (linear_solver == LsLDLT)) { apply_rates(); }
  
  if (verbose > 0)
    { cout << "linear_solver_init: calling " << routine << "..." << endl; }
  
  if ((me_type == Band) && (linear_solver == LsCHO))
    { LAPACK(dpbtrf)(uplo, N, kd, U, ldab, info); }
  else if ((me_type == Dens) && (linear_solver == LsCHO))
    { LAPACK(dpotrf)(uplo, N, U, N, info); }
  else if ((me_type == Dens) && (linear_solver == LsLU))
    { LAPACK(dgetrf)(N, N, U, N, ipiv, info); }
  else if ((me_type == Dens) && (linear_solver == LsLDLT)) {
    // check work array size
    LAPACK(dsytrf)(uplo, N, U, N, ipiv, &optimal_lwork, lwork, info);
    if (info == 0) {
      lwork = (lapack_int)optimal_lwork;
      if (verbose > 1) {
        cout << "linear_solver_init: LWORK = " << lwork << " ("
          << fmt_size(sizeof(double)*lwork) << ")" << endl;
      }
      if ((res = allocate_array(work, (int64_t)lwork)) < 0) { return res; }
      // factorize with optimal lwork
      LAPACK(dsytrf)(uplo, N, U, N, ipiv, work, lwork, info);
      delete_array(work, (int64_t)lwork);
    }
  }
  else if (linear_solver == LsCG) {
    if ((res = cg_init()) != 0) { return res; }
  }
  if (info != 0) {
    cout << "ERROR: linear_solver_init: " << routine << " info = " << info << endl;
    return ErrLapack;
  }
  return 0;
}


int MESolver::linear_solver_solve(double *b) {
  lapack_int N=(lapack_int)nsiz, kd=(lapack_int)bwidth, ldab=kd+1;
  lapack_int nrhs=1, info=0;
  int res=0;
  char uplo='U', trans='N';
  std::string routine;
  
  if ((me_type == Band) && (linear_solver == LsCHO)) { routine = "DPBTRS"; }
  else if ((me_type == Dens) && (linear_solver == LsCHO)) { routine = "DPOTRS"; }
  else if ((me_type == Dens) && (linear_solver == LsLU)) { routine = "DGETRS"; }
  else if ((me_type == Dens) && (linear_solver == LsLDLT)) { routine = "DSYTRS"; }
  else if ((me_type == Band) && (linear_solver == LsCG)) { routine = "CG_SOLVE(BAND)"; }
  else if ((me_type == Dens) && (linear_solver == LsCG)) { routine = "CG_SOLVE"; }
  else if ((me_type == MW) && (linear_solver == LsCG)) { routine = "CG_SOLVE(MW)"; }
  else {
    cout << "ERROR: linear_solver_solve: (me_type,linear_solver) = ("
      << me_type << "," << linear_solver << ") not implemented." << endl;
    return ErrSolver;
  }
  
  if (verbose > 0)
    { cout << "linear_solver_solve: calling " << routine << "..." << endl; }
  
  if ((me_type == Band) && (linear_solver == LsCHO))
    { LAPACK(dpbtrs)(uplo, N, kd, nrhs, U, ldab, b, N, info); }
  else if ((me_type == Dens) && (linear_solver == LsCHO))
    { LAPACK(dpotrs)(uplo, N, nrhs, U, N, b, N, info); }
  else if ((me_type == Dens) && (linear_solver == LsLU))
    { LAPACK(dgetrs)(trans, N, nrhs, U, N, ipiv, b, N, info); }
  else if ((me_type == Dens) && (linear_solver == LsLDLT))
    { LAPACK(dsytrs)(uplo, N, nrhs, U, N, ipiv, b, N, info); }
  else if (linear_solver == LsCG) {
    if ((res = cg_solve(b)) != 0) { return res; }
  }
  
  if (info != 0) {
    cout << "ERROR: linear_solver_solve: " << routine << " info = " << info << endl;
    return (double)ErrLapack;
  }
  return 0;
}



void MESolver::linear_solver_clear(){
  delete_array(ipiv, ipiv_size); ipiv_size = 0;
  delete_array(dwork_cg, dwork_cg_size); dwork_cg_size = 0;
  delete_array(qwork_cg, qwork_cg_size); qwork_cg_size = 0;
  delete_array(dwork_pc, dwork_pc_size); dwork_pc_size = 0;
  return;
}



// ---- preconditioner ----

int MESolver::pc_init() {
  int res;
  int64_t wsiz=nsiz;
  lapack_int N=(lapack_int)nsiz, kd, ldab, info;
  char uplo='U';

  if (preconditioner == PcDiag) {
    if (verbose > 0) { cout << "pc_init: diagonal preconditioner" << endl; }
  }
  else if (preconditioner == PcBand) {
    Pc_bwidth = get_reduced_bwidth(Pc_RelB);
    wsiz += nsiz * Pc_bwidth;
    if (verbose > 0) {
      cout << "pc_init: band preconditioner, size=("
        << nsiz << "," << Pc_bwidth + 1 << "), "
          << fmt_size(sizeof(double)*wsiz) << endl;
    }
  }
  
  delete_array(dwork_pc, dwork_pc_size);
  if ((res = allocate_array(dwork_pc, wsiz)) < 0) { return res; }
  dwork_pc_size = wsiz;
  
  if (preconditioner == PcDiag) { diagonal_U(dwork_pc); }
  else if (preconditioner == PcBand) {
    reduced_banded_U(dwork_pc, Pc_bwidth);
    kd = (lapack_int)Pc_bwidth;
    ldab = kd + 1;
    if (verbose > 0) { cout << "pc_init: calling DPBTRF..." << endl; }
    LAPACK(dpbtrf)(uplo, N, kd, dwork_pc, ldab, info);
    if (info != 0) {
      cout << "ERROR: pc_init: DPBTRF info = " << info << endl;
      return ErrLapack;
    }
  }
  return 0;
}


int MESolver::pc_solve(double *b) {
  int64_t i;
  lapack_int N=(lapack_int)nsiz, kd, ldab, nrhs=1, info;
  char uplo='U';
  
  if (preconditioner == PcDiag)
    { for (i=0; i<nsiz; i++) { b[i] /= dwork_pc[i]; } }
  else if (preconditioner == PcBand) {
    kd = (lapack_int)Pc_bwidth;
    ldab = kd + 1;
    LAPACK(dpbtrs)(uplo, N, kd, nrhs, dwork_pc, ldab, b, N, info);
    if (info != 0) {
      cout << "ERROR: pc_solve: Band solver, DPBTRF info = " << info << endl;
      return ErrLapack;
    }
  }
  return 0;
}



// ---- conjugate gradient ----

int MESolver::cg_init() {
  int res;

  delete_array(dwork_cg, dwork_cg_size);
  if ((res = allocate_array(dwork_cg, nsiz*5)) < 0) { return res; }
  dwork_cg_size = nsiz*5;

  if (CG_RTol_Semiquad >= 0.) {
    check_quad();
    delete_array(qwork_cg, qwork_cg_size);
    if ((res = allocate_array(qwork_cg, nsiz*6)) < 0) { return res; }
    qwork_cg_size = nsiz*6;
  }
  
  if ((res = pc_init()) < 0) { return res; }
  CG_guess_stored = false;
  return 0;
}


int MESolver::cg_solve(double *b) {
  // preconditioned conjugate gradient method.
  // ref: R.Barrett et.al., "Templates for the Solution of Linear Systems: Building
  //      Blocks for Iterative Methods, 2nd Edition"; SIAM, Philadelphia, 1994.
  //      <http://www.netlib.org/linalg/html_templates/Templates.html>
  int64_t i;
  int res, it=0, converged=false;
  double rtol, relnorm, bnorm, alpha, beta, rho, rho_old=1.;
  double *x = dwork_cg; // guess/solution
  double *r = dwork_cg + nsiz;
  double *z = dwork_cg + 2*nsiz;
  double *p = dwork_cg + 3*nsiz;
  double *q = dwork_cg + 4*nsiz;
  
  // use cg_solve_semiquad if CG_RTol < 0
  if (CG_RTol < 0.) {
    if (CG_RTol_Semiquad >= 0.) { return cg_solve_semiquad(b, false); }
    else {
      cout << "ERROR: cg_solve: both CG_RTol and CG_RTol_Semiquad are negative." << endl;
      return ErrSolver;
    }
  }

  if (CG_RTol > std::numeric_limits<double>::epsilon()) { rtol = CG_RTol; }
  else { rtol = std::numeric_limits<double>::epsilon(); }
  if (verbose > 0) { cout << "cg_solve: rtol = " << rtol << endl; }
  
  bnorm = nrm2(b);
  if (bnorm == 0.) { bnorm = 1.; }
  
  if (!CG_guess_stored) { // set initial guess vector to x
    for (i=0; i<nsiz; i++) { x[i] = b[i]; }
    if ((res = pc_solve(x)) < 0) { return res; } // x = M^-1 * b
    // scale magnitude of x if ||U*x|| >> ||b|| or ||U*x|| << ||b||
    multiply_U(x, r);
    relnorm = bnorm / nrm2(r);
    if ((relnorm < 0.01) || (relnorm > 100.))
      { for (i=0; i<nsiz; i++) { x[i] *= relnorm; } }
    CG_guess_stored = true;
  }
  
  multiply_U(x, r);  // r = U*x
  for (i=0; i<nsiz; i++) { r[i] = b[i] - r[i]; }  // r = b - U*x
  relnorm = nrm2(r) / bnorm;
  if (verbose > 0) { cout << "cg_solve: initial ||b-U*x||/||b|| = " << relnorm << endl; }
  if (relnorm < rtol) {
    if (verbose > 0) { cout << "cg_solve: converged with initial vector." << endl; }
    for (i=0; i<nsiz; i++) { b[i] = x[i]; }
    return 0;
  }
  
  it = 0;
  while (true) {
    it++;
    if (it > CG_MaxIter) { break; }
    for (i=0; i<nsiz; i++) { z[i] = r[i]; }
    if ((res = pc_solve(z)) < 0) { return res; } // z = M^-1 r
    rho = dot(r, z);
    if (it == 1) { for (i=0; i<nsiz; i++) { p[i] = z[i]; } }
    else {
      beta = rho / rho_old;
      for (i=0; i<nsiz; i++) { p[i] = z[i] + beta * p[i]; }
    }
    multiply_U(p, q);  // q = U*p
    alpha = rho / dot(p, q);
    for (i=0; i<nsiz; i++) { x[i] += alpha * p[i]; r[i] -= alpha * q[i]; }
    relnorm = nrm2(r) / bnorm;
    
    if (relnorm < rtol) { converged = true; }
    if (((verbose > 0) && (converged || (it % 50 == 0) || (it == CG_MaxIter)))
        || ((verbose > 1) && (it % 10 == 0)) || (verbose > 2)) {
      cout << "cg_solve: it = " << it << ", ||b-U*x||/||b|| = " << relnorm;
      if (converged) { cout << ", converged."; }
      cout << endl;
    }
    if (converged) { break; }
    rho_old = rho;
  }
  if (!converged) {
    cout << "ERROR: cg_solve: max iterations exceeded." << endl;
    return ErrCGMaxIter;
  }
  
  // calculate actual relnorm
  multiply_U(x, r);  // r = U*x
  for (i=0; i<nsiz; i++) { r[i] = b[i] - r[i]; }  // r = b - U*x
  relnorm = nrm2(r) / bnorm;

  if ((CG_RTol_Semiquad >= 0.) && (relnorm > CG_RTol_Semiquad)) {
    // refine solution
    if (verbose > 0) {
      cout << "cg_solve: actual ||b-U*x||/||b|| = "
        << relnorm << ", refine with CG_Semiquad."<< endl;
    }
    if ((res = cg_solve_semiquad(b, true)) < 0) { return res; }
  }
  else {
    if (verbose > 0) { cout << "cg_solve: actual  ||b-U*x||/||b|| = " << relnorm << endl; }
    for (i=0; i<nsiz; i++) { b[i] = x[i]; }
  }
  return 0;
}


int MESolver::cg_solve_semiquad(double *b, int x_as_guess) {
  int64_t i;
  int res, it=0, converged=false;
  double rtol, relnorm, bnorm;
  quad_float qbnorm, qalpha, qbeta, qrho, qrho_old=1.;
  quad_float *qx = qwork_cg; // guess/solution
  quad_float *qr = qwork_cg + nsiz;
  quad_float *qz = qwork_cg + 2*nsiz;
  quad_float *qp = qwork_cg + 3*nsiz;
  quad_float *qq = qwork_cg + 4*nsiz;
  quad_float *qb = qwork_cg + 5*nsiz;
  double *x = dwork_cg;
  double *r = dwork_cg + nsiz;
  double *z = dwork_cg + 2*nsiz;
  
  if (CG_RTol_Semiquad > std::numeric_limits<double>::epsilon()) { rtol = CG_RTol_Semiquad; }
  else { rtol = std::numeric_limits<double>::epsilon(); }
  if (verbose > 0) { cout << "cg_solve_semiquad: rtol = " << rtol << endl; }
  
  for (i=0; i<nsiz; i++) { qb[i] = b[i]; }
  bnorm = nrm2(b);
  if (bnorm == 0.) { bnorm = 1.; }
  qbnorm = nrm2(qb);
  if (qbnorm == 0.) { qbnorm = 1.; }

  if (x_as_guess) {
    for (i=0; i<nsiz; i++) { qx[i] = x[i]; }
  }
  else if (!CG_guess_stored) { // set initial guess vector to qx
    for (i=0; i<nsiz; i++) { x[i] = b[i]; }
    if ((res = pc_solve(x)) < 0) { return res; } // x = M^-1 * b
    for (i=0; i<nsiz; i++) { qx[i] = x[i]; }
    // scale magnitude of x if ||U*x|| >> ||b|| or ||U*x|| << ||b||
    multiply_U_semiquad(qx, qr);
    relnorm = (double)(qbnorm / nrm2(qr));
    if ((relnorm < 0.01) || (relnorm > 100.))
      { for (i=0; i<nsiz; i++) { qx[i] *= relnorm; } }
    CG_guess_stored = true;
  }
  
  multiply_U_semiquad(qx, qr);  // r = U*x
  for (i=0; i<nsiz; i++) { qr[i] = qb[i] - qr[i]; }  // r = b - U*x
  relnorm = (double)(nrm2(qr) / qbnorm);
  if (verbose > 0) { cout << "cg_solve_semiquad: initial ||b-U*x||/||b|| = " << relnorm << endl; }
  if (relnorm < rtol) {
    if (verbose > 0) { cout << "cg_solve_semiquad: converged with initial vector." << endl; }
    for (i=0; i<nsiz; i++) { qb[i] = qx[i]; }
    for (i=0; i<nsiz; i++) { b[i] = (double)(qx[i]); }
    return 0;
  }
  
  it = 0;
  while (true) {
    it++;
    if (it > CG_MaxIter) { break; }
    // use double float for preconditioning
    for (i=0; i<nsiz; i++) { z[i] = (double)(qr[i]); }
    if ((res = pc_solve(z)) < 0) { return res; } // z = M^-1 r
    for (i=0; i<nsiz; i++) { qz[i] = z[i]; }
    qrho = dot(qr, qz);
    if (it == 1) { for (i=0; i<nsiz; i++) { qp[i] = qz[i]; } }
    else {
      qbeta = qrho / qrho_old;
      for (i=0; i<nsiz; i++) { qp[i] = qz[i] + qbeta * qp[i]; }
    }
    multiply_U_semiquad(qp, qq);  // q = U*p
    qalpha = qrho / dot(qp, qq);
    for (i=0; i<nsiz; i++) { qx[i] += qalpha * qp[i]; qr[i] -= qalpha * qq[i]; }
    relnorm = (double)(nrm2(qr) / qbnorm);
    
    if (relnorm < rtol) { converged = true; }
    if (((verbose > 0) && (converged || (it % 50 == 0) || (it == CG_MaxIter)))
        || ((verbose > 1) && (it % 10 == 0)) || (verbose > 2)) {
      cout << "cg_solve_semiquad: it = " << it << ", ||b-U*x||/||b|| = " << relnorm;
      if (converged) { cout << ", converged."; }
      cout << endl;
    }
    if (converged) { break; }
    qrho_old = qrho;
  }
  if (!converged) {
    cout << "ERROR: cg_solve_semiquad: max iterations exceeded." << endl;
    return ErrCGMaxIter;
  }
  
  // calculate actual relnorm (quad)
  multiply_U_semiquad(qx, qr);  // r = U*x
  for (i=0; i<nsiz; i++) { qr[i] = qb[i] - qr[i]; }  // r = b - U*x
  relnorm = (double)(nrm2(qr) / qbnorm);
  if (verbose > 0)
    { cout << "cg_solve_semiquad: actual  ||b-U*x||/||b|| (quad)   = " << relnorm << endl; }
  
  // calculate actual relnorm (double)
  for (i=0; i<nsiz; i++) { x[i] = (double)(qx[i]); }
  multiply_U(x, r);  // r = U*x
  for (i=0; i<nsiz; i++) { r[i] = b[i] - r[i]; }  // r = b - U*x
  relnorm = nrm2(r) / bnorm;
  if (verbose > 0)
    { cout << "cg_solve_semiquad: actual  ||b-U*x||/||b|| (double) = " << relnorm << endl; }
  
  for (i=0; i<nsiz; i++) { qb[i] = qx[i]; }
  for (i=0; i<nsiz; i++) { b[i] = x[i]; }
  return 0;
}



// ---- direct eigen solver ----

int MESolver::eigen_solve(int64_t neig, double *vals, double *z) {
  if ((me_type == Dens) && (eigen_solver == EsDSYEVR))
    { return eigen_solve_dsyevr(neig, vals, z); }
  
  cout << "ERROR: eigen_solve: (me_type,eigen_solver) = ("
    << me_type << "," << eigen_solver << ") not implemented." << endl;
  return ErrSolver;
}


int MESolver::eigen_solve_dsyevr(int64_t neig, double *vals, double *z) {
  int res;
  int64_t i;
  lapack_int N=(lapack_int)nsiz, il=1, iu=(lapack_int)neig, m;
  lapack_int lwork=-1, liwork=-1, optimal_liwork, info;
  double vl=0., vu=0., abstol=0., optimal_lwork;
  char jobz='V', range='I', uplo='U';
  double *w = nullptr;
  double *work = nullptr;
  lapack_int *isuppz = nullptr;
  lapack_int *iwork = nullptr;
  
  invert_sign(); // invert sign so that vals[0] is the least-negative eigenvalue
  apply_rates(); // apply rate constants to U for direct solver
  
  if (verbose > 0)
    { cout << "eigen_solve_dsyevr: calling DSYEVR..." << endl; }

  // check work array sizes
  LAPACK(dsyevr)(jobz, range, uplo, N, U, N, vl, vu, il, iu, abstol, m, w, z,
                 N, isuppz, &optimal_lwork, lwork, &optimal_liwork, liwork, info);
  if (info == 0) {
    lwork = (lapack_int)(optimal_lwork) + 6*N;
    liwork = optimal_liwork;
    if (verbose > 1) {
      cout << "eigen_solve_dsyevr: LWORK = " << lwork << ", LIWORK = " << liwork
        << " (" << fmt_size(sizeof(double)*lwork + sizeof(lapack_int)*(liwork))
          << ")" << endl;
    }
    if ((res = allocate_array(w, nsiz)) < 0) { return res; }
    if ((res = allocate_array(work, (int64_t)lwork)) < 0) { return res; }
    if ((res = allocate_array(isuppz, 2*neig)) < 0) { return res; }
    if ((res = allocate_array(iwork, (int64_t)liwork)) < 0) { return res; }
    
    // solve eigen problem with optimal lwork and liwork
    LAPACK(dsyevr)(jobz, range, uplo, N, U, N, vl, vu, il, iu, abstol,
                   m, w, z, N, isuppz, work, lwork, iwork, liwork, info);
    for (i=0; i<neig; i++) { vals[i] = -w[i]; }
    delete_array(w, nsiz);
    delete_array(work, (int64_t)lwork);
    delete_array(isuppz, 2*neig);
    delete_array(iwork, (int64_t)liwork);
  }
  
  if (info != 0) {
    cout << "ERROR: eigen_solve_dsyevr: DSYEVR info = " << info << endl;
    return ErrLapack;
  }
  for (i=0; i<neig; i++) { desymmetrize_vector(z+i*nsiz); normalize_vector(z+i*nsiz); }
  invert_sign();
  if (verbose > 0) {
    cout << "eigen_solve_dsyevr: solved, neig = " << m << ", least-negative eig. = "
      << vals[0] << endl;
  }
  return 0;
}



