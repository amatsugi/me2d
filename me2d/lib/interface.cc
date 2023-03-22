#include "interface.h"
#include "funcs.h"
#include "mesolver.h"
#include "mesingle.h"
#include "memulti.h"
#include <iostream>
#include <string>
using std::cout;
using std::endl;

/*
interface for ME solvers
  [solve1d, solve2d] 1D/2D ME for single well system
    int64_t nsiz: problem size
    int64_t neig: No. of eigenpairs to be calculated
    double *vals: [neig] (output) calculated eigenvalues
    double *z: [neig*nsiz] (input for iterative solvers) guess,
               (input for LinEq) flux vector, (output) eigenvectors
    double *Ea: [nsiz] total energies
    double *ea: [nsiz] (only in 2D) active energies (used in P(E,J;E',J'))
    double *Ja: [nsiz] (only in 2D) total angular momenta
    double *rhoa: [nsiz] density of states
    double *ka: [nsiz] microscopic rate coefficients
    double y_e: parameter y_e for P(E;E') and P(E,J;E',J')
    double y_J: (only in 2D) parameter y_J for P(E,J;E',J')
    double *ainv_ea: [nsiz] parameters 1/alpha_e for P(E;E') and P(E,J;E',J')
    double *ainv_Ja: [nsiz] (only in 2D) parameters 1/alpha_J for P(E,J;E',J')
    int64_t ptype: type of probability function:  0=symmetric, -1=downward, 1=upward
    double bandpcrit: criteria for banded matrix; negative to use dense matrix
    double ZM: Z[M]
    double kbT: kbT
    char *solver: solver specification (see MESolver::show_solvers)
    char *chkfn: file name for storing matrix data (used with the solver options save/load)
    int64_t verbose: verbose level
  
  [solve1d_mw, solve2d_mw] 1D/2D ME for multiple well system
    int64_t nwell: No. of wells
    int64_t *nsiz: [nwell] well sizes
    int64_t neig: No. of eigenpairs to be calculated
    double *vals: [neig] (output) calculated eigenvalues
    double *z: [neig*sum(nsiz)] (input for iterative solvers) guess,
               (input for LinEq) flux vector, (output) eigenvectors
    double *Ea: [sum(nsiz)] total energies
    double *ea: [sum(nsiz)] (only in 2D) active energies (used in P(E,J;E',J'))
    double *Ja: [sum(nsiz)] (only in 2D) total angular momenta
    double *rhoa: [sum(nsiz)] density of states
    double *ka: [sum(nsiz)] microscopic rate coefficients
    double *y_e: [nwell] parameter y_e for P(E;E') and P(E,J;E',J')
    double *y_J: [nwell] (only in 2D) parameter y_J for P(E,J;E',J')
    double *ainv_ea: [sum(nsiz)] parameters 1/alpha_e for P(E;E') and P(E,J;E',J')
    double *ainv_Ja: [sum(nsiz)] (only in 2D) parameters 1/alpha_J for P(E,J;E',J')
    int64_t *ptype: [nwell] type of probability function:  0=symmetric, -1=downward, 1=upward
    int64_t nkisom: No. of microscopic rate coefficients for isomerization
    double *kisom:  [nkisom] microscopic rate coefficients for isomerization
    int64_t *kisom_i: [nkisom] index i for kisom
    int64_t *kisom_j: [nkisom] index j for kisom
    double *bandpcrit: [nwell] criteria for banded matrix; negative to use dense matrix
    double *ZM: [nwell] Z[M]
    double kbT: kbT
    char *solver: solver specification (see MESolver::show_solvers)
    int64_t reactant: index of reactant well (only for InvIter solver)
    char *chkfn: file name for storing matrix data (used with the solver options save/load)
    int64_t verbose: verbose level

 */

void set_me_maxmem_GB(double memGB) {
  MESolver::max_array_size = (int64_t)(memGB*1024.*1024.*1024.);
}

double get_me_maxmem_GB() {
  return (double)MESolver::max_array_size /1024./1024./1024.;
}

void show_solvers() {
  MESolver::show_solvers();
}

int64_t solve1d(int64_t nsiz, int64_t neig, double *vals, double *z,
                double *Ea, double *rhoa, double *ka,
                double y_e, double *ainv_ea, int64_t ptype,
                double bandpcrit, double ZM, double kbT, char *solver, char *chkfn, int64_t verbose) {
  double *ea = Ea;  // ea is used in the probability function.
  double *Ja = nullptr; // set Ja = nullptr for 1D
  double *ainv_Ja = nullptr; // not referenced if Ja == nullptr
  double y_J = 1.; // dummy
  if (verbose > 0) { cout << "Solve 1D master equation." << endl; }
  return (int64_t)solve(nsiz, neig, vals, z, Ea, ea, Ja, rhoa, ka,
                        y_e, y_J, ainv_ea, ainv_Ja, ptype,
                        bandpcrit, ZM, kbT, solver, chkfn, verbose);
}

int64_t solve2d(int64_t nsiz, int64_t neig, double *vals, double *z,
            double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
            double y_e, double y_J, double *ainv_ea, double *ainv_Ja, int64_t ptype,
            double bandpcrit, double ZM, double kbT, char *solver, char *chkfn, int64_t verbose) {
  if (verbose > 0) { cout << "Solve 2D master equation." << endl; }
  return (int64_t)solve(nsiz, neig, vals, z, Ea, ea, Ja, rhoa, ka,
                        y_e, y_J, ainv_ea, ainv_Ja, ptype,
                        bandpcrit, ZM, kbT, solver, chkfn, verbose);
}


int64_t solve1d_mw(int64_t nwell, int64_t *nsiz, int64_t neig, double *vals, double *z,
                   double *Ea, double *rhoa, double *ka,
                   double *y_e, double *ainv_ea, int64_t *ptype,
                   int64_t nkisom, double *kisom, int64_t *kisom_i, int64_t *kisom_j,
                   double *bandpcrit, double *ZM, double kbT,
                   char *solver, int64_t reactant, char *chkfn, int64_t verbose){
  double *ea = Ea;  // ea is used in the probability function.
  double *Ja = nullptr; // set Ja = nullptr for 1D
  double *ainv_Ja = nullptr; // not referenced if Ja == nullptr
  double *y_J = nullptr; // not referenced if Ja == nullptr
  if (verbose > 0) { cout << "Solve multiple-well 1D master equation." << endl; }
  return (int64_t)solve_mw(nwell, nsiz, neig, vals, z, Ea, ea, Ja, rhoa, ka,
                           y_e, y_J, ainv_ea, ainv_Ja, ptype,
                           nkisom, kisom, kisom_i, kisom_j,
                           bandpcrit, ZM, kbT, solver, reactant, chkfn, verbose);
}

int64_t solve2d_mw(int64_t nwell, int64_t *nsiz, int64_t neig, double *vals, double *z,
                   double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
                   double *y_e, double *y_J, double *ainv_ea, double *ainv_Ja, int64_t *ptype,
                   int64_t nkisom, double *kisom, int64_t *kisom_i, int64_t *kisom_j,
                   double *bandpcrit, double *ZM, double kbT,
                   char *solver, int64_t reactant, char *chkfn, int64_t verbose){
  if (verbose > 0) { cout << "Solve multiple-well 2D master equation." << endl; }
  return (int64_t)solve_mw(nwell, nsiz, neig, vals, z, Ea, ea, Ja, rhoa, ka,
                           y_e, y_J, ainv_ea, ainv_Ja, ptype,
                           nkisom, kisom, kisom_i, kisom_j,
                           bandpcrit, ZM, kbT, solver, reactant, chkfn, verbose);
}


int solve(int64_t nsiz, int64_t neig, double *vals, double *z,
          double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
          double y_e, double y_J, double *ainv_ea, double *ainv_Ja, int64_t ptype,
          double bandpcrit, double ZM, double kbT, char *solver, char *chkfn, int64_t verbose) {
  int res;
  double t0;
  std::string solver_str = solver;
  std::string chkfn_str = chkfn;
  t0 = get_wtime();
  
  MESingle *me = new MESingle();
  me->verbose = verbose;
  
  if (bandpcrit < 0.) { res = me->init_dens(nsiz); }
  else { res = me->init_band(nsiz, ea, Ja, y_e, y_J, ainv_ea, ainv_Ja, bandpcrit); }
  if (res < 0) { delete me; return res; }
  if (verbose > 0) { cout << "solve: initialized. " << fmt_elapsed(t0) << endl; }
  
  if ((res = me->set_solver(solver_str)) < 0) { delete me; return res; }
  if ((res = me->set_chkfn(chkfn_str)) < 0) { delete me; return res; }
  
  me->set_ka(ka);
  me->set_sa(rhoa, Ea, kbT);
  me->set_prob(ea, Ja, y_e, y_J, ainv_ea, ainv_Ja, ptype);
  if (verbose > 0) { cout << "solve: probability set. " << fmt_elapsed(t0) << endl; }
  
  if ((res = me->normalize()) < 0) { delete me; return res; }
  if (verbose > 0) { cout << "solve: normalized. " << fmt_elapsed(t0) << endl; }
  
  if ((res = me->symmetrize()) < 0) { delete me; return res; }
  if (verbose > 0) { cout << "solve: symmetrized. " << fmt_elapsed(t0) << endl; }
  
  me->apply_collision_rate(ZM);
  
  if (verbose > 0) { cout << "solve: call solver. " << fmt_elapsed(t0) << endl; }
  
  if ((res = me->solve(neig, vals, z)) < 0) { delete me; return res; }
  if (verbose > 0) {
    cout << "solve: solved. " << fmt_size(me->total_array_size_used) << 
      " of memory used. " << fmt_elapsed(t0) << endl;
  }
  delete me;
  
  if (verbose > 0) {
    if (MESolver::total_array_size > 0) {
      cout << "WARNING: solve: total_array_size = "
        << MESolver::total_array_size << " remains to be released." << endl;
    }
  }
  MESolver::total_array_size_used = MESolver::total_array_size;
  return 0;
}


int solve_mw(int64_t nwell, int64_t *nsiz, int64_t neig, double *vals, double *z,
             double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
             double *y_e, double *y_J, double *ainv_ea, double *ainv_Ja, int64_t *ptype,
             int64_t nkisom, double *kisom, int64_t *kisom_i, int64_t *kisom_j,
             double *bandpcrit, double *ZM, double kbT,
             char *solver, int64_t reactant, char *chkfn, int64_t verbose) {
  int res;
  int64_t iwell, pos=0;
  MESingle *me1;
  double t0;
  std::string solver_str = solver;
  std::string chkfn_str = chkfn;
  t0 = get_wtime();
  
  MEMulti *multi = new MEMulti();
  multi->verbose = verbose;
  
  if ((res = multi->init(nwell)) < 0) { delete multi; return res; }
  if ((res = multi->set_solver(solver_str)) < 0) { delete multi; return res; }
  if ((res = multi->set_chkfn(chkfn_str)) < 0) { delete multi; return res; }

  for (iwell=0; iwell<nwell; iwell++) {
    me1 = &(multi->wells[iwell]);
    if (iwell == 0) { pos = 0; }
    else { pos += nsiz[iwell-1]; }
    
    if (bandpcrit[iwell] < 0.) { res = me1->init_dens(nsiz[iwell]); }
    else if (Ja == nullptr) {
      res = me1->init_band(nsiz[iwell], ea+pos, nullptr, y_e[iwell], 1.,
                           ainv_ea+pos, nullptr, bandpcrit[iwell]);
    }
    else {
      res = me1->init_band(nsiz[iwell], ea+pos, Ja+pos, y_e[iwell], y_J[iwell],
                           ainv_ea+pos, ainv_Ja+pos, bandpcrit[iwell]);
    }
    if (res < 0) { delete multi; return res; }
    if (verbose > 0) {
      cout << "solve_mw: well #" << iwell+1 << " initialized. " << fmt_elapsed(t0) << endl;
    }
    
    me1->set_ka(ka+pos);
    me1->set_sa(rhoa+pos, Ea+pos, kbT);
    if (Ja == nullptr) {
      me1->set_prob(ea+pos, nullptr, y_e[iwell], 1., ainv_ea+pos, nullptr, ptype[iwell]);
    }
    else {
      me1->set_prob(ea+pos, Ja+pos, y_e[iwell], y_J[iwell], ainv_ea+pos,
                   ainv_Ja+pos, ptype[iwell]);
    }
    if (verbose > 0) {
      cout << "solve_mw: well #" << iwell+1 << " probability set. " << fmt_elapsed(t0) << endl;
    }
  
    if ((res = me1->normalize()) < 0) { delete multi; return res; }
    if (verbose > 0) {
      cout << "solve_mw: well #" << iwell+1 << " normalized. " << fmt_elapsed(t0) << endl;
    }
  
    if ((res = me1->symmetrize()) < 0) { delete multi; return res; }
    if (verbose > 0) {
      cout << "solve_mw: well #" << iwell+1 << " symmetrized. " << fmt_elapsed(t0) << endl;
    }
    
    me1->apply_collision_rate(ZM[iwell]);
  }
  
  if ((res = multi->post_set_wells()) < 0) { delete multi; return res; }
  if ((res = multi->set_kisom(nkisom, kisom, kisom_i, kisom_j)) < 0) { delete multi; return res; }
  if (multi->me_type == multi->Dens) {
    if ((res = multi->set_full_matrix()) < 0) { delete multi; return res; }
  }
  
  if ((reactant >= 0) && (reactant < nwell)) {
    if ((res = multi->set_reactant(reactant)) < 0) { delete multi; return res; }
  }
  
  if (verbose > 0) { cout << "solve_mw: call solver. " << fmt_elapsed(t0) << endl; }
  
  if ((res = multi->solve(neig, vals, z)) < 0) { delete multi; return res; }
  if (verbose > 0) {
    cout << "solve_mw: solved. " << fmt_size(multi->total_array_size_used) << 
      " of memory used. " << fmt_elapsed(t0) << endl;
  }
  
  delete multi;
  
  if (verbose > 0) {
    if (MESolver::total_array_size > 0) {
      cout << "WARNING: solve_mw: total_array_size = "
        << MESolver::total_array_size << " remains to be released." << endl;
    }
  }
  MESolver::total_array_size_used = MESolver::total_array_size;
  return 0;
}
