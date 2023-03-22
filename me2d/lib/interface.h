#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "defs.h"

extern "C" {
void set_me_maxmem_GB(double memGB);
double get_me_maxmem_GB();
void show_solvers();

/* use standard integer types for exported functions */

int64_t solve1d(int64_t nsiz, int64_t neig, double *vals, double *z,
            double *Ea, double *rhoa, double *ka,
            double y_e, double *ainv_ea, int64_t ptype,
            double bandpcrit, double ZM, double kbT, char *solver, char *chkfn, int64_t verbose);
int64_t solve2d(int64_t nsiz, int64_t neig, double *vals, double *z,
            double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
            double y_e, double y_J, double *ainv_ea, double *ainv_Ja, int64_t ptype,
            double bandpcrit, double ZM, double kbT, char *solver, char *chkfn, int64_t verbose);

int64_t solve1d_mw(int64_t nwell, int64_t *nsiz, int64_t neig, double *vals, double *z,
               double *Ea, double *rhoa, double *ka,
               double *y_e, double *ainv_ea, int64_t *ptype,
               int64_t nkisom, double *kisom, int64_t *kisom_i, int64_t *kisom_j,
               double *bandpcrit, double *ZM, double kbT,
               char *solver, int64_t reactant, char *chkfn, int64_t verbose);
int64_t solve2d_mw(int64_t nwell, int64_t *nsiz, int64_t neig, double *vals, double *z,
               double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
               double *y_e, double *y_J, double *ainv_ea, double *ainv_Ja, int64_t *ptype,
               int64_t nkisom, double *kisom, int64_t *kisom_i, int64_t *kisom_j,
               double *bandpcrit, double *ZM, double kbT,
               char *solver, int64_t reactant, char *chkfn, int64_t verbose);
}

int solve(int64_t nsiz, int64_t neig, double *vals, double *z,
          double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
          double y_e, double y_J, double *ainv_ea, double *ainv_Ja, int64_t ptype,
          double bandpcrit, double ZM, double kbT, char *solver, char *chkfn, int64_t verbose);

int solve_mw(int64_t nwell, int64_t *nsiz, int64_t neig, double *vals, double *z,
             double *Ea, double *ea, double *Ja, double *rhoa, double *ka,
             double *y_e, double *y_J, double *ainv_ea, double *ainv_Ja, int64_t *ptype,
             int64_t nkisom, double *kisom, int64_t *kisom_i, int64_t *kisom_j,
             double *bandpcrit, double *ZM, double kbT,
             char *solver, int64_t reactant, char *chkfn, int64_t verbose);


#endif
