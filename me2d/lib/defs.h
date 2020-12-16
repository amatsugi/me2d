#ifndef _DEFS_H_
#define _DEFS_H_

#include <cmath>
#include <cstdint>
using std::int64_t;

// ---- quad ----

#ifndef NO_QUAD
#define NO_QUAD 0
#endif

#if NO_QUAD
typedef long double quad_float;
#elif defined(__GNUC__)
extern "C" {
#include <quadmath.h>
}
typedef __float128 quad_float;
namespace std {
  inline quad_float sqrt(quad_float x) { return sqrtq(x); }
  inline quad_float fabs(quad_float x) { return fabsq(x); }
}
#elif defined(__INTEL_COMPILER)
typedef _Quad quad_float;
extern "C" {
_Quad __sqrtq(_Quad);
_Quad __fabsq(_Quad);
}
namespace std {
  inline quad_float sqrt(quad_float x) { return __sqrtq(x); }
  inline quad_float fabs(quad_float x) { return __fabsq(x); }
}
#else
#error "quadruple-precision float not enabled; use GCC or Intel, or define NO_QUAD=1"
#endif


// ---- external libs. ----

#ifndef NO_ARPACK
#define NO_ARPACK 0
#endif

#ifndef INTERFACE64
#define INTERFACE64 0
#endif

#if INTERFACE64
typedef int64_t blas_int;
typedef int64_t lapack_int;
typedef int64_t arpack_int;
#else
typedef int blas_int;
typedef int lapack_int;
typedef int arpack_int;
#endif

#define BLAS(NAME) NAME##_
#define LAPACK(NAME) NAME##_
#define ARPACK(NAME) NAME##_



#endif

