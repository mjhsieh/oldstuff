#ifndef _gbmacros_
#define _gbmacros_

// These are random useful macros and definitions used in scallop.
// Actually, some may not be useful anymore.  At some point I'll go
// through and remove the useless stuff.

#include "def3.h"

#ifdef DIM2
#define FortranRegion FortranRegion2
#define FORTRAN_REGION FORTRAN_REGION2
#define FORTRAN_ARGS FORTRAN_ARGS2
#else
#define FortranRegion FortranRegion3
#define FORTRAN_REGION FORTRAN_REGION3
#define FORTRAN_ARGS FORTRAN_ARGS3
#endif


#define FORT_ARRAY_ARG(A, I) A(I).data(), FORT_REGION(A(I).domain)
#define FORT_ARRAY_ARG_PTR(A, I) A->ptr(I)->data(), \
FORT_REGION(A->ptr(I)->domain)
#define FORT_GRID_PTR(A) A->data(), FORT_REGION(A->domain)
#define FORT_ARRAY_BOUND_PTR(A, I) FORT_REGION(A->ptr(I)->domain)
#define FORT_ARRAY_BOUND(A, I) FORT_REGION(A(I).domain)
#define FORT_GRID_BOUND_PTR(A) FORT_REGION(A->domain)

#ifdef DIM2
#define FAA_TYPE double*, const int*, const int*, const int*, const int*
#define FR_TYPE const int*, const int*, const int*, const int*
#define FORT_REGION(R) &R.l0, &R.l1, &R.u0, &R.u1
#else
#define FAA_TYPE double*, const int*, const int*, const int*, const int*, \
const int*, const int*
#define FR_TYPE const int*, const int*, const int*, const int*, const int*, \
const int*
#define FORT_REGION(R) &R.l0, &R.l1, &R.l2, &R.u0, &R.u1, &R.u2
#endif


#define LEFT    0
#define RIGHT   1
#define FRONT   2
#define BACK    3
#define BOTTOM  4
#define TOP     5

#define ECBC 0
#define CCBC 1

#define CENTERED 0
#define PERIODIC 1
#define FROMFILE 2
#define OFFSET 3
#define MPERIODIC 4	
#define RANDOMISH 5
#define ONE 6

#endif
