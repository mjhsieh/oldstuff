/*  report used CPU time to a Fortran program (Cray Fortran) */
#include <sys/time.h>
#include <sys/resource.h>
#include <fortran.h>

GETCPU(e) float *e; {
  clock_t x;
   x = clock();
   *e = x/CLOCKS_PER_SEC;
}


