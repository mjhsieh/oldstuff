
#include <sys/time.h>
#include <sys/resource.h>

double cputime(start) double start; {
  double e;
  int w = 0;
  struct rusage o;
  getrusage( w, &o );
  e = ( (double)o.ru_utime.tv_sec + (double)o.ru_utime.tv_usec*1.0e-6 )
      - start;
  return e;
}

getcpu_( start, cpu ) float *start, *cpu; {
  *cpu = (float)(cputime((double)(*start)));
}

