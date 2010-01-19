#include <time.h>
#include <sys/resource.h>

getcpu_(e) double *e; {
   *e = (double)clock()/(double)CLOCKS_PER_SEC;
}
