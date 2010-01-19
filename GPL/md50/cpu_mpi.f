*   This procedure get wall clock cpu time using "MPI_Wtime" function
*   (from the MPI library)
*
       real*8 function cputime(was)
       include 'mpif.h'
       real*8 now,was
       now=MPI_Wtime()
       cputime = now-was
       end
