        subroutine sfft3d(f,N)

        implicit none
#include "fftw_f77.i"

c Computes 3D sine transform on a cubic grid using FFTW. 
        integer N
        real*8 f(0:N,0:N,0:N)

        real*8 in1d(0:N-1)
        real*8 sincoef(0:N-1)
	    real*8 out1d(0:N-1)

	    integer i,j,k,ll
	    integer plan
        real*8 pi
c N must be an even number.
#if 0
        if ((N/2)*2 .ne. N) then
          write(6,*) "sfft3d: N is odd, N = ",N
          call abort
        endif
#endif
c       
c       Set up preliminary stuff.

        call rfftw_f77_create_plan(plan,N,FFTW_REAL_TO_COMPLEX,
     &          FFTW_ESTIMATE)

        pi = 4*datan(1.0d0)

        do ll = 0,N-1

            sincoef(ll) = sin(ll*pi/N)

        enddo

c Sine transform in the x direction.

        do k = 0,N-1
          do j = 0,N-1
            in1D(0) = 0.d0
            do i = 1,N-1

              in1D(i) = sincoef(i)*(f(i,j,k) + f(N-i,j,k)) +
     &                  (f(i,j,k) - f(N-i,j,k))/2
            enddo

	        call rfftw_f77_one(plan,in1d,out1d)
           
            f(1,j,k) = out1d(0)/2
            f(0,j,k) = 0.d0
            
            do i = 1,N/2-1

              f(2*i,j,k) = -out1d(N - i)
              f(2*i+1,j,k) = out1d(i) + f(2*i-1,j,k)
            
            enddo
          enddo
        enddo

c Sine transform in the y direction.

        do k = 0,N-1
          do j = 0,N-1
            in1D(0) = 0.d0
            do i = 1,N-1

              in1D(i) = sincoef(i)*(f(j,i,k) + f(j,N-i,k)) +
     &                  (f(j,i,k) - f(j,N-i,k))/2
            enddo

	        call rfftw_f77_one(plan,in1d,out1d)
           
            f(j,1,k) = out1d(0)/2
            f(j,0,k) = 0.d0
            
            do i = 1,N/2-1

              f(j,2*i,k) = -out1d(N - i)
              f(j,2*i+1,k) = out1d(i) + f(j,2*i-1,k)
            
            enddo
          enddo
        enddo

c Sine transform in the z direction.

        do k = 0,N-1
          do j = 0,N-1
            in1D(0) = 0.d0
            do i = 1,N-1

              in1D(i) = sincoef(i)*(f(j,k,i) + f(j,k,N-i)) +
     &                  (f(j,k,i) - f(j,k,N-i))/2
            enddo

	        call rfftw_f77_one(plan,in1d,out1d)
           
            f(j,k,0) = 0.d0
            f(j,k,1) = out1d(0)/2
            
            do i = 1,N/2-1

              f(j,k,2*i) = -out1d(N - i)
              f(j,k,2*i+1) = out1d(i) + f(j,k,2*i-1)
            
            enddo
          enddo
        enddo
	    call rfftw_f77_destroy_plan(plan)
        return
        end
