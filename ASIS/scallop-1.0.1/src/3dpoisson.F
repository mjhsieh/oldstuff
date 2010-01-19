        subroutine poisson3d(f,N,H,Op)

c Nodal-point Poisson solver with Dirichlet boundary conditions on a cube.
c inputs:
c N : number of unknowns in each direction, + 1. n should be a highly
c     composite number.
c f(0:N,0:N,0:N): rhs on input, solution to Poisson's equation on
c                 output. Rhs is defined on interior nodes [1:n-1]^3. 
c H : mesh spacing.
c Op: operator type (0=Spectral, 1=Mehrstellen, 2=7point)

c The algorithm uses sfft3d, a 3d sine transform that uses
c the real-to-complex 1D FFT's in the FFTW library. This program must be
c linked with librfftw.a and libfftw.a . 

c The structure of the program is as follows:
c (1) Compute the 3d sine transform of the inputs - sfft updates f in
c     place.
c (2) Multiply each Fourier component by the symbol of the operator 
c     evaluated at that component. These are easily derived using 
c     the addition formulas for sin(A+B). In the case of the spectral 
c     approximation, multiply the Fourier coefficients by the symbol 
c     of the Laplacian. We also scale the coefficients to account for
c     the fact that the sine trasform applied twice is (n/2)^3 times
c     the identity operator.
c (3) Apply the sine transform again, which is its own inverse.

        implicit none
!Args
        integer N,Op
        real*8  f(0:n,0:n,0:n)
        real*8  H
!Locals
        real*8 coscoef(N-1),pi,ampl
        integer i,j,k,ll,nh

c Apply the sine transform.

        call sfft3d(f,n)

c Compute some setup information.

        nh = n/2

        pi = ATAN(1.d0)*4
        do ll = 1,n-1
            coscoef(ll) = cos(ll*pi/n)
        enddo

c Multiply each Fourier coefficient by the symbol of the operator.

        if( op .EQ. 0 )then
            ! Spectral discretization.
            do k = 1,n-1
                do j = 1,n-1
                    do i = 1,n-1
                        ampl = -(i*i + j*j + k*k)*pi*pi/(n*h)/(n*h)
                        ! divide by the symbol, and scale by 1/(n/2)^3.
                        f(i,j,k) = f(i,j,k)/ampl/nh/nh/nh
                    enddo
                enddo
            enddo
        elseif( op .EQ. 1 )then
            ! Mehrstellen operator.
            do k = 1,n-1
                do j = 1,n-1
                    do i = 1,n-1
                        ampl = ( coscoef(i)*coscoef(j)*coscoef(k)*8
     &                         + 12*( coscoef(i)*coscoef(j)
     &                              + coscoef(j)*coscoef(k)
     &                              + coscoef(i)*coscoef(k) )
     &                         + 28*( coscoef(i)+coscoef(j)+coscoef(k) )
     &                         - 128.d0 ) /30/h/h
                        ! divide by the symbol, and scale by 1/(n/2)^3.
                        f(i,j,k) = f(i,j,k)*(1.d0 + h*h*ampl/12)
     &                  /ampl/nh/nh/nh
                    enddo
                enddo
            enddo
        elseif( op .EQ. 2 )then
            ! Standard 7-point operator.
            do k = 1,n-1
                do j = 1,n-1
                    do i = 1,n-1
                        ampl = (coscoef(i)*2 + coscoef(j)*2
     &                         + coscoef(k)*2 - 6.d0)/h/h
                        ! divide by the symbol, and scale by 1/(n/2)^3.
                        f(i,j,k) = f(i,j,k)/ampl/nh/nh/nh
                    enddo
                enddo
            enddo
        else
            write(6,*) 'error: poisson3d(): unknown op argument: ',op
            stop 'poisson3'
        endif

c Compute inverse sine transform.        

        call sfft3d(f,n)

        return
        end
