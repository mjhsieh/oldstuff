c LAMMPS 99 - Molecular Dynamics Simulator
c Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
c Steve Plimpton, sjplimp@sandia.gov
c
c Copyright (2003) Sandia Corporation.  Under the terms of Contract
c DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
c certain rights in this software.  This software is distributed under 
c the GNU General Public License.
c
c See the README file in the top-level LAMMPS directory.

c -------------------------------------------------------------------------
c end of timestepping

      subroutine finish
      include "lammps.h"
      include "mpif.h"

      integer ihisto(10),ihistotmp(10)

 900  format(' ',a,f15.6,f13.4)
 901  format(' ',a,f13.4,a,f13.4,a,f13.4,a)
 902  format(' ',a,10i5)

c timings stats

      time_other = time_loop - (time_nonbond+time_long+time_bond+
     $     time_angle+time_dihedral+time_improper+time_neigh1+
     $     time_neigh2+time_exch+time_comm+time_fcomm+time_io)

      rtmp = time_loop
      call mpi_allreduce(rtmp,time_loop,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      time_loop = time_loop/nprocs
      if (node.eq.0) then
        write (6,*) 'Loop time:',time_loop,
     $       ' on',nprocs,' procs for',natoms,' atoms'
        write (1,*) 'Loop time:',time_loop,
     $       ' on',nprocs,' procs for',natoms,' atoms'
      endif

      if (time_loop.eq.0.0) time_loop = 1.0

      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      call mpi_allreduce(time_nonbond,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Nbond time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Nbond time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_long,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Long  time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Long  time/%:',rtmp,rtmp/time_loop*100
      endif
      
      call mpi_allreduce(time_bond,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Bond  time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Bond  time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_angle,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Angle time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Angle time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_dihedral,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Dihed time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Dihed time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_improper,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Impro time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Impro time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_neigh1,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Nay-1 time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Nay-1 time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_neigh2,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Nay-2 time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Nay-2 time/%:',rtmp,rtmp/time_loop*100
      endif
      
      call mpi_allreduce(time_exch,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Exch  time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Exch  time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_comm,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Comm  time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Comm  time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_fcomm,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Fcomm time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Fcomm time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_io,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'I/O   time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'I/O   time/%:',rtmp,rtmp/time_loop*100
      endif
      
      call mpi_allreduce(time_other,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Other time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Other time/%:',rtmp,rtmp/time_loop*100
      endif
      
      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      call stats(time_nonbond,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nbond time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nbond time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_long,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Long  time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Long  time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      call stats(time_bond,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Bond  time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Bond  time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_angle,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Angle time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Angle time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_dihedral,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Dihed time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Dihed time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_improper,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Impro time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Impro time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_neigh1,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nay-1 time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nay-1 time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_neigh2,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nay-2 time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nay-2 time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      call stats(time_exch,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Exch  time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Exch  time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_comm,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Comm  time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Comm  time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_fcomm,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Fcomm time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Fcomm time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_io,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'I/O   time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'I/O   time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif
      
      call stats(time_other,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Other time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Other time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

c timing breakdown for Ewald

      if (coulstyle.eq.3) then

        if (node.eq.0) then
          write (6,*)
          write (6,*) 'Ewald timing info:'
          write (1,*)
          write (1,*) 'Ewald timing info:'
        endif

        do i = 1,kcount
          sfacrl(i) = 0.0
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,10
          call mpi_allreduce(sfacrl,sfacrl_all,kcount,
     $         mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        call mpi_allreduce(time_long,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        ave = rtmp/nprocs
        if (ave.eq.0.0) ave = 1.0
        rtmp = (time2-time1)/10 * 2

        if (node.eq.0) then
          write (6,*) '  Allreduce time per timestep =',rtmp
          write (6,*) '  Total allreduce time =',rtmp*nsteps
          write (6,*) '  Allreduce % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  Allreduce time per timestep =',rtmp
          write (1,*) '  Total allreduce time =',rtmp*nsteps
          write (1,*) '  Allreduce % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        if (ensemble.eq.3) then

          call mpi_barrier(mpi_comm_world,ierror)
          time1 = mpi_wtime()
          do i = 1,10
            call ewald_coeff
          enddo
          call mpi_barrier(mpi_comm_world,ierror)
          time2 = mpi_wtime()

          call mpi_allreduce(time_other,rtmp,1,
     $         mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
          ave = rtmp/nprocs
          if (ave.eq.0.0) ave = 1.0
          rtmp = (time2-time1)/10 * 2

          if (node.eq.0) then
            write (6,*) '  Setup time per timestep =',rtmp
            write (6,*) '  Total setup time =',rtmp*nsteps
            write (6,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
            write (1,*) '  Setup time per timestep =',rtmp
            write (1,*) '  Total setup time =',rtmp*nsteps
            write (1,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
          endif

        endif

      endif

c timing breakdown for PPPM

      if (coulstyle.eq.4) then

        if (node.eq.0) then
          write (6,*)
          write (6,*) 'PPPM timing info:'
          write (1,*)
          write (1,*) 'PPPM timing info:'
        endif

        call mpi_allreduce(time_long,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        ave = rtmp/nprocs
        if (ave.eq.0.0) ave = 1.0

        call mpi_allreduce(time_rho,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs

        if (node.eq.0) then
          write (6,*) '  Make_rho time =',rtmp
          write (6,*) '  Make_rho % of long time =',
     $         rtmp/ave*100.0
          write (1,*) '  Make_rho time =',rtmp
          write (1,*) '  Make_rho % of long time =',
     $         rtmp/ave*100.0
        endif

        call mpi_allreduce(time_poiss,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs

        if (node.eq.0) then
          write (6,*) '  Poisson time =',rtmp
          write (6,*) '  Poisson % of long time =',
     $         rtmp/ave*100.0
          write (1,*) '  Poisson time =',rtmp
          write (1,*) '  Poisson % of long time =',
     $         rtmp/ave*100.0
        endif

        call mpi_allreduce(time_field,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs

        if (node.eq.0) then
          write (6,*) '  Electric_field time =',rtmp
          write (6,*) '  Electric_field % of long time =',
     $         rtmp/ave*100.0
          write (1,*) '  Electric_field time =',rtmp
          write (1,*) '  Electric_field % of long time =',
     $         rtmp/ave*100.0
        endif

        do indx = 1,(nxhi_out-nxlo_out+1)*(nyhi_out-nylo_out+1)*
     $       (nzhi_out-nzlo_out+1)
          density_brick(indx) = 0.0
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,10
          call brick2fft(density_brick,density_fft)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        rtmp = (time2-time1)/10

        if (node.eq.0) then
          write (6,*) '  Brick2fft time per timestep =',rtmp
          write (6,*) '  Total Brick2fft time =',rtmp*nsteps
          write (6,*) '  Brick2fft % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  Brick2fft time per timestep =',rtmp
          write (1,*) '  Total brick2fft time =',rtmp*nsteps
          write (1,*) '  Brick2fft % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        do indx = 1,(nxhi_out-nxlo_out+1)*(nyhi_out-nylo_out+1)*
     $       (nzhi_out-nzlo_out+1)
          vdx_brick(indx) = 0.0
          vdy_brick(indx) = 0.0
          vdz_brick(indx) = 0.0
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,10
          call fillbrick(vdx_brick,vdy_brick,vdz_brick)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        rtmp = (time2-time1)/10

        if (node.eq.0) then
          write (6,*) '  Fillbrick time per timestep =',rtmp
          write (6,*) '  Total fillbrick time =',rtmp*nsteps
          write (6,*) '  Fillbrick % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  Fillbrick time per timestep =',rtmp
          write (1,*) '  Total fillbrick time =',rtmp*nsteps
          write (1,*) '  Fillbrick % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        do i = 1,nfft
          workvec1(i) = cmplx(0.0,0.0)
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,5
          call fft_3d(workvec1,workvec1,1,plan1_fft)
          call fft_3d(workvec1,workvec1,-1,plan2_fft)
          call fft_3d(workvec1,workvec1,-1,plan2_fft)
          call fft_3d(workvec1,workvec1,-1,plan2_fft)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        rtmp = (time2-time1)/5

        if (node.eq.0) then
          write (6,*) '  FFT time per timestep =',rtmp
          write (6,*) '  Total FFT time =',rtmp*nsteps
          write (6,*) '  FFT % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  FFT time per timestep =',rtmp
          write (1,*) '  Total FFT time =',rtmp*nsteps
          write (1,*) '  FFT % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        if (ensemble.eq.3) then

          call mpi_barrier(mpi_comm_world,ierror)
          time1 = mpi_wtime()
          do i = 1,10
            call pppm_coeff(1)
          enddo
          call mpi_barrier(mpi_comm_world,ierror)
          time2 = mpi_wtime()

          call mpi_allreduce(time_other,rtmp,1,
     $         mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
          ave = rtmp/nprocs
          if (ave.eq.0.0) ave = 1.0
          rtmp = (time2-time1)/10 * 2

          if (node.eq.0) then
            write (6,*) '  Setup time per timestep =',rtmp
            write (6,*) '  Total setup time =',rtmp*nsteps
            write (6,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
            write (1,*) '  Setup time per timestep =',rtmp
            write (1,*) '  Total setup time =',rtmp*nsteps
            write (1,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
          endif

        endif

      endif

c load-balance info

      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      r8tmp = nlocal
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nlocal:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nlocal:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = nother
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nother:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nother:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = nbondlocal
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nbonds:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nbonds:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = nanglelocal
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nangle:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nangle:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = ndihedlocal
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Ndihed:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Ndihed:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = nimprolocal
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nimpro:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nimpro:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = neightop
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Neighs:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = nslist(nswap+1) - 1
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nswaps:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nswaps:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call mpi_allreduce(neightop,itmp1,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)

      itmp2 = 0
      i = atompnt
      do ii = 1,nlocal
        itmp2 = itmp2 + num3bond(i)
        i = list(i)
      enddo

      r8tmp = itmp2
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)

      call mpi_allreduce(itmp2,itmp,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      itmp2 = itmp

      if (node.eq.0) then
        write (6,901) 'Nspecs:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nspecs:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      if (node.eq.0) then
        write (6,*)
        write (6,*) 'Ave neighs/atom =',float(itmp1)/natoms
        write (6,*) 'Ave nspecs/atom =',float(itmp2)/natoms
        write (6,*) 'Number of reneighborings =',numneigh
        write (6,*) 'Dangerous reneighborings =',ndanger
        write (6,*)
        write (1,*)
        write (1,*) 'Ave neighs/atom =',float(itmp1)/natoms
        write (1,*) 'Ave nspecs/atom =',float(itmp2)/natoms
        write (1,*) 'Number of reneighborings =',numneigh
        write (1,*) 'Dangerous reneighborings =',ndanger
        write (1,*)
      endif

      itmp = max_nlocal
      call mpi_allreduce(itmp,max_nlocal,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_bond
      call mpi_allreduce(itmp,max_bond,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_angle
      call mpi_allreduce(itmp,max_angle,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_dihed
      call mpi_allreduce(itmp,max_dihed,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_impro
      call mpi_allreduce(itmp,max_impro,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_neigh
      call mpi_allreduce(itmp,max_neigh,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_exch
      call mpi_allreduce(itmp,max_exch,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_slist
      call mpi_allreduce(itmp,max_slist,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_swap
      call mpi_allreduce(itmp,max_swap,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_exch
      call mpi_allreduce(itmp,max_exch,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = mbinx*mbiny*mbinz
      call mpi_allreduce(itmp,max_bin,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      
      max_bondlocal = 0
      max_anglelocal = 0
      max_dihedlocal = 0
      max_improlocal = 0
      i = atompnt
      do ii = 1,nlocal
        if (numbond(i).gt.max_bondlocal) max_bondlocal = numbond(i)
        if (numangle(i).gt.max_anglelocal) max_anglelocal = numangle(i)
        if (numdihed(i).gt.max_dihedlocal) max_dihedlocal = numdihed(i)
        if (numimpro(i).gt.max_improlocal) max_improlocal = numimpro(i)
        i = list(i)
      enddo

      itmp = max_bondlocal
      call mpi_allreduce(itmp,max_bondlocal,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_anglelocal
      call mpi_allreduce(itmp,max_anglelocal,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_dihedlocal
      call mpi_allreduce(itmp,max_dihedlocal,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_improlocal
      call mpi_allreduce(itmp,max_improlocal,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (node.eq.0) then
        
        write (1,*) 'Max # of local atoms =',max_nlocal,
     $       ' out of',maxown
        write (1,*) 'Max # of other atoms =',max_nother,
     $       ' out of',maxother
        write (1,*) 'Max # of bonds =',max_bond,
     $       ' out of',maxbond
        write (1,*) 'Max # of angles =',max_angle,
     $       ' out of',maxangle
        write (1,*) 'Max # of dihedrals =',max_dihed,
     $       ' out of',maxdihed
        write (1,*) 'Max # of impropers =',max_impro,
     $       ' out of',maximpro
        write (1,*) 'Max # of bonds/atom =',max_bondlocal,
     $       ' out of',maxbondper
        write (1,*) 'Max # of angles/atom =',max_anglelocal,
     $       ' out of',maxangleper
        write (1,*) 'Max # of dihedrals/atom =',max_dihedlocal,
     $       ' out of',maxdihedper
        write (1,*) 'Max # of impropers/atom =',max_improlocal,
     $       ' out of',maximproper
        write (1,*) 'Max # of neighbors =',max_neigh,
     $       ' out of',maxneigh*maxown
        write (1,*) 'Max used in exchange buffer =',max_exch,
     $       ' out of',maxexch
        write (1,*) 'Max sent in all swaps =',max_slist,
     $       ' out of',maxsend
        write (1,*) 'Max sent in one swap =',max_swap,
     $       ' out of',maxsendone
        if (max_bin.gt.0) write (1,*) 'Max # of bins =',max_bin,
     $       ' out of',maxbin
        write (1,*)
        write (1,*) '# of swaps =',nswap,
     $       ' Needs =',need(1),need(2),need(3)
        write (1,*) 'Cutneigh =',cutneigh,' Cut/Box =',
     $       cutneigh*pgrid(1)/xprd,cutneigh*pgrid(2)/yprd,
     $       cutneigh*pgrid(3)/zprd
        
      endif

c end-of-run clean-up

c reset special neighs after Ewald/PPPM

      if (coulstyle.ge.3) then
        special(1) = special(1) - 2.0
        special(2) = special(2) - 2.0
        special(3) = special(3) - 2.0
      endif

c free up PPPM plans

      if (coulstyle.eq.4) then
        call fft_3d_destroy_plan(plan1_fft)
        call fft_3d_destroy_plan(plan2_fft)
        call remap_3d_destroy_plan(plan_remap)
      endif

      return
      end
