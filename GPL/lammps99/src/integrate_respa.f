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
c 4-level timestep integrator using RESPA
c outer loop (long) is conceptually equivalent to regular timestep
c  in standard integrator
c 3 inner-nested timestep loops (short, intra, stretch)
c  conceptually replace nve_x and reneighboring in standard integrator
c communicate routine is only called once - in inner loop which is the
c  only place positions are ever updated
c reverse communication occurs whenever forces are computed, within each
c  of the 4 force routines
c neighbor check only needed at nonbond short-range level
c  prior to force_short computation since bonded forces are correct 
c  no matter how far atoms move
c virial flag is only turned on for last inner loop iterations
c NVT integrator is coupled to nonbond short-range forces


      subroutine integrate_respa
      include "lammps.h"
      include "mpif.h"

      integer vflag,v1flag
      dimension p_omega_dt(3)

      call mpi_barrier(mpi_comm_world,ierror)
      time_loop = mpi_wtime()

c long-range loop
c 1st half of long-range update
c inner iteration on short,intra,bond
c long-range forces
c 2nd half of long-range update

      do itime = 1,nsteps

        ntimestep = ntimestep + 1

c should virial be computed this timestep

        vflag = 0
        if (ensemble.eq.3) then
          vflag = 1
        else if (nthermo_next.eq.ntimestep) then
          vflag = 1
        endif

        if (ensemble.eq.1) then
          call nve_v(f_long,dthalf_long)
        else if (ensemble.eq.2) then
          call nvt_eta(dthalf_long)
          call nvt_v(f_long,dthalf_long,1)
        else if (ensemble.eq.3) then
          call nve_v(f_long,dthalf_long)
        endif

c nonbond loop
c 1st half of nonbond update
c inner iteration on intra,bond
c neighboring
c nonbond forces
c 2nd half of nonbond update

        do jshort = 1,nshort 

          if (ensemble.eq.1) then
             call nve_v(f_short,dthalf_short)
          else if (ensemble.eq.2) then
             call nve_v(f_short,dthalf_short)
          else if (ensemble.eq.3) then
             call nvt_eta(dthalf_short)
             call npt_v(f_short,dthalf_short,1)
             call npt_omegadot(dthalf_short,p_omega_dt)
             call npt_x(dthalf_short,2,p_omega_dt)
          endif

c 3/4-body loop
c 1st half of 3/4-body update
c inner iteration on bond
c 3/4-body forces
c 2nd half of 3/4-body update

          do jintra = 1,nintra

            if (ensemble.eq.1) then
               call nve_v(f_intra,dthalf_intra)
            else if (ensemble.eq.2) then
               call nve_v(f_intra,dthalf_intra)
            else if (ensemble.eq.3) then
               call nve_v(f_intra,dthalf_intra)
            end if

c bond loop
c 1st half of bond update
c communicate new positions
c bond forces
c 2nd half of bond update

            do jstretch = 1,nstretch

               if (ensemble.eq.1) then
                  call nve_v(f_stretch,dthalf)
               else if (ensemble.eq.2) then
                  call nve_v(f_stretch,dthalf)
               else if (ensemble.eq.3) then
                  call nve_v(f_stretch,dthalf)
               end if
               call nve_x

#ifdef SYNC
               call mpi_barrier(mpi_comm_world,ierror)
#endif
               time1 = mpi_wtime()
               call communicate
               time2 = mpi_wtime()
               time_comm = time_comm + time2-time1

               v1flag = vflag
               if (jstretch.lt.nstretch) v1flag = 0
               call force_stretch(v1flag)

               if (ensemble.eq.1) then
                  call nve_v(f_stretch,dthalf)
               else if (ensemble.eq.2) then
                  call nve_v(f_stretch,dthalf)
               else if (ensemble.eq.3) then
                  call nve_v(f_stretch,dthalf)
               end if

            enddo

            v1flag = vflag
            if (jintra.lt.nintra) v1flag = 0
            call force_intra(v1flag)

            if (ensemble.eq.1) then
               call nve_v(f_intra,dthalf_intra)
            else if (ensemble.eq.2) then
               call nve_v(f_intra,dthalf_intra)
            else if (ensemble.eq.3) then
               call nve_v(f_intra,dthalf_intra)
            end if

          enddo

          call check_neighbor(nflag)

          if (nflag.eq.1) then
#ifdef SYNC
            call mpi_barrier(mpi_comm_world,ierror)
#endif
            time1 = mpi_wtime()
            if (perflagx+perflagy+perflagz.gt.0.or.ensemble.eq.3)
     $           call setup_box
            call exchange
            time2 = mpi_wtime()
            call borders
            time3 = mpi_wtime()
            call neighbor
            time4 = mpi_wtime()
            if (nbonds.gt.0) call neighbor_bond
            if (nangles.gt.0) call neighbor_angle
            if (ndihedrals.gt.0) call neighbor_dihedral
            if (nimpropers.gt.0) call neighbor_improper
            time5 = mpi_wtime()
            time_exch = time_exch + time2-time1
            time_comm = time_comm + time3-time2
            time_neigh1 = time_neigh1 + time4-time3
            time_neigh2 = time_neigh2 + time5-time4
          endif

          v1flag = vflag
          if (jshort.lt.nshort) v1flag = 0
          call force_short(v1flag)

          if (ensemble.eq.1) then
             call nve_v(f_short,dthalf_short)
          else if (ensemble.eq.2) then
             call nve_v(f_short,dthalf_short)
          else if (ensemble.eq.3) then
             call npt_v(f_short,dthalf_short,2)
             call npt_x(dthalf_short,2,p_omega_dt)
             call pressure
             call npt_omegadot(dthalf_short,p_omega_dt)
             call temperature
             call nvt_eta(dthalf_short)
          endif

        enddo

        call force_long(vflag)

        if (ensemble.eq.2) then
           call nvt_v(f_long,dthalf_long,2)
           call temperature
           call nvt_eta(dthalf_long)
        else if (ensemble.eq.3) then
           call nve_v(f_long,dthalf_long)
        else
           call nve_v(f_long,dthalf_long)
        end if
      
c global Langevin temperature control

c        if (tempflag.eq.3) call langevin
        
c constraints

c        if (nfixes.gt.0) call fix_apply

c global temperature control (rescale and replace)

        if (tempflag.eq.1) then
          if (mod(itime,t_every).eq.0) then
            call temperature
            target = t_start + float(itime)/nsteps*(t_stop-t_start)
            if (abs(target-t_current).gt.t_window)
     $           call temp_rescale(target)
          endif
        else if (tempflag.eq.2) then
          if (mod(itime,t_every).eq.0) then
            target = t_start + float(itime)/nsteps*(t_stop-t_start)
            call temp_replace(target)
          endif
        endif

c outputs

        if (ntimestep.eq.noutput_next) then
          if (ntimestep.eq.nthermo_next) call thermo(0)
          if (ntimestep.eq.ndumpatom_next) call dump_atom
          if (ntimestep.eq.ndumpvel_next) call dump_vel
          if (ntimestep.eq.ndumpforce_next) call dump_force
          if (ntimestep.eq.nrestart_next) call write_restart
          if (ntimestep.eq.ndiag_next) call diagnostic
          noutput_next = nthermo_next
          noutput_next = min(noutput_next,ndumpatom_next)
          noutput_next = min(noutput_next,ndumpvel_next)
          noutput_next = min(noutput_next,ndumpforce_next)
          noutput_next = min(noutput_next,nrestart_next)
          noutput_next = min(noutput_next,ndiag_next)
        endif

      enddo

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()
      time_loop = time1-time_loop

      return
      end


c -------------------------------------------------------------------------

      subroutine nvt_f(f_partial,dt_partial,nvt_flag)
      include "lammps.h"
      dimension f_partial(3,*)

c nvt_flag = 1 for 1st part of 'Verlet' velocity calculations
c nvt_flag = 2 for 2nd part of 'Verlet' velocity calculations

      factor = exp(-dthalf_long*eta_dot)

      if (nvt_flag.eq.1) then
         factor = 1.0 + factor
      else if (nvt_flag.eq.2) then
         factor = (1.0 + factor)/factor
      end if

      i = atompnt
      do ii = 1,nlocal
        dt_4 = dt_partial/(2.0*mass(type(i)))*factor
        v(1,i) = v(1,i) + dt_4*f_partial(1,i)
        v(2,i) = v(2,i) + dt_4*f_partial(2,i)
        v(3,i) = v(3,i) + dt_4*f_partial(3,i)
        i = list(i)
      enddo

      return
      end
