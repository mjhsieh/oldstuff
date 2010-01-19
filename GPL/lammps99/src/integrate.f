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
c timestep integrator

      subroutine integrate
      include "lammps.h"
      include "mpif.h"

      integer vflag
      dimension p_omega_dt(3)
      
      call mpi_barrier(mpi_comm_world,ierror)
      time_loop = mpi_wtime()

      do itime = 1,nsteps

        ntimestep = ntimestep + 1

c should virial be computed this timestep

        vflag = 0
        if (ensemble.eq.3) then
          vflag = 1
        else if (nthermo_next.eq.ntimestep) then
          vflag = 1
        endif

c 1st half of update
        
        if (ensemble.eq.1) then
          call nve_v(f,dthalf)
        else if (ensemble.eq.2) then
          call nvt_eta(dthalf)
          call nvt_v(f,dthalf,1)
        else if (ensemble.eq.3) then
          call nvt_eta(dthalf)
          call npt_omegadot(dthalf,p_omega_dt)
          call npt_v(f,dthalf,1)
          call npt_x(dthalf,2,p_omega_dt)
        endif

        call nve_x

c neighboring

        call check_neighbor(nflag)

        if (nflag.eq.0) then
#ifdef SYNC
          call mpi_barrier(mpi_comm_world,ierror)
#endif
          time1 = mpi_wtime()
          call communicate
          time2 = mpi_wtime()
          time_comm = time_comm + time2-time1
        else
#ifdef SYNC
          call mpi_barrier(mpi_comm_world,ierror)
#endif
          time1 = mpi_wtime()
          if (perflagx+perflagy+perflagz.gt.0.or.ensemble.eq.3)
     $         call setup_box
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

c zero forces and virial

        call zero_force

c bonded forces

#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()

        if (nbonds.gt.0) then
          if (bondstyle.eq.1) then
            call bond_harmonic(vflag)
          else if (bondstyle.eq.2) then
            call bond_fene_standard(vflag)
          else if (bondstyle.eq.3) then
            call bond_fene_shift(vflag)
          else if (bondstyle.eq.4) then
            call bond_nonlinear(vflag)
          else
            call bond_class2(vflag)
          endif
        endif

        time2 = mpi_wtime()

        if (nangles.gt.0) then
          if (anglestyle.eq.1) then
            call angle_harmonic(vflag)
          else if (anglestyle.eq.2) then
            call angle_class2(vflag)
          endif
        endif

        time3 = mpi_wtime()

        if (ndihedrals.gt.0) then
          if (dihedstyle.eq.1) then
            call dihedral_harmonic(vflag)
          else if (dihedstyle.eq.2) then
            call dihedral_class2(vflag)
          endif
        endif

        time4 = mpi_wtime()

        if (nimpropers.gt.0) then
          if (improstyle.eq.1) then
            call improper_harmonic(vflag)
          else if (improstyle.eq.2) then
            call improper_cvff(vflag)
          else if (improstyle.eq.3) then
            call improper_class2(vflag)
            call angleangle_class2(vflag)
          endif
        endif

#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time5 = mpi_wtime()

c nonbond forces

        if (nonstyle.eq.0) then
          if (coulstyle.eq.1) then
            call lj_cut_coul_cut(vflag)
          else if (coulstyle.eq.2) then
            call lj_cut_coul_smooth(vflag)
          else if (coulstyle.ge.3) then
            call lj_cut_coul_long(vflag)
          endif
        else if (nonstyle.eq.1) then
          if (coulstyle.eq.0) then
            call lj_cut(vflag)
          else if (coulstyle.eq.1) then
            call lj_cut_coul_cut(vflag)
          else if (coulstyle.eq.2) then
            call lj_cut_coul_smooth(vflag)
          else if (coulstyle.ge.3) then
            call lj_cut_coul_long(vflag)
          endif
        else if (nonstyle.eq.2) then
          if (coulstyle.eq.0) then
            call lj_smooth(vflag)
          else if (coulstyle.eq.1) then
            call lj_smooth_coul_cut(vflag)
          else if (coulstyle.eq.2) then
            call lj_smooth_coul_smooth(vflag)
          else if (coulstyle.ge.3) then
            call lj_smooth_coul_long(vflag)
          endif
        else if (nonstyle.eq.3) then
          call lj_shift(vflag)
        else if (nonstyle.eq.4) then
          call soft(vflag)
        else if (nonstyle.eq.5) then
          if (coulstyle.eq.0) then
            call ljclass2_cut(vflag)
          else if (coulstyle.eq.1) then
            call ljclass2_cut_coul_cut(vflag)
          else if (coulstyle.ge.3) then
            call ljclass2_cut_coul_long(vflag)
          endif
        endif

#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time6 = mpi_wtime()

c long-range forces

        if (coulstyle.eq.3) then
          call ewald(vflag)
        else if (coulstyle.eq.4) then
          call pppm(vflag)
        endif
        
#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time7 = mpi_wtime()

c timers

        time_bond = time_bond + time2-time1
        time_angle = time_angle + time3-time2
        time_dihedral = time_dihedral + time4-time3
        time_improper = time_improper + time5-time4
        time_nonbond = time_nonbond + time6-time5
        time_long = time_long + time7-time6
        
c reverse force communication

        if (newton.ge.1) then
          time1 = mpi_wtime()
          call reverse_comm
          time2 = mpi_wtime()
          time_fcomm = time_fcomm + time2-time1
        endif

c global Langevin temperature control

        if (tempflag.eq.3) call langevin

c constraints

        if (nfixes.gt.0) call fix_apply

c 2nd half of update

        if (ensemble.eq.1) then
          call nve_v(f,dthalf)
        else if (ensemble.eq.2) then
          call nvt_v(f,dthalf,2)
          call temperature
          call nvt_eta(dthalf)
        else if (ensemble.eq.3) then
          call npt_v(f,dthalf,2)
          call npt_x(dthalf,1,p_omega_dt)
          call temperature
          call nvt_eta(dthalf)
          call pressure
          call npt_omegadot(dthalf,p_omega_dt)
        endif

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
c constant NVE ensemble

      subroutine nve_v(f_partial,dt_partial)
      include "lammps.h"
      dimension f_partial(3,*)

      i = atompnt
      do ii = 1,nlocal
        dtmass = dt_partial/mass(type(i))
        v(1,i) = v(1,i) + dtmass*f_partial(1,i)
        v(2,i) = v(2,i) + dtmass*f_partial(2,i)
        v(3,i) = v(3,i) + dtmass*f_partial(3,i)
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
      subroutine nve_x
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        x(1,i) = x(1,i) + dt*v(1,i)
        if (perflagx.eq.0) then
          if (x(1,i).ge.xboundhi) then
            x(1,i) = x(1,i) - xprd
            true(i) = true(i) + 1
          endif
          if (x(1,i).lt.xboundlo) then
            x(1,i) = x(1,i) + xprd
            true(i) = true(i) - 1
          endif
        endif
        x(2,i) = x(2,i) + dt*v(2,i)
        if (perflagy.eq.0) then
          if (x(2,i).ge.yboundhi) then
            x(2,i) = x(2,i) - yprd
            true(i) = true(i) + 1000
          endif
          if (x(2,i).lt.yboundlo) then
            x(2,i) = x(2,i) + yprd
            true(i) = true(i) - 1000
          endif
        endif
        x(3,i) = x(3,i) + dt*v(3,i)
        if (perflagz.eq.0) then
          if (x(3,i).ge.zboundhi) then
            x(3,i) = x(3,i) - zprd
            true(i) = true(i) + 1000000
          endif
          if (x(3,i).lt.zboundlo) then
            x(3,i) = x(3,i) + zprd
            true(i) = true(i) - 1000000
          endif
        endif
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c constant NVT ensemble

      subroutine nvt_eta(dt_partial)
      include "lammps.h"

      t_target = t_start + float(itime)/nsteps*(t_stop-t_start)
      f_eta = t_freq*t_freq * (t_current/t_target - 1.0)
      eta_dot = eta_dot + f_eta*dt_partial

      return
      end


c -------------------------------------------------------------------------

      subroutine nvt_v(f_partial,dt_partial,nvt_flag)
      include "lammps.h"
      dimension f_partial(3,*)

c nvt_flag = 1 for 1st half of 'Verlet' velocity calculations
c nvt_flag = 2 for 2nd half of 'Verlet' velocity calculations

      factor = exp(-dt_partial*eta_dot)

      if (nvt_flag.eq.1) then
         factor1 = 1.0
      else if (nvt_flag.eq.2) then
         factor1 = factor
      end if

      i = atompnt
      do ii = 1,nlocal
	dt_4 = dt_partial/(mass(type(i)))*factor1
        v(1,i) = v(1,i)*factor + dt_4*f_partial(1,i)
        v(2,i) = v(2,i)*factor + dt_4*f_partial(2,i)
        v(3,i) = v(3,i)*factor + dt_4*f_partial(3,i)
        i = list(i)
      enddo

      return
      end

c -------------------------------------------------------------------------
c constant NPT ensemble

      subroutine npt_omegadot(dt_partial,p_omega_dt)
      include "lammps.h"

      dimension p_omega_dt(3),f_omega(3)

      dtime = itime/nsteps
      denskt = (natoms*boltz*t_target) / (xprd*yprd*zprd)

c advance omega_dot in time

      p_target(1) = p_start(1) + dtime*(p_stop(1)-p_start(1))
      f_omega(1) = p_freq2(1)* (p_current(1)-p_target(1))/denskt 
      omega_dot(1) = omega_dot(1) + f_omega(1)*dt_partial
      p_omega_dt(1) = dt_partial*omega_dot(1)

      p_target(2) = p_start(2) + dtime*(p_stop(2)-p_start(2))
      f_omega(2) = p_freq2(2)* (p_current(2)-p_target(2))/denskt 
      omega_dot(2) = omega_dot(2) + f_omega(2)*dt_partial
      p_omega_dt(2) = dt_partial*omega_dot(2)

      p_target(3) = p_start(3) + dtime*(p_stop(3)-p_start(3))
      f_omega(3) = p_freq2(3)* (p_current(3)-p_target(3))/denskt 
      omega_dot(3) = omega_dot(3) + f_omega(3)*dt_partial
      p_omega_dt(3) = dt_partial*omega_dot(3)

      return
      end


c -------------------------------------------------------------------------

      subroutine npt_v(f_partial,dt_partial,npt_flag)
      include "lammps.h"

      dimension f_partial(3,*)
      dimension factor(3),factor1(3)

c nvt_flag = 1 for 1st part of 'Verlet' velocity calculations
c nvt_flag = 2 for 2nd part of 'Verlet' velocity calculations

      factor(1) = exp(-dt_partial*(eta_dot+omega_dot(1)))
      factor(2) = exp(-dt_partial*(eta_dot+omega_dot(2)))
      factor(3) = exp(-dt_partial*(eta_dot+omega_dot(3)))

      if (npt_flag.eq.1) then
         factor1(1) = 1.0
         factor1(2) = 1.0
         factor1(3) = 1.0
      else if (npt_flag.eq.2) then
         factor1(1) = factor(1)
         factor1(2) = factor(2)
         factor1(3) = factor(3)
      end if

      i = atompnt
      do ii = 1,nlocal
        dt_4 = dt_partial/(mass(type(i)))
        v(1,i) = v(1,i)*factor(1) + dt_4*f_partial(1,i)*factor1(1)
        v(2,i) = v(2,i)*factor(2) + dt_4*f_partial(2,i)*factor1(2)
        v(3,i) = v(3,i)*factor(3) + dt_4*f_partial(3,i)*factor1(3)
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------

      subroutine npt_x(dt_partial,npt_flag,p_omega_dt)
      include "lammps.h"
      include "mpif.h"

      real*8 cm(3),cmall(3)
      dimension p_omega_dt(3),dilation(3)

c propagate omega
c dilate box & particles: it is important use the same dilation
c   for the box & the particles.

c change the volume (~ omega)

      omega(1) = omega(1) + p_omega_dt(1)
      dilation(1) = exp(p_omega_dt(1))
      omega(2) = omega(2) + p_omega_dt(2)
      dilation(2) = exp(p_omega_dt(2))
      omega(3) = omega(3) + p_omega_dt(3)
      dilation(3) = exp(p_omega_dt(3))
      volume = exp(omega(1)+omega(2)+omega(3))

c center of mass USING TRUE POSITIONS

      cm(1) = 0.0
      cm(2) = 0.0
      cm(3) = 0.0
      i = atompnt
      do ii = 1,nlocal
        rmass = mass(type(i))
        ix = mod(true(i),1000) - 500
        iy = mod(true(i)/1000,1000) - 500
        iz = true(i)/1000000 - 500
        cm(1) = cm(1) + rmass*(x(1,i) + ix * xprd)
        cm(2) = cm(2) + rmass*(x(2,i) + iy * yprd)
        cm(3) = cm(3) + rmass*(x(3,i) + iz * zprd)
        i = list(i)
      enddo
      call mpi_allreduce(cm,cmall,3,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      cm(1) = cmall(1)/masssum
      cm(2) = cmall(2)/masssum
      cm(3) = cmall(3)/masssum

c dilate the cell about the center of mass 

      xboundlo = (xboundlo-cm(1))*dilation(1) + cm(1)
      xboundhi = (xboundhi-cm(1))*dilation(1) + cm(1)
      yboundlo = (yboundlo-cm(2))*dilation(2) + cm(2)
      yboundhi = (yboundhi-cm(2))*dilation(2) + cm(2)
      zboundlo = (zboundlo-cm(3))*dilation(3) + cm(3)
      zboundhi = (zboundhi-cm(3))*dilation(3) + cm(3)
      xprd = xboundhi - xboundlo
      yprd = yboundhi - yboundlo
      zprd = zboundhi - zboundlo
      xmc = xprd - cutforce
      ymc = yprd - cutforce
      zmc = zprd - cutforce

      i = atompnt
      do ii = 1,nlocal
        x(1,i) = cm(1) + (x(1,i)-cm(1))*dilation(1)
        x(2,i) = cm(2) + (x(2,i)-cm(2))*dilation(2)
        x(3,i) = cm(3) + (x(3,i)-cm(3))*dilation(3)
        i = list(i)
      enddo

c redo long-range coeffs because volume has changed
c TODO: do not do every call, if possible

      if (npt_flag.eq.2) then
        if (coulstyle.eq.3) then
          call ewald_coeff
        else if (coulstyle.eq.4) then
          call pppm_coeff(1)
        endif
      endif

      return
      end


c -------------------------------------------------------------------------
c Langevin thermostat 
 
      subroutine langevin
      include "lammps.h"
      
      target = t_start + float(itime)/nsteps*(t_stop-t_start)
      gconst = sqrt(24.0*boltz*target*t_freq/dt)

      i = atompnt
      do ii = 1,nlocal
        rmass = mass(type(i))
        gamma1 = -t_freq*rmass
        gamma2 = gconst*sqrt(rmass)
        f(1,i) = f(1,i) + gamma1*v(1,i) + gamma2*(ranmars(0)-0.5)
        f(2,i) = f(2,i) + gamma1*v(2,i) + gamma2*(ranmars(0)-0.5)
        f(3,i) = f(3,i) + gamma1*v(3,i) + gamma2*(ranmars(0)-0.5)
        i = list(i)
      enddo
 
      return
      end
