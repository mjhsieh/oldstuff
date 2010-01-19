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
c RESPA partial force routines

      subroutine force_stretch(vflag)
      include "lammps.h"
      include "mpif.h"
      integer vflag
      double precision time1,time2

      call zero_force

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
      time_bond = time_bond + time2-time1

      if (newton_bond.eq.1) then
          time1 = mpi_wtime()
          call reverse_comm
          time2 = mpi_wtime()
          time_fcomm = time_fcomm + time2-time1
      endif

      call copy_force(f_stretch,vir_stretch)

      return
      end

      
      subroutine force_intra(vflag)
      include "lammps.h"
      include "mpif.h"
      integer vflag
      double precision time1,time2,time3,time4

      call zero_force

      time1 = mpi_wtime()

      if (nangles.gt.0) then
        if (anglestyle.eq.1) then
          call angle_harmonic(vflag)
        else if (anglestyle.eq.2) then
          call angle_class2(vflag)
        endif
      endif
      
      time2 = mpi_wtime()
      
      if (ndihedrals.gt.0) then
        if (dihedstyle.eq.1) then
          call dihedral_harmonic(vflag)
        else if (dihedstyle.eq.2) then
          call dihedral_class2(vflag)
        endif
      endif

      time3 = mpi_wtime()

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

      time4 = mpi_wtime()

      time_angle = time_angle + time2-time1
      time_dihedral = time_dihedral + time3-time2
      time_improper = time_improper + time4-time3

      if (newton_bond.eq.1) then
        time1 = mpi_wtime()
        call reverse_comm
        time2 = mpi_wtime()
        time_fcomm = time_fcomm + time2-time1
      endif

      call copy_force(f_intra,vir_intra)

      return
      end


      subroutine force_short(vflag)
      include "lammps.h"
      include "mpif.h"
      integer vflag
      double precision time1,time2

      call zero_force

      time1 = mpi_wtime()

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

      time2 = mpi_wtime()
      time_nonbond = time_nonbond + time2-time1

      if (newton_nonbond.eq.1) then
        time1 = mpi_wtime()
        call reverse_comm
        time2 = mpi_wtime()
        time_fcomm = time_fcomm + time2-time1
      endif

      call copy_force(f_short,vir_short)

      return
      end


      subroutine force_long(vflag)
      include "lammps.h"
      include "mpif.h"
      integer vflag
      double precision time1,time2

      call zero_force

      time1 = mpi_wtime()

      if (coulstyle.eq.3) then
        call ewald(vflag)
      else if (coulstyle.eq.4) then
        call pppm(vflag)
      endif

      time2 = mpi_wtime()
      time_long = time_long + time2-time1

      call copy_force(f_long,vir_long)

      return
      end 
