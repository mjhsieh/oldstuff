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
c Optimization routines authored by Todd Plantenga, Sandia National Labs
c

C**********************************************************************
C
C   FILE min_support.f
C
C   Externally callable routines:
C      minimize
C      Compute_f_g
C      Get_Smallest_Box_Dim
C      Get_Atom_Positions
C      Put_Atom_Positions
C      Compute_Inf_Norm_PAR
C      Compute_L2_Norm_PAR
C      Compute_ddot_PAR
C      Normalize_ddot
C      Find_Biggest_Jump
C      Repartitn_and_Reneighbor
C
C   Local routines (in order)
C      Write_Energy_Vals
C      Write_Opt_Timer
C
C   Code development:
C      21 Oct 96 - Originated by T. Plantenga, Sandia National Labs.
C**********************************************************************


C**********************************************************************
      SUBROUTINE minimize
C
C       orig:  T. Plantenga   4 Aug 97
C
C
C   This routine is called from "lammps.f" to find a molecular
C   configuration that gives a local minimum of the potential energy.
C
C   The following LAMMPS input command lines pass in relevant parameters:
C   "neighbor      real int int int int" --> skin, neighstyle, neighfreq,
C                                             neighdelay, neightrigger
C   "restart    int ..."                 --> nrestart
C   "min file      string"               --> opt_outfile, optfileflag
C   "min style     string"               --> optstyle (always 1 for hftn)
C   "minimize  real int int"             --> opt_stop_tol, opt_max_iters,
C                                             opt_max_fns
C
C   The passed algorithm control parameters have the following effects:
C        neighfreq     = exchange & reneighbor after this many acc steps
C        neighdelay    = delay reneighboring until after this many acc steps
C        nrestart      = write restart file after this many evals (0=none)
C        opt_outfile   = output filename
C        optfileflag   = 0 if no output, 1 if output
C        opt_stop_tol  = applied to ||grad||_inf
C        opt_max_iters = maximum number of outer iterations
C        opt_max_fns   = maximum number of force/energy evaluations
C
C   At the end of execution, LAMMPS dumps some timing information.
C   For minimization, these mean:
C        time_loop      - total time in minimize
C        time_bond      - computing forces due to covalent bonds
C        time_angle     - computing forces due to angles
C        time_dihedral  - computing forces due to dihedral angles
C        time_improper  - computing forces due to improper angles
C        time_nonbond   - computing forces due to long-range interactions
C        time_fcomm     - comm for updating force calculations
C        time_comm      - comm for updating atom positions
C        time_exch      - comm for shifting atoms to new processors
C        time_neigh1    - computing new neighbor lists
C        time_neigh2    - verifying new neighbor lists
C        opt_time1      - total time in Compute_f_g
C        opt_time2      - Compute_f_g doing nonbonded energy
C        opt_time3      - synchronous parallel BLAS
C
C   External subroutines called:
C        error                   "misc.f"
C        dump_force              "output.f"
C**********************************************************************

      INCLUDE "lammps.h"
      INCLUDE "mpif.h"

C     ---- LOCAL VARIABLES.
      INTEGER           outUnit
      INTEGER           procNum
      INTEGER           totalN, localN
      INTEGER           num_fns
      DOUBLE PRECISION  endTime
      REAL*8            rTmp
      DOUBLE PRECISION  avgTime, maxTime
      INTEGER           iflag

C----------------------------------------------------------------------

      procNum = node

C     ---- OPEN FILE FOR OPTIMIZATION OUTPUT (ONLY PROCESSOR ZERO).
      IF ((procNum .EQ. 0) .AND. optfileflag.eq.1) THEN
         outUnit = 16
         OPEN (UNIT=outUnit, FILE=opt_outfile, ERR=900,
     $        status='unknown')
         close (outUnit,status='delete')
         OPEN (UNIT=outUnit, FILE=opt_outfile, ERR=900,
     $        status='new')
      ELSE
         outUnit = 0
      ENDIF

C     ---- DISPLAY MINIMIZATION PARAMETERS.
      IF (procNum .EQ. 0) THEN
         WRITE (6,*)
         WRITE (6,*)
         WRITE (6,600) nprocs, pgrid(1), pgrid(2), pgrid(3)
         IF (optstyle .EQ. 1)  WRITE (6,601)
         WRITE (6,602)
         WRITE (6,603) opt_stop_tol
         WRITE (6,604) opt_max_fns
         WRITE (6,605) opt_max_iters
         WRITE (6,606) neighfreq
         WRITE (6,607) neighdelay

         IF (outUnit .GT. 0) THEN
            WRITE (outUnit,600) nprocs, pgrid(1), pgrid(2), pgrid(3)
            IF (optstyle .EQ. 1)  WRITE (outUnit,601)
            WRITE (outUnit,602)
            WRITE (outUnit,603) opt_stop_tol
            WRITE (outUnit,604) opt_max_fns
            WRITE (outUnit,605) opt_max_iters
            WRITE (outUnit,606) neighfreq
            WRITE (outUnit,607) neighdelay
         ENDIF
      ENDIF

C     ---- CHECK THE SMOOTHNESS OF THE LONG RANGE FORCE FIELDS.
      IF (coulstyle .EQ. 2) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,610) cutcoulint, cutcoul
            IF (outUnit .GT. 0)  WRITE (outUnit,610) cutcoulint, cutcoul
         ENDIF
      ELSE
         IF (procNum .EQ. 0) THEN
            WRITE (6,611)
            IF (outUnit .GT. 0)  WRITE (outUnit,611)
         ENDIF
      ENDIF
      IF (nonstyle .EQ. 1) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,612) cutlj
            IF (outUnit .GT. 0)
     +         WRITE (outUnit,612) cutlj
         ENDIF
      ELSE IF (nonstyle .EQ. 2) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,613) cutljinterior, cutlj
            IF (outUnit .GT. 0)
     +         WRITE (outUnit,613) cutljinterior, cutlj
         ENDIF
      ELSE
         IF (procNum .EQ. 0) THEN
            WRITE (6,614)
            IF (outUnit .GT. 0)  WRITE (outUnit,614)
         ENDIF
      ENDIF

C     ---- DUMP THE INITIAL POTENTIAL ENERGY (COMPUTED IN "start.f").
      CALL Write_Energy_Vals (outUnit, procNum)

C     ---- CREATE VARIABLES FOR THE MINIMIZATION ALGORITHM.
      totalN = 3 * natoms
      localN = 3 * nlocal

C     ---- CALL THE MINIMIZATION ALGORITHM.
      opt_time1 = 0.0D0
      opt_time2 = 0.0D0
      opt_time3 = 0.0D0
      time_loop = mpi_wtime()
      IF (optstyle .EQ. 1) THEN
         CALL Min_HFTN (procNum, outUnit, totalN, localN,
     +                  opt_stop_tol, opt_max_fns, opt_max_iters,
     +                  neighfreq, neighdelay, num_fns, nrestart)
      ENDIF
      endTime = mpi_wtime()
      time_loop = endTime - time_loop

C     ---- PRINT OUT MINIMIZATION TIMES.
      IF (procNum .EQ. 0) THEN
         WRITE (6,620)
         IF (outUnit .GT. 0)  WRITE (outUnit,620)
      ENDIF
      CALL Write_Opt_Timer (time_loop, nprocs, num_fns,
     +                      outUnit, procNum)
      IF (procNum .EQ. 0) THEN
         WRITE (6,621)
         IF (outUnit .GT. 0)  WRITE (outUnit,621)
      ENDIF
      CALL Write_Opt_Timer (opt_time1, nprocs, num_fns,
     +                      outUnit, procNum)
      IF (procNum .EQ. 0) THEN
         WRITE (6,622)
         IF (outUnit .GT. 0)  WRITE (outUnit,622)
      ENDIF
      CALL Write_Opt_Timer (opt_time2, nprocs, num_fns,
     +                      outUnit, procNum)
      IF (procNum .EQ. 0) THEN
         WRITE (6,623)
         IF (outUnit .GT. 0)  WRITE (outUnit,623)
      ENDIF
      CALL Write_Opt_Timer (opt_time3, nprocs, num_fns,
     +                      outUnit, procNum)

C     ---- DUMP THE FINAL POTENTIAL ENERGY.
      CALL Write_Energy_Vals (outUnit, procNum)

C     ---- DUMP THE FINAL FORCE VECTOR IF REQUESTED IN THE COMMAND FILE.
      IF (ndumpforce .GT. 1)  CALL dump_force

C     ---- CLOSE THE OUTPUT FILE.
      IF (outUnit .GT. 0)  CLOSE (outUnit)

c     ---- write final restart file if requested
      if (nrestart.gt.0) call write_restart

      RETURN

C----------------------------------------------------------------------
C   THE REMAINING LINES OF CODE ARE I/O STATEMENTS REFERENCED EARLIER.
C----------------------------------------------------------------------

 900  CONTINUE
      CALL error ('Could not open min file for output')

 600  FORMAT (1X, I4, ' procs gridded as x = ', I4,
     +        ', y = ', I4, ', z = ', I4 /)

 601  FORMAT ('Find local energy minimum using truncated Newton CG')
 602  FORMAT ('  Stop minimization if:')
 603  FORMAT ('       ||grad||_inf .LE. ', 1PE9.3, '            OR')
 604  FORMAT ('       max number force calculations > ', I6, ' OR')
 605  FORMAT ('       max number iterations         > ', I6)
 606  FORMAT ('  Exchange & reneighbor after ', I6, ' accepted steps' /)
 607  FORMAT ('  Exch & reneigh delay = ', I6, ' accepted steps' /)

 610  FORMAT ('Coulomb forces smooth, cut off between ',
     +        F8.2, ' and ', F8.2, ' Angstroms')
 611  FORMAT ('Warning - smooth Coulomb cutoffs not used')
 612  FORMAT ('L-J     forces sharply cut off at      ',
     +        F8.2, ' Angstroms')
 613  FORMAT ('L-J     forces decay to zero between   ', F8.2, ' and ',
     +        F8.2, ' Angstroms -- discontinuous at inner edge')
 614  FORMAT ('Warning - L-J cutoffs not used')

 620  FORMAT ('CPU time for minimization, per proc (secs):')
 621  FORMAT ('CPU time for Compute_f_g, per proc (secs):')
 622  FORMAT ('CPU time for energy in Compute_f_g, per proc (secs):')
 623  FORMAT ('CPU time for synchronous BLAS, per proc (secs):')

      END


C**********************************************************************
      SUBROUTINE Compute_f_g
     +   (obj_f, g_vec_size, g_vec)
C
C       orig: 22 Oct 96  T. Plantenga
C       mod:   8 May 97  (TDP)  ignore constrained forces
C       mod:  27 Jul 97  (TDP)  adapted for LAMMPS V4.0
C
C   Input/output variable definitions:
C       obj_f            O  potential energy for minimization
C       g_vec_size      I
C       g_vec            O  gradient vector
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
      INCLUDE "mpif.h"
C
      DOUBLE PRECISION  obj_f
      INTEGER           g_vec_size
      DOUBLE PRECISION  g_vec(g_vec_size)
C
C
C   This routine consists of calls extracted from "integrate" and
C   "thermo" that compute interatomic forces and potential energy.
C   Calculations are done at the current atom positions stored in
C   the global variable "x".
C
C   LAMMPS V4.0 combines the computation of energy and forces for
C   covalent bond terms, but not for the non-bonded interactions.
C   This subroutine should support all bonded and non-bonded model
C   "styles", but only the following have been tested:
C        harmonic (bond, angle, dihedral, improper)
C        lj/cutoff, lj/smooth
C        coulomb/cutoff, coul/smooth
C
C   Force computations are done in parallel, but without interprocessor
C   communication until the very end of this routine.
C   
C   External subroutines called:
C        zero_force              "misc.f"
C        error                   "misc.f"
C        bond_harmonic           "force_bond.f"
C        bond_fene_standard      "force_bond.f"
C        bond_fene_shift         "force_bond.f"
C        bond_nonlinear          "force_bond.f"
C        bond_class2             "force_bond.f"
C        angle_harmonic          "force_many.f"
C        angle_class2            "force_class2.f"
C        dihedral_harmonic       "force_many.f"
C        dihedral_class2         "force_class2.f"
C        improper_harmonic       "force_many.f"
C        improper_class2         "force_class2.f"
C        angleangle_class2       "force_class2.f"
C        lj_cut_coul_cut         "force.f"
C        lj_cut_coul_long        "force.f"
C        lj_cut                  "force.f"
C        lj_smooth               "force.f"
C        lj_smooth_coul_cut      "force.f"
C        lj_smooth_coul_long     "force.f"
C        lj_shift                "force.f"
C        soft                    "force.f"
C        ljclass2_cut            "force.f"
C        ljclass2_cut_coul_cut   "force.f"
C        ljclass2_cut_coul_long  "force.f"
C        lj_cut_coul_smooth      "force.f"
C        lj_smooth_coul_smooth   "force.f"
C        ewald                   "ewald.f"
C        pppm                    "pppm.f"
C        reverse_comm            "communicate.f"
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER           i, j, k
      INTEGER           iflag
      DOUBLE PRECISION  time1, time2, time3, time4, time5, time6, time7
      DOUBLE PRECISION  time_top
      DOUBLE PRECISION  e_potential, rnorm

C----------------------------------------------------------------------
c        i = atompnt
c        do k = 1, nlocal
c         write (6,900) node, k, x(1,i), x(2,i), x(3,i)
c 900     format ('TBD node ', I1, ' has atom ', I3,
c     +           ' at [', F10.5, 1X, F10.5, 1X, F10.5, ']')
c         i = list(i)
c        enddo

C     ---- SET IFLAG TO ZERO TO DISABLE VIRIAL CALCULATIONS.
      iflag = 0

C     ---- ZERO THE GLOBAL FORCE VECTOR (SEE "integrate").
      CALL zero_force ()

C----------------------------------------------------------------------
C   COMPUTE COVALENT BOND FORCES AND ENERGIES SIMULTANEOUSLY (TAKEN
C   FROM "integrate").  RESULTS ARE STORED AS SIDE EFFECTS IN THE
C   GLOBAL VARIABLES:
C        FORCES   - f(1..3,1..maxatom)
C        ENERGIES - e_bond
C                   e_angle
C                   e_dihedral
C                   e_improper
C----------------------------------------------------------------------

C     ---- THE barrier CALL GIVES MORE ACCURATE TIMING INFORMATION.
      CALL mpi_barrier(mpi_comm_world,ierror)

      time1 = mpi_wtime()
      time_top = time1

C     ---- BOND LENGTH TERMS.
      IF (nbonds .GT. 0) THEN
         IF (bondstyle .EQ. 1) THEN
            CALL bond_harmonic (iflag)
         ELSE IF (bondstyle .EQ. 2) THEN
            CALL bond_fene_standard (iflag)
         ELSE IF (bondstyle .EQ. 3) THEN
            CALL bond_fene_shift (iflag)
         ELSE IF (bondstyle .EQ. 4) THEN
            CALL bond_nonlinear (iflag)
         ELSE
            CALL bond_class2 (iflag)
         ENDIF
      ENDIF

      time2 = mpi_wtime()

C     ---- BOND ANGLE TERMS.
      IF (nangles .GT. 0) THEN
         IF (anglestyle .EQ. 1) THEN
            CALL angle_harmonic (iflag)
         ELSE IF (anglestyle .EQ. 2) THEN
            CALL angle_class2 (iflag)
         ENDIF
      ENDIF

      time3 = mpi_wtime()

C     ---- DIHEDRAL ANGLE TERMS.
      IF (ndihedrals .GT. 0) THEN
         IF (dihedstyle .EQ. 1) THEN
            CALL dihedral_harmonic (iflag)
         ELSE IF (dihedstyle .EQ. 2) THEN
            CALL dihedral_class2 (iflag)
         ENDIF
      ENDIF

      time4 = mpi_wtime()

C     ---- IMPROPER ANGLE TERMS.
      IF (nimpropers .GT. 0) THEN
         IF (improstyle .EQ. 1) THEN
            CALL improper_harmonic (iflag)
         ELSE IF (improstyle .EQ. 2) then
            call improper_cvff(iflag)
         ELSE IF (improstyle .EQ. 3) THEN
            CALL improper_class2 (iflag)
            CALL angleangle_class2 (iflag)
         ENDIF
      ENDIF

C----------------------------------------------------------------------
C   COMPUTE NON-BOND FORCES (TAKEN FROM "integrate").  THESE ROUTINES
C   DO NOT SIMULTANEOUSLY CALCULATE ENERGY.  RESULTS APPEAR IN THE
C   GLOBAL VARIABLE "f(1..3,1..maxatom)".
C----------------------------------------------------------------------

C     ---- THE barrier CALL GIVES MORE ACCURATE TIMING INFORMATION.
      CALL mpi_barrier(mpi_comm_world,ierror)

      time5 = mpi_wtime()

C     ---- LENNARD-JONES (AKA VAN DER WAALS) FORCES.
      IF (nonstyle .EQ. 0) THEN
         IF (coulstyle .EQ. 1) THEN
            CALL lj_cut_coul_cut (iflag)
         ELSE IF (coulstyle .EQ. 2) THEN
            CALL lj_cut_coul_smooth (iflag)
         ELSE IF (coulstyle .GE. 3) THEN
            CALL lj_cut_coul_long (iflag)
         ENDIF
      ELSE IF (nonstyle .EQ. 1) THEN
         IF (coulstyle .EQ. 0) THEN
            CALL lj_cut (iflag)
         ELSE IF (coulstyle .EQ. 1) THEN
            CALL lj_cut_coul_cut (iflag)
         ELSE IF (coulstyle .EQ. 2) THEN
            CALL lj_cut_coul_smooth (iflag)
         ELSE IF (coulstyle .GE. 3) THEN
            CALL lj_cut_coul_long (iflag)
         ENDIF
      ELSE IF (nonstyle .EQ. 2) THEN
         IF (coulstyle .EQ. 0) THEN
            CALL lj_smooth (iflag)
         ELSE IF (coulstyle .EQ. 1) THEN
            CALL lj_smooth_coul_cut (iflag)
         ELSE IF (coulstyle .EQ. 2) THEN
            CALL lj_smooth_coul_smooth (iflag)
         ELSE IF (coulstyle .GE. 3) THEN
            CALL lj_smooth_coul_long (iflag)
         ENDIF
      ELSE IF (nonstyle .EQ. 3) THEN
         CALL lj_shift (iflag)
      ELSE IF (nonstyle .EQ. 4) THEN
         CALL soft (iflag)
      ELSE IF (nonstyle .EQ. 5) THEN
         IF (coulstyle .EQ. 0) THEN
            CALL ljclass2_cut (iflag)
         ELSE IF (coulstyle .EQ. 1) THEN
            CALL ljclass2_cut_coul_cut (iflag)
         ELSE IF (coulstyle .GE. 3) THEN
            CALL ljclass2_cut_coul_long (iflag)
         ENDIF
      ENDIF

C     ---- THE barrier CALL GIVES MORE ACCURATE TIMING INFORMATION.
c      CALL mpi_barrier(mpi_comm_world,ierror)

      time6 = mpi_wtime()

C     ---- COULOMB FORCES.
      IF (coulstyle .EQ. 3) THEN
         CALL ewald (iflag)
      ELSE IF (coulstyle .EQ. 4) THEN
         CALL pppm (iflag)
      ENDIF

C     ---- THE barrier CALL GIVES MORE ACCURATE TIMING INFORMATION.
      CALL mpi_barrier(mpi_comm_world,ierror)

      time7 = mpi_wtime()

C     ---- UPDATE GLOBAL TIMERS.
      time_bond     = time_bond     + (time2 - time1)
      time_angle    = time_angle    + (time3 - time2)
      time_dihedral = time_dihedral + (time4 - time3)
      time_improper = time_improper + (time5 - time4)
      time_nonbond  = time_nonbond  + (time6 - time5)
      time_long     = time_long     + (time7 - time6)

C     ---- SWAP FORCE INFORMATION BETWEEN PROCESSORS AND GET THE
C     ---- TOTAL FORCE ON EACH LOCAL ATOM.
      IF (newton .GE. 1) THEN
         time1 = mpi_wtime()
         CALL reverse_comm ()
         time2 = mpi_wtime()
         time_fcomm = time_fcomm + (time2 - time1)
      ENDIF

C     ---- COMPUTE NON-BONDED ENERGIES (SEE "energy" IN "thermo.f").
C     ---- RESULTS APPEAR IN THE GLOBAL VARIABLES "e_vdwl" AND "e_coul".
      time7 = mpi_wtime()
      CALL energy ()


C----------------------------------------------------------------------
C   EXTRACT INFORMATION FROM GLOBAL VARIABLES TO GET f AND g.
C----------------------------------------------------------------------

C     ---- SWAP AND ACCUMULATE EACH ENERGY COMPONENT (THESE CALLS
C     ---- ARE COPIED FROM SUBROUTINE "thermo").  AT THE END, EACH
C     ---- PROCESSOR KNOWS THE TOTAL ENERGIES.

      tmp = e_vdwl
      call mpi_allreduce(tmp,e_vdwl,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_coul
      call mpi_allreduce(tmp,e_coul,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      tmp = e_bond
      call mpi_allreduce(tmp,e_bond,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_angle
      call mpi_allreduce(tmp,e_angle,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_dihedral
      call mpi_allreduce(tmp,e_dihedral,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_improper
      call mpi_allreduce(tmp,e_improper,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      if (anglestyle.eq.2) then
        tmp = e_bondbond
        call mpi_allreduce(tmp,e_bondbond,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_bondangle
        call mpi_allreduce(tmp,e_bondangle,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif
      if (dihedstyle.eq.2) then
        tmp = e_midbondtorsion
        call mpi_allreduce(tmp,e_midbondtorsion,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_endbondtorsion
        call mpi_allreduce(tmp,e_endbondtorsion,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_angletorsion
        call mpi_allreduce(tmp,e_angletorsion,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_angleangletorsion
        call mpi_allreduce(tmp,e_angleangletorsion,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_bondbond13
        call mpi_allreduce(tmp,e_bondbond13,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif
      if (improstyle.eq.3) then
        tmp = e_angleangle
        call mpi_allreduce(tmp,e_angleangle,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif

C     ---- LUMP IN CLASS 2 TERMS SOMEHOW (THEY'RE ALREADY MERGED).
      IF (thermostyle .LE. 1) THEN
         IF (anglestyle .EQ. 2)
     +      e_angle = e_angle + e_bondbond + e_bondangle
         IF (dihedstyle .EQ. 2)
     +      e_dihedral = e_dihedral + e_midbondtorsion +
     +                   e_endbondtorsion + e_angletorsion +
     +                   e_angleangletorsion
         IF (improstyle .EQ. 3)
     +      e_improper = e_improper + e_angleangle
      ENDIF

C     ---- ADD EVERYTHING TO GET THE POTENTIAL ENERGY.
      IF (thermostyle .LE. 1) THEN
         e_potential = e_long + e_vdwl + e_coul +
     +                 e_bond + e_angle + e_dihedral + e_improper
      ELSE
         e_potential = e_long + e_vdwl + e_coul +
     +                 e_bond + e_angle + e_dihedral + e_improper +
     +                 e_bondbond + e_bondangle + e_midbondtorsion +
     +                 e_endbondtorsion + e_angletorsion +
     +                 e_angleangletorsion + e_angleangle
      ENDIF

C     ---- NORMALIZE TO THE NUMBER OF ATOMS IF USER REQUESTED IT.
      IF (units .EQ. 1) THEN
         rnorm = 1.0D0 / (DBLE (natoms))
         e_long = e_long*rnorm
         e_vdwl = e_vdwl*rnorm
         e_coul = e_coul*rnorm
         e_bond = e_bond*rnorm
         e_angle = e_angle*rnorm
         e_dihedral = e_dihedral*rnorm
         e_improper = e_improper*rnorm
         IF (anglestyle .EQ. 2) THEN
            e_bondbond = e_bondbond*rnorm
            e_bondangle = e_bondangle*rnorm
         ENDIF
         IF (dihedstyle .EQ. 2) THEN
            e_midbondtorsion = e_midbondtorsion*rnorm
            e_endbondtorsion = e_endbondtorsion*rnorm
            e_angletorsion = e_angletorsion*rnorm
            e_angleangletorsion = e_angleangletorsion*rnorm
         ENDIF
         IF (improstyle .EQ. 3) e_angleangle = e_angleangle*rnorm
         e_potential = e_potential*rnorm
         e_total = e_total*rnorm
      ENDIF


C----------------------------------------------------------------------
C   COPY INFORMATION FROM GLOBAL VARIABLES.
C----------------------------------------------------------------------

      obj_f = e_potential

C     ---- COPY THE FORCE VECTOR.  BY CONVENTION, THE LAMMPS FORCE IS
C     ---- STEEPEST DESCENT, SO NEGATE IT INTO A GRADIENT.

      IF (g_vec_size .LT. (3 * nlocal)) THEN
         CALL error ('g_vec_size in min_support must be larger')
      ENDIF

      i = atompnt
      j = 0
      DO k = 1, nlocal
         g_vec(j+1) = (-1.0D0) * f(1,i)
         g_vec(j+2) = (-1.0D0) * f(2,i)
         g_vec(j+3) = (-1.0D0) * f(3,i)
         j = j + 3
         i = list(i)
      ENDDO
c        i = atompnt
c        do k = 1, nlocal
c         write (6,901) node, k, f(1,i), f(2,i), f(3,i)
c 901     format ('TBD node ', I1, ' has total force on ', I3,
c     +           ' = [', F10.5, 1X, F10.5, 1X, F10.5, ']')
c         i = list(i)
c        enddo

      time1 = mpi_wtime()
      opt_time1 = opt_time1 + (time1 - time_top)
      opt_time2 = opt_time2 + (time1 - time7)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Get_Smallest_Box_Dim
     +   (smallest_dim)
C
C       orig: 30 Jan 97  T. Plantenga
C
C   Input/output variable definitions:
C       pgrid(3)         I  (global variable)
C       xboundlo         I  (global variable)
C       xboundhi         I  (global variable)
C       yboundlo         I  (global variable)
C       yboundhi         I  (global variable)
C       zboundlo         I  (global variable)
C       zboundhi         I  (global variable)
C       smallest_dim      O MIN (x length, y length, z length)
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
C
      DOUBLE PRECISION  smallest_dim
C
C
C   This routine accesses global variables to compute the dimensions
C   of the box size used on each processor.
C**********************************************************************

C     ---- LOCAL VARIABLES
      DOUBLE PRECISION  x_dim, y_dim, z_dim

C----------------------------------------------------------------------

C     ---- PROCESSOR BOX SIZE IS DETERMINED BY SIMPLY DIVIDING THE
C     ---- TOTAL EXTENT OF THE MOLECULE IN EACH DIRECTION BY THE NUMBER
C     ---- OF PROCESSORS IN THE DIRECTION.  LAMMPS DETERMINES THE
C     ---- PROCESSOR LAYOUT ONCE AT STARTUP, STORING THE VALUES IN THE
C     ---- ARRAY pgrid.  THE EXTENT OF THE MOLECULE IS STORED IN
C     ---- xboundlo, xboundhi, yboundlo, ..., ETC.  THESE ARE COMPUTED
C     ---- IN SUBROUTINE "setup_bounds" (FROM FILE "setup.f"), WHICH IS
C     ---- CALLED BY SUBROUTINE "setup_box".

      x_dim = xboundhi - xboundlo
      x_dim = x_dim / (DBLE (pgrid(1)) )

      y_dim = yboundhi - yboundlo
      y_dim = y_dim / (DBLE (pgrid(2)) )

      z_dim = zboundhi - zboundlo
      z_dim = z_dim / (DBLE (pgrid(3)) )

      smallest_dim = MIN (x_dim, y_dim)
      smallest_dim = MIN (z_dim, smallest_dim)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Get_Atom_Positions
     +   (localN, x_vec)
C
C       orig: 25 Oct 96  T. Plantenga
C
C   Input/output variable definitions:
C       nlocal           I  number local atoms (global variable)
C       localN            O 3 * nlocal
C       x                I  atom positions (global variable)
C       x_vec             O atom positions
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
C
      INTEGER           localN
      DOUBLE PRECISION  x_vec(3*maxown)
C
C
C   This routine copies the atom positions vector from the global
C   LAMMPS array to the local x_vec variable.  It also updates n
C   from nlocal, which could change after Repartitn_and_Reneighbor.
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER  i, j, k

C----------------------------------------------------------------------

      i = atompnt
      j = 0
      DO k = 1, nlocal
         x_vec(j+1) = x(1,i)
         x_vec(j+2) = x(2,i)
         x_vec(j+3) = x(3,i)
         j = j + 3
         i = list(i)
      ENDDO

      localN = 3 * nlocal

      RETURN
      END


C**********************************************************************
      SUBROUTINE Put_Atom_Positions
     +   (n, x_vec, update_per_trues)
C
C       orig: 25 Oct 96  T. Plantenga
C
C   Input/output variable definitions:
C       x_vec            I  atom positions
C       update_per_trues I  boolean
C       perflagx         I  (global variable)
C       perflagy         I  (global variable)
C       perflagz         I  (global variable)
C       xboundlo         I  (global variable)
C       xboundhi         I  (global variable)
C       yboundlo         I  (global variable)
C       yboundhi         I  (global variable)
C       zboundlo         I  (global variable)
C       zboundhi         I  (global variable)
C       xprd             I  (global variable)
C       yprd             I  (global variable)
C       zprd             I  (global variable)
C       true             I  (global variable)
C       x                 O atom positions (global variable)
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
      INCLUDE "mpif.h"
C
      INTEGER           n
      DOUBLE PRECISION  x_vec(n)
      INTEGER           update_per_trues
C
C
C   This routine copies the atom positions x_vec into the global
C   LAMMPS array.  Position information is then swapped to other
C   processors that need to know non-resident atom positions.
C
C   Before storing positions, they are adjusted if necessary to stay
C   within a periodic problem boundary.  This adjustment is not visible
C   to the calling routine, so it is important to recover atom positions
C   using "Get_Atom_Positions".  If "update_per_trues" is nonzero,
C   the "true" flags are updated to track absolute atom positions.
C   
C   External subroutines called:
C        communicate         "communicate.f"
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER           i, j, k
      DOUBLE PRECISION  time1, time2

C----------------------------------------------------------------------

      i = atompnt
      j = 0
      DO k = 1, nlocal
         x(1,i) = x_vec(j+1)
         x(2,i) = x_vec(j+2)
         x(3,i) = x_vec(j+3)

C        ---- STAY WITHIN PERIODIC BOUNDARY CONDITIONS
C        ---- (SEE SUBROUTINE nve_x IN "integrate.f" FOR AN EXAMPLE).
C        ---- "true" FLAGS COUNT THE CUMULATIVE NONPERIODIC DISTANCE
C        ---- MOVED BY EACH ATOM; USED ONLY AT THE END OF A LAMMPS RUN.
         IF (perflagx .EQ. 0) THEN
            IF (x(1,i) .GE. xboundhi) THEN
               x(1,i) = x(1,i) - xprd
               IF (update_per_trues .NE. 0)
     $              true(i) = true(i) + 1
            ENDIF
            IF (x(1,i) .LT. xboundlo) THEN
               x(1,i) = x(1,i) + xprd
               IF (update_per_trues .NE. 0)
     $              true(i) = true(i) - 1
            ENDIF
         ENDIF
         IF (perflagy .EQ. 0) THEN
            IF (x(2,i) .GE. yboundhi) THEN
               x(2,i) = x(2,i) - yprd
               IF (update_per_trues .NE. 0)
     $              true(i) = true(i) + 1000
            ENDIF
            IF (x(2,i) .LT. yboundlo) THEN
               x(2,i) = x(2,i) + yprd
               IF (update_per_trues .NE. 0)
     $              true(i) = true(i) - 1000
            ENDIF
         ENDIF
         IF (perflagz .EQ. 0) THEN
            IF (x(3,i) .GE. zboundhi) THEN
               x(3,i) = x(3,i) - zprd
               IF (update_per_trues .NE. 0)
     $              true(i) = true(i) + 1000000
            ENDIF
            IF (x(3,i) .LT. zboundlo) THEN
               x(3,i) = x(3,i) + zprd
               IF (update_per_trues .NE. 0)
     $              true(i) = true(i) - 1000000
            ENDIF
         ENDIF

         j = j + 3
         i = list(i)
      ENDDO

      time1 = mpi_wtime()
      CALL communicate ()
      time2 = mpi_wtime()

C     ---- UPDATE GLOBAL TIMER.
      time_comm = time_comm + (time2 - time1)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Compute_Inf_Norm_PAR
     +   (n, x_vec, inf_norm)
C
C       orig: 23 Oct 96  T. Plantenga
C       mod:  21 May 97  TDP - changed from fn to sub for C
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
      INCLUDE "mpif.h"
C
      INTEGER           n
      DOUBLE PRECISION  x_vec(n)
      DOUBLE PRECISION  inf_norm
C
C
C   This routine returns the infinity norm of x_vec, which is a vector
C   strung out across all processors.  Each processor computes the
C   infinity norm of its locally stored portion (in parallel), then
C   merges this with all the other processor norms (requiring
C   interprocessor communication).  At the end, all processors know
C   ||x_vec||_inf.
C
C**********************************************************************

C     ---- LOCAL VARIABLES
C     ---- (VARIABLES COMMUNICATED VIA MPI MUST USE REAL*8).
      INTEGER           i
      DOUBLE PRECISION  dd1, local_max
      DOUBLE PRECISION  time1, time2

C----------------------------------------------------------------------

      time1 = mpi_wtime()

      local_max = 0.0D0
      DO i = 1, n
         dd1 = DABS (x_vec(i))
         IF (dd1 .GT. local_max)  local_max = dd1
      ENDDO

      call mpi_allreduce(local_max,inf_norm,1,mpi_double_precision,
     $     mpi_max,mpi_comm_world,ierror)

      time2 = mpi_wtime()
      opt_time3 = opt_time3 + (time2 - time1)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Compute_L2_Norm_PAR
     +   (n, x_vec, L2_norm)
C
C       orig: 25 Oct 96  T. Plantenga
C       mod:  21 May 97  TDP - changed from fn to sub for C
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
      INCLUDE "mpif.h"
C
      INTEGER           n
      DOUBLE PRECISION  x_vec(n)
      DOUBLE PRECISION  L2_norm
C
C
C   This routine returns the L2 norm of a distributed vector x_vec.
C   Each processor computes the L2 norm of its locally stored portion
C   using BLAS (in parallel), then merges this with all the other
C   processor norms (requiring interprocessor communication).  At the
C   end, all processors know ||x_vec||_2.
C
C   External subroutines called:
C        dnrm2               "blas/dnrm2.f"
C**********************************************************************

      DOUBLE PRECISION  dnrm2
      EXTERNAL          dnrm2

C     ---- LOCAL VARIABLES
C     ---- (VARIABLES COMMUNICATED VIA MPI MUST USE REAL*8).
      DOUBLE PRECISION  x_2norm
      DOUBLE PRECISION  time1, time2
      REAL*8            global_real

C----------------------------------------------------------------------

      time1 = mpi_wtime()

      x_2norm = dnrm2 (n, x_vec, 1)

      tmp = x_2norm * x_2norm
      call mpi_allreduce(tmp,global_real,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      L2_norm = SQRT (global_real)

      time2 = mpi_wtime()
      opt_time3 = opt_time3 + (time2 - time1)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Compute_ddot_PAR
     +   (n, x1_vec, x2_vec, dot_product)
C
C       orig: 31 Oct 96  T. Plantenga
C       mod:  21 May 97  TDP - changed from fn to sub for C
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
      INCLUDE "mpif.h"
C
      INTEGER           n
      DOUBLE PRECISION  x1_vec(n), x2_vec(n)
      DOUBLE PRECISION  dot_product
C
C
C   This routine returns the dot product of two distributed vectors.
C   Each processor computes the dot product of its locally stored
C   portions using BLAS (in parallel), then adds this with all the
C   other processor norms (requiring interprocessor communication).
C   At the end, all processors know x1_vec'x1_vec
C
C   External subroutines called:
C        ddot                "blas/ddot.f"
C**********************************************************************

      DOUBLE PRECISION  ddot
      EXTERNAL          ddot

C     ---- LOCAL VARIABLES
C     ---- (VARIABLES COMMUNICATED VIA MPI MUST USE REAL*8).
      DOUBLE PRECISION  x1_dot_x2
      DOUBLE PRECISION  time1, time2
      REAL*8            global_real

C----------------------------------------------------------------------

      time1 = mpi_wtime()

      x1_dot_x2 = ddot (n, x1_vec, 1, x2_vec, 1)

      call mpi_allreduce(x1_dot_x2,dot_product,1,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

      time2 = mpi_wtime()
      opt_time3 = opt_time3 + (time2 - time1)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Normalize_ddot(dot_product,n)
      include "lammps.h"
C
C     orig: 15 Jan 01  S. Plimpton
C
      DOUBLE PRECISION dot_product
      INTEGER n
C
C   This routine normalizes a previously computed dot-product
C   by the length of the vector n
C   This is necessary if units are LJ for energy estimation performed
C   by the minimizer via Approximate_Bv

C**********************************************************************

      if (units .eq. 1) dot_product = dot_product/n

      RETURN
      END


C**********************************************************************
      SUBROUTINE Find_Biggest_Jump
     +   (boxes_jumped, grid_num_procs, periodic_flag)
C
C       orig: 31 Jan 97  T. Plantenga
C
C   Input/output variable definitions:
C       x                I  proposed atom positions (global variable)
C       pgrid            I  num procs in each direction (global variable)
C       xboundlo         I  (global variable)
C       xboundhi         I  (global variable)
C       yboundlo         I  (global variable)
C       yboundhi         I  (global variable)
C       zboundlo         I  (global variable)
C       zboundhi         I  (global variable)
C       perflagx         I  (global variable)
C       perflagy         I  (global variable)
C       perflagz         I  (global variable)
C       border           I  borders of local processor (global variable)
C       boxes_jumped      O biggest jump in terms of boxes
C       grid_num_procs    O number processors in direction of biggest jump
C       periodic_flag     O zero if all boundaries are periodic
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
      INCLUDE "mpif.h"
C
      DOUBLE PRECISION  boxes_jumped
      INTEGER           grid_num_procs, periodic_flag
C
C   This routine has each processor check the current position of their
C   local atoms to see if any have left the local "box" domain.  The
C   number of boxes away (in any one direction) that an atom has moved
C   is computed.  The single largest jump made by any atom in the whole
C   set of processors is returned, along with the number of processors
C   in that direction.
C
C   Parallel communication is used at the end.
C   
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER           i, j, k
      REAL*8            dd1, dd2, dd3
      REAL*8            box_dim(3), prob_edges(3,2)
      REAL*8            border_lo, border_hi
      INTEGER           per_flags(3)
      DOUBLE PRECISION  biggest_jump(3)

C----------------------------------------------------------------------

C     ---- CALCULATE THE SPATIAL DIMENSIONS OF A GENERIC PROCESSOR BOX.
      box_dim(1) = (xboundhi - xboundlo) / ( FLOAT (pgrid(1)) )
      box_dim(2) = (yboundhi - yboundlo) / ( FLOAT (pgrid(2)) )
      box_dim(3) = (zboundhi - zboundlo) / ( FLOAT (pgrid(3)) )

C     ---- PUT THE OVERALL EDGES OF THE (PERIODIC) PROBLEM IN AN ARRAY
C     ---- FOR EASIER ACCESS.
      prob_edges(1,1) = xboundlo
      prob_edges(1,2) = xboundhi
      prob_edges(2,1) = yboundlo
      prob_edges(2,2) = yboundhi
      prob_edges(3,1) = zboundlo
      prob_edges(3,2) = zboundhi

C     ---- LOOK AT PERIODICITY OF THE BOUNDARIES.
      per_flags(1) = perflagx
      per_flags(2) = perflagy
      per_flags(3) = perflagz
      IF ( (per_flags(1) .EQ. 0) .AND.
     +     (per_flags(2) .EQ. 0) .AND.
     +     (per_flags(3) .EQ. 0)       ) THEN
         periodic_flag = 0
      ELSE
         periodic_flag = 1
      ENDIF

C     ---- THIS LOOP IS ADAPTED FROM SUBROUTINE "exchange".
      DO k = 1, 3
         dd1 = 0.0
         IF (pgrid(k) .GT. 1) THEN

            border_lo = border(1,k)
            border_hi = border(2,k)

            i = atompnt
            DO j = 1, nlocal

C              ---- CALCULATE HOW MANY BOXES AWAY THE ATOM HAS MOVED.
C              ---- COMPENSATE FOR PERIODIC BOUNDARIES IF NECESSARY.

               IF (x(k,i) .LT. border_lo) THEN
                  dd2 = (border_lo - x(k,i)) / box_dim(k)
                  IF (per_flags(k) .EQ. 0) THEN
                     dd3 = (x(k,i) - prob_edges(k,1)) / box_dim(k)
                     dd2 = MIN (dd2, dd3)
                  ENDIF
                  dd1 = MAX (dd1, dd2)
               ENDIF

               IF (x(k,i) .GE. border_hi) THEN
                  dd2 = (x(k,i) - border_hi) / box_dim(k)
                  IF (per_flags(k) .EQ. 0) THEN
                     dd3 = (prob_edges(k,2) - x(k,i)) / box_dim(k)
                     dd2 = MIN (dd2, dd3)
                  ENDIF
                  dd1 = MAX (dd1, dd2)
               ENDIF

               i = list(i)
            ENDDO

C           ---- MERGE EACH PROCESSOR'S BIGGEST JUMP TO FIND OUT THE
C           ---- BIGGEST JUMP MADE BY ANY ATOM IN THE PROBLEM.

            tmp = dd1
            call mpi_allreduce(tmp,dd1,1,
     $           mpi_double_precision,mpi_max,mpi_comm_world,ierror)

         ENDIF
         biggest_jump(k) = dd1
      ENDDO

C     ---- FIND THE BIGGEST JUMP OVER THE X, Y, AND Z DIMENSIONS.
      boxes_jumped   = 0.0D0
      grid_num_procs = 1
      DO k = 1, 3
        IF (biggest_jump(k) .GT. boxes_jumped) THEN
           boxes_jumped   = biggest_jump(k)
           grid_num_procs = pgrid(k)
        ENDIF
      ENDDO

      RETURN
      END


C**********************************************************************
      SUBROUTINE Repartitn_and_Reneighbor
C
C       orig: 20 Jan 97  T. Plantenga
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
      INCLUDE "mpif.h"
C
C
C   This routine consists of calls extracted from "integrate" that
C   redefine each processor's domain, move atoms into the correct domains,
C   and then build new neighbor lists for nonbonded interactions.
C
C   Operations are performed in parallel, with substantial interprocessor
C   communication.
C   
C   External subroutines called:
C        setup_box           "setup.f"
C        exchange            "communicate.f"
C        borders             "communicate.f"
C        neighbor            "neighbor.f"
C        neighbor_bond       "neighbor.f"
C        neighbor_angle      "neighbor.f"
C        neighbor_dihedral   "neighbor.f"
C        neighbor_improper   "neighbor.f"
C**********************************************************************

C     ---- LOCAL VARIABLES.
      DOUBLE PRECISION  time1, time2

C----------------------------------------------------------------------

C     ---- UNLESS THE PROBLEM IS COMPLETELY PERIODIC, RECOMPUTE THE
C     ---- GEOMETRIC REGION EACH PROCESSOR IS RESPONSIBLE FOR.
      time1 = mpi_wtime()
      IF ( ((perflagx + perflagy + perflagz) .GT. 0) .OR.
     +     (ensemble .EQ. 3) ) THEN
         CALL setup_box ()
      ENDIF

C     ---- ERASE THE "OTHER" ATOM LIST FOR NONBONDED INTERACTIONS AND
C     ---- SHIFT ANY ATOMS NO LONGER IN THEIR LOCAL PROCESSOR.  (AN ATOM
C     ---- CAN BE SHIFTED TO AT MOST THE NEAREST OF THE 26 BOX NEIGHBORS;
C     ---- IF IT'S NEW PROCESSOR HOME IS FARTHER AWAY, IT WILL BE LOST.)
C     ---- AS A SIDE EFFECT, exchange CALLS setup_comm TO GENERATE A
C     ---- SWAPPING SCHEME FOR UPDATING THE NEW "OTHER" ATOM LISTS.
      CALL exchange ()
      time2 = mpi_wtime()
      time_exch = time_exch + (time2 - time1)

C     ---- SWAP ATOM INFORMATION SO THAT EACH LOCAL PROCESSOR HAS ALL
C     ---- THE "OTHER" ATOMS WITHIN THE MAXIMUM CUTOFF DISTANCE OF THE
C     ---- NONBONDED FORCES.
      time1 = mpi_wtime()
      CALL borders ()
      time2 = mpi_wtime()
      time_comm = time_comm + (time2 - time1)

C     ---- MAKE NEW NEIGHBOR LISTS FOR EACH LOCAL ATOM, COMPRISING ALL
C     ---- ATOMS WITHIN THE CUTOFF RADIUS.
      time1 = mpi_wtime()
      CALL neighbor ()
      time2 = mpi_wtime()
      time_neigh1 = time_neigh1 + (time2 - time1)

C     ---- VERIFY THAT THE NEW NEIGHBOR LIST IS LEGAL.
      time1 = mpi_wtime()
      IF (nbonds     .GT. 0)  CALL neighbor_bond ()
      IF (nangles    .GT. 0)  CALL neighbor_angle ()
      IF (ndihedrals .GT. 0)  CALL neighbor_dihedral ()
      IF (nimpropers .GT. 0)  CALL neighbor_improper ()
      time2 = mpi_wtime()
      time_neigh2 = time_neigh2 + (time2 - time1)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Write_Energy_Vals
     +   (outUnit, procNum)
C
C       orig: 14 Nov 96  T. Plantenga
C
C     ---- MUST DO THIS BEFORE DECLARING ANY VARIABLES BECAUSE OF
C     ---- FORTRAN implicit STATEMENT BURIED INSIDE.
      INCLUDE "lammps.h"
C
      INTEGER  outUnit, procNum
C
C
C   This routine obtains and writes out the many components of
C   potential energy.  It is assumed energy has been previously
C   calculated by a call to "Compute_f" or "Compute_f_g".
C**********************************************************************

C     ---- LOCAL VARIABLES.
      DOUBLE PRECISION  e_potential

C----------------------------------------------------------------------

      IF (procNum .NE. 0)  RETURN

C     ---- CALCULATION OF e_potential MATCHES CODE IN SUBROUTINES
C     ---- "Compute_f_g" AND "thermo".
      IF (thermostyle .LE. 1) THEN
         e_potential = e_long + e_vdwl + e_coul +
     +                 e_bond + e_angle + e_dihedral + e_improper
      ELSE
         e_potential = e_long + e_vdwl + e_coul +
     +                 e_bond + e_angle + e_dihedral + e_improper +
     +                 e_bondbond + e_bondangle + e_midbondtorsion +
     +                 e_endbondtorsion + e_angletorsion +
     +                 e_angleangletorsion + e_angleangle
      ENDIF

      WRITE (6,600)
      WRITE (6,601) e_potential
      WRITE (6,602) e_bond, e_coul
      WRITE (6,603) e_angle, e_vdwl
      WRITE (6,604) e_dihedral, e_long
      WRITE (6,605) e_improper
      WRITE (6,606)
      IF (outUnit .GT. 0) THEN
         WRITE (outUnit,600)
         WRITE (outUnit,601) e_potential
         WRITE (outUnit,602) e_bond, e_coul
         WRITE (outUnit,603) e_angle, e_vdwl
         WRITE (outUnit,604) e_dihedral, e_long
         WRITE (outUnit,605) e_improper
         WRITE (outUnit,606)
      ENDIF

      RETURN

C----------------------------------------------------------------------
C   THE REMAINING LINES OF CODE ARE I/O STATEMENTS REFERENCED EARLIER.
C----------------------------------------------------------------------

 600  FORMAT (/'======================================================',
     +        '==============')
 601  FORMAT ('  Total PE  = ', F14.4)
 602  FORMAT ('  E (bond)  = ', F14.4, 10X, 'E (coulomb) = ', F14.4)
 603  FORMAT ('  E (angle) = ', F14.4, 10X, 'E (VDW)     = ', F14.4)
 604  FORMAT ('  E (dih)   = ', F14.4, 10X, 'E (long)    = ', F14.4)
 605  FORMAT ('  E (imprp) = ', F14.4)
 606  FORMAT ('======================================================',
     +        '=============='/)

      END


C**********************************************************************
      SUBROUTINE Write_Opt_Timer
     +   (localTimerValue, totalNumProcs, numFnEvals,
     +    outUnit, procNum)
      INCLUDE "mpif.h"
C
C       orig:  4 Sep 97  T. Plantenga
C
      DOUBLE PRECISION  localTimerValue
      INTEGER           totalNumProcs, numFnEvals
      INTEGER           outUnit, procNum
C
C
C   This routine uses interprocessor communication to calculate the
C   average and maximum timer values over all processors.  Results
C   are written out by processor zero.
C   
C**********************************************************************

C     ---- LOCAL VARIABLES.
      DOUBLE PRECISION  avgTime, maxTime

C----------------------------------------------------------------------

      call mpi_allreduce(localTimerValue,avgTime,1,
     $     mpi_double_precision,mpi_max,mpi_comm_world,ierror)
      avgTime = avgTime / totalNumProcs

      call mpi_allreduce(localTimerValue,maxTime,1,
     $     mpi_double_precision,mpi_max,mpi_comm_world,ierror)

      IF (procNum .EQ. 0) THEN
         WRITE (6,600) avgTime, (avgTime / (DBLE (numFnEvals)))
         WRITE (6,601) maxTime
         IF (outUnit .GT. 0) THEN
            WRITE (outUnit,600) avgTime, (avgTime / (DBLE (numFnEvals)))
            WRITE (outUnit,601) maxTime
         ENDIF
      ENDIF

      RETURN

C----------------------------------------------------------------------
C   THE REMAINING LINES OF CODE ARE I/O STATEMENTS REFERENCED EARLIER.
C----------------------------------------------------------------------

 600  FORMAT ('  (raw) avg = ', F10.2, '  (per f&g eval) avg = ', F8.3)
 601  FORMAT ('        max = ', F10.2)

      END


C***********************End of source file*****************************
