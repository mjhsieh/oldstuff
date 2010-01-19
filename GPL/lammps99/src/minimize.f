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
C   FILE minimize.f
C
C   Externally callable routines:
C      Min_HFTN
C
C   Local routines (in order)
C      Approximate_Bv
C      Compute_HFTN_Step
C      Extend_to_TR
C      Update_TR
C      Display_Opt_Step
C
C   Code development:
C      16 Nov 94 - Originated by T. Plantenga, Sandia National Labs.
C      21 Oct 96 - (TDP) adapted for LAMMPS version 3.0.
C      27 Jul 97 - (TDP) adapted for LAMMPS version 4.0.
C**********************************************************************


C**********************************************************************
      SUBROUTINE Min_HFTN
     +   (procNum, outUnit, totalN, localN,
     +    stop_tol, max_num_f_g, max_iters, repartitn_freq,
     +    repartitn_delay, num_f_g, checkpoint_flag)
C
C       orig:  T. Plantenga  21 Oct 96
C
C   Input/output variable definitions:
C       procNum          I  processor number (starting from zero)
C       outUnit          I  if nonzero, output file unit
C       localN           I  number of unknowns on local processor
C       totalN           I  total number of unknowns (3 * num atoms)
C       stop_tol         I  applied to ||grad||_inf
C       max_num_f_g      I  max number calls to Compute_f_g
C       max_iters        I  max number outer loops of algorithm
C       repartitn_freq   I  number of accepted steps between repartitions
C       repartitn_delay  I  delay repartition for this many accepted steps
C       num_f_g          O  number of energy/force evaluations made
C       checkpoint_flag  I  how often to write restart file (0=never)
C
      INTEGER           procNum, outUnit
      INTEGER           totalN, localN
      DOUBLE PRECISION  stop_tol
      INTEGER           max_num_f_g, max_iters
      INTEGER           repartitn_freq,repartitn_delay
      INTEGER           num_f_g
      INTEGER           checkpoint_flag
C
C
C   This routine is called from "lammps.f" to find a molecular
C   configuration that gives a local minimum of the potential energy.
C   The outline of the optmization algorithm is as follows:
C
C        (the start point is already selected)
C        Initialize parameters
C        Compute f and g at the start point
C        LOOP
C          IF converged OR max_num_f_g exceeded THEN RETURN
C          Compute a step d (thru an inner CG iteration)
C          Compute the predicted reduction of the quadratic model
C          Compute the actual potential energy reduction
C          IF ared/pred > eta
C            THEN x = x + d
C                 TR_size >= TR_size
C                 Possibly repartition and reneighbor; if so, recompute f
C                 Compute the gradient g at the new point
C            ELSE TR_size < ||d||
C        END LOOP
C
C   Convergence is defined as either
C        ||g||_inf <= stop_tol  or  |f(x_k) - f(x_k + d)| < DBL_epsilon.
C   The last test is a safeguard since DBL_epsilon is so small.
C   In the first test the gradient should be divided by max{1, |f|} to
C   provide scale invariance; however, this makes it hard to compare with
C   constrained optimization algorithms.
C
C   This algorithm runs in parallel on each processor.  All processors
C   execute the same code (except for file I/O by processor zero).
C
C   External subroutines called:
C        Compute_f_g               "min_support.f"
C        Compute_Inf_Norm_PAR      "min_support.f"
C        Compute_ddot_PAR          "min_support.f"
C        Get_Smallest_Box_Dim      "min_support.f"
C        Get_Atom_Positions        "min_support.f"
C        Put_Atom_Positions        "min_support.f"
C        Find_Biggest_Jump         "min_support.f"
C        Repartitn_and_Reneighbor  "min_support.f"
C        write_restart             "output.f"
C        dcopy                     "blas/dcopy.f"
C**********************************************************************

C     ---- MAXIMUM PROBLEM SIZES FOR THE LOCAL PROCESSOR ARE DEFINED HERE.
      INCLUDE 'param.h'

      DOUBLE PRECISION  DBL_epsilon
      PARAMETER (DBL_epsilon = 2.2204460492503131D-16)

C     ---- LOCAL VARIABLES.
C        f            - objective value (potential energy)
C        x_vec        - <x,y,z> atom positions for the local atoms
C        g_vec        - gradient vector (negative of force vector)
C        d_vec        - trial step
C        TR_size      - trust region radius
C        TR_eta       - trust region step acceptance threshold
C        ared         - actual reduction in f made by the trial step
C        pred         - predicted reduction in f made by the trial step
C        findiff_flag - 'f' for forward, 'c' for central differences

      INTEGER           num_evals
      INTEGER           iteration
      INTEGER           steps_since_repartitn, num_repartitns
      INTEGER           fevals_since_chkpnt
      INTEGER           HFTN_return_code, CG_iters
      CHARACTER         findiff_flag
      DOUBLE PRECISION  f, f_new
      DOUBLE PRECISION  x_vec(3*maxown)
      DOUBLE PRECISION  g_vec(3*maxown), g_Inorm
      DOUBLE PRECISION  d_vec(3*maxown), d_2norm
      DOUBLE PRECISION  TR_size, TR_eta
      DOUBLE PRECISION  ared, pred

      INTEGER           num_procs_in_jump, periodic_bnds
      DOUBLE PRECISION  boxes_jumped, jump_threshold

      DOUBLE PRECISION  dd1, dd2
      DOUBLE PRECISION  tmp_vec1(3*maxown), tmp_vec2(3*maxown)

C----------------------------------------------------------------------

C     ---- DISPLAY INITIAL HEADER (ONLY PROCESSOR ZERO).
      IF (procNum .EQ. 0) THEN
         WRITE (6,600) totalN
         WRITE (6,601)
         WRITE (6,602)
         IF (outUnit .GT. 0) THEN
            WRITE (outUnit,600) totalN
            WRITE (outUnit,601)
            WRITE (outUnit,602)
         ENDIF
      ENDIF

C     ---- INITIALIZE ALGORITHM PARAMETERS.
      TR_eta                = 0.3D0
      iteration             = 0
      steps_since_repartitn = 0
      num_repartitns        = 0
      fevals_since_chkpnt   = 1

C     ---- INITIALIZE THE TRUST REGION TO ALLOW THE AVERAGE ATOM TO
C     ---- MOVE 0.1 ANGSTROMS.  RESTRICT IF NECESSARY TO KEEP THE AVERAGE
C     ---- FROM MOVING MORE THAN ONE BOX PER ITERATION.
      CALL Get_Smallest_Box_Dim (dd1)
      TR_size = MIN (0.1D0, dd1) * SQRT (DBLE (totalN))


C======================================================================
C   ---- Evaluating the function and gradient in LAMMPS ----
C
C   LAMMPS evaluates at the atom positions stored in global variable "x".
C   The "x" positions should be extracted and overwritten using the
C   routines "Get_Atom_Positions" and "Put_Atom_Positions".
C   LAMMPS calculates the potential energy and forces with some
C   interprocessor communication (see "Compute_f_g").  It returns the
C   gradient vector; i.e., direction of steepest ascent.  The flag
C   "no_nonbonded_E" (first argument in the call) must be set to zero
C   if the returned energy is to be meaningful; setting it to one
C   save execution time.
C
C   Periodic boundary conditions cause atoms to wrap around if they
C   leave the box, but LAMMPS also tracks absolute position.  To keep
C   these consistent, "Put_Atom_Positions" and "Get_Atom_Positions"
C   must be used to move to a trial position, and to move back.
C======================================================================

C     ---- COMPUTE THE INITIAL f AND g.  LAMMPS HAS ALREADY CALLED
C     ---- "start", WHICH SETS UP NEIGHBOR LISTS AND CALCULATES THE
C     ---- FORCES.  THUS, THIS CALL DUPLICATES THE FORCE CALCULATION.
      CALL Compute_f_g (f, 3*maxown, g_vec)
      num_f_g = 1
      CALL Compute_Inf_Norm_PAR (localN, g_vec, g_Inorm)

C     ---- DISPLAY START POINT DATA (ONLY PROCESSOR ZERO).
      IF (procNum .EQ. 0) THEN
         WRITE (6,610) iteration, f, g_Inorm, TR_size, num_f_g
         IF (outUnit .GT. 0) THEN
            WRITE (outUnit,610) iteration, f, g_Inorm, TR_size,
     +                          num_f_g
         ENDIF
      ENDIF

C     ---- INITIALIZE FINITE DIFFERENCE SCHEME.
      findiff_flag = 'f'
      IF (procNum .EQ. 0) THEN
         WRITE (6,*) '  <forward diffs>'
         IF (outUnit .GT. 0)  WRITE (outUnit,*) '  <forward diffs>'
      ENDIF


C======================================================================
C   ---- Trust region minimization ----
C
C   The minimization algorithm computes iterative steps that move all
C   atoms toward a lower potential energy state.  Each iteration, a
C   trial step is computed that reduces a quadratic model of the
C   objective energy.  The model is composed of the analytic gradient
C   plus finite difference Hessian approximations in a given direction.
C
C   The trial step is then evaluated using the full potential energy.
C   If (actual reduction) / (predicted reduction) exceeds TR_eta, then
C   the step is accepted; otherwise, it is rejected and a new step
C   computed using a smaller trust region.  Neighbor lists may be
C   updated after accepting a step, which could in principle cause
C   the energy to increase.
C======================================================================

C     ---- GET THE INITIAL ATOM POSITIONS INTO x_vec.
      CALL Get_Atom_Positions (localN, x_vec)

C----------------------------------------------------------------------
C   BEGIN MAIN OPTIMIZATION LOOP.
C----------------------------------------------------------------------

 100  CONTINUE

C     ---- DECIDE WHETHER TO STOP.

      IF (g_Inorm .LE. stop_tol) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,*) '+++ ||g||_inf less than tolerance'
            IF (outUnit .GT. 0)
     +         WRITE (outUnit,*) '+++ ||g||_inf less than tolerance'
         ENDIF
         GOTO 200
      ENDIF

      IF (num_f_g .GE. max_num_f_g) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,*) '+++ Max number of force evals reached'
            IF (outUnit .GT. 0)
     +         WRITE (outUnit,*) '+++ Max number of force evals reached'
         ENDIF
         GOTO 200
      ENDIF

      iteration = iteration + 1
      IF (iteration .GT. max_iters) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,*) '+++ Max number iterations reached'
            IF (outUnit .GT. 0)
     +         WRITE (outUnit,*) '+++ Max number iterations reached'
         ENDIF
         GOTO 200
      ENDIF


C     ---- CHOOSE THE FINITE DIFFERENCING SCHEME.
      IF (g_Inorm .LE. 1.0D-3) THEN
         IF (findiff_flag .EQ. 'f') THEN
            IF (procNum .EQ. 0) THEN
               WRITE (6,*) '  <central diffs>'
               IF (outUnit .GT. 0) WRITE (outUnit,*) '  <central diffs>'
            ENDIF
         ENDIF
         findiff_flag = 'c'
      ELSE
         IF (findiff_flag .EQ. 'c') THEN
            IF (procNum .EQ. 0) THEN
               WRITE (6,*) '  <forward diffs>'
               IF (outUnit .GT. 0) WRITE (outUnit,*) '  <forward diffs>'
            ENDIF
         ENDIF
         findiff_flag = 'f'
      ENDIF

C     ---- COMPUTE THE TRIAL STEP d_vec WITH LENGTH d_2norm.
      CALL Compute_HFTN_Step
     +   (localN, totalN, x_vec, g_vec, TR_size,
     +    iteration, max_num_f_g, num_f_g, findiff_flag,
     +    d_vec, d_2norm, CG_iters, num_evals, HFTN_return_code)
      num_f_g             = num_f_g             + num_evals
      fevals_since_chkpnt = fevals_since_chkpnt + num_evals

C     ---- COMPUTE THE PREDICTED REDUCTION MADE BY THE STEP IN tmp_vec1.
      CALL Approximate_Bv
     +   (findiff_flag, localN, x_vec, g_vec, d_vec, d_2norm,
     +    tmp_vec2, num_evals, tmp_vec1)
      num_f_g             = num_f_g             + num_evals
      fevals_since_chkpnt = fevals_since_chkpnt + num_evals
      CALL Compute_ddot_PAR (localN, d_vec, g_vec,    dd1)
      CALL Compute_ddot_PAR (localN, d_vec, tmp_vec1, dd2)
      pred = -dd1 - (0.5D0 * dd2)
      CALL Normalize_ddot(pred,localN)

C     ---- LOAD THE TRIAL POINT INTO tmp_vec1 AND COMPUTE THE ACTUAL
C     ---- REDUCTION MADE. SAVE THE NEW GRADIENT IN tmp_vec2.
      CALL dcopy (localN, x_vec, 1, tmp_vec1, 1)
      CALL daxpy (localN, 1.0D0, d_vec, 1, tmp_vec1, 1)
      CALL Put_Atom_Positions (localN, tmp_vec1, 1)
      CALL Compute_f_g (f_new, 3*maxown, tmp_vec2)
      num_f_g             = num_f_g             + 1
      fevals_since_chkpnt = fevals_since_chkpnt + 1
      ared = f - f_new

C     ---- EXIT IF THE NEW STEP MAKES NO CHANGE IN THE ENERGY.
      IF (DABS (ared) .LT. DBL_epsilon) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,*) '+++ Last step made no progress'
            IF (outUnit .GT. 0)
     +         WRITE (outUnit,*) '+++ Last step made no progress'
         ENDIF
         CALL Display_Opt_Step
     +      ('x', procNum, outUnit, iteration, f_new, 0.0D0, 0.0D0,
     +       d_2norm, num_f_g, ared, pred, CG_iters, HFTN_return_code)
         GOTO 200
      ENDIF

C     ---- FIND THE MAXIMUM NUMBER OF BOXES ANY ATOM MOVED, INCLUSIVE
C     ---- OVER EACH DIMENSION (BOXES ARE NOT CUBES).
      CALL Find_Biggest_Jump
     +   (boxes_jumped, num_procs_in_jump, periodic_bnds)

C     ---- THE Repartitn_and_Reneighbor SUBROUTINE RELOCATES ATOMS
C     ---- THAT LEAVE THEIR LOCAL PROCESSOR'S SPATIAL DOMAIN; HOWEVER,
C     ---- IT CAN ONLY RELOCATE TO ONE OF THE 8 NEAREST BOX NEIGHBORS.
C     ---- IN ADDITION, THE SPATIAL EXTENT OF EACH PROCESSOR MAY BE
C     ---- REDEFINED BY Repartitn_and_Reneighbor (UNLESS THE PROBLEM HAS
C     ---- PERIODIC BOUNDARIES).  THE WORST CASE OCCURS IF THE TOTAL
C     ---- GEOMETRY SHRINKS IN ONE DIRECTION WHILE AN ATOM MOVES IN THE
C     ---- OPPOSITE DIRECTION.  AS AN ABSOLUTE SAFEGUARD, REJECT THE STEP
C     ---- IF SOME ATOM MOVES FAR ENOUGH THAT A TWO-BOX SHIFT IS POSSIBLE.
      IF (periodic_bnds .EQ. 0) THEN
         jump_threshold = 1.0D0
      ELSE
         jump_threshold = DBLE (num_procs_in_jump)
     +                    / DBLE (1 + num_procs_in_jump)
      ENDIF
      IF (boxes_jumped .GT. jump_threshold) THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,*) '--- Step rejected, atom moved too far'
            IF (outUnit .GT. 0)
     +         WRITE (outUnit,*) '--- Step rejected, atom moved too far'
         ENDIF
      ENDIF


C----------------------------------------------------------------------
C   TEST FOR AGREEMENT BETWEEN ared AND pred.
C----------------------------------------------------------------------

C     --------------------------------------------
      IF ( (pred .GE. DBL_epsilon)     .AND.
     +     ((ared / pred) .GE. TR_eta) .AND.
     +     (boxes_jumped .LE. jump_threshold) ) THEN
C     --------------------------------------------

C        ---- THE STEP IS ACCEPTED.

         CALL Update_TR ('a', pred, ared, TR_eta, d_2norm, TR_size)

C        ---- RESTRICT THE TRUST REGION TO KEEP THE AVERAGE ATOM
C        ---- FROM MOVING MORE THAN ONE BOX PER ITERATION.
         CALL Get_Smallest_Box_Dim (dd1)
         TR_size = MIN ( TR_size, (dd1 * SQRT (DBLE (totalN))) )

C        ---- MAKE THE TRIAL POINT THE CURRENT POSITION.  GET IT FROM
C        ---- LAMMPS IN CASE ADJUSTMENTS WERE MADE FOR PERIODICITY.
         CALL Get_Atom_Positions (localN, x_vec)
         f = f_new
         CALL dcopy (localN, tmp_vec2, 1, g_vec, 1)
         CALL Compute_Inf_Norm_PAR (localN, g_vec, g_Inorm)

         CALL Display_Opt_Step
     +      ('a', procNum, outUnit, iteration, f, g_Inorm, TR_size,
     +       d_2norm, num_f_g, ared, pred, CG_iters, HFTN_return_code)

C        ---- REPARTITION EVERY K-TH STEP.  (AN ALTERNATIVE STRATEGY
C        ---- MIGHT BE BASED ON ACCUMULATED boxes_jumped NUMBERS.)
C        ---- ALSO, FORCE A REPARTITION IF A CHECK-POINT FILE IS TO BE
C        ---- PRODUCED, SO THAT THE RESTART WILL BE EXACT.

         steps_since_repartitn = steps_since_repartitn + 1

c        ---- decide whether to reneighbor
c        ---- yes, if restart files are requested and its been long enough
c        ---- yes, if neighbor delay and neighbor freq are exceeded

         nflag = 0
         if ( (checkpoint_flag .GT. 0) .AND.
     $        (fevals_since_chkpnt .GE. checkpoint_flag) ) nflag = 1
         if ( (steps_since_repartitn .GE. repartitn_delay) .AND.
     $        (steps_since_repartitn .GE. repartitn_freq) ) nflag = 1

         if (nflag .EQ. 1) then

            steps_since_repartitn = 0

C           ---- REDEFINE EACH PROCESSOR'S GEOGRAPHIC DOMAIN, MOVE ATOMS
C           ---- TO THOSE DOMAINS, AND CONSTRUCT NEW NEIGHBOR LISTS FOR
C           ---- THE NON-BONDED INTERACTIONS.
            CALL Repartitn_and_Reneighbor ()
            num_repartitns = num_repartitns + 1

C           ---- IF IT'S TIME, WRITE OUT A RESTART FILE.
            IF (checkpoint_flag .GT. 0. .and.
     $           fevals_since_chkpnt .GE. checkpoint_flag) THEN
               fevals_since_chkpnt = 0
               CALL write_restart ()
            ENDIF

C           ---- GET THE NEW SET OF LOCAL ATOM POSITIONS
C           ---- (localN MAY ALSO CHANGE).
            CALL Get_Atom_Positions (localN, x_vec)

C           ---- RECOMPUTE THE ENERGY AND FORCES BECAUSE NEIGHBOR LISTS
C           ---- MAY CHANGE.  IT IS POSSIBLE FOR f TO INCREASE, BUT THIS
C           ---- IMPLEMENTATION ACCEPTS THE STEP ANYWAY.
            CALL Compute_f_g (f, 3*maxown, g_vec)
            num_f_g             = num_f_g             + 1
            fevals_since_chkpnt = fevals_since_chkpnt + 1
            CALL Compute_Inf_Norm_PAR (localN, g_vec, g_Inorm)

C           ---- DISPLAY NEW ENERGY (ONLY PROCESSOR ZERO).
            IF (procNum .EQ. 0) THEN
               WRITE (6,611) f, g_Inorm, num_f_g
               IF (outUnit .GT. 0)
     +            WRITE (outUnit,611) f, g_Inorm, num_f_g
            ENDIF

         ENDIF


C     ----
      ELSE
C     ----

C        ---- THE STEP IS REJECTED.

         CALL Update_TR ('r', pred, ared, TR_eta, d_2norm, TR_size)

         CALL Display_Opt_Step
     +      ('r', procNum, outUnit, iteration, f_new, 0.0D0, TR_size,
     +       d_2norm, num_f_g, ared, pred, CG_iters, HFTN_return_code)

C        ---- EXIT IF THE TRUST REGION IS TOO SMALL.
         CALL Compute_Inf_Norm_PAR (localN, x_vec, dd1)
         dd1 = MAX (1.0D0, dd1)
         IF (TR_size .LT. (10.0D0 * DBL_epsilon * dd1)) THEN
            IF (procNum .EQ. 0) THEN
               WRITE (6,*) '*** Trust region too small'
               IF (outUnit .GT. 0)
     +            WRITE (outUnit,*) '*** Trust region too small'
            ENDIF
            GOTO 200
         ENDIF

C        ---- SUBTRACT THE TRIAL STEP TO GET BACK TO THE LAST GOOD
C        ---- CONFIGURATION.  USE Get AND Put TO DO THIS IN CASE LAMMPS
C        ---- MADE ADJUSTMENTS FOR PERIODIC BOUNDARY CONDITIONS.
         CALL Get_Atom_Positions (localN, x_vec)
         CALL daxpy (localN, -1.0D0, d_vec, 1, x_vec, 1)
         CALL Put_Atom_Positions (localN, x_vec, 1)
         CALL Get_Atom_Positions (localN, x_vec)

C     -----
      ENDIF
C     -----

      GOTO 100

C----------------------------------------------------------------------
C   END OF MAIN OPTIMIZATION LOOP.
C----------------------------------------------------------------------


 200  CONTINUE

      IF (procNum .EQ. 0) THEN
         WRITE (6,620) num_repartitns
         IF (outUnit .GT. 0)  WRITE (outUnit,620) num_repartitns
      ENDIF

      RETURN

C----------------------------------------------------------------------
C   THE REMAINING LINES OF CODE ARE I/O STATEMENTS REFERENCED EARLIER.
C----------------------------------------------------------------------

 600  FORMAT (/ 'Truncated Newton CG with trust regions on ', I6,
     +        ' unknowns' /)
 601  FORMAT ('     Iter      f(x)    ||grad||_inf   Delta',
     +        '     ||step||  f evals     ared        pred',
     +        '       CG iters')
 602  FORMAT ('    -----  -----------  ----------  ---------',
     +        '  ---------  -------  ----------  ----------',
     +        '    --------')
 610  FORMAT (4X,     I5, 2X, 1PE11.4, 2X, 1PE10.4, 2X, 1PE9.3,
     +        13X,            I7)
 611  FORMAT (3X, '(       ', 1PE11.4, 2X, 1PE10.4, 2X,     9X,
     +        13X,            I7, '      after reneighboring )')

 620  FORMAT (/ '  Number of repartitionings = ', I8 /)

      END


C**********************************************************************
      SUBROUTINE Approximate_Bv
     +   (method, n, x_vec, g_vec, v_vec, v_2norm, tmp_vec,
     +    num_evals, Bv_vec)
C
C       orig:  T. Plantenga  28 Oct 96
C
C   Input/output variable definitions:
C       method           I  'f' for forward, 'c' for central differences
C       n                I  unknowns on this processor
C       x_vec            I  current atom positions
C       g_vec            I  steepest ascent direction at x (used in FD)
C       v_vec            I  direction vector
C       v_2norm          I
C       tmp_vec          IO scratch space of length 3*maxown
C       num_evals         O number f&g evals use
C       Bv_vec            O Hessian times v_vec (approximate)
C
      CHARACTER         method
      INTEGER           n
      DOUBLE PRECISION  x_vec(n), g_vec(n)
      DOUBLE PRECISION  v_vec(n), v_2norm
      DOUBLE PRECISION  tmp_vec(n)
      INTEGER           num_evals
      DOUBLE PRECISION  Bv_vec(n)
C
C
C   This routine approximates the Hessian B times a vector v_vec by a
C   finite difference calculation of the directional derivative.
C   The calling routines specifies forward or central differences.
C
C   External subroutines called:
C        Compute_L2_Norm_PAR       "min_support.f"
C        Put_Atom_Positions        "min_support.f"
C        Compute_f_and_minusg      "min_support.f"
C        daxpy                     "blas/daxpy.f"
C        dcopy                     "blas/dcopy.f"
C        dscal                     "blas/dscal.f"
C**********************************************************************

C-- THESE TWO ACCURACY MEASUREMENTS DETERMINE THE FINITE DIFFERENCE
C-- SPACING.  grad_inacc REFLECTS HOW MANY DIGITS SEEM TO BE LOST IN
C-- COMPUTING THE FORWARD DIFFERENCE [grad(x + hv) - grad(x)] WHEN
C-- DIFFERENT NUMBERS OF PROCESSORS ARE USED.
      DOUBLE PRECISION  DBL_epsilon, grad_inacc
      PARAMETER (DBL_epsilon = 2.2204460492503131D-16)
      PARAMETER (grad_inacc  = 1.0D3)

C     ---- LOCAL VARIABLES.
      DOUBLE PRECISION  h, f

C----------------------------------------------------------------------

C     -------------------------
      IF (method .EQ. 'f') THEN
C     -------------------------

C        ---- FORWARD DIFFERENCE DIRECTIONAL DERIVATIVE.

         h = 2.0D0 * SQRT (grad_inacc * DBL_epsilon) / v_2norm

         CALL dcopy (n, x_vec, 1, tmp_vec, 1)
         CALL daxpy (n, h, v_vec, 1, tmp_vec, 1)
         CALL Put_Atom_Positions (n, tmp_vec, 0)

         CALL Compute_f_g (f, n, Bv_vec)
         num_evals = 1

C        ---- THE FORWARD DIFFERENCE FROM TAYLOR SERIES IS
C        ----    [grad(x + hv) - grad(x)] / h.
         CALL daxpy (n, -1.0D0, g_vec, 1, Bv_vec, 1)
         CALL dscal (n, (1.0D0 / h), Bv_vec, 1)

C     ----
      ELSE
C     ----

C        ---- CENTRAL DIFFERENCE DIRECTIONAL DERIVATIVE.

         h = (3.0D0 * grad_inacc * DBL_epsilon)**(1.0D0 / 3.0D0)
     +       / v_2norm

         CALL dcopy (n, x_vec, 1, tmp_vec, 1)
         CALL daxpy (n, h, v_vec, 1, tmp_vec, 1)
         CALL Put_Atom_Positions (n, tmp_vec, 0)

         CALL Compute_f_g (f, n, Bv_vec)

         CALL dcopy (n, x_vec, 1, tmp_vec, 1)
         CALL daxpy (n, -h, v_vec, 1, tmp_vec, 1)
         CALL Put_Atom_Positions (n, tmp_vec, 0)

         CALL Compute_f_g (f, n, tmp_vec)
         num_evals = 2

C        ---- THE CENTRAL DIFFERENCE FROM TAYLOR SERIES IS
C        ----    [grad(x + hv) - grad(x - hv)] / 2h.
         CALL daxpy (n, -1.0D0, tmp_vec, 1, Bv_vec, 1)
         CALL dscal (n, (0.5D0 / h), Bv_vec, 1)

C     -----
      ENDIF
C     -----

C     ---- RESTORE THE GLOBAL VARIABLE x TO THE ORIGINAL ATOM POSITIONS.
C     ---- NOTE THAT THE THIRD ARGUMENT IS ZERO IN ALL CALLS MADE BY
C     ---- THIS SUBROUTINE.  THIS PREVENTS THE LAMMPS "true" FLAGS FROM
C     ---- CHANGING, AND IS REASONABLE ONLY BECAUSE ORIGINAL POSITIONS
C     ---- ARE EXACTLY RESTORED.
      CALL Put_Atom_Positions (n, x_vec, 0)

      RETURN
      END


C**********************************************************************
      SUBROUTINE Compute_HFTN_Step
     +   (n, totalN, x_vec, g_vec, TR_size,
     +    outer_iteration, max_num_f_g, num_f_g, findiff_flag,
     +    d, d_2norm, inner_iterations, num_evals, return_code)
C
C       orig:  T. Plantenga  25 Oct 96
C
C   Input/output variable definitions:
C       n                I  unknowns on this processor
C       totalN           I  total number of unknowns
C       x_vec            I  current atom positions
C       g_vec            I  gradient at x_vec
C       TR_size          I
C       outer_iteration  I  used to set the stop tolerance
C       max_num_f        I
C       num_f_g          I  current total, used to set stop tolerance
C       findiff_flag     I
C       d                 O trial step
C       d_2norm           O
C       inner_iterations  O
C       num_evals         O number f&g evals to get the step
C       return_code       O character string
C
      INTEGER           n, totalN
      DOUBLE PRECISION  x_vec(n), g_vec(n)
      DOUBLE PRECISION  TR_size
      INTEGER           outer_iteration, inner_iterations
      INTEGER           max_num_f_g, num_f_g, num_evals
      CHARACTER         findiff_flag
      DOUBLE PRECISION  d(n), d_2norm
      INTEGER           return_code
C
C
C   This routine solves the trust region subproblem by Steihaug's method.
C   It performs CG iterations until either
C      1 - "ng" a search direction of negative curvature is found,
C      2 - "TR" the trust region bound is reached,
C      3 - "Nw" approximate convergence is achieved,
C      4 - "it" the maximum number of iterations is reached, or
C      5 - "FD" a CG step causes the model objective to increase.
C   The string above is returned in return_code.
C
C   External subroutines called:
C        Compute_L2_Norm_PAR       "min_support.f"
C        Compute_ddot_PAR          "min_support.f"
C        daxpy                     "blas/daxpy.f"
C        dcopy                     "blas/dcopy.f"
C        dscal                     "blas/dscal.f"
C**********************************************************************

C-- MAXIMUM PROBLEM SIZES FOR ONE PROCESSOR ARE DEFINED HERE.
      INCLUDE 'param.h'

      DOUBLE PRECISION  DBL_epsilon
      PARAMETER (DBL_epsilon = 2.2204460492503131D-16)

C-- USE THIS TO SET THE MAXIMUM NUMBER OF INNER ITERATIONS.
      INTEGER  CG_ITER_LIMIT
      PARAMETER (CG_ITER_LIMIT = 2)

C     ---- LOCAL VARIABLES.
      INTEGER           i
      INTEGER           max_iters
      DOUBLE PRECISION  r(3*maxown), p(3*maxown), Bp(3*maxown)
      DOUBLE PRECISION  pBp
      DOUBLE PRECISION  g_2norm, p_2norm
      DOUBLE PRECISION  rr_old, rr_new, stop_tol
      DOUBLE PRECISION  alpha, beta
      DOUBLE PRECISION  dd1, dd2, tmp_vec(3*maxown)

C----------------------------------------------------------------------

C     ---- INITIALIZE r = -g, p = r, rr_old = r'r.
      CALL dcopy (n, g_vec, 1, r, 1)
      CALL dscal (n, -1.0D0, r, 1)
      CALL dcopy (n, r, 1, p, 1)
      CALL Compute_L2_Norm_PAR (n, g_vec, g_2norm)
      rr_old = g_2norm * g_2norm
      p_2norm = g_2norm

C     ---- INITIALIZE d = 0.
      DO i = 1, n
         d(i) = 0.0D0
      ENDDO
      d_2norm = 0.0D0

C     ---- SET THE STOP TOLERANCE (FOR QUADRATIC RATE OF CONVERGENCE).
      stop_tol = 0.1D0 / DBLE (outer_iteration)
      stop_tol = MIN (stop_tol, g_2norm)

      max_iters = CG_ITER_LIMIT * totalN

      inner_iterations = 0
      num_evals        = 0
      return_code      = 0

C======================================================================
C   ---- Conjugate gradient ----
C
C   LOOP
C      Compute Bp by finite difference
C      IF p'Bp < small positive number
C         THEN Extend to trust region along p and RETURN
C      alpha = r'r / p'Bp
C      d_new = d + alpha*p
C      IF ||d_new||_2 > delta
C         THEN Extend to trust region along p and RETURN
C      r_new = r - alpha*Bp
C      IF ||r_new||_2 < stop_tol*||r_0||_2
C         THEN RETURN
C      beta = r_new'r_new / r'r
C      p_new = r_new + beta*p
C   CONTINUE
C======================================================================

C----------------------------------------------------------------------
C   BEGIN CG INNER OPTIMIZATION LOOP.
C----------------------------------------------------------------------

 100  CONTINUE

C     ---- STOP IF THE ITERATION OR FUNCTION COUNT IS TOO LARGE.
      IF ( ((num_f_g + num_evals) .GE. max_num_f_g) .OR.
     +     (inner_iterations .GE. max_iters)             ) THEN
         return_code = 4
         GOTO 200
      ENDIF

      inner_iterations = inner_iterations + 1

      CALL Approximate_Bv
     +   (findiff_flag, n, x_vec, g_vec, p, p_2norm, tmp_vec, i, Bp)
      num_evals = num_evals + i

C     ---- FORM p'Bp AND CHECK THAT IT IS SUFFICIENTLY POSITIVE.
C     ---- IF NOT, EXTEND THE STEP TO THE TRUST REGION BOUNDARY AND EXIT.
C     ---- THE RAYLEIGH QUOTIENT CHECK I USE IS QUITE SIMILAR TO SCHLICK
C     ---- AND OVERTON (J. COMP. CHEM., 1987, 8:1025).  NOTE THAT
C     ---- Extend_to_TR WILL LOOK IN BOTH DIRECTIONS ALONG p.
      CALL Compute_ddot_PAR (n, p, Bp, pBp)
      dd1 = p_2norm * p_2norm
      IF (pBp .LT. (dd1 * SQRT (DBL_epsilon))) THEN
         CALL Compute_ddot_PAR (n, p, g_vec, dd1)
         CALL Compute_ddot_PAR (n, d, Bp, dd2)
         CALL Extend_to_TR (n, p, TR_size, dd1, pBp, dd2, d)
         CALL Compute_L2_Norm_PAR (n, d, d_2norm)
         return_code = 1
         GOTO 200
      ENDIF

C     ---- COMPUTE alpha, THE DISTANCE TO THE QUADRATIC MINIMUM ALONG p.
      alpha = rr_old / pBp

C     ---- UPDATE THE STEP AS tmp_vec = d + alpha*p.
      CALL dcopy (n, d, 1, tmp_vec, 1)
      CALL daxpy (n, alpha, p, 1, tmp_vec, 1)

C     ---- CHECK THAT THE NEW STEP IS STILL INSIDE THE TRUST REGION.
C     ---- IF NOT, EXTEND THE STEP TO THE TRUST REGION BOUNDARY AND EXIT.
      CALL Compute_L2_Norm_PAR (n, tmp_vec, d_2norm)
      IF (d_2norm .GE. TR_size) THEN
         CALL Compute_ddot_PAR (n, p, g_vec, dd1)
         CALL Compute_ddot_PAR (n, d, Bp, dd2)
         CALL Extend_to_TR (n, p, TR_size, dd1, pBp, dd2, d)
         CALL Compute_L2_Norm_PAR (n, d, d_2norm)
         return_code = 2
         GOTO 200
      ENDIF

C     ---- ACCEPT THE STEP, OVERWRITING d WITH tmp_vec.
      CALL dcopy (n, tmp_vec, 1, d, 1)

C     ---- COMPUTE THE RESIDUAL AT THE NEW POINT.
      CALL daxpy (n, -alpha, Bp, 1, r, 1)
      CALL Compute_ddot_PAR (n, r, r, rr_new)


C     ---- CHECK FOR CONVERGENCE USING THE DEMBO, EISENSTAT, AND STEIHAUG
C     ---- (SIAM J. NUMER. ANAL, 1982, 19:400) TEST THAT ALLOWS A QUADRATIC
C     ---- RATE OF CONVERGENCE OF THE OVERALL OPTIMIZATION ALGORITHM.

C     ---- THEORETICAL (EXACT ARITHMETIC) CG STOP TEST.
      IF (rr_new .EQ. 0.0D0) THEN
         return_code = 3
         GOTO 200
      ENDIF

      IF (SQRT (rr_new) .LT. (stop_tol * g_2norm)) THEN
         return_code = 3
         GOTO 200
      ENDIF

C     ---- COMPUTE THE NEXT CONJUGATE SEARCH DIRECTION, p = r + beta*p.
      beta = rr_new / rr_old
      rr_old = rr_new
      CALL dscal (n, beta, p, 1)
      CALL daxpy (n, 1.0D0, r, 1, p, 1)
      CALL Compute_L2_Norm_PAR (n, p, p_2norm)


      GOTO 100

C----------------------------------------------------------------------
C   END OF CG INNER OPTIMIZATION LOOP.
C----------------------------------------------------------------------

 200  CONTINUE
      
      RETURN
      END


C**********************************************************************
      SUBROUTINE Extend_to_TR
     +   (n, p, TR_size, g_dot_p, pBp, dBp, d)
C
C       orig:  T. Plantenga   4 Nov 96
C
      INTEGER           n
      DOUBLE PRECISION  p(n)
      DOUBLE PRECISION  TR_size, g_dot_p, pBp, dBp
      DOUBLE PRECISION  d(n)
C
C
C   This routine calculates a point along p that intersects a spherical
C   trust region of radius TR_size.  Two such points exist; this routine
C   returns the one which minimizes a quadratic objective with gradient
C   g and Hessian B.
C
C   The routine is adapted from "c_libs/trsubs.c".  Cannot use my
C   library routine because it has nonparallel ddot calls.
C
C   External subroutines called:
C        Compute_ddot_PAR          "min_support.f"
C        daxpy                     "blas/daxpy.f"
C**********************************************************************

      DOUBLE PRECISION  DBL_epsilon
      PARAMETER (DBL_epsilon = 2.2204460492503131D-16)

C     ---- LOCAL VARIABLES
      DOUBLE PRECISION  tau_neg, tau_pos
      DOUBLE PRECISION  p_dot_p, d_dot_p, d_dot_d
      DOUBLE PRECISION  dd1, discr
      DOUBLE PRECISION  obj_neg, obj_pos

C----------------------------------------------------------------------

      CALL Compute_ddot_PAR (n, p, p, p_dot_p)
      CALL Compute_ddot_PAR (n, d, p, d_dot_p)
      CALL Compute_ddot_PAR (n, d, d, d_dot_d)

C     ---- CHECK IF p IS A ZERO VECTOR.
      IF (p_dot_p .EQ. 0.0D0)  RETURN


C     ---- THE else STATEMENT BELOW SAFEGUARDS THE CALCULATION OF tau FOR
C     ---- THE CASE THAT TR_size NEARLY EQUALS x'x  (CF. ACCURACY AND
C     ---- STABILITY OF NUMERICAL ALGORITHMS, N.J. HIGHAM, P. 11).

      dd1 = p_dot_p * ((TR_size * TR_size) - d_dot_d)

      IF (dd1 .GT. ((d_dot_p * d_dot_p) * SQRT (DBL_epsilon))) THEN
         discr = SQRT (dd1 + (d_dot_p * d_dot_p))
         tau_neg = (-d_dot_p - discr) / p_dot_p
         tau_pos = (-d_dot_p + discr) / p_dot_p

      ELSE
         discr = SQRT (dd1 + (d_dot_p * d_dot_p))
         IF (d_dot_p .GE. 0.0D0) THEN
            tau_neg = (-d_dot_p - discr) / p_dot_p
            tau_pos = ((d_dot_d - (TR_size * TR_size)) / p_dot_p)
     +                / tau_neg
            if (tau_pos .LT. 0.0D0)  tau_pos = 0.0D0
         ELSE
            tau_pos = (-d_dot_p + discr) / p_dot_p
            tau_neg = ((d_dot_d - (TR_size * TR_size)) / p_dot_p)
     +                / tau_pos
            if (tau_neg .LT. 0.0D0)  tau_neg = 0.0D0
         ENDIF
      ENDIF

C     ---- EVALUATE THE QUADRATIC OBJECTIVE FOR tau_pos AND tau_neg.

      obj_neg = (tau_neg * (g_dot_p + dBp))
     +          + (0.5D0 * tau_neg * tau_neg * pBp)
      obj_pos = (tau_pos * (g_dot_p + dBp))
     +          + (0.5D0 * tau_pos * tau_pos * pBp)

      IF (obj_neg .LE. obj_pos) THEN
         CALL daxpy (n, tau_neg, p, 1, d, 1)
      ELSE
         CALL daxpy (n, tau_pos, p, 1, d, 1)
      ENDIF


      RETURN
      END


C**********************************************************************
      SUBROUTINE Update_TR
     +   (step_flag, pred, ared, eta, d_2norm, TR_size)
C
C       orig:  T. Plantenga  25 Oct 96
C
C   Input/output variable definitions:
C       step_flag        I  'a' if step accepted, 'r' if rejected
C       pred             I  predicted reduction
C       ared             I  actual reduction
C       eta              I  trust region acceptance threshold
C       d_2norm          I
C       TR_size           O
C
      CHARACTER         step_flag
      DOUBLE PRECISION  pred, ared, eta, d_2norm
      INTEGER           n
      DOUBLE PRECISION  TR_size
C
C
C   This routine adjusts the trust region size.  The size is increased
C   if the step was accepted, decreased if rejected.
C**********************************************************************

C     ---- LOCAL VARIABLES
      DOUBLE PRECISION  dd1

C----------------------------------------------------------------------

      IF (step_flag .EQ. 'a') THEN

C        ---- INCREASE THE TRUST REGION SIZE.

         IF ((ared / pred) .GE. 0.9D0) THEN
            dd1 = 2.0D0 * d_2norm
         ELSE
            dd1 = 1.0D0 * d_2norm
         ENDIF

         IF (dd1 .GT. TR_size)  TR_size = dd1

      ELSE

C        ---- DECREASE THE TRUST REGION SIZE USING LINEAR INTERPOLATION.

         dd1 = (1.0D0 - eta) / (1.0D0 - (ared / pred))
         IF (dd1 .LT. 0.1D0)  dd1 = 0.1D0
         IF (dd1 .GT. 0.5D0)  dd1 = 0.5D0
         TR_size = dd1 * d_2norm

      ENDIF

      RETURN
      END


C**********************************************************************
      SUBROUTINE Display_Opt_Step
     +   (step_flag, procNum, outUnit, iteration, f, g_Inorm, TR_size,
     +    d_2norm, num_f_g, ared, pred, CG_iters, HFTN_return_code)
C
C       orig:  T. Plantenga  31 Oct 96
C
      CHARACTER         step_flag
      INTEGER           procNum, outUnit
      INTEGER           iteration
      DOUBLE PRECISION  f, g_Inorm, TR_size, d_2norm
      INTEGER           num_f_g
      DOUBLE PRECISION  ared, pred
      INTEGER           CG_iters, HFTN_return_code
C
C
C   This routine writes step information to the screen and/or the
C   output file.  The variable "step_flag" specifies whether the step
C   was accepted ('a'), rejected ('r'), or not evaluated ('x').
C**********************************************************************

C     ---- LOCAL VARIABLES
      CHARACTER  HFTN_status*2

C----------------------------------------------------------------------

      IF (HFTN_return_code .EQ. 1) THEN
         HFTN_status = 'ng'
      ELSE IF (HFTN_return_code .EQ. 2) THEN
         HFTN_status = 'TR'
      ELSE IF (HFTN_return_code .EQ. 3) THEN
         HFTN_status = 'Nw'
      ELSE IF (HFTN_return_code .EQ. 4) THEN
         HFTN_status = 'it'
      ELSE IF (HFTN_return_code .EQ. 5) THEN
         HFTN_status = 'FD'
      ELSE
         HFTN_status = '??'
      ENDIF


      IF (step_flag .EQ. 'a') THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,600)
     +              iteration, f, g_Inorm, TR_size, d_2norm,
     +              num_f_g, ared, pred, CG_iters, HFTN_status
            IF (outUnit .GT. 0)
     +           WRITE (outUnit,600)
     +              iteration, f, g_Inorm, TR_size, d_2norm,
     +              num_f_g, ared, pred, CG_iters, HFTN_status
         ENDIF

      ELSE IF (step_flag .EQ. 'r') THEN
         IF (procNum .EQ. 0) THEN
            WRITE (6,601)
     +              iteration, f,          TR_size, d_2norm,
     +              num_f_g, ared, pred, CG_iters, HFTN_status
            IF (outUnit .GT. 0)
     +           WRITE (outUnit,601)
     +              iteration, f,          TR_size, d_2norm,
     +              num_f_g, ared, pred, CG_iters, HFTN_status
         ENDIF

      ELSE
         IF (procNum .EQ. 0) THEN
            WRITE (6,602)
     +              iteration, f,                   d_2norm,
     +              num_f_g, ared, pred, CG_iters, HFTN_status
            IF (outUnit .GT. 0)
     +           WRITE (outUnit,602)
     +              iteration, f,                   d_2norm,
     +              num_f_g, ared, pred, CG_iters, HFTN_status
         ENDIF

      ENDIF

      RETURN


 600  FORMAT (3X,    1X, I5, 2X, 1PE11.4, 2X, 1PE10.4, 2X, 1PE9.3,
     +        2X, 1PE9.3, 2X, I7, 2X, 1PE10.3, 2X, 1PE10.3,
     +        2X, I6, 2X, A2)
 601  FORMAT ('rej', 1X, I5, 2X, 1PE11.4, 2X, 10X,     2X, 1PE9.3,
     +        2X, 1PE9.3, 2X, I7, 2X, 1PE10.3, 2X, 1PE10.3,
     +        2X, I6, 2X, A2)
 602  FORMAT (' - ', 1X, I5, 2X, 1PE11.4, 2X, 10X,     2X,     9X,
     +        2X, 1PE9.3, 2X, I7, 2X, 1PE10.3, 2X, 1PE10.3,
     +        2X, I6, 2X, A2)

      END


C***********************End of source file*****************************
