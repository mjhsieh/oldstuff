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

c user-specified diagnostic routines
c user needs to
c   (1) modify diag_init routine in 2 places (see ***)
c   (2) modify diagnostic routine in 1 place (see ***)
c   (3) add own routines for computing diagnostics (see samples at ***)
c       can put them in this file or in separate file(s)
c       if in separate file(s), have to add file(s) to Makefile

c ----------------------------------------------------------------------
c initializization of total # of diagnostic routines and their names

      subroutine diag_init
      include "lammps.h"

c *** start of user-modified section
c set numdiag = # of diagnostic routines defined
c numdiag must be <= maxdiag setting in param.h

      numdiag = 0     

c      numdiag = 2

c *** end of user-modified section

c error check and initialization

      if (numdiag.gt.maxdiag)
     $     call error('Too many diagnostic routines - boost maxdiag')

      do idiag = 1,numdiag
        ndiag(idiag) = 0
        diagfileflag(idiag) = 0
      enddo

c *** start of user-modified section
c assign name of each diagnostic routine from 1 to numdiag
c these names are then used in diagnostic commands in input script
c each name must be a lowercase word with no spaces, 16 characters max

c      diagnames(1) = 'diffusion'
c      diagnames(2) = 'energy'

c *** end of user-modified section

      return
      end


c ----------------------------------------------------------------------
c driver that calls user diagnostic routines

      subroutine diagnostic
      include "lammps.h"

c determine which user routines should be called this timestep
c set diagcall = 0 if not called, 1 if called
c only call if this diagnostic is currently enabled
c only call if diagnext = current timestep
c allow routine to be called even if previously called on this timestep,
c   this allows user to do setup at beginning of a new run,
c   diagnostic routine itself can choose to exit on this condition
c update diagnext to the next time this routine should be called

      do idiag = 1,numdiag
        diagcall(idiag) = 0
        if (ndiag(idiag).gt.0) then
          if (ntimestep.eq.diagnext(idiag)) then
            diagcall(idiag) = 1
            diagnext(idiag) = min(ntimestep+ndiag(idiag),ntime_last)
          endif
        endif
      enddo

c *** start of user-modified section
c call each user diagnostic routine from 1 to numdiag
c argument of called routine is same as index of diagcall
c routine names can be arbitrary
c routines must be provided by user

      if (diagcall(1).eq.1) call diag_diffusion(1)
      if (diagcall(2).eq.1) call diag_energy(2)

c      if (diagcall(3).eq.1) call diag_xxx(3)
c      if (diagcall(4).eq.1) call diag_xxx(4)
c      if (diagcall(5).eq.1) call diag_xxx(5)

c *** end of user-modified section

c ndiag_next = earliest timestep any diagnostic routine should next be called
c only check enabled diagnostics

      ndiag_next = ntime_last + 1
      do idiag = 1,numdiag
        if (ndiag(idiag).gt.0)
     $       ndiag_next = min(ndiag_next,diagnext(idiag))
      enddo

      return
      end


c ----------------------------------------------------------------------
c *** start of user-modified section

c sample user routines, called from "diagnostic" routine up above
c one routine for each diagnostic from 1 to numdiag
c idiag argument can be used to index diagnostic params and file info
c   diagnparams(idiag) = # of user-specified parameters for this diagnostic
c   diagparam(n,idiag) = 1 to n parameters for this diagnostic
c   diagfileflag(idiag) = 0,1 if diagnostic file is closed/open
c also, itime in "lammps.h" tells when in the run this routine is called
c   itime = 0, before run (from start routine)
c   itime < nsteps, during run (from integrate routine)
c   itime = nsteps, last call of run (from integrate routine)

c 1st example routine

      subroutine diag_diffusion(idiag)
      include "lammps.h"
      integer idiag

c check if this routine was called previously on same timestep
c could have occurred in a previous run
c if so, may want to just exit

      if (ntimestep.eq.diagprev(idiag)) return

c ifile = handle for file opened when the "diagnostic" command was issued

      ifile = 20 + idiag

c set diagprev 

      diagprev(idiag) = ntimestep

c itime = 0 on 1st call to diagnostic in each run
c can do setup here if wish, open own files, etc

c      if (itime.eq.0) then
c        setup code
c      endif

c -------------------------------------------
c user code to compute the desired diagnostic
c this is just an example of simple code

c write file header if this is 1st call to this diagnostic

      if (diagfileflag(idiag).eq.1) then
        diagfileflag(idiag) = 2
        if (node.eq.0)
     $       write (ifile,*) 'Dummy diagnostic for diffusion:'
      endif

c do a computation, using params, write result to file

      if (node.eq.0) write (ifile,*) ntimestep,
     $     (diagparam(i,idiag),i=1,diagnparams(idiag))

c end of user code for the diagnostic
c -------------------------------------------

c itime = nsteps on final call to diagnostic in each run
c can do clean-up here if wish
c note that this call occurs within integrate,
c   so you may wish to have also computed the diagnostic

c      if (itime.eq.nsteps)
c        clean-up code
c      endif

      return
      end


c 2nd example routine

      subroutine diag_energy(idiag)
      include "lammps.h"
      integer idiag

c check if this routine was called previously on same timestep
c could have occurred in a previous run
c if so, may want to just exit

      if (ntimestep.eq.diagprev(idiag)) return

c ifile = handle for file opened when the "diagnostic" command was issued

      ifile = 20 + idiag

c set diagprev 

      diagprev(idiag) = ntimestep

c itime = 0 on 1st call to diagnostic in each run
c can do setup here if wish, open own files, etc

c      if (itime.eq.0) then
c        setup code
c      endif

c -------------------------------------------
c user code to compute the desired diagnostic
c this is just an example of simple code

c write file header if this is 1st call to this diagnostic

      ifile = 20 + idiag

c write file header if this is 1st call to this diagnostic

      if (diagfileflag(idiag).eq.1) then
        diagfileflag(idiag) = 2
        if (node.eq.0)
     $       write (ifile,*) 'Dummy diagnostic for energy:'
      endif

c do a computation, using params, write result to file

      if (node.eq.0) write (ifile,*) ntimestep,
     $     (diagparam(i,idiag),i=1,diagnparams(idiag))

c end of user code for the diagnostic
c -------------------------------------------

c itime = nsteps on final call to diagnostic in each run
c can do clean-up here if wish
c note that this call occurs within integrate,
c   so you may wish to have also computed the diagnostic

c      if (itime.eq.nsteps)
c        clean-up code
c      endif

      return
      end

c *** end of user-modified section

c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c user should not modify routines from here to end-of-file

c ----------------------------------------------------------------------
c open diagnostic file

      subroutine diag_open(idiag)
      include "lammps.h"

      ifile = 20 + idiag

      if (node.eq.0) then
        open (ifile,file=diagfile(idiag),status='unknown')
        close (ifile,status='delete')
        open (ifile,file=diagfile(idiag),status='new')
      endif

c no diagnostics have been written into new file
      
      diagprev(idiag) = -1

      return
      end


c ----------------------------------------------------------------------
c close diagnostic file

      subroutine diag_close(idiag)
      include "lammps.h"

      ifile = 20 + idiag
      if (node.eq.0) close (ifile)

      return
      end


c ----------------------------------------------------------------------
c lookup the name of a diagnostic and return an index
c error return of 0 means name not found

      subroutine diag_lookup(name,idiag)
      include "lammps.h"
      character*(*) name
      integer idiag

      do idiag = 1,numdiag
        if (name.eq.diagnames(idiag)) return
      enddo

      idiag = 0

      return
      end
