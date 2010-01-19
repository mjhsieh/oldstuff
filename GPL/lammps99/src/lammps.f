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
      program lammps
      include "lammps.h"
      
      call initialize
      
      do infinity = 1,100000000
        
        call input(iopflag)

        if (iopflag.eq.1) then
          call read_data
        else if (iopflag.eq.2) then
          call temp_create
        else if (iopflag.eq.3) then
          call read_restart
        else if (iopflag.eq.4) then
          call set_fix
        else if (iopflag.eq.5) then
          call start
          if (nrespa.eq.0) then
            call integrate
          else
            call integrate_respa
          endif
          call finish
        else if (iopflag.eq.6) then
          call start
          call minimize
          call finish
        endif
        
      enddo
      
      end
