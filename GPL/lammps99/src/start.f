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
c beginning of timestep

      subroutine start
      include "lammps.h"
      include "mpif.h"
      integer iopflag

      if (node.eq.0) write (6,*) 'Setting up run ...'

      if (readflag.eq.0)
     $     call error('Data file must be read before run')

c other errors/warnings

      if (peratom_comm.gt.maxperatom)
     $     call error('Too many atom quantities - boost maxperatom')

      if ((neighfreq.lt.neighdelay.and.
     $     mod(neighdelay,neighfreq).ne.0).or.
     $     (neighdelay.gt.0.and.neighdelay.lt.neighfreq)) then
        if (node.eq.0) write (6,*) 'WARNING: Mismatch between ',
     $       ' neighbor frequency and delay factors'
      endif

      if ((tempflag.eq.2.or.tempflag.eq.3).and.nrestart.gt.0) then
        if (node.eq.0) write (6,*)
     $       'WARNING: Temp control incompatible with restart files'
      endif

      if ((coulstyle.ge.3).and.
     $     (perflagx.eq.1.or.perflagy.eq.1.or.perflagz.eq.1))
     $     call error('Ewald/PPPM requires periodic BC')

      if (coulstyle.ge.3.and.idimension.eq.2)
     $     call error('Ewald/PPPM requires 3-d')

      if ((pressflag.eq.1.or.
     $     xpressflag.eq.1.or.ypressflag.eq.1.or.zpressflag.eq.1).and.
     $     (perflagx.eq.1.or.perflagy.eq.1.or.perflagz.eq.1))
     $     call error('Constant pressure requires periodic BC')

      if ((nonstyle.eq.3.or.nonstyle.eq.4).and.coulstyle.ne.0)
     $     call error('Mismatch between nonbond and coulomb styles')

      if (nonstyle.eq.5.and.coulstyle.eq.2)
     $     call error('Cannot combine class2 with coul/smooth style')

      if (nrespa.eq.1.and.tempflag.eq.3) call error(
     $     'RESPA does not support Langevin temp control')

      if (nrespa.eq.1.and.nfixes.gt.0)
     $     call error('RESPA does not support constraints')

c safety check that length of nlist is sufficient to use it as work
c  array for typechecking - in start.f

      if (ntypes.gt.maxown.or.
     $     nbondtypes.gt.maxown.or.nangletypes.gt.maxown.or.
     $     ndihedtypes.gt.maxown.or.nimprotypes.gt.maxown)
     $     call error('Too many types for typechecking - boost maxown')

c safety check that length of nlist is sufficient to use it as work
c  array for typechecking - in start.f
c also serves as safety check that length of buf4 is sufficient to use it 
c  as work array for constraint summations - in fix.f

      if (nfixes.gt.maxown)
     $     call error('Too many fixes for check/sum - boost maxown')

c check for inefficient use of Coulomb style with no charges

      if (coulstyle.gt.0) then

        iflag = 0
        i = atompnt
        do ii = 1,nlocal
          if (q(i).ne.0.0) iflag = 1
          i = list(i)
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.eq.0)
     $       call error('Coulomb style set but all charges = 0')

      endif

c check for skipped constraint number

      if (nfixes.gt.0) then
        do i = 1,nfixes
          if (fixstyle(i).eq.0) call error('Missing fix')
        enddo
      endif

c check no atom has been assigned to missing constraint
c check that every constraint has at least one atom assigned to it
c  use nlist and nliststart as work arrays

      do i = 1,nfixes
        nlist(i) = 0
      enddo

      iflag = 0
      i = atompnt
      do ii = 1,nlocal
        itmp = fix(i)
 40     if (itmp.gt.0) then
          iwhich = mod(itmp,maxfix)
          if (iwhich.gt.nfixes) then
            iflag = 1
          else
            nlist(iwhich) = 1
          endif
          itmp = itmp/maxfix
          goto 40
        endif
        i = list(i)
      enddo

      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) call error('Atom assigned to missing fix')

      call mpi_allreduce(nlist,nliststart,nfixes,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      do i = 1,nfixes
        if (nliststart(i).eq.0)
     $       call error('No atom assigned to a fix')
      enddo

c do not allow thermostatting fixes in conjunction with global temp control

      if (tempflag.gt.0) then
        iflag = 0
        do i = 1,nfixes
          if (fixstyle(i).eq.4.or.fixstyle(i).eq.5.or.
     $         fixstyle(i).eq.6) iflag = 1
        enddo
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
        if (iflag.gt.0) call error(
     $       'Temperature control incompatible with temperature fix')
      endif

c check consistency of atom types and nonbond type pairs
c  use nlist and nliststart as work arrays

      do i = 1,ntypes
        nlist(i) = 0
      enddo

      iflag = 0
      i = atompnt
      do ii = 1,nlocal
        if (type(i).le.0.or.type(i).gt.ntypes) then
          iflag = 1
        else
          nlist(type(i)) = 1
        endif
        i = list(i)
      enddo
      
      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) call error('Illegal atom type')

      call mpi_allreduce(nlist,nliststart,ntypes,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      do i = 1,ntypes
        if (nliststart(i).eq.0.and.node.eq.0)
     $       write (6,*) 'WARNING: Non-existent atom type',i
      enddo

      iflag = 0
      do i = 1,ntypes
        do j = i,ntypes
          if (nontypeflag(i,j).ne.nonstyle) iflag = 1
        enddo
      enddo

      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) call error('Undefined nonbond type pair')

c check consistency of bond types

      if (nbonds.gt.0) then

        do i = 1,nbondtypes
          nlist(i) = 0
        enddo

        iflag = 0
        i = atompnt
        do ii = 1,nlocal
          do j = 1,numbond(i)
            if (bondtype(j,i).le.0.or.
     $           bondtype(j,i).gt.nbondtypes) then
              iflag = 1
            else
              nlist(bondtype(j,i)) = 1
            endif
          enddo
          i = list(i)
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Illegal bond type')

        call mpi_allreduce(nlist,nliststart,nbondtypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,nbondtypes
          if (nliststart(i).eq.0.and.node.eq.0)
     $         write (6,*) 'WARNING: Non-existent bond type',i
        enddo

        iflag = 0
        do i = 1,nbondtypes
          if (bondtypeflag(i).ne.bondstyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Undefined bond type')

      endif

c check consistency of angle types

      if (nangles.gt.0) then

        do i = 1,nangletypes
          nlist(i) = 0
        enddo

        iflag = 0
        i = atompnt
        do ii = 1,nlocal
          do j = 1,numangle(i)
            if (angletype(j,i).le.0.or.
     $           angletype(j,i).gt.nangletypes) then
              iflag = 1
            else
              nlist(angletype(j,i)) = 1
            endif
          enddo
          i = list(i)
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Illegal angle type')

        call mpi_allreduce(nlist,nliststart,nangletypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,nangletypes
          if (nliststart(i).eq.0.and.node.eq.0)
     $         write (6,*) 'WARNING: Non-existent angle type',i
        enddo

        iflag = 0
        do i = 1,nangletypes
          if (angletypeflag(i).ne.anglestyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Undefined angle type')

      endif

c check consistency of dihedral types

      if (ndihedrals.gt.0) then

        do i = 1,ndihedtypes
          nlist(i) = 0
        enddo

        iflag = 0
        i = atompnt
        do ii = 1,nlocal
          do j = 1,numdihed(i)
            if (dihedtype(j,i).le.0.or.
     $           dihedtype(j,i).gt.ndihedtypes) then
              iflag = 1
            else
              nlist(dihedtype(j,i)) = 1
            endif
          enddo
          i = list(i)
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Illegal dihed type')

        call mpi_allreduce(nlist,nliststart,ndihedtypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,ndihedtypes
          if (nliststart(i).eq.0.and.node.eq.0)
     $         write (6,*) 'WARNING: Non-existent dihedral type',i
        enddo

        iflag = 0
        do i = 1,ndihedtypes
          if (dihedtypeflag(i).ne.dihedstyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Undefined dihedral type')

      endif

c check consistency of improper types

      if (nimpropers.gt.0) then

        do i = 1,nimprotypes
          nlist(i) = 0
        enddo

        iflag = 0
        i = atompnt
        do ii = 1,nlocal
          do j = 1,numimpro(i)
            if (improtype(j,i).le.0.or.
     $           improtype(j,i).gt.nimprotypes) then
              iflag = 1
            else
              nlist(improtype(j,i)) = 1
            endif
          enddo
          i = list(i)
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Illegal impro type')

        call mpi_allreduce(nlist,nliststart,nimprotypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,nimprotypes
          if (nliststart(i).eq.0.and.node.eq.0)
     $         write (6,*) 'WARNING: Non-existent improper type',i
        enddo

        iflag = 0
        do i = 1,nimprotypes
          if (improtypeflag(i).ne.improstyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) call error('Undefined improper type')

      endif

c pre-compute coeffs and cutoff distances for all nonbond/Coulombic styles
c  assumes noncoeffs are set for all type pairs with i <= j

      cutcoulsq = cutcoul*cutcoul
      if (coulstyle.eq.2) cutcoulintsq = cutcoulint*cutcoulint
      triggersq = (skin/2.0)*(skin/2.0)
      cutforce = 0.0

      do i = 1,ntypes
        do j = i,ntypes

          if (nonstyle.eq.0) then
            cutljsq(i,j) = 0.0
            cut = cutcoul
          else if (nonstyle.eq.1) then
            lj1(i,j) = 48.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj2(i,j) = 24.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj4(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            cutljsq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            if (offsetflag.eq.0.or.noncoeff3(i,j).eq.0.0) then
              offset(i,j) = 0.0
            else
              ratio = noncoeff2(i,j)/noncoeff3(i,j)
              offset(i,j) = 4.0*noncoeff1(i,j)*(ratio**12 - ratio**6)
            endif
            cut = max(cutcoul,noncoeff3(i,j))
          else if (nonstyle.eq.2) then
            lj1(i,j) = 48.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj2(i,j) = 24.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj4(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            cutljinner(i,j) = noncoeff3(i,j)
            cutljinnersq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            cutljsq(i,j) = noncoeff4(i,j)*noncoeff4(i,j)
            if (noncoeff3(i,j).ne.0.0) then
              r6inv = 1.0/noncoeff3(i,j)**6
              t = noncoeff4(i,j) - noncoeff3(i,j)
              tsq = t*t
              ratio = noncoeff2(i,j)/noncoeff3(i,j)
              ljsw0(i,j) = 4.0*noncoeff1(i,j)*(ratio**12 - ratio**6)
              ljsw1(i,j) = r6inv*(lj1(i,j)*r6inv-lj2(i,j)) /
     $             noncoeff3(i,j)
              ljsw2(i,j) = -r6inv *
     $             (13.0*lj1(i,j)*r6inv - 7.0*lj2(i,j)) /
     $             cutljinnersq(i,j)
              ljsw3(i,j) = -(3.0/tsq) *
     $             (ljsw1(i,j) + 2.0/3.0*ljsw2(i,j)*t)
              ljsw4(i,j) = -1.0/(3.0*tsq) *
     $             (ljsw2(i,j) + 2.0*ljsw3(i,j)*t)
              offset(i,j) = ljsw0(i,j) - ljsw1(i,j)*t -
     $             ljsw2(i,j)*tsq/2.0 - ljsw3(i,j)*tsq*t/3.0 -
     $             ljsw4(i,j)*tsq*tsq/4.0
            else
              ljsw0(i,j) = 0.0
              ljsw1(i,j) = 0.0
              ljsw2(i,j) = 0.0
              ljsw3(i,j) = 0.0
              ljsw4(i,j) = 0.0
              offset(i,j) = 0.0
            endif
            cut = max(cutcoul,noncoeff4(i,j))
          else if (nonstyle.eq.3) then
            lj1(i,j) = 48.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj2(i,j) = 24.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj4(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj5(i,j) = noncoeff3(i,j)
            cutljsq(i,j) = (noncoeff4(i,j)+noncoeff3(i,j))*
     $           (noncoeff4(i,j)+noncoeff3(i,j))
            if (offsetflag.eq.0.or.noncoeff4(i,j).eq.0.0) then
              offset(i,j) = 0.0
            else
              ratio = noncoeff2(i,j)/noncoeff4(i,j)
              offset(i,j) = 4.0*noncoeff1(i,j)*(ratio**12 - ratio**6)
            endif
            cut = max(cutcoul,noncoeff4(i,j)+noncoeff3(i,j))
          else if (nonstyle.eq.4) then
            lj1(i,j) = noncoeff1(i,j)
            lj2(i,j) = noncoeff2(i,j)
            lj3(i,j) = noncoeff3(i,j)
            cutljsq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            cut = noncoeff3(i,j)
          else if (nonstyle.eq.5) then
            lj1(i,j) = 18.0*noncoeff1(i,j)*noncoeff2(i,j)**9
            lj2(i,j) = 18.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 2.0*noncoeff1(i,j)*noncoeff2(i,j)**9
            lj4(i,j) = 3.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            cutljsq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            if (offsetflag.eq.0.or.noncoeff3(i,j).eq.0.0) then
              offset(i,j) = 0.0
            else
              ratio = noncoeff2(i,j)/noncoeff3(i,j)
              offset(i,j) = noncoeff1(i,j) *
     $             (2.0*ratio**9 - 3.0*ratio**6)
            endif
            cut = max(cutcoul,noncoeff3(i,j))
          endif

          cutforcesq(i,j) = cut*cut
          cutneighsq(i,j) = (cut+skin)*(cut+skin)
          cutforce = max(cutforce,cut)

        enddo
      enddo

      cutneigh = cutforce + skin

c zero out all potentials where user-specified cutoff = 0.0
c  setting values to -1.0 means no interaction will ever be computed
c   in force or neighbor routines

      do i = 1,ntypes
        do j = i,ntypes
          if (cutljsq(i,j).eq.0.0) cutljsq(i,j) = -1.0
          if (nonstyle.eq.2.and.cutljinnersq(i,j).eq.0.0)
     $         cutljinnersq(i,j) = -1.0
          if (cutforcesq(i,j).eq.0.0) then
            cutforcesq(i,j) = -1.0
            cutneighsq(i,j) = -1.0
          endif
        enddo
      enddo

      if (cutforce.eq.0.0) cutneigh = 0.0

c symmetrize nonbond j-i interaction to be equivalent to i-j
c  do not need noncoeff values for force/energy routines

      do i = 1,ntypes
        do j = 1,i-1

          if (nonstyle.eq.0) then
            cutljsq(i,j) = cutljsq(j,i)
          else if (nonstyle.eq.1) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            offset(i,j) = offset(j,i)
          else if (nonstyle.eq.2) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            cutljinner(i,j) = cutljinner(j,i)
            cutljinnersq(i,j) = cutljinnersq(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            ljsw0(i,j) = ljsw0(j,i)
            ljsw1(i,j) = ljsw1(j,i)
            ljsw2(i,j) = ljsw2(j,i)
            ljsw3(i,j) = ljsw3(j,i)
            ljsw4(i,j) = ljsw4(j,i)
            offset(i,j) = offset(j,i)
          else if (nonstyle.eq.3) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            lj5(i,j) = lj5(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            offset(i,j) = offset(j,i)
          else if (nonstyle.eq.4) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            cutljsq(i,j) = cutljsq(j,i)
          else if (nonstyle.eq.5) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            offset(i,j) = offset(j,i)
          endif

          cutforcesq(i,j) = cutforcesq(j,i)
          cutneighsq(i,j) = cutneighsq(j,i)

        enddo
      enddo

c check bond lengths against max neighbor cutoff

      if (nbonds.gt.0) then

        if (bondstyle.eq.1) then
          do i = 1,nbondtypes
            if (cutneigh.lt.bondcoeff(2,i)) then
              if (node.eq.0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < average length of bond type',i
            endif
          enddo
        else if (bondstyle.eq.2) then
          do i = 1,nbondtypes
            if (cutneigh.lt.bondcoeff(2,i)) then
              if (node.eq.0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < maximum length of bond type',i
            endif
          enddo
        else if (bondstyle.eq.3) then
          do i = 1,nbondtypes
            if (cutneigh.lt.bondcoeff(2,i)+bondcoeff(5,i)) then
              if (node.eq.0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < maximum length of bond type',i
            endif
          enddo
        else if (bondstyle.eq.4) then
          do i = 1,nbondtypes
            if (cutneigh.lt.bondcoeff(2,i)+bondcoeff(3,i)) then
              if (node.eq.0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < maximum length of bond type',i
            endif
          enddo
        else if (bondstyle.eq.5) then
          do i = 1,nbondtypes
            if (cutneigh.lt.bondcoeff(1,i)) then
              if (node.eq.0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < average length of bond type',i
            endif
          enddo
        endif

      endif

c include dielectric constant in Coulomb prefactor

      coulpre = efactor/dielectric

c change special coeffs for Ewald/PPPM

      if (coulstyle.ge.3) then
        special(1) = special(1) + 2.0
        special(2) = special(2) + 2.0
        special(3) = special(3) + 2.0
      endif

c zero out velocity flags for subsequent runs

      i = atompnt
      do ii = 1,nlocal
        velflag(i) = 0
        i = list(i)
      enddo

c timestep parameters

      dthalf = dt/2.0
      dthalf_intra = dthalf * nstretch
      dthalf_short = dthalf * nintra*nstretch
      dthalf_long = dthalf * nshort*nintra*nstretch

c ensemble type: (1) NVE, (2) NVT, or (3) NPT

      ensemble = 1
      if (tempflag.eq.4) ensemble = 2
      if (pressflag.eq.1) ensemble = 3
      if (xpressflag.eq.1) ensemble = 3
      if (ypressflag.eq.1) ensemble = 3
      if (zpressflag.eq.1) ensemble = 3

      if (ensemble.eq.3.and.tempflag.ne.4)
     $     call error('Must use Nose/Hoover temp control if NPT')

c initialized time counter
c  used by soft potential, fixes, and diagnostics

      itime = 0

c initialize before 1st communication and neighboring and I/O

      max_nlocal = 0
      max_nother = 0
      max_neigh = 0
      max_slist = 0
      max_exch = 0
      max_swap = 0
      max_bond = 0

      numneigh = 0
      time_io = 0.0

c initialize all energies
c  in case force routines never called so thermo routine has good values

      e_long = 0.0
      e_vdwl = 0.0
      e_coul = 0.0

      e_bond = 0.0
      e_angle = 0.0
      e_dihedral = 0.0
      e_improper = 0.0

      e_bondbond = 0.0
      e_bondangle = 0.0
      e_midbondtorsion = 0.0
      e_endbondtorsion = 0.0
      e_angletorsion = 0.0
      e_angleangletorsion = 0.0
      e_bondbond13 = 0.0
      e_angleangle = 0.0

c initialize pointers into fix vectors

      if (nfixes.gt.0) call fix_setup

c setup global box and sub-domains
c  forces setup_comm and setup_neigh to be called 1st time thru
c  if this run follows non-PBC run, insures box is current
c  if PBC, sets PBC test distances with current cutoffs

      call setup_box

c acquire nearby atoms

      call exchange
      call borders

c create special and initial neighbor lists

      call setup_special
      call neighbor

      call neighbor_bond
      call neighbor_angle
      call neighbor_dihedral
      call neighbor_improper

c print run header

      if (node.eq.0) then
        write (6,*)
        write (1,*)
        if (thermostyle.eq.1) then
          write (6,*) 'Step Temp E_nbond E_bond E_long ',
     $         'E_total Pressure Volume'
          write (1,*) 'Step Temp E_nbond E_bond E_long ',
     $         'E_total Pressure Volume'
        endif
      endif

c first force call to initialize velocity Verlet integrator

      call zero_force

      if (nbonds.gt.0) then
        if (bondstyle.eq.1) then
          call bond_harmonic(1)
        else if (bondstyle.eq.2) then
          call bond_fene_standard(1)
        else if (bondstyle.eq.3) then
          call bond_fene_shift(1)
        else if (bondstyle.eq.4) then
          call bond_nonlinear(1)
        else if (bondstyle.eq.5) then
          call bond_class2(1)
        endif
      endif

      if (nrespa.eq.1) then
        if (newton_bond.eq.1) call reverse_comm
        call copy_force(f_stretch,vir_stretch)
        call zero_force
      endif

      if (nangles.gt.0) then
        if (anglestyle.eq.1) then
          call angle_harmonic(1)
        else if (anglestyle.eq.2) then
          call angle_class2(1)
        endif
      endif

      if (ndihedrals.gt.0) then
        if (dihedstyle.eq.1) then
          call dihedral_harmonic(1)
        else if (dihedstyle.eq.2) then
          call dihedral_class2(1)
        endif
      endif

      if (nimpropers.gt.0) then
        if (improstyle.eq.1) then
          call improper_harmonic(1)
        else if (improstyle.eq.2) then
          call improper_cvff(1)
        else if (improstyle.eq.3) then
          call improper_class2(1)
          call angleangle_class2(1)
        endif
      endif
      
      if (nrespa.eq.1) then
        if (newton_bond.eq.1) call reverse_comm
        call copy_force(f_intra,vir_intra)
        call zero_force
      endif

      if (nonstyle.eq.0) then
        if (coulstyle.eq.0) then
          continue
        else if (coulstyle.eq.1) then
          call lj_cut_coul_cut(1)
        else if (coulstyle.eq.2) then
          call lj_cut_coul_smooth(1)
        else if (coulstyle.eq.3) then
          call ewald_coeff
          call lj_cut_coul_long(1)
        else if (coulstyle.eq.4) then
          call pppm_coeff(0)
          call lj_cut_coul_long(1)
        endif
      else if (nonstyle.eq.1) then
        if (coulstyle.eq.0) then
          call lj_cut(1)
        else if (coulstyle.eq.1) then
          call lj_cut_coul_cut(1)
        else if (coulstyle.eq.2) then
          call lj_cut_coul_smooth(1)
        else if (coulstyle.eq.3) then
          call ewald_coeff
          call lj_cut_coul_long(1)
        else if (coulstyle.eq.4) then
          call pppm_coeff(0)
          call lj_cut_coul_long(1)
        endif
      else if (nonstyle.eq.2) then
        if (coulstyle.eq.0) then
          call lj_smooth(1)
        else if (coulstyle.eq.1) then
          call lj_smooth_coul_cut(1)
        else if (coulstyle.eq.2) then
          call lj_smooth_coul_smooth(1)
        else if (coulstyle.eq.3) then
          call ewald_coeff
          call lj_smooth_coul_long(1)
        else if (coulstyle.eq.4) then
          call pppm_coeff(0)
          call lj_smooth_coul_long(1)
        endif
      else if (nonstyle.eq.3) then
        call lj_shift(1)
      else if (nonstyle.eq.4) then
        call soft(1)
      else if (nonstyle.eq.5) then
        if (coulstyle.eq.0) then
          call ljclass2_cut(1)
        else if (coulstyle.eq.1) then
          call ljclass2_cut_coul_cut(1)
        else if (coulstyle.eq.3) then
          call ewald_coeff
          call ljclass2_cut_coul_long(1)
        else if (coulstyle.eq.4) then
          call pppm_coeff(0)
          call ljclass2_cut_coul_long(1)
        endif
      endif

      if (nrespa.eq.1) then
        if (newton_nonbond.eq.1) call reverse_comm
        call copy_force(f_short,vir_short)
        call zero_force
      endif

      if (coulstyle.eq.3) then
        call ewald(1)
      else if (coulstyle.eq.4) then
        time_rho = 0.0
        time_poiss = 0.0
        time_field = 0.0
        call pppm(1)
      endif

      if (nrespa.eq.0) then
        if (newton.ge.1) call reverse_comm
      else if (nrespa.eq.1) then
        call copy_force(f_long,vir_long)
      endif

      if (tempflag.eq.3) call langevin
      if (nfixes.gt.0) call fix_apply

c initial thermodynamics

      ntime_last = ntimestep + nsteps

      call thermo(-1)
      if (nthermo.eq.0) nthermo_next = ntime_last

c initial dumps and restart file

      if (ndumpatom.gt.0) call dump_atom
      if (ndumpatom.eq.0) ndumpatom_next = ntime_last + 1
      if (ndumpvel.gt.0) call dump_vel
      if (ndumpvel.eq.0) ndumpvel_next = ntime_last + 1
      if (ndumpforce.gt.0) call dump_force
      if (ndumpforce.eq.0) ndumpforce_next = ntime_last + 1

      if (nrestart.gt.0)
     $     nrestart_next = min(ntimestep+nrestart,ntime_last)
      if (nrestart.eq.0) nrestart_next = ntime_last + 1

c initial diagnostics
c set all diagnext to this timestep to trigger calling of user routines

      do idiag = 1,numdiag
        diagnext(idiag) = ntimestep
      enddo
      if (numdiag.gt.0) call diagnostic
      if (numdiag.eq.0) ndiag_next = ntime_last + 1

      noutput_next = nthermo_next
      noutput_next = min(noutput_next,ndumpatom_next)
      noutput_next = min(noutput_next,ndumpvel_next)
      noutput_next = min(noutput_next,ndumpforce_next)
      noutput_next = min(noutput_next,nrestart_next)
      noutput_next = min(noutput_next,ndiag_next)
      
c zero out timers and counters

      time_nonbond = 0.0
      time_long = 0.0
      time_rho = 0.0
      time_poiss = 0.0
      time_field = 0.0

      time_bond = 0.0
      time_angle = 0.0
      time_dihedral = 0.0
      time_improper = 0.0

      time_neigh1 = 0.0
      time_neigh2 = 0.0
      time_exch = 0.0

      time_comm = 0.0
      time_fcomm = 0.0

      time_io = 0.0

      numneigh = 0
      ndanger = 0

      return
      end
