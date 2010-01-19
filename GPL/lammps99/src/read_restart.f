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
c read a restart file

      subroutine read_restart
      include "lammps.h"
      include "mpif.h"
      
      call mpi_barrier(mpi_comm_world,ierror)
      if (node.eq.0) write (6,*) 'Reading restart file ...'
      
      readflag = 1

c check restart file version #
c read but ignore nbyte,ibyte (version 5)

      if (node.eq.0) then
        open (unit=2,file=restart_in,status='old',form='unformatted')
        read (2) iversion_tmp,nbyte_tmp,ibyte_tmp
      endif

      call mpi_bcast(iversion_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      if (iversion_tmp.ne.6.and.iversion_tmp.ne.5)
     $     call error('Version mis-match - cannot read file')

c read remainder of restart file header
c read but ignore eta (version 5)

      nbondtypes = 0
      nangletypes = 0
      ndihedtypes = 0
      nimprotypes = 0

      if (node.eq.0) then
        read (2) ntimestep,nprocs_tmp,igridx,igridy,igridz
        read (2) idim,iflagx,iflagy,iflagz
        read (2) iunits,inewton
        read (2) natoms,nbonds,nangles,ndihedrals,nimpropers
        read (2) ntypes
        if (nbonds.gt.0) read (2) nbondtypes
        if (nangles.gt.0) read (2) nangletypes
        if (ndihedrals.gt.0) read (2) ndihedtypes
        if (nimpropers.gt.0) read (2) nimprotypes
        read (2) max_bondper,max_angleper,max_dihedper,max_improper
        read (2) xboundlo,xboundhi,yboundlo,yboundhi,
     $       zboundlo,zboundhi
        read (2) eta,eta_dot,omega_dot(1),omega_dot(2),omega_dot(3)
      endif

c broadcast info

      call mpi_bcast(ntimestep,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nprocs_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(igridx,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(igridy,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(igridz,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(idim,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(iflagx,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(iflagy,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(iflagz,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(iunits,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(inewton,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(natoms,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nbonds,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nangles,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(ndihedrals,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nimpropers,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(ntypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nbondtypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nangletypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(ndihedtypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nimprotypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(max_bondper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(max_angleper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(max_dihedper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(max_improper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(xboundlo,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(xboundhi,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(yboundlo,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(yboundhi,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(zboundlo,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(zboundhi,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(eta_dot,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(omega_dot,3,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

c error check
      
      if (nprocs_tmp.ne.nprocs) then
        if (node.eq.0) write (6,*) 'WARNING: Restarting on ',
     $       'different number of processors'
      endif

      if (idim.ne.idimension) then
        if (node.eq.0) write (6,*) 'WARNING: Old dimension ',
     $       'does not match current dimension'
      endif

      if (iflagx.ne.perflagx.or.iflagy.ne.perflagy.or.
     $     iflagz.ne.perflagz) then
        if (node.eq.0) write (6,*) 'WARNING: Old periodicity ',
     $       'does not match current periodicity'
      endif

      if (iunits.ne.units) then
        if (node.eq.0) write (6,*) 'WARNING: Old units ',
     $       'does not match current units'
      endif

      if (inewton.ne.newton) then
        if (node.eq.0) write (6,*) 'WARNING: Old newton flag ',
     $       'does not match current newton flag'
      endif

      if (natoms.gt.maxtotal)
     $     call error('Too many atoms - boost maxtotal')
      if (ntypes.gt.maxtype)
     $     call error('Too many atom types - boost maxtype')
      if (nbondtypes.gt.maxbondtype)
     $     call error('Too many bond types - boost maxtype')
      if (nangletypes.gt.maxangletype)
     $     call error('Too many angle types - boost maxangletype')
      if (ndihedtypes.gt.maxdihedtype)
     $     call error('Too many dihedral types - boost maxdihedtype')
      if (nimprotypes.gt.maximprotype)
     $     call error('Too many improper types - boost maximprotype')

      if (max_bondper.gt.maxbondper)
     $     call error('Too many bonds/atom - boost maxbondper')
      if (max_angleper.gt.maxangleper)
     $     call error('Too many angles/atom - boost maxangleper')
      if (max_dihedper.gt.maxdihedper)
     $     call error('Too many dihedrals/atom - boost maxdihedper')
      if (max_improper.gt.maximproper)
     $     call error('Too many impropers/atom - boost maximproper')

      if (idimension.eq.2.and.pgrid(1).gt.0) then
        if (pgrid(3).ne.1)
     $       call error('Can only use one z-processor for 2-d')
      endif

c define box
c for non-PBC, assumes bounds are already augmented by "small" when
c  restart file was written

      xprd = xboundhi - xboundlo
      yprd = yboundhi - yboundlo
      zprd = zboundhi - zboundlo

c define mesh topology and error check against restart file values

      call mesh_3d(node,nprocs,idimension,xprd,yprd,zprd,
     $     pgrid,me,mpart)

      if (node.eq.0) then
        write (6,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
        write (1,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
      endif

      if (igridx.ne.pgrid(1).or.igridy.ne.pgrid(2).or.
     $     igridz.ne.pgrid(3)) then
        if (node.eq.0) write (6,*) 'WARNING: Restarting on ',
     $       'different 3-d grid of processors'
      endif

c define sub-domain boundaries

      border(1,1) = xboundlo + float(me(1))/pgrid(1) * xprd
      border(2,1) = xboundlo + float(me(1)+1)/pgrid(1) * xprd
      border(1,2) = yboundlo + float(me(2))/pgrid(2) * yprd
      border(2,2) = yboundlo + float(me(2)+1)/pgrid(2) * yprd
      border(1,3) = zboundlo + float(me(3))/pgrid(3) * zprd
      border(2,3) = zboundlo + float(me(3)+1)/pgrid(3) * zprd

c initialize my atom stack

      freepnt = 1
      do i = 1,maxown
        list(i) = i + 1
      enddo
      list(maxown) = 0
      nlocal = 0
      nother = 0
      atompnt = maxown + 1
      
c initialize global -> local pointers

      do i = 1,natoms
        localptr(i) = 0
      enddo

c read all coefficients

      if (node.eq.0) then

        write (6,*) '  Coefficients ...'

        read (2) (mass(i),i=1,ntypes)

        read (2) nonstyle
        if (nonstyle.eq.0) then
          write (6,*) 'Nonbond style set to none'
        else if (nonstyle.eq.1) then
          write (6,*) 'Nonbond style set to lj/cutoff'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.2) then
          write (6,*) 'Nonbond style set to lj/smooth'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.3) then
          write (6,*) 'Nonbond style set to lj/shift'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.4) then
          write (6,*) 'Nonbond style set to soft'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.5) then
          write (6,*) 'Nonbond style set to class2/cutoff'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        endif

        if (nbonds.gt.0) then
          read (2) bondstyle
          if (bondstyle.eq.0) then
            write (6,*) 'Bond style set to none'
          else if (bondstyle.eq.1) then
            write (6,*) 'Bond style set to harmonic'
            read (2) ((bondcoeff(i,j),i=1,2),j=1,nbondtypes)
          else if (bondstyle.eq.2) then
            write (6,*) 'Bond style set to fene/standard'
            read (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          else if (bondstyle.eq.3) then
            write (6,*) 'Bond style set to fene/shift'
            read (2) ((bondcoeff(i,j),i=1,5),j=1,nbondtypes)
          else if (bondstyle.eq.4) then
            write (6,*) 'Bond style set to nonlinear'
            read (2) ((bondcoeff(i,j),i=1,3),j=1,nbondtypes)
          else if (bondstyle.eq.5) then
            write (6,*) 'Bond style set to class2'
            read (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          endif
        endif

        if (nangles.gt.0) then
          read (2) anglestyle
          if (anglestyle.eq.0) then
            write (6,*) 'Angle style set to none'
          else if (anglestyle.eq.1) then
            write (6,*) 'Angle style set to harmonic'
            read (2) ((anglecoeff(i,j),i=1,2),j=1,nangletypes)
          else if (anglestyle.eq.2) then
            write (6,*) 'Angle style set to class2'
            read (2) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
            read (2) ((bondbondcoeff(i,j),i=1,3),j=1,nangletypes)
            read (2) ((bondanglecoeff(i,j),i=1,4),j=1,nangletypes)
          endif
        endif

        if (ndihedrals.gt.0) then
          read (2) dihedstyle
          if (dihedstyle.eq.0) then
            write (6,*) 'Dihedral style set to none'
          else if (dihedstyle.eq.1) then
            write (6,*) 'Dihedral style set to harmonic'
            read (2) ((dihedcoeff(i,j),i=1,3),j=1,ndihedtypes)
          else if (dihedstyle.eq.2) then
            write (6,*) 'Dihedral style set to class2'
            read (2) ((dihedcoeff(i,j),i=1,6),j=1,ndihedtypes)
            read (2) ((midbondtorsioncoeff(i,j),i=1,4),
     $           j=1,ndihedtypes)
            read (2) ((endbondtorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            read (2) ((angletorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            read (2) ((angleangletorsioncoeff(i,j),i=1,3),
     $           j=1,ndihedtypes)
          endif
        endif

        if (nimpropers.gt.0) then
          read (2) improstyle
          if (improstyle.eq.0) then
            write (6,*) 'Improper style set to none'
          else if (improstyle.eq.1) then
            write (6,*) 'Improper style set to harmonic'
            read (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
          else if (improstyle.eq.2) then
            write (6,*) 'Improper style set to cvff'
            read (2) ((improcoeff(i,j),i=1,3),j=1,nimprotypes)
          else if (improstyle.eq.3) then
            write (6,*) 'Improper style set to class2'
            read (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
            read (2) ((angleanglecoeff(i,j),i=1,6),j=1,nimprotypes)
          endif
        endif

      endif

c broadcast all coefficients

      call mpi_bcast(mass,ntypes,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(nonstyle,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      if (nonstyle.eq.1) then
        call mpi_bcast(noncoeff1,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.2) then
        call mpi_bcast(noncoeff1,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff4,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.3) then
        call mpi_bcast(noncoeff1,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff4,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.4) then
        call mpi_bcast(noncoeff1,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.5) then
        call mpi_bcast(noncoeff1,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*maxtype,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      endif

      do i = 1,ntypes
        do j = i,ntypes
          nontypeflag(i,j) = nonstyle
        enddo
      enddo

      if (nbonds.gt.0) then
        call mpi_bcast(bondstyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (bondstyle.gt.0)
     $       call mpi_bcast(bondcoeff,5*nbondtypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        do i = 1,nbondtypes
          bondtypeflag(i) = bondstyle
        enddo
      endif

      if (nangles.gt.0) then
        call mpi_bcast(anglestyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (anglestyle.gt.0)
     $       call mpi_bcast(anglecoeff,4*nangletypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        if (anglestyle.eq.2) then
          call mpi_bcast(bondbondcoeff,3*nangletypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(bondanglecoeff,4*nangletypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
        endif
        do i = 1,nangletypes
          angletypeflag(i) = anglestyle
        enddo
      endif

      if (ndihedrals.gt.0) then
        call mpi_bcast(dihedstyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (dihedstyle.gt.0)
     $       call mpi_bcast(dihedcoeff,6*ndihedtypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        if (dihedstyle.eq.2) then
          call mpi_bcast(midbondtorsioncoeff,4*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(endbondtorsioncoeff,8*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(angletorsioncoeff,8*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(angleangletorsioncoeff,3*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
        endif
        do i = 1,ndihedtypes
          dihedtypeflag(i) = dihedstyle
        enddo
      endif

      if (nimpropers.gt.0) then
        call mpi_bcast(improstyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (improstyle.gt.0)
     $       call mpi_bcast(improcoeff,3*nimprotypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        if (improstyle.eq.3) then
          call mpi_bcast(angleanglecoeff,6*nimprotypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
        endif
        do i = 1,nimprotypes
          improtypeflag(i) = improstyle
        enddo
      endif

c read and broadcast atoms
c  add atom to my list if in my sub-box
c  do not check for PBC or remap since assume restart file is "correct"
c  only read 1000 at a time to avoid comm buffer overflow

      if (node.eq.0) write (6,*) '  Reading',natoms,' atoms ...'

      do i = 1,natoms

        if (node.eq.0) then
          read (2) (buf5(k),k=1,peratom_io)
          j = peratom_io
          numbond_tmp = nint(buf5(j-3))
          numangle_tmp = nint(buf5(j-2))
          numdihed_tmp = nint(buf5(j-1))
          numimpro_tmp = nint(buf5(j))
          do jj = 1,numbond_tmp
            read (2) buf5(j+1),buf5(j+2),buf5(j+3)
            j = j + 3
          enddo
          do jj = 1,numangle_tmp
            read (2) buf5(j+1),buf5(j+2),buf5(j+3),buf5(j+4)
            j = j + 4
          enddo
          do jj = 1,numdihed_tmp
            read (2) buf5(j+1),buf5(j+2),buf5(j+3),
     $           buf5(j+4),buf5(j+5)
            j = j + 5
          enddo
          do jj = 1,numimpro_tmp
            read (2) buf5(j+1),buf5(j+2),buf5(j+3),
     $           buf5(j+4),buf5(j+5)
            j = j + 5
          enddo
        endif
        call mpi_bcast(numbond_tmp,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        call mpi_bcast(numangle_tmp,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        call mpi_bcast(numdihed_tmp,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        call mpi_bcast(numimpro_tmp,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        isize = (peratom_io + 3*numbond_tmp + 4*numangle_tmp +
     $       5*numdihed_tmp + 5*numimpro_tmp)
        call mpi_bcast(buf5,isize,mpi_double_precision,0,
     $       mpi_comm_world,ierror)
        
        if (buf5(1).ge.border(1,1).and.buf5(1).lt.border(2,1).and.
     $       buf5(2).ge.border(1,2).and.buf5(2).lt.border(2,2).and.
     $       buf5(3).ge.border(1,3).and.buf5(3).lt.border(2,3)) then
          nlocal = nlocal + 1
          if (freepnt.ne.0) then
            jtmp = atompnt
            atompnt = freepnt
            freepnt = list(freepnt)
            list(atompnt) = jtmp
            x(1,atompnt) = buf5(1)
            x(2,atompnt) = buf5(2)
            x(3,atompnt) = buf5(3)
            itag = nint(buf5(4))
            localptr(itag) = atompnt
            tag(atompnt) = itag
            molecule(atompnt) = nint(buf5(5))
            type(atompnt) = nint(buf5(6))
            q(atompnt) = buf5(7)
            v(1,atompnt) = buf5(8)
            v(2,atompnt) = buf5(9)
            v(3,atompnt) = buf5(10)
            true(atompnt) = nint(buf5(11))
            fix(atompnt) = nint(buf5(12))
            numbond(atompnt) = nint(buf5(13))
            numangle(atompnt) = nint(buf5(14))
            numdihed(atompnt) = nint(buf5(15))
            numimpro(atompnt) = nint(buf5(16))
            velflag(atompnt) = 1
            j = peratom_io
            do jj = 1,numbond(atompnt)
              bondatom1(jj,atompnt) = nint(buf5(j+1))
              bondatom2(jj,atompnt) = nint(buf5(j+2))
              bondtype(jj,atompnt) = nint(buf5(j+3))
              j = j + 3
            enddo
            do jj = 1,numangle(atompnt)
              angleatom1(jj,atompnt) = nint(buf5(j+1))
              angleatom2(jj,atompnt) = nint(buf5(j+2))
              angleatom3(jj,atompnt) = nint(buf5(j+3))
              angletype(jj,atompnt) = nint(buf5(j+4))
              j = j + 4
            enddo
            do jj = 1,numdihed(atompnt)
              dihedatom1(jj,atompnt) = nint(buf5(j+1))
              dihedatom2(jj,atompnt) = nint(buf5(j+2))
              dihedatom3(jj,atompnt) = nint(buf5(j+3))
              dihedatom4(jj,atompnt) = nint(buf5(j+4))
              dihedtype(jj,atompnt) = nint(buf5(j+5))
              j = j + 5
            enddo
            do jj = 1,numimpro(atompnt)
              improatom1(jj,atompnt) = nint(buf5(j+1))
              improatom2(jj,atompnt) = nint(buf5(j+2))
              improatom3(jj,atompnt) = nint(buf5(j+3))
              improatom4(jj,atompnt) = nint(buf5(j+4))
              improtype(jj,atompnt) = nint(buf5(j+5))
              j = j + 5
            enddo
          endif
        endif

        if (mod(i,1000).eq.0) then
          call mpi_barrier(mpi_comm_world,ierror)
          if (node.eq.0) write (6,*) i
        endif

      enddo

c check if too many atoms on any proc

      call mpi_allreduce(nlocal,iflagmax,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (iflagmax.gt.maxown) then
        if (node.eq.0) write (6,*) 'Nlocal,maxown =',iflagmax,maxown
        call error('Too many atoms/processor - boost maxown')
      endif

c check that all atoms were assigned uniquely to a processor

      call mpi_allreduce(nlocal,natoms_tmp,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (natoms.ne.natoms_tmp)
     $     call error('Did not read all atoms correctly')

c close restart file

      if (node.eq.0) close (2)

c compute total mass and charge

      masssum = 0.0
      qsum = 0.0
      qsqsum = 0.0
      i = atompnt
      do ii = 1,nlocal
        masssum = masssum + mass(type(i))
        qsum = qsum + q(i)
        qsqsum = qsqsum + q(i)*q(i)
        i = list(i)
      enddo

      tmp = masssum
      call mpi_allreduce(tmp,masssum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      tmp = qsum
      call mpi_allreduce(tmp,qsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      tmp = qsqsum
      call mpi_allreduce(tmp,qsqsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)

      call mpi_barrier(mpi_comm_world,ierror)

      return
      end
