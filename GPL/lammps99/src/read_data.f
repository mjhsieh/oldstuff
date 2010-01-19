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

      subroutine read_data
      include "lammps.h"
      include "mpif.h"

      parameter (small=1.0E-4)
      logical match
      real*8 tmp(10)
      integer itmp(5)
      character*80 str

 900  format (a)

      call mpi_barrier(mpi_comm_world,ierror)
      if (node.eq.0) write (6,*) 'Reading data file ...'

      readflag = 1
      ntimestep = 0

      nbondtypes = 0
      nangletypes = 0
      ndihedtypes = 0
      nimprotypes = 0

c read data file header - fixed format

      if (node.eq.0) then
        open (unit=2,file=datafile,status='old')
        read (2,*)
        read (2,*)
        read (2,*) natoms
        read (2,*) nbonds
        read (2,*) nangles
        read (2,*) ndihedrals
        read (2,*) nimpropers
        read (2,*)
        read (2,*) ntypes
        if (nbonds.gt.0) read (2,*) nbondtypes
        if (nangles.gt.0) read (2,*) nangletypes
        if (ndihedrals.gt.0) read (2,*) ndihedtypes
        if (nimpropers.gt.0) read (2,*) nimprotypes
        read (2,*)
        read (2,*) xboundlo,xboundhi
        read (2,*) yboundlo,yboundhi
        if (idimension.eq.3) read (2,*) zboundlo,zboundhi
      endif

c broadcast values

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
      call mpi_bcast(xboundlo,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(xboundhi,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(yboundlo,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(yboundhi,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      if (idimension.eq.2) then
        zboundlo = -0.1
        zboundhi = 0.1
      else
        call mpi_bcast(zboundlo,1,mpi_double_precision,0,
     $       mpi_comm_world,ierror)
        call mpi_bcast(zboundhi,1,mpi_double_precision,0,
     $       mpi_comm_world,ierror)
      endif

c initialize ensemble variables

      eta_dot = 0.0
      omega_dot(1) = 0.0
      omega_dot(2) = 0.0
      omega_dot(3) = 0.0

c error check

      if (natoms.le.0) call error('Atoms <= 0')
      if (natoms.gt.maxtotal)
     $     call error('Too many atoms - boost maxtotal')

      if (nbonds.lt.0) call error('Bonds < 0')
      if (nangles.lt.0) call error('Angles < 0')
      if (ndihedrals.lt.0) call error('Dihedrals < 0')
      if (nimpropers.lt.0) call error('Impropers < 0')

      if (ntypes.le.0) call error('Atom types <= 0')
      if (ntypes.gt.maxtype)
     $     call error('Too many atom types - boost maxtype')

      if (nbonds.gt.0.and.nbondtypes.le.0)
     $     call error('Bond types <= 0')
      if (nbondtypes.gt.maxbondtype)
     $     call error('Too many bond types - boost maxbondtype')

      if (nangles.gt.0.and.nangletypes.le.0)
     $     call error('Angle types <= 0')
      if (nangletypes.gt.maxangletype)
     $     call error('Too many angle types - boost maxangletype')

      if (ndihedrals.gt.0.and.ndihedtypes.le.0)
     $     call error('Dihedral types <= 0')
      if (ndihedtypes.gt.maxdihedtype)
     $     call error('Too many dihedral types - boost maxdihedtype')

      if (nimpropers.gt.0.and.nimprotypes.le.0)
     $     call error('Improper types <= 0')
      if (nimprotypes.gt.maximprotype)
     $     call error('Too many improper types - boost maximprotype')

      if (xboundlo.ge.xboundhi.or.yboundlo.ge.yboundhi.or.
     $     zboundlo.ge.zboundhi) call error('Illegal box dimensions')

      if (idimension.eq.2.and.pgrid(1).gt.0) then
        if (pgrid(3).ne.1)
     $       call error('Can only use one z-processor for 2-d')
      endif

c if non-PBC add small to bounds
c  assumes input bounds are exact
c  must add small to insure every atom is inside some procs sub-domain

      if (perflagx.eq.1) then
        xboundlo = xboundlo - small
        xboundhi = xboundhi + small
      endif
      if (perflagy.eq.1) then
        yboundlo = yboundlo - small
        yboundhi = yboundhi + small
      endif
      if (perflagz.eq.1) then
        zboundlo = zboundlo - small
        zboundhi = zboundhi + small
      endif

c define box

      xprd = xboundhi - xboundlo
      yprd = yboundhi - yboundlo
      zprd = zboundhi - zboundlo

c define mesh topology

      call mesh_3d(node,nprocs,idimension,xprd,yprd,zprd,
     $     pgrid,me,mpart)

      if (node.eq.0) then
        write (6,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
        write (1,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
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

c read remainder of data file - free format

      massflag = 0
      natomflag = 0
      nbondflag = 0
      nangleflag = 0
      ndihedflag = 0
      nimproflag = 0

      noncoeffflag = 0
      nbondcoeffflag = 0
      nanglecoeffflag = 0
      ndihedcoeffflag = 0
      nimprocoeffflag = 0

      nbondbondflag = 0
      nbondangleflag = 0
      nmidbondtorsionflag = 0
      nendbondtorsionflag = 0
      nangletorsionflag = 0
      nangleangletorsionflag = 0
      nbondbond13flag = 0
      nangleangleflag = 0

      do infinity = 1,100000000

c read identifier string

        if (node.eq.0) then
          read (2,*,end=999,err=999)
          read (2,900,end=999,err=999) str
          read (2,*,end=999,err=999)
        endif
        call mpi_bcast(str,80,mpi_character,0,mpi_comm_world,ierror)

c end of file for all procs except 0

        if (match('All Done',str,m)) then

          goto 999

c read and broadcast atom masses

        else if (match('Masses',str,m)) then

          if (node.eq.0) write (6,*) '  Masses ...'
          massflag = 1
          if (node.eq.0) then
            do i = 1,ntypes
              read (2,*) jtmp,mass(i)
            enddo
          endif

          call mpi_bcast(mass,ntypes,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

c read and broadcast atoms
c  add atom to my list if in my sub-box
c  put them in periodic box no matter how far away they are
c  only read 1000 at a time to avoid comm buffer overflow

        else if (match('Atoms',str,m)) then

          if (node.eq.0) write (6,*) '  Atoms ...'
          natomflag = 1

          inum = 7
          if (trueflag.gt.1) inum = 10
          isize = inum

          do i = 1,natoms

            if (node.eq.0) read (2,*) (tmp(j),j=1,inum)
            call mpi_bcast(tmp,isize,mpi_double_precision,0,
     $           mpi_comm_world,ierror)

            if (trueflag.le.1) then
              itrue = 500500500
            else
              itrue = 500+nint(tmp(8)) + (500+nint(tmp(9)))*1000 +
     $             (500+nint(tmp(10)))*1000000
            endif

            call remap(tmp(5),tmp(6),tmp(7),itrue)

            if (tmp(5).ge.border(1,1).and.tmp(5).lt.border(2,1).and.
     $           tmp(6).ge.border(1,2).and.tmp(6).lt.border(2,2).and.
     $           tmp(7).ge.border(1,3).and.tmp(7).lt.border(2,3)) then
              nlocal = nlocal + 1
              if (freepnt.ne.0) then
                jtmp = atompnt
                atompnt = freepnt
                freepnt = list(freepnt)
                list(atompnt) = jtmp
                itag = nint(tmp(1))
                localptr(itag) = atompnt
                tag(atompnt) = itag
                molecule(atompnt) = nint(tmp(2))
                type(atompnt) = nint(tmp(3))
                q(atompnt) = tmp(4)
                x(1,atompnt) = tmp(5)
                x(2,atompnt) = tmp(6)
                x(3,atompnt) = tmp(7)
                true(atompnt) = itrue
                fix(atompnt) = 0
                numbond(atompnt) = 0
                numangle(atompnt) = 0
                numdihed(atompnt) = 0
                numimpro(atompnt) = 0
                velflag(atompnt) = 0
                v(1,atompnt) = 0.0
                v(2,atompnt) = 0.0
                v(3,atompnt) = 0.0
              endif
            endif

            if (mod(i,1000).eq.0) then
              call mpi_barrier(mpi_comm_world,ierror)
              if (node.eq.0) write (6,*) i
            endif

          enddo

c check if too many atoms on any proc

          call mpi_allreduce(nlocal,iflagmax,1,mpi_integer,mpi_max,
     $         mpi_comm_world,ierror)
          if (iflagmax.gt.maxown) then
            if (node.eq.0)
     $           write (6,*) 'Nlocal,maxown =',iflagmax,maxown
            call error('Too many atoms/processor - boost maxown')
          endif

c check that all atoms were assigned uniquely to a processor

          call mpi_allreduce(nlocal,natoms_tmp,1,mpi_integer,
     $         mpi_sum,mpi_comm_world,ierror)
          if (natoms.ne.natoms_tmp)
     $         call error('Did not read all atoms correctly')

c read and broadcast bonds
c  only read 1000 at a time to avoid comm buffer overflow

        else if (match('Bonds',str,m)) then

          if (node.eq.0) write (6,*) '  Bonds ...'
          if (nbonds.eq.0.or.nbondtypes.eq.0)
     $         call error('Should be no Bonds entry')

          nbondflag = 1
          isize = 3
          do i = 1,nbonds

            if (node.eq.0) read (2,*) jtmp,itmp(1),itmp(2),itmp(3)
            call mpi_bcast(itmp,isize,mpi_integer,0,
     $           mpi_comm_world,ierror)

            ipnt = localptr(itmp(2))
            if (ipnt.ne.0) then
              numbond(ipnt) = numbond(ipnt) + 1
              if (numbond(ipnt).le.maxbondper) then
                bondtype(numbond(ipnt),ipnt) = itmp(1)
                bondatom1(numbond(ipnt),ipnt) = itmp(2)
                bondatom2(numbond(ipnt),ipnt) = itmp(3)
              endif
            endif

            if (newton_bond.eq.0) then
              ipnt = localptr(itmp(3))
              if (ipnt.ne.0) then
                numbond(ipnt) = numbond(ipnt) + 1
                if (numbond(ipnt).le.maxbondper) then
                  bondtype(numbond(ipnt),ipnt) = itmp(1)
                  bondatom1(numbond(ipnt),ipnt) = itmp(2)
                  bondatom2(numbond(ipnt),ipnt) = itmp(3)
                endif
              endif
            endif

            if (mod(i,1000).eq.0) then
              call mpi_barrier(mpi_comm_world,ierror)
              if (node.eq.0) write (6,*) i
            endif

          enddo

c check for too many bonds/atom

          max_bondper = 0
          i = atompnt
          do ii = 1,nlocal
            if (numbond(i).gt.max_bondper) max_bondper = numbond(i)
            i = list(i)
          enddo

          jtmp = max_bondper
          call mpi_allreduce(jtmp,max_bondper,1,mpi_integer,
     $         mpi_max,mpi_comm_world,ierror)
          if (max_bondper.gt.maxbondper) then
            if (node.eq.0) write (6,*) 'Max bonds/atom =',max_bondper
            call error('Too many bonds/atom - boost maxbondper')
          endif

c check that all bonds were assigned uniquely to a processor

          nbonds_tmp = 0
          i = atompnt
          do ii = 1,nlocal
            nbonds_tmp = nbonds_tmp + numbond(i)
            i = list(i)
          enddo
          
          jtmp = nbonds_tmp
          call mpi_allreduce(jtmp,nbonds_tmp,1,mpi_integer,
     $         mpi_sum,mpi_comm_world,ierror)
          ifactor = 1
          if (newton_bond.eq.0) ifactor = 2
          if (nbonds_tmp.ne.ifactor*nbonds)
     $         call error('Bonds assigned incorrectly')

c read and broadcast angles
c  only read 1000 at a time to avoid comm buffer overflow

        else if (match('Angles',str,m)) then

          if (node.eq.0) write (6,*) '  Angles ...'
          if (nangles.eq.0.or.nangletypes.eq.0)
     $         call error('Should be no Angles entry')

          nangleflag = 1
          isize = 4
          do i = 1,nangles

            if (node.eq.0) read (2,*) jtmp,
     $           itmp(1),itmp(2),itmp(3),itmp(4)
            call mpi_bcast(itmp,isize,mpi_integer,0,
     $           mpi_comm_world,ierror)

            ipnt = localptr(itmp(3))
            if (ipnt.ne.0) then
              numangle(ipnt) = numangle(ipnt) + 1
              if (numangle(ipnt).le.maxangleper) then
                angletype(numangle(ipnt),ipnt) = itmp(1)
                angleatom1(numangle(ipnt),ipnt) = itmp(2)
                angleatom2(numangle(ipnt),ipnt) = itmp(3)
                angleatom3(numangle(ipnt),ipnt) = itmp(4)
              endif
            endif

            if (newton_bond.eq.0) then
              ipnt = localptr(itmp(2))
              if (ipnt.ne.0) then
                numangle(ipnt) = numangle(ipnt) + 1
                if (numangle(ipnt).le.maxangleper) then
                  angletype(numangle(ipnt),ipnt) = itmp(1)
                  angleatom1(numangle(ipnt),ipnt) = itmp(2)
                  angleatom2(numangle(ipnt),ipnt) = itmp(3)
                  angleatom3(numangle(ipnt),ipnt) = itmp(4)
                endif
              endif
              
              ipnt = localptr(itmp(4))
              if (ipnt.ne.0) then
                numangle(ipnt) = numangle(ipnt) + 1
                if (numangle(ipnt).le.maxangleper) then
                  angletype(numangle(ipnt),ipnt) = itmp(1)
                  angleatom1(numangle(ipnt),ipnt) = itmp(2)
                  angleatom2(numangle(ipnt),ipnt) = itmp(3)
                  angleatom3(numangle(ipnt),ipnt) = itmp(4)
                endif
              endif
            endif

            if (mod(i,1000).eq.0) then
              call mpi_barrier(mpi_comm_world,ierror)
              if (node.eq.0) write (6,*) i
            endif

          enddo

c check for too many angles/atom

          max_angleper = 0
          i = atompnt
          do ii = 1,nlocal
            if (numangle(i).gt.max_angleper) max_angleper = numangle(i)
            i = list(i)
          enddo

          jtmp = max_angleper
          call mpi_allreduce(jtmp,max_angleper,1,mpi_integer,
     $         mpi_max,mpi_comm_world,ierror)
          if (max_angleper.gt.maxangleper) then
            if (node.eq.0) write (6,*) 'Max angles/atom =',max_angleper
            call error('Too many angles/atom - boost maxangleper')
          endif

c check that all angles were assigned uniquely to a processor

          nangles_tmp = 0
          i = atompnt
          do ii = 1,nlocal
            nangles_tmp = nangles_tmp + numangle(i)
            i = list(i)
          enddo
          
          jtmp = nangles_tmp
          call mpi_allreduce(jtmp,nangles_tmp,1,mpi_integer,
     $         mpi_sum,mpi_comm_world,ierror)
          ifactor = 1
          if (newton_bond.eq.0) ifactor = 3
          if (nangles_tmp.ne.ifactor*nangles)
     $         call error('Angles assigned incorrectly')

c read and broadcast dihedrals
c  only read 1000 at a time to avoid comm buffer overflow

        else if (match('Dihedrals',str,m)) then

          if (node.eq.0) write (6,*) '  Dihedrals ...'
          if (ndihedrals.eq.0.or.ndihedtypes.eq.0)
     $         call error('Should be no Dihedrals entry')

          ndihedflag = 1
          isize = 5
          do i = 1,ndihedrals

            if (node.eq.0) read (2,*) jtmp,
     $           itmp(1),itmp(2),itmp(3),itmp(4),itmp(5)
            call mpi_bcast(itmp,isize,mpi_integer,0,
     $           mpi_comm_world,ierror)

            ipnt = localptr(itmp(3))
            if (ipnt.ne.0) then
              numdihed(ipnt) = numdihed(ipnt) + 1
              if (numdihed(ipnt).le.maxdihedper) then
                dihedtype(numdihed(ipnt),ipnt) = itmp(1)
                dihedatom1(numdihed(ipnt),ipnt) = itmp(2)
                dihedatom2(numdihed(ipnt),ipnt) = itmp(3)
                dihedatom3(numdihed(ipnt),ipnt) = itmp(4)
                dihedatom4(numdihed(ipnt),ipnt) = itmp(5)
              endif
            endif

            if (newton_bond.eq.0) then
              ipnt = localptr(itmp(2))
              if (ipnt.ne.0) then
                numdihed(ipnt) = numdihed(ipnt) + 1
                if (numdihed(ipnt).le.maxdihedper) then
                  dihedtype(numdihed(ipnt),ipnt) = itmp(1)
                  dihedatom1(numdihed(ipnt),ipnt) = itmp(2)
                  dihedatom2(numdihed(ipnt),ipnt) = itmp(3)
                  dihedatom3(numdihed(ipnt),ipnt) = itmp(4)
                  dihedatom4(numdihed(ipnt),ipnt) = itmp(5)
                endif
              endif

              ipnt = localptr(itmp(4))
              if (ipnt.ne.0) then
                numdihed(ipnt) = numdihed(ipnt) + 1
                if (numdihed(ipnt).le.maxdihedper) then
                  dihedtype(numdihed(ipnt),ipnt) = itmp(1)
                  dihedatom1(numdihed(ipnt),ipnt) = itmp(2)
                  dihedatom2(numdihed(ipnt),ipnt) = itmp(3)
                  dihedatom3(numdihed(ipnt),ipnt) = itmp(4)
                  dihedatom4(numdihed(ipnt),ipnt) = itmp(5)
                endif
              endif
              
              ipnt = localptr(itmp(5))
              if (ipnt.ne.0) then
                numdihed(ipnt) = numdihed(ipnt) + 1
                if (numdihed(ipnt).le.maxdihedper) then
                  dihedtype(numdihed(ipnt),ipnt) = itmp(1)
                  dihedatom1(numdihed(ipnt),ipnt) = itmp(2)
                  dihedatom2(numdihed(ipnt),ipnt) = itmp(3)
                  dihedatom3(numdihed(ipnt),ipnt) = itmp(4)
                  dihedatom4(numdihed(ipnt),ipnt) = itmp(5)
                endif
              endif
            endif

            if (mod(i,1000).eq.0) then
              call mpi_barrier(mpi_comm_world,ierror)
              if (node.eq.0) write (6,*) i
            endif

          enddo

c check for too many dihedrals/atom

          max_dihedper = 0
          i = atompnt
          do ii = 1,nlocal
            if (numdihed(i).gt.max_dihedper) max_dihedper = numdihed(i)
            i = list(i)
          enddo

          jtmp = max_dihedper
          call mpi_allreduce(jtmp,max_dihedper,1,mpi_integer,
     $         mpi_max,mpi_comm_world,ierror)
          if (max_dihedper.gt.maxdihedper) then
            if (node.eq.0) write (6,*) 'Max diheds/atom =',max_dihedper
            call error('Too many diheds/atom - boost maxdihedper')
          endif

c check that all dihedrals were assigned uniquely to a processor

          ndiheds_tmp = 0
          i = atompnt
          do ii = 1,nlocal
            ndiheds_tmp = ndiheds_tmp + numdihed(i)
            i = list(i)
          enddo
          
          jtmp = ndiheds_tmp
          call mpi_allreduce(jtmp,ndiheds_tmp,1,mpi_integer,
     $         mpi_sum,mpi_comm_world,ierror)
          ifactor = 1
          if (newton_bond.eq.0) ifactor = 4
          if (ndiheds_tmp.ne.ifactor*ndihedrals)
     $         call error('Dihedrals assigned incorrectly')

c read and broadcast impropers
c  only read 1000 at a time to avoid comm buffer overflow

        else if (match('Impropers',str,m)) then

          if (node.eq.0) write (6,*) '  Impropers ...'
          if (nimpropers.eq.0.or.nimprotypes.eq.0)
     $         call error('Should be no Impropers entry')

          nimproflag = 1
          isize = 5
          do i = 1,nimpropers

            if (node.eq.0) read (2,*) jtmp,
     $           itmp(1),itmp(2),itmp(3),itmp(4),itmp(5)
            call mpi_bcast(itmp,isize,mpi_integer,0,
     $           mpi_comm_world,ierror)

            ipnt = localptr(itmp(3))
            if (ipnt.ne.0) then
              numimpro(ipnt) = numimpro(ipnt) + 1
              if (numimpro(ipnt).le.maximproper) then
                improtype(numimpro(ipnt),ipnt) = itmp(1)
                improatom1(numimpro(ipnt),ipnt) = itmp(2)
                improatom2(numimpro(ipnt),ipnt) = itmp(3)
                improatom3(numimpro(ipnt),ipnt) = itmp(4)
                improatom4(numimpro(ipnt),ipnt) = itmp(5)
              endif
            endif

            if (newton_bond.eq.0) then
              ipnt = localptr(itmp(2))
              if (ipnt.ne.0) then
                numimpro(ipnt) = numimpro(ipnt) + 1
                if (numimpro(ipnt).le.maximproper) then
                  improtype(numimpro(ipnt),ipnt) = itmp(1)
                  improatom1(numimpro(ipnt),ipnt) = itmp(2)
                  improatom2(numimpro(ipnt),ipnt) = itmp(3)
                  improatom3(numimpro(ipnt),ipnt) = itmp(4)
                  improatom4(numimpro(ipnt),ipnt) = itmp(5)
                endif
              endif

              ipnt = localptr(itmp(4))
              if (ipnt.ne.0) then
                numimpro(ipnt) = numimpro(ipnt) + 1
                if (numimpro(ipnt).le.maximproper) then
                  improtype(numimpro(ipnt),ipnt) = itmp(1)
                  improatom1(numimpro(ipnt),ipnt) = itmp(2)
                  improatom2(numimpro(ipnt),ipnt) = itmp(3)
                  improatom3(numimpro(ipnt),ipnt) = itmp(4)
                  improatom4(numimpro(ipnt),ipnt) = itmp(5)
                endif
              endif
              
              ipnt = localptr(itmp(5))
              if (ipnt.ne.0) then
                numimpro(ipnt) = numimpro(ipnt) + 1
                if (numimpro(ipnt).le.maximproper) then
                  improtype(numimpro(ipnt),ipnt) = itmp(1)
                  improatom1(numimpro(ipnt),ipnt) = itmp(2)
                  improatom2(numimpro(ipnt),ipnt) = itmp(3)
                  improatom3(numimpro(ipnt),ipnt) = itmp(4)
                  improatom4(numimpro(ipnt),ipnt) = itmp(5)
                endif
              endif
            endif

            if (mod(i,1000).eq.0) then
              call mpi_barrier(mpi_comm_world,ierror)
              if (node.eq.0) write (6,*) i
            endif

          enddo

c check for too many impropers/atom

          max_improper = 0
          i = atompnt
          do ii = 1,nlocal
            if (numimpro(i).gt.max_improper) max_improper = numimpro(i)
            i = list(i)
          enddo

          jtmp = max_improper
          call mpi_allreduce(jtmp,max_improper,1,mpi_integer,
     $         mpi_max,mpi_comm_world,ierror)
          if (max_improper.gt.maximproper) then
            if (node.eq.0) write (6,*) 'Max impropers/atom =',
     $           max_improper
            call error('Too many impropers/atom - boost maximproper')
          endif

c check that all impropers were assigned uniquely to a processor

          nimpros_tmp = 0
          i = atompnt
          do ii = 1,nlocal
            nimpros_tmp = nimpros_tmp + numimpro(i)
            i = list(i)
          enddo
          
          jtmp = nimpros_tmp
          call mpi_allreduce(jtmp,nimpros_tmp,1,mpi_integer,
     $         mpi_sum,mpi_comm_world,ierror)
          ifactor = 1
          if (newton_bond.eq.0) ifactor = 4
          if (nimpros_tmp.ne.ifactor*nimpropers)
     $         call error('Impropers assigned incorrectly')

c read and broadcast nonbond coeffs

        else if (match('Nonbond Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  Nonbond Coeffs ...'
          if (nonstyle.eq.0)
     $         call error('Nonbond Coeffs entry not allowed')
          noncoeffflag = 1

          if (nonstyle.eq.1) then
            do i = 1,ntypes
              if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
              call mpi_bcast(tmp,2,mpi_double_precision,0,
     $             mpi_comm_world,ierror)
              noncoeff1(i,i) = tmp(1)
              noncoeff2(i,i) = tmp(2)
              noncoeff3(i,i) = cutlj
            enddo
          else if (nonstyle.eq.2) then
            do i = 1,ntypes
              if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
              call mpi_bcast(tmp,2,mpi_double_precision,0,
     $             mpi_comm_world,ierror)
              noncoeff1(i,i) = tmp(1)
              noncoeff2(i,i) = tmp(2)
              noncoeff3(i,i) = cutljinterior
              noncoeff4(i,i) = cutlj
            enddo
          else if (nonstyle.eq.3) then
            do i = 1,ntypes
              if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,3)
                call mpi_bcast(tmp,3,mpi_double_precision,0,
     $             mpi_comm_world,ierror)
                noncoeff1(i,i) = tmp(1)
                noncoeff2(i,i) = tmp(2)
                noncoeff3(i,i) = tmp(3)
                noncoeff4(i,i) = cutlj
            enddo
          else if (nonstyle.eq.4) then
            do i = 1,ntypes
              if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
              call mpi_bcast(tmp,2,mpi_double_precision,0,
     $             mpi_comm_world,ierror)
              noncoeff1(i,i) = tmp(1)
              noncoeff2(i,i) = tmp(2)
              noncoeff3(i,i) = cutlj
            enddo
          else if (nonstyle.eq.5) then
            do i = 1,ntypes
              if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
              call mpi_bcast(tmp,2,mpi_double_precision,0,
     $             mpi_comm_world,ierror)
              noncoeff1(i,i) = tmp(1)
              noncoeff2(i,i) = tmp(2)
              noncoeff3(i,i) = cutlj
            enddo
          endif

c read and broadcast bond coeffs

        else if (match('Bond Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  Bond Coeffs ...'
          if (nbonds.eq.0.or.bondstyle.eq.0)
     $         call error('Bond Coeffs entry not allowed')
          nbondcoeffflag = 1

          if (bondstyle.eq.1) then
            if (node.eq.0) then
              do i = 1,nbondtypes
                read (2,*) jtmp,(bondcoeff(j,i),j=1,2)
              enddo
            endif
          else if (bondstyle.eq.2) then
            if (node.eq.0) then
              do i = 1,nbondtypes
                read (2,*) jtmp,(bondcoeff(j,i),j=1,4)
              enddo
            endif
          else if (bondstyle.eq.3) then
            if (node.eq.0) then
              do i = 1,nbondtypes
                read (2,*) jtmp,(bondcoeff(j,i),j=1,5)
              enddo
            endif
          else if (bondstyle.eq.4) then
            if (node.eq.0) then
              do i = 1,nbondtypes
                read (2,*) jtmp,(bondcoeff(j,i),j=1,3)
              enddo
            endif
          else if (bondstyle.eq.5) then
            if (node.eq.0) then
              do i = 1,nbondtypes
                read (2,*) jtmp,(bondcoeff(j,i),j=1,4)
              enddo
            endif
          endif
          
          isize = nbondtypes*5
          call mpi_bcast(bondcoeff,isize,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

          do i = 1,nbondtypes
            bondtypeflag(i) = bondstyle
          enddo

c read and broadcast angle coeffs

        else if (match('Angle Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  Angle Coeffs ...'
          if (nangles.eq.0.or.anglestyle.eq.0)
     $         call error('Angle Coeffs entry not allowed')
          nanglecoeffflag = 1

          if (anglestyle.eq.1) then
            if (node.eq.0) then
              do i = 1,nangletypes
                read (2,*) jtmp,(anglecoeff(j,i),j=1,2)
              enddo
            endif
          else if (anglestyle.eq.2) then
            if (node.eq.0) then
              do i = 1,nangletypes
                read (2,*) jtmp,(anglecoeff(j,i),j=1,4)
              enddo
            endif
          endif

          isize = nangletypes*4
          call mpi_bcast(anglecoeff,isize,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

          do i = 1,nangletypes
            angletypeflag(i) = anglestyle
          enddo

c read and broadcast dihedral coeffs

        else if (match('Dihedral Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  Dihedral Coeffs ...'
          if (ndihedrals.eq.0.or.dihedstyle.eq.0)
     $         call error('Dihedral Coeffs entry not allowed')
          ndihedcoeffflag = 1

          if (dihedstyle.eq.1) then
            if (node.eq.0) then
              do i = 1,ndihedtypes
                read (2,*) jtmp,(dihedcoeff(j,i),j=1,3)
              enddo
            endif
          else if (dihedstyle.eq.2) then
            if (node.eq.0) then
              do i = 1,ndihedtypes
                read (2,*) jtmp,(dihedcoeff(j,i),j=1,6)
              enddo
            endif
          endif

          isize = ndihedtypes*6
          call mpi_bcast(dihedcoeff,isize,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

          do i = 1,ndihedtypes
            dihedtypeflag(i) = dihedstyle
          enddo

c read and broadcast improper coeffs

        else if (match('Improper Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  Improper Coeffs ...'
          if (nimpropers.eq.0.or.improstyle.eq.0)
     $         call error('Improper Coeffs entry not allowed')
          nimprocoeffflag = 1

          if (improstyle.eq.1) then
            if (node.eq.0) then
              do i = 1,nimprotypes
                read (2,*) jtmp,(improcoeff(j,i),j=1,2)
              enddo
            endif
          else if (improstyle.eq.2) then
            if (node.eq.0) then
              do i = 1,nimprotypes
                read (2,*) jtmp,(improcoeff(j,i),j=1,3)
              enddo
            endif
          else if (improstyle.eq.3) then
            if (node.eq.0) then
              do i = 1,nimprotypes
                read (2,*) jtmp,(improcoeff(j,i),j=1,2)
              enddo
            endif
          endif

          isize = nimprotypes*3
          call mpi_bcast(improcoeff,isize,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

          do i = 1,nimprotypes
            improtypeflag(i) = improstyle
          enddo

c read and broadcast bond-bond coeffs for class 2 force field

        else if (match('BondBond Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  BondBond Coeffs ...'
          if (nangles.eq.0.or.anglestyle.ne.2)
     $         call error('BondBond Coeffs entry not allowed')
          nbondbondflag = 1

          if (node.eq.0) then
            do i = 1,nangletypes
              read (2,*) jtmp,(bondbondcoeff(j,i),j=1,3)
            enddo
          endif

          isize = nangletypes*3
          call mpi_bcast(bondbondcoeff,isize,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

c read and broadcast bond-angle coeffs for class 2 force field

        else if (match('BondAngle Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  BondAngle Coeffs ...'
          if (nangles.eq.0.or.anglestyle.ne.2)
     $         call error('BondAngle Coeffs entry not allowed')
          nbondangleflag = 1

          if (node.eq.0) then
            do i = 1,nangletypes
              read (2,*) jtmp,(bondanglecoeff(j,i),j=1,4)
            enddo
          endif

          isize = nangletypes*4
          call mpi_bcast(bondanglecoeff,isize,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

c read and broadcast class 2 middle-bond-torsion coeffs

        else if (match('MiddleBondTorsion Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  MiddleBondTorsion Coeffs ...'
          if (ndihedrals.eq.0.or.dihedstyle.ne.2) call error(
     $         'MiddleBondTorsion Coeffs entry not allowed')
          nmidbondtorsionflag = 1

          if (node.eq.0) then
            do i = 1,ndihedtypes
              read (2,*) jtmp,(midbondtorsioncoeff(j,i),j=1,4)
            enddo
          endif

          isize = ndihedtypes*4
          call mpi_bcast(midbondtorsioncoeff,
     $         isize,mpi_double_precision,0,
     $         mpi_comm_world,ierror)

c read and broadcast class 2 end-bond-torsion coeffs

        else if (match('EndBondTorsion Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  EndBondTorsion Coeffs ...'
          if (ndihedrals.eq.0.or.dihedstyle.ne.2) call error(
     $         'EndBondTorsion Coeffs entry not allowed')
          nendbondtorsionflag = 1

          if (node.eq.0) then
            do i = 1,ndihedtypes
              read (2,*) jtmp,(endbondtorsioncoeff(j,i),j=1,8)
            enddo
          endif

          isize = ndihedtypes*8
          call mpi_bcast(endbondtorsioncoeff,
     $         isize,mpi_double_precision,0,mpi_comm_world,ierror)

c read and broadcast class 2 angle-torsion coeffs

        else if (match('AngleTorsion Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  AngleTorsion Coeffs ...'
          if (ndihedrals.eq.0.or.dihedstyle.ne.2) call error(
     $         'AngleTorsion Coeffs entry not allowed')
          nangletorsionflag = 1

          if (node.eq.0) then
            do i = 1,ndihedtypes
              read (2,*) jtmp,(angletorsioncoeff(j,i),j=1,8)
            enddo
          endif

          isize = ndihedtypes*8
          call mpi_bcast(angletorsioncoeff,
     $         isize,mpi_double_precision,0,mpi_comm_world,ierror)

c read and broadcast class 2 angle-angle-torsion coeffs

        else if (match('AngleAngleTorsion Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  AngleAngleTorsion Coeffs ...'
          if (ndihedrals.eq.0.or.dihedstyle.ne.2) call error(
     $         'AngleAngleTorsion Coeffs entry not allowed')
          nangleangletorsionflag = 1

          if (node.eq.0) then
            do i = 1,ndihedtypes
              read (2,*) jtmp,(angleangletorsioncoeff(j,i),j=1,3)
            enddo
          endif

          isize = ndihedtypes*3
          call mpi_bcast(angleangletorsioncoeff,
     $         isize,mpi_double_precision,0,mpi_comm_world,ierror)

c read and broadcast class 2 bondbond13 coeffs

        else if (match('BondBond13 Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  BondBond13 Coeffs ...'
          if (ndihedrals.eq.0.or.dihedstyle.ne.2) call error(
     $         'BondBond13 Coeffs entry not allowed')
          nbondbond13flag = 1

          if (node.eq.0) then
            do i = 1,ndihedtypes
              read (2,*) jtmp,(bondbond13coeff(j,i),j=1,3)
            enddo
          endif

          isize = ndihedtypes*3
          call mpi_bcast(bondbond13coeff,
     $         isize,mpi_double_precision,0,mpi_comm_world,ierror)

c read and broadcast class 2 angle-angle coeffs

        else if (match('AngleAngle Coeffs',str,m)) then

          if (node.eq.0) write (6,*) '  AngleAngle Coeffs ...'
          if (nimpropers.eq.0.or.improstyle.ne.3) call error(
     $         'AngleAngle Coeffs entry not allowed')
          nangleangleflag = 1

          if (node.eq.0) then
            do i = 1,nimprotypes
              read (2,*) jtmp,(angleanglecoeff(j,i),j=1,6)
            enddo
          endif

          isize = nimprotypes*6
          call mpi_bcast(angleanglecoeff,
     $         isize,mpi_double_precision,0,mpi_comm_world,ierror)

c unknown identifier

        else
          call error('Unknown identifier in data file')
        endif

      enddo

c close data file

 999  continue

      if (node.eq.0) then
        str = 'All Done'
        call mpi_bcast(str,80,mpi_character,0,mpi_comm_world,ierror)
        close (2)
      endif

c check that all required options were read in from data file

      if (massflag.eq.0) call error('No Masses in data file')
      if (natomflag.eq.0) call error('No Atoms in data file')
      if (nbonds.gt.0.and.nbondflag.eq.0)
     $     call error('No Bonds in data file')
      if (nangles.gt.0.and.nangleflag.eq.0)
     $     call error('No Angles in data file')
      if (ndihedrals.gt.0.and.ndihedflag.eq.0)
     $     call error('No Dihedrals in data file')
      if (nimpropers.gt.0.and.nimproflag.eq.0)
     $     call error('No Impropers in data file')

      if (nbonds.gt.0.and.bondstyle.eq.5.and.nbondcoeffflag.eq.0)
     $     call error('No Bond Coeffs in data file')
      if (nangles.gt.0.and.nanglecoeffflag.eq.0)
     $     call error('No Angle Coeffs in data file')
      if (ndihedrals.gt.0.and.ndihedcoeffflag.eq.0)
     $     call error('No Dihedral Coeffs in data file')
      if (nimpropers.gt.0.and.nimprocoeffflag.eq.0)
     $     call error('No Improper Coeffs in data file')

      if (nangles.gt.0.and.anglestyle.eq.2.and.nbondbondflag.eq.0)
     $     call error('No BondBond Coeffs in data file')
      if (nangles.gt.0.and.anglestyle.eq.2.and.nbondangleflag.eq.0)
     $     call error('No BondAngle Coeffs in data file')

      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nmidbondtorsionflag.eq.0)
     $     call error('No MiddleBondTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nendbondtorsionflag.eq.0)
     $     call error('No EndBondTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nangletorsionflag.eq.0)
     $     call error('No AngleTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nangleangletorsionflag.eq.0)
     $     call error('No AngleAngleTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nbondbond13flag.eq.0)
     $     call error('No BondBond13 Coeffs in data file')

      if (nimpropers.gt.0.and.improstyle.eq.3.and.
     $     nangleangleflag.eq.0)
     $     call error('No AngleAngle Coeffs in data file')

c compute nonbond mixing rules (geometric,arithmetic,sixthpower)
c  only if nonbond coeffs were input

      if (noncoeffflag.eq.1) then

        do i = 1,ntypes-1
          do j = i+1,ntypes
            if (mixstyle.eq.1.or.mixstyle.eq.2.or.nonstyle.eq.4) then
              noncoeff1(i,j) = sqrt(noncoeff1(i,i)*noncoeff1(j,j))
            else if (mixstyle.eq.3) then
              noncoeff1(i,j) = 2.0 *
     $             sqrt(noncoeff1(i,i)*noncoeff1(j,j)) *
     $             noncoeff2(i,i)**3 * noncoeff2(j,j)**3 /
     $             (noncoeff2(i,i)**6 + noncoeff2(j,j)**6)
            endif
            if (mixstyle.eq.1.or.nonstyle.eq.4) then
              noncoeff2(i,j) = sqrt(noncoeff2(i,i)*noncoeff2(j,j))
            else if (mixstyle.eq.2) then
              noncoeff2(i,j) = (noncoeff2(i,i)+noncoeff2(j,j))/2.0
            else if (mixstyle.eq.3) then
              noncoeff2(i,j) =
     $             ((noncoeff2(i,i)**6 + noncoeff2(j,j)**6) / 2.0) **
     $             (1.0/6.0)
            endif
            if (nonstyle.eq.3) noncoeff3(i,j) =
     $           (noncoeff3(i,i)+noncoeff3(j,j))/2.0
          enddo
        enddo

        do i = 1,ntypes
          do j = i,ntypes
            nontypeflag(i,j) = nonstyle
          enddo
        enddo
        
      endif

c perform unit conversions where necessary
c  convert angle theta from degrees to radians
c  convert improper omega from degrees to radians
c  convert dihedral phi from degrees to radians
c  convert class 2 angles from degrees to radians

      if (nanglecoeffflag.eq.1) then
        if (anglestyle.eq.1) then
          do i = 1,nangletypes
            anglecoeff(2,i) = anglecoeff(2,i)/180.0 * 3.1415926
          enddo
        else if (anglestyle.eq.2) then
          do i = 1,nangletypes
            anglecoeff(1,i) = anglecoeff(1,i)/180.0 * 3.1415926
          enddo
        endif
      endif

      if (improstyle.ne.2.and.nimprocoeffflag.eq.1) then
        do i = 1,nimprotypes
          improcoeff(2,i) = improcoeff(2,i)/180.0 * 3.1415926
        enddo
      endif

      if (dihedstyle.eq.2.and.ndihedcoeffflag.eq.1) then
        do i = 1,ndihedtypes
          dihedcoeff(2,i) = dihedcoeff(2,i)/180.0 * 3.1415926
          dihedcoeff(4,i) = dihedcoeff(4,i)/180.0 * 3.1415926
          dihedcoeff(6,i) = dihedcoeff(6,i)/180.0 * 3.1415926
        enddo
      endif

      if (dihedstyle.eq.2.and.nangletorsionflag.eq.1) then
        do i = 1,ndihedtypes
          angletorsioncoeff(7,i) =
     $         angletorsioncoeff(7,i)/180.0 * 3.1415926
          angletorsioncoeff(8,i) =
     $         angletorsioncoeff(8,i)/180.0 * 3.1415926
        enddo
      endif

      if (dihedstyle.eq.2.and.nangleangletorsionflag.eq.1) then
        do i = 1,ndihedtypes
          angleangletorsioncoeff(2,i) =
     $         angleangletorsioncoeff(2,i)/180.0 * 3.1415926
          angleangletorsioncoeff(3,i) =
     $         angleangletorsioncoeff(3,i)/180.0 * 3.1415926
        enddo
      endif

      if (improstyle.eq.3.and.nangleangleflag.eq.1) then
        do i = 1,nimprotypes
          angleanglecoeff(4,i) =
     $         angleanglecoeff(4,i)/180.0 * 3.1415926
          angleanglecoeff(5,i) =
     $         angleanglecoeff(5,i)/180.0 * 3.1415926
          angleanglecoeff(6,i) =
     $         angleanglecoeff(6,i)/180.0 * 3.1415926
        enddo
      endif

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

      rtmp = masssum
      call mpi_allreduce(rtmp,masssum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      rtmp = qsum
      call mpi_allreduce(rtmp,qsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      rtmp = qsqsum
      call mpi_allreduce(rtmp,qsqsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)

      call mpi_barrier(mpi_comm_world,ierror)

      return
      end
