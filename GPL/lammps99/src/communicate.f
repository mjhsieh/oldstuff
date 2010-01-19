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

c -----------------------------------------------------------------------
c interprocessor communication via 6-way stencil

c -------------------------------------------------------------------------
c swap slabs of atoms with other processors in all 3 directions
c done every timestep
c send/receive message always to another proc - never to self
      
      subroutine communicate
      include "lammps.h"
      include "mpif.h"
      
      integer istatus(mpi_status_size)

      do k = 1,nswap

        j = 0
        do ii = nslist(k),nslist(k+1)-1
          i = slist(ii)
          buf3(j+1) = x(1,i)
          buf3(j+2) = x(2,i)
          buf3(j+3) = x(3,i)
          j = j + 3
        enddo

        isize_send = nslist(k+1) - nslist(k)
        isize_recv = nrlist(k+1) - nrlist(k)
        call mpi_irecv(x(1,nrlist(k)),3*isize_recv,
     $       mpi_double_precision,rpart(k),0,
     $       mpi_comm_world,irequest,ierror)
        call mpi_send(buf3,3*isize_send,
     $       mpi_double_precision,spart(k),0,
     $       mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)

      enddo
      
      return
      end
      
      
c -------------------------------------------------------------------------
c reverse swap slabs of forces with other processors in all 3 directions
c sum to existing forces
c done every timestep
c send/receive message always to another proc - never to self
      
      subroutine reverse_comm
      include "lammps.h"
      include "mpif.h"

      integer istatus(mpi_status_size)

      do k = nswap,1,-1

        isize_send = nrlist(k+1) - nrlist(k)
        isize_recv = nslist(k+1) - nslist(k)
        call mpi_irecv(buf3,3*isize_recv,
     $       mpi_double_precision,spart(k),0,
     $       mpi_comm_world,irequest,ierror)
        call mpi_send(f(1,nrlist(k)),3*isize_send,
     $       mpi_double_precision,rpart(k),0,
     $       mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)

        j = 0
        do ii = nslist(k),nslist(k+1)-1
          i = slist(ii)
          f(1,i) = f(1,i) + buf3(j+1)
          f(2,i) = f(2,i) + buf3(j+2)
          f(3,i) = f(3,i) + buf3(j+3)
          j = j + 3
        enddo

      enddo
      
      return
      end
      

c -------------------------------------------------------------------------
c send out atoms that have left my box, receive ones entering my box
c  since last reneighboring
c done in all 3 directions
c all info associated with an atom migrates with the atom
c called right before each reneighboring
      
      subroutine exchange
      include "lammps.h"
      include "mpif.h"

      integer istatus(mpi_status_size)
      
c setup communication volumes if box has changed

      if (boxflag.eq.1) call setup_comm

c clear out other atoms from global table before exchanged atoms migrate

      do i = maxown+1,maxown+nother
        localptr(tag(i)) = 0
      enddo
      nother = 0

c loop over all 3 dimensions

      do k = 1,3
        
        if (pgrid(k).gt.1) then
          
          blo = border(1,k)
          bhi = border(2,k)
          ndelete = 0
          
c fill buffer with atoms leaving my box, update local list
          
          iprev = 0
          j = 0
          i = atompnt
          do ii = 1,nlocal
            if (x(k,i).lt.blo.or.x(k,i).ge.bhi) then
              ndelete = ndelete + 1
              if (ndelete.le.maxexch) then
                buf1(j+1) = x(1,i)
                buf1(j+2) = x(2,i)
                buf1(j+3) = x(3,i)
                buf1(j+4) = tag(i)
                buf1(j+5) = molecule(i)
                buf1(j+6) = type(i)
                buf1(j+7) = q(i)
                buf1(j+8) = v(1,i)
                buf1(j+9) = v(2,i)
                buf1(j+10) = v(3,i)
                buf1(j+11) = true(i)
                buf1(j+12) = fix(i)
                if (nrespa.eq.0) then
                  buf1(j+13) = numbond(i)
                  buf1(j+14) = numangle(i)
                  buf1(j+15) = numdihed(i)
                  buf1(j+16) = numimpro(i)
                  buf1(j+17) = num1bond(i)
                  buf1(j+18) = num2bond(i)
                  buf1(j+19) = num3bond(i)
                else
                  buf1(j+13) = f_stretch(1,i)
                  buf1(j+14) = f_stretch(2,i)
                  buf1(j+15) = f_stretch(3,i)
                  buf1(j+16) = f_intra(1,i)
                  buf1(j+17) = f_intra(2,i)
                  buf1(j+18) = f_intra(3,i)
                  buf1(j+19) = numbond(i)
                  buf1(j+20) = numangle(i)
                  buf1(j+21) = numdihed(i)
                  buf1(j+22) = numimpro(i)
                  buf1(j+23) = num1bond(i)
                  buf1(j+24) = num2bond(i)
                  buf1(j+25) = num3bond(i)
                endif
                j = j + peratom_comm
                do jj = 1,numbond(i)
                  buf1(j+1) = bondatom1(jj,i)
                  buf1(j+2) = bondatom2(jj,i)
                  buf1(j+3) = bondtype(jj,i)
                  j = j + 3
                enddo
                do jj = 1,numangle(i)
                  buf1(j+1) = angleatom1(jj,i)
                  buf1(j+2) = angleatom2(jj,i)
                  buf1(j+3) = angleatom3(jj,i)
                  buf1(j+4) = angletype(jj,i)
                  j = j + 4
                enddo
                do jj = 1,numdihed(i)
                  buf1(j+1) = dihedatom1(jj,i)
                  buf1(j+2) = dihedatom2(jj,i)
                  buf1(j+3) = dihedatom3(jj,i)
                  buf1(j+4) = dihedatom4(jj,i)
                  buf1(j+5) = dihedtype(jj,i)
                  j = j + 5
                enddo
                do jj = 1,numimpro(i)
                  buf1(j+1) = improatom1(jj,i)
                  buf1(j+2) = improatom2(jj,i)
                  buf1(j+3) = improatom3(jj,i)
                  buf1(j+4) = improatom4(jj,i)
                  buf1(j+5) = improtype(jj,i)
                  j = j + 5
                enddo
                do jj = 1,num3bond(i)
                  j = j + 1
                  buf1(j) = specbond(jj,i)
                enddo
              endif
              localptr(tag(i)) = 0
              if (iprev.eq.0) then
                atompnt = list(i)
              else
                list(iprev) = list(i)
              endif
              itmp = list(i)
              list(i) = freepnt
              freepnt = i
              i = itmp
            else
              iprev = i
              i = list(i)
            endif
          enddo
          
          nlocal = nlocal - ndelete
          max_exch = max(max_exch,ndelete)
          
          if (ndelete.gt.maxexch) then
            write (6,*)
     $           'Sending too many exchange atoms - boost maxexch:',
     $           node,ntimestep,ndelete,k
            call exit(0)
          endif
          
c send them out in both directions (if neighboring nodes are different)
          
          call mpi_irecv(buf2,maxexchtot,
     $         mpi_double_precision,mpart(2,k),0,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf1,j,mpi_double_precision,mpart(1,k),0,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_double_precision,
     $         icnt,ierror)

          if (pgrid(k).gt.2) then
            call mpi_irecv(buf2(icnt+1),2*maxexchtot-icnt,
     $           mpi_double_precision,mpart(1,k),0,
     $           mpi_comm_world,irequest,ierror)
            call mpi_send(buf1,j,mpi_double_precision,mpart(2,k),0,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
            call mpi_get_count(istatus,mpi_double_precision,
     $           itmp,ierror)
            icnt = icnt + itmp
          endif

c check incoming atoms to see if they are in my box (could be in node's
c  box on other side of the sender)
          
          j = 0
 10       if (j.lt.icnt) then
            tmp = buf2(j+k)
            if (tmp.ge.blo.and.tmp.lt.bhi) then
              nlocal = nlocal + 1
              if (freepnt.ne.0) then
                itmp = atompnt
                atompnt = freepnt
                freepnt = list(freepnt)
                list(atompnt) = itmp
                x(1,atompnt) = buf2(j+1)
                x(2,atompnt) = buf2(j+2)
                x(3,atompnt) = buf2(j+3)
                tag(atompnt) = nint(buf2(j+4))
                molecule(atompnt) = nint(buf2(j+5))
                type(atompnt) = nint(buf2(j+6))
                q(atompnt) = buf2(j+7)
                v(1,atompnt) = buf2(j+8)
                v(2,atompnt) = buf2(j+9)
                v(3,atompnt) = buf2(j+10)
                true(atompnt) = nint(buf2(j+11))
                fix(atompnt) = nint(buf2(j+12))
                if (nrespa.eq.0) then
                  numbond(atompnt) = nint(buf2(j+13))
                  numangle(atompnt) = nint(buf2(j+14))
                  numdihed(atompnt) = nint(buf2(j+15))
                  numimpro(atompnt) = nint(buf2(j+16))
                  num1bond(atompnt) = nint(buf2(j+17))
                  num2bond(atompnt) = nint(buf2(j+18))
                  num3bond(atompnt) = nint(buf2(j+19))
                else
                  f_stretch(1,atompnt) = buf2(j+13)
                  f_stretch(2,atompnt) = buf2(j+14)
                  f_stretch(3,atompnt) = buf2(j+15)
                  f_intra(1,atompnt) = buf2(j+16)
                  f_intra(2,atompnt) = buf2(j+17)
                  f_intra(3,atompnt) = buf2(j+18)
                  numbond(atompnt) = nint(buf2(j+19))
                  numangle(atompnt) = nint(buf2(j+20))
                  numdihed(atompnt) = nint(buf2(j+21))
                  numimpro(atompnt) = nint(buf2(j+22))
                  num1bond(atompnt) = nint(buf2(j+23))
                  num2bond(atompnt) = nint(buf2(j+24))
                  num3bond(atompnt) = nint(buf2(j+25))
                endif
                j = j + peratom_comm
                do jj = 1,numbond(atompnt)
                  bondatom1(jj,atompnt) = nint(buf2(j+1))
                  bondatom2(jj,atompnt) = nint(buf2(j+2))
                  bondtype(jj,atompnt) = nint(buf2(j+3))
                  j = j + 3
                enddo
                do jj = 1,numangle(atompnt)
                  angleatom1(jj,atompnt) = nint(buf2(j+1))
                  angleatom2(jj,atompnt) = nint(buf2(j+2))
                  angleatom3(jj,atompnt) = nint(buf2(j+3))
                  angletype(jj,atompnt) = nint(buf2(j+4))
                  j = j + 4
                enddo
                do jj = 1,numdihed(atompnt)
                  dihedatom1(jj,atompnt) = nint(buf2(j+1))
                  dihedatom2(jj,atompnt) = nint(buf2(j+2))
                  dihedatom3(jj,atompnt) = nint(buf2(j+3))
                  dihedatom4(jj,atompnt) = nint(buf2(j+4))
                  dihedtype(jj,atompnt) = nint(buf2(j+5))
                  j = j + 5
                enddo
                do jj = 1,numimpro(atompnt)
                  improatom1(jj,atompnt) = nint(buf2(j+1))
                  improatom2(jj,atompnt) = nint(buf2(j+2))
                  improatom3(jj,atompnt) = nint(buf2(j+3))
                  improatom4(jj,atompnt) = nint(buf2(j+4))
                  improtype(jj,atompnt) = nint(buf2(j+5))
                  j = j + 5
                enddo
                do jj = 1,num3bond(atompnt)
                  j = j + 1
                  specbond(jj,atompnt) = nint(buf2(j))
                enddo
                localptr(tag(atompnt)) = atompnt
              else
                j = j + peratom_comm +
     $               3*nint(buf2(j+peratom_comm-6)) +
     $               4*nint(buf2(j+peratom_comm-5)) +
     $               5*nint(buf2(j+peratom_comm-4)) +
     $               5*nint(buf2(j+peratom_comm-3)) +
     $               nint(buf2(j+peratom_comm))
              endif
            else
                j = j + peratom_comm +
     $               3*nint(buf2(j+peratom_comm-6)) +
     $               4*nint(buf2(j+peratom_comm-5)) +
     $               5*nint(buf2(j+peratom_comm-4)) +
     $               5*nint(buf2(j+peratom_comm-3)) +
     $               nint(buf2(j+peratom_comm))
            endif
            goto 10
          endif
          
          if (nlocal.gt.maxown) then
            write (6,*)
     $           'Recvd too many exchange atoms - boost maxown:',
     $           node,ntimestep,nlocal,k
            call exit(0)
          endif
          
        endif
        
      enddo

      return
      end
      
      
c -------------------------------------------------------------------------
c make lists of nearby atoms to send to neighboring nodes at every timestep
c one list made for every swap that will be made
c as list is made, actually do swaps
c this does equivalent of a communicate (so don't need to explicitly
c  call communicate routine on reneighboring timestep)
c tag, type, 3 atom positions are what is swapped
c  since they are needed for force computation
c called right before each reneighboring
      
      subroutine borders
      include "lammps.h"
      include "mpif.h"
      
      integer irequest(4),istatus(mpi_status_size,4)

      nswap = 0
      npnt = 1
      
      do k = 1,3
        do kk = 0,2*need(k)-1
          
c buffer up all atoms I own (plus those previously received) that are
c  inside slab boundaries
c also store pointers to those atoms in slist for communicate routine
c  to use in future timesteps
c ship 3 positions, tag, and type for each atom
          
          nswap = nswap + 1
          nslist(nswap) = npnt
          nrlist(nswap) = maxown + nother + 1
          blo = boundlo(nswap)
          bhi = boundhi(nswap)
          
          nsend = 0
          j = 0
          i = atompnt
          do ii = 1,nlocal+nother
            if (x(k,i).ge.blo.and.x(k,i).lt.bhi) then
              if (npnt.le.maxsend) slist(npnt) = i
              npnt = npnt + 1
              nsend = nsend + 1
              if (nsend.le.maxsendone) then
                buf3(j+1) = x(1,i)
                buf3(j+2) = x(2,i)
                buf3(j+3) = x(3,i)
                ibuf1(nsend) = tag(i)
                ibuf2(nsend) = type(i)
                buf6(nsend) = q(i)
              endif
              j = j + 3
            endif
            if (ii.le.nlocal) then
              i = list(i)
            else
              i = i + 1
            endif
          enddo
          
          max_swap = max(max_swap,nsend)
          
          if (npnt.gt.maxsend) then
            write (6,*)
     $           'Too many atoms in border list - boost maxsend:',
     $           node,ntimestep,npnt,k,kk
            call exit(0)
          endif
          
          if (nsend.gt.maxsendone) then
            write (6,*)
     $           'Sending too many border atoms - boost maxsendone:',
     $           node,ntimestep,nsend,k,kk
            call exit(0)
          endif
          
c swap atoms
c  put incoming ones at end of my position, charge, tag, type arrays
c irecvmax = max # I can recv
c check on it by handshaking first since MPI may otherwise crash
          
          istart = nrlist(nswap)
          irecvmax = maxother - nother

          call mpi_send(nsend,1,mpi_integer,spart(nswap),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nrecv,1,mpi_integer,rpart(nswap),0,
     $         mpi_comm_world,istatus,ierror)

          if (nrecv.gt.irecvmax) then
            write (6,*)
     $           'Recvd too many border atoms - boost maxother:',
     $           node,ntimestep,nother+nrecv,k,kk
            call exit(0)
          endif

          call mpi_irecv(x(1,istart),3*nrecv,
     $         mpi_double_precision,rpart(nswap),0,
     $         mpi_comm_world,irequest(1),ierror)
          call mpi_irecv(q(istart),nrecv,
     $         mpi_double_precision,rpart(nswap),0,
     $         mpi_comm_world,irequest(2),ierror)
          call mpi_irecv(tag(istart),nrecv,
     $         mpi_integer,rpart(nswap),0,
     $         mpi_comm_world,irequest(3),ierror)
          call mpi_irecv(type(istart),nrecv,
     $         mpi_integer,rpart(nswap),0,
     $         mpi_comm_world,irequest(4),ierror)

          call mpi_send(buf3,j,
     $         mpi_double_precision,spart(nswap),0,
     $         mpi_comm_world,ierror)
          call mpi_send(buf6,nsend,
     $         mpi_double_precision,spart(nswap),0,
     $         mpi_comm_world,ierror)
          call mpi_send(ibuf1,nsend,
     $         mpi_integer,spart(nswap),0,
     $         mpi_comm_world,ierror)
          call mpi_send(ibuf2,nsend,
     $         mpi_integer,spart(nswap),0,
     $         mpi_comm_world,ierror)

          call mpi_waitall(4,irequest,istatus,ierror)
          nother = nother + nrecv

        enddo
      enddo
      
      nslist(nswap+1) = npnt
      nrlist(nswap+1) = maxown + nother + 1
      
c fill in global ptrs to new other atoms

      do i = maxown+1,maxown+nother
        localptr(tag(i)) = i
      enddo

      return
      end
