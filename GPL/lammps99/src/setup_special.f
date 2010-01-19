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
c create all 1st,2nd,3rd neighbor lists

      subroutine setup_special
      include "lammps.h"
      include "mpif.h"

      integer istatus(mpi_status_size)

      call mpi_barrier(mpi_comm_world,ierror)

c use list of bonds to find 1st nearest neighbors of all my atoms

      iorig = atompnt
      do i = 1,nlocal
        itag = tag(iorig)
        num1bond(iorig) = 0
        do j = 1,numbond(iorig)
          if (itag.eq.bondatom1(j,iorig)) newtag = bondatom2(j,iorig)
          if (itag.eq.bondatom2(j,iorig)) newtag = bondatom1(j,iorig)
          num1bond(iorig) = num1bond(iorig) + 1
          if (num1bond(iorig).le.maxsneigh)
     $         specbond(num1bond(iorig),iorig) = newtag
        enddo
        iorig = list(iorig)
      enddo

c check for overflow of specbond list

      iflag = 0
      i = atompnt
      do ii = 1,nlocal
        if (num1bond(i).gt.max(maxsneigh,iflag)) iflag = num1bond(i)
        i = list(i)
      enddo

      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) then
        if (node.eq.0) write (6,*) 'Max special neighs =',jflag
        call error('Too many 1st neighs - boost maxsneigh')
      endif

c if newton flag set, this only acquired 1/2 of all bonds - get other 1/2

      if (newton_bond.eq.1) then

c fill buffer with my bond pairs before checking own partners

        ipnt = 0
        iorig = atompnt
        do i = 1,nlocal
          itag = tag(iorig)
          do j = 1,num1bond(iorig)
            ibuf3(ipnt+1) = tag(iorig)
            ibuf3(ipnt+2) = specbond(j,iorig)
            ipnt = ipnt + 2
          enddo
          iorig = list(iorig)
        enddo

c check my atoms (out of buffer) to see if I own 1st bond partners
c if so, fill special list of bond partner with original atom

        do i = 2,ipnt,2
          inew = localptr(ibuf3(i))
          if (inew.ne.0.and.inew.le.maxown) then
            num1bond(inew) = num1bond(inew) + 1
            if (num1bond(inew).le.maxsneigh)
     $           specbond(num1bond(inew),inew) = ibuf3(i-1)
          endif
        enddo

c check for overflow of specbond list

        iflag = 0
        i = atompnt
        do ii = 1,nlocal
          if (num1bond(i).gt.max(maxsneigh,iflag)) iflag = num1bond(i)
          i = list(i)
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) then
          if (node.eq.0) write (6,*) 'Max special neighs =',jflag
          call error('Too many 1st neighs - boost maxsneigh')
        endif

c cycle previously filled buffer around to all nodes in handshaked fashion
c when receive buffer, scan bond partner tags for atoms I own
c when find one, add original atom to bond partners list

        inext = node + 1
        iprev = node - 1
        if (inext.eq.nprocs) inext = 0
        if (iprev.lt.0) iprev = nprocs - 1

        imesssize = ipnt
        imesstag = 2

        do iwhich = 1,nprocs-1
          call mpi_irecv(ibuf4,maxspec,mpi_integer,iprev,imesstag,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(ibuf3,imesssize,mpi_integer,inext,imesstag,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_integer,imesssize,ierror)
          do i = 1,imesssize
            ibuf3(i) = ibuf4(i)
          enddo
          do i = 2,imesssize,2
            inew = localptr(ibuf3(i))
            if (inew.ne.0.and.inew.le.maxown) then
              num1bond(inew) = num1bond(inew) + 1
              if (num1bond(inew).le.maxsneigh)
     $             specbond(num1bond(inew),inew) = ibuf3(i-1)
            endif
          enddo
        enddo

c check for overflow of specbond list

        iflag = 0
        i = atompnt
        do ii = 1,nlocal
          if (num1bond(i).gt.max(maxsneigh,iflag)) iflag = num1bond(i)
          i = list(i)
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)
        if (jflag.gt.0) then
          if (node.eq.0) write (6,*) 'Max special neighs =',jflag
          call error('Too many 1st neighs - boost maxsneigh')
        endif

      endif

c find 2nd nearest neighbors of all my atoms

c check my atoms to see if I own 1st bond partners
c if so, fill special list of original atom with 2nd nearest neighbors
c for each new neighbor, check that it's not self or already in list

      iorig = atompnt
      do i = 1,nlocal
        itag = tag(iorig)
        num2bond(iorig) = num1bond(iorig)
        do j = 1,num1bond(iorig)
          inew = localptr(specbond(j,iorig))
          if (inew.ne.0.and.inew.le.maxown) then
            do k = 1,num1bond(inew)
              newtag = specbond(k,inew)
              iflag = 0
              do m = 1,num2bond(iorig)
                if (newtag.eq.specbond(m,iorig).or.newtag.eq.itag)
     $               iflag = 1
              enddo
              if (iflag.eq.0) then
                num2bond(iorig) = num2bond(iorig) + 1
                if (num2bond(iorig).le.maxsneigh)
     $               specbond(num2bond(iorig),iorig) = newtag
              endif
            enddo
          endif
        enddo
        iorig = list(iorig)
      enddo

c check for overflow of specbond list

      iflag = 0
      i = atompnt
      do ii = 1,nlocal
        if (num2bond(i).gt.max(maxsneigh,iflag)) iflag = num2bond(i)
        i = list(i)
      enddo

      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) then
        if (node.eq.0) write (6,*) 'Max special neighs =',jflag
        call error('Too many 2nd neighs - boost maxsneigh')
      endif
      
c identify 1st bond partners I don't own
c put original tag and tag of bond partner in buffer
c nreply = # of messages I will eventually receive back

      ipnt = 0
      iorig = atompnt
      do i = 1,nlocal
        do j = 1,num1bond(iorig)
          inew = localptr(specbond(j,iorig))
          if (inew.eq.0.or.inew.gt.maxown) then
            ibuf3(ipnt+1) = tag(iorig)
            ibuf3(ipnt+2) = specbond(j,iorig)
            ipnt = ipnt + 2
          endif
        enddo
        iorig = list(iorig)
      enddo

      nreply = ipnt/2

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan bond partner tags for atoms I own
c for each owned atom send reply message back to requesting node
c reply message contains tags of all neighbors of the bond partner 

      inext = node + 1
      iprev = node - 1
      if (inext.eq.nprocs) inext = 0
      if (iprev.lt.0) iprev = nprocs - 1

      imesstag = 3
      imesssize = ipnt
      ireplytag = 4

      do iwhich = 1,nprocs-1
        call mpi_irecv(ibuf4,maxspec,mpi_integer,iprev,imesstag,
     $       mpi_comm_world,irequest,ierror)
        call mpi_send(ibuf3,imesssize,mpi_integer,inext,imesstag,
     $       mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
        call mpi_get_count(istatus,mpi_integer,imesssize,ierror)
        do i = 1,imesssize
          ibuf3(i) = ibuf4(i)
        enddo
        do i = 1,imesssize,2
          inew = localptr(ibuf3(i+1))
          if (inew.ne.0.and.inew.le.maxown) then
            ibuf4(1) = ibuf3(i)
            do j = 1,num1bond(inew)
              ibuf4(j+1) = specbond(j,inew)
            enddo
            irequester = node - iwhich
            if (irequester.lt.0) irequester = irequester + nprocs
            ireplysize = num1bond(inew) + 1
            call mpi_send(ibuf4,ireplysize,mpi_integer,irequester,
     $           ireplytag,mpi_comm_world,ierror)
          endif
        enddo
      enddo

c read all incoming replies
c fill special list of original atom with received 2nd nearest neighbors
c for each new neighbor, check that it's not self or already in list

      ireplytag = 4
      do i = 1,nreply
        call mpi_recv(ibuf4,maxspec,mpi_integer,mpi_any_source,
     $       ireplytag,mpi_comm_world,istatus,ierror)
        call mpi_get_count(istatus,mpi_integer,ireplysize,ierror)
        iorig = localptr(ibuf4(1))
        itag = tag(iorig)
        do j = 2,ireplysize
          newtag = ibuf4(j)
          iflag = 0
          do k = 1,num2bond(iorig)
            if (newtag.eq.specbond(k,iorig).or.newtag.eq.itag)
     $           iflag = 1
          enddo
          if (iflag.eq.0) then
            num2bond(iorig) = num2bond(iorig) + 1
            if (num2bond(iorig).le.maxsneigh)
     $           specbond(num2bond(iorig),iorig) = newtag
          endif
        enddo
      enddo

c check for overflow of specbond list

      iflag = 0
      i = atompnt
      do ii = 1,nlocal
        if (num2bond(i).gt.max(maxsneigh,iflag)) iflag = num2bond(i)
        i = list(i)
      enddo

      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) then
        if (node.eq.0) write (6,*) 'Max special neighs =',jflag
        call error('Too many 2nd neighs - boost maxsneigh')
      endif

c repeat process for 3rd nearest neighbors

c check my atoms to see if I own 2nd bond partners
c if so, fill special list of original atom with 3rd nearest neighbors
c for each new neighbor, check that it's not self or already in list

      iorig = atompnt
      do i = 1,nlocal
        itag = tag(iorig)
        num3bond(iorig) = num2bond(iorig)
        do j = num1bond(iorig)+1,num2bond(iorig)
          inew = localptr(specbond(j,iorig))
          if (inew.ne.0.and.inew.le.maxown) then
            do k = 1,num1bond(inew)
              newtag = specbond(k,inew)
              iflag = 0
              do m = 1,num3bond(iorig)
                if (newtag.eq.specbond(m,iorig).or.newtag.eq.itag)
     $               iflag = 1
              enddo
              if (iflag.eq.0) then
                num3bond(iorig) = num3bond(iorig) + 1
                if (num3bond(iorig).le.maxsneigh)
     $               specbond(num3bond(iorig),iorig) = newtag
              endif
            enddo
          endif
        enddo
        iorig = list(iorig)
      enddo

c check for overflow of specbond list

      iflag = 0
      i = atompnt
      do ii = 1,nlocal
        if (num3bond(i).gt.max(maxsneigh,iflag)) iflag = num3bond(i)
        i = list(i)
      enddo

      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) then
        if (node.eq.0) write (6,*) 'Max special neighs =',jflag
        call error('Too many 3rd neighs - boost maxsneigh')
      endif
      
c identify 2nd bond partners I don't own
c put original tag and tag of bond partner in buffer
c nreply = # of messages I will eventually receive back

      ipnt = 0
      iorig = atompnt
      do i = 1,nlocal
        do j = num1bond(iorig)+1,num2bond(iorig)
          inew = localptr(specbond(j,iorig))
          if (inew.eq.0.or.inew.gt.maxown) then
            ibuf3(ipnt+1) = tag(iorig)
            ibuf3(ipnt+2) = specbond(j,iorig)
            ipnt = ipnt + 2
          endif
        enddo
        iorig = list(iorig)
      enddo

      nreply = ipnt/2

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan bond partner tags for atoms I own
c for each owned atom send reply message back to requesting node
c reply message contains tags of all neighbors of the bond partner 

      inext = node + 1
      iprev = node - 1
      if (inext.eq.nprocs) inext = 0
      if (iprev.lt.0) iprev = nprocs - 1

      imesstag = 5
      imesssize = ipnt
      ireplytag = 6

      do iwhich = 1,nprocs-1
        call mpi_irecv(ibuf4,maxspec,mpi_integer,iprev,imesstag,
     $       mpi_comm_world,irequest,ierror)
        call mpi_send(ibuf3,imesssize,mpi_integer,inext,imesstag,
     $       mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
        call mpi_get_count(istatus,mpi_integer,imesssize,ierror)
        do i = 1,imesssize
          ibuf3(i) = ibuf4(i)
        enddo
        do i = 1,imesssize,2
          inew = localptr(ibuf3(i+1))
          if (inew.ne.0.and.inew.le.maxown) then
            ibuf4(1) = ibuf3(i)
            do j = 1,num1bond(inew)
              ibuf4(j+1) = specbond(j,inew)
            enddo
            irequester = node - iwhich
            if (irequester.lt.0) irequester = irequester + nprocs
            ireplysize = num1bond(inew) + 1
            call mpi_send(ibuf4,ireplysize,mpi_integer,irequester,
     $           ireplytag,mpi_comm_world,ierror)
          endif
        enddo
      enddo

c read all incoming replies
c fill special list of original atom with received 3rd nearest neighbors
c for each new neighbor, check that it's not self or already in list

      ireplytag = 6
      do i = 1,nreply
        call mpi_recv(ibuf4,maxspec,mpi_integer,mpi_any_source,
     $       ireplytag,mpi_comm_world,istatus,ierror)
        call mpi_get_count(istatus,mpi_integer,ireplysize,ierror)
        iorig = localptr(ibuf4(1))
        itag = tag(iorig)
        do j = 2,ireplysize
          newtag = ibuf4(j)
          iflag = 0
          do k = 1,num3bond(iorig)
            if (newtag.eq.specbond(k,iorig).or.newtag.eq.itag)
     $           iflag = 1
          enddo
          if (iflag.eq.0) then
            num3bond(iorig) = num3bond(iorig) + 1
            if (num3bond(iorig).le.maxsneigh)
     $           specbond(num3bond(iorig),iorig) = newtag
          endif
        enddo
      enddo

c check for overflow of specbond list

      iflag = 0
      i = atompnt
      do ii = 1,nlocal
        if (num3bond(i).gt.max(maxsneigh,iflag)) iflag = num3bond(i)
        i = list(i)
      enddo

      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) then
        if (node.eq.0) write (6,*) 'Max special neighs =',jflag
        call error('Too many 3rd neighs - boost maxsneigh')
      endif

      call mpi_barrier(mpi_comm_world,ierror)

      return
      end
