c Create a collection of polymer chains of various lengths and bead sizes
c  write out as a LAMMPS data file (version 99)

c syntax:  setup_chain < chain.def > data.file

c chain.def is an input file that specifies the chains
c data.file in an output file that will be an input dataset for LAMMPS

      program setup_chain
      parameter (maxatom=100000,maxset=5)

      integer nchain(maxset),nmonomer(maxset)
      integer ntype(maxset),nbondtype(maxset)
      real x(maxatom),y(maxatom),z(maxatom)
      real bondlength(maxset),restrict(maxset)
      common xprd,yprd,zprd,xboundlo,xboundhi,
     $     yboundlo,yboundhi,zboundlo,zboundhi

c read chain definitions

      read (5,*)
      read (5,*)
      read (5,*) rhostar
      read (5,*) iseed
      read (5,*) nsets

      if (nsets.gt.maxset) then
        write (6,*) 'Too many chain sets'
        stop
      endif
      
      do iset = 1,nsets
        read (5,*)
        read (5,*) nchain(iset)
        read (5,*) nmonomer(iset)
        read (5,*) ntype(iset)
        read (5,*) nbondtype(iset)
        read (5,*) bondlength(iset)
        read (5,*) restrict(iset)
      enddo

c sum total atoms and error check

      natoms = 0
      do iset = 1,nsets
        natoms = natoms + nchain(iset)*nmonomer(iset)
      enddo

      if (natoms.gt.maxatom) then
        write (6,*) 'Too many total atoms'
        stop
      endif

c setup box size (sigma = 1.0)

      volume = natoms/rhostar
      xprd = volume**(1.0/3.0)
      yprd = xprd
      zprd = xprd

      xboundlo = -xprd/2.0
      xboundhi = -xboundlo
      yboundlo = xboundlo
      yboundhi = xboundhi
      zboundlo = xboundlo
      zboundhi = xboundhi

c generate random chains
c  loop over sets and chains in each set
      
      n = 0
      do iset = 1,nsets
        do ichain = 1,nchain(iset)

c random starting point for the chain in the box

          x1 = 0.0
          y1 = 0.0
          z1 = 0.0
          x2 = xboundlo + random(iseed)*xprd
          y2 = yboundlo + random(iseed)*yprd
          z2 = zboundlo + random(iseed)*zprd
          call pbc(x2,y2,z2)
          n = n + 1
          x(n) = x2
          y(n) = y2
          z(n) = z2

c generate rest of monomers in this chain

          do imonomer = 2,nmonomer(iset)

            x0 = x1
            y0 = y1
            z0 = z1
            x1 = x2
            y1 = y2
            z1 = z2

c random point inside sphere of unit radius

 10         xinner = 2.0*random(iseed) - 1.0
            yinner = 2.0*random(iseed) - 1.0
            zinner = 2.0*random(iseed) - 1.0
            rsq = xinner*xinner + yinner*yinner + zinner*zinner
            if (rsq.gt.1.0) goto 10

c project point to surface of sphere of unit radius

            r = sqrt(rsq)
            xsurf = xinner/r
            ysurf = yinner/r
            zsurf = zinner/r

c create new point by scaling unit offsets by bondlength (sigma = 1.0)

            x2 = x1 + xsurf*bondlength(iset)
            y2 = y1 + ysurf*bondlength(iset)
            z2 = z1 + zsurf*bondlength(iset)

c check that new point meets restriction requirement
c  only for 3rd monomer and beyone

            dx = x2 - x0
            dy = y2 - y0
            dz = z2 - z0
            r = sqrt(dx*dx + dy*dy + dz*dz)

            if (imonomer.gt.2.and.r.le.restrict(iset)) goto 10

c store new point

            call pbc(x2,y2,z2)
            n = n + 1
            x(n) = x2
            y(n) = y2
            z(n) = z2

          enddo

        enddo
      enddo

c compute quantities needed for LAMMPS file

      nbonds = 0
      ntypes = 0
      nbondtypes = 0
      do iset = 1,nsets
        nbonds = nbonds + nchain(iset)*(nmonomer(iset)-1)
        if (ntype(iset).gt.ntypes) ntypes = ntype(iset)
        if (nbondtype(iset).gt.nbondtypes)
     $       nbondtypes = nbondtype(iset)
      enddo

c write out LAMMPS file

      call prints('LAMMPS Description\n')
      call printn
      call printsi('%d atoms\n',natoms)
      call printsi('%d bonds\n',nbonds)
      call prints('0 angles\n')
      call prints('0 dihedrals\n')
      call prints('0 impropers\n')
      call printn
      call printsi('%d atom types\n',ntypes)
      if (nbonds.gt.0) call printsi('%d bond types\n',nbondtypes)
      call printn
      call printf(xboundlo)
      call printf(xboundhi)
      call prints('xlo xhi\n')
      call printf(yboundlo)
      call printf(yboundhi)
      call prints('ylo yhi\n')
      call printf(zboundlo)
      call printf(zboundhi)
      call prints('zlo zhi\n')

      call printn
      call prints('Masses\n')
      call printn

      do i = 1,ntypes
        call printi(i)
        call printf(1.0)
        call printn
      enddo
      
      call printn
      call prints('Atoms\n')
      call printn

      n = 0
      molecule = 0
      do iset = 1,nsets
        do ichain = 1,nchain(iset)
          molecule = molecule + 1
          do imonomer = 1,nmonomer(iset)
            n = n + 1
            call printi(n)
            call printi(molecule)
            call printi(ntype(iset))
            call printf(0.0)
            call printf(x(n))
            call printf(y(n))
            call printf(z(n))
            call printn
          enddo
        enddo
      enddo

      if (nbonds.gt.0) then

        call printn
        call prints('Bonds\n')
        call printn

        n = 0
        m = 0
        do iset = 1,nsets
          do ichain = 1,nchain(iset)
            do imonomer = 1,nmonomer(iset)
              n = n + 1
              if (imonomer.ne.nmonomer(iset)) then
                m = m + 1
                call printi(m)
                call printi(nbondtype(iset))
                call printi(n)
                call printi(n+1)
                call printn
              endif
            enddo
          enddo
        enddo

      endif

      end

        
c ************
c Subroutines
c ************

c periodic boundary conditions - map atom back into periodic box

      subroutine pbc(x,y,z)
      common xprd,yprd,zprd,xboundlo,xboundhi,
     $     yboundlo,yboundhi,zboundlo,zboundhi

      if (x.lt.xboundlo) x = x + xprd
      if (x.ge.xboundhi) x = x - xprd
      if (y.lt.yboundlo) y = y + yprd
      if (y.ge.yboundhi) y = y - yprd
      if (z.lt.zboundlo) z = z + zprd
      if (z.ge.zboundhi) z = z - zprd

      return
      end


c RNG - compute in double precision, return single
      
      real*4 function random(iseed)
      real*8 aa,mm,sseed
      parameter (aa=16807.0D0,mm=2147483647.0D0)
      
      sseed = iseed
      sseed = mod(aa*sseed,mm)
      random = sseed/mm
      iseed = sseed

      return
      end
