c Create LAMMPS data file (version 99)
c  for 3-d LJ simulation of monomers
c  can be used to setup mixture runs - but will take a long time if
c  other componenets are very high percentage of total

c syntax: setup_lj data.file
c
c prompts for cubic lattice size: total atoms = 4 * nx * ny * nz

      program setup_lj
      parameter (maxtotal=100000,maxtype=3)
      real*4 x(3,maxtotal)
      integer target(maxtype),type(maxtotal)
      character*80 polyfile

c liquid Argon Lennard-Jones values

      real*4 epsilon,sigma,mass
c      parameter (epsilon=120.0,sigma=3.4,mass=39.948)
      parameter (epsilon=1.0,sigma=1.0,mass=1.0)

c prompts

      if (iargc().ne.1) then
        write (6,*) 'Syntax: setup_lj data.file'
        call exit(0)
      else
        call getarg(1,polyfile)
      endif

      write (6,*) '(0) FCC Lattice or (1) Number of atoms'
      read (5,*) isetup
      if (isetup.eq.0) then
        write (6,*) 'Lattice: nx,ny,nz'
        read (5,*) nx,ny,nz
        natoms = 4*nx*ny*nz
      else
        write (6,*) 'Number of atoms: n'
        read (5,*) natoms
      endif

      if (natoms.gt.maxtotal) then
        write (6,*) 'Too many atoms'
        call exit(0)
      endif

      write (6,*) 'Rhostar'
      read (5,*) rhostar

      write (6,*) 'Random # seed (8 digit max)'
      read (5,*) iseed

      write (6,*) 'Number of atom types'
      read (5,*) ntypes

      if (ntypes.gt.maxtype) then
        write (6,*) 'Too many atom types'
        call exit(0)
      endif

      if (ntypes.gt.1) then
        nsum = 0
        do i = 1,ntypes
          write (6,*) '# of type',i
          read (5,*) target(i)
          nsum = nsum + target(i)
        enddo
        if (nsum.ne.natoms) then
          write (6,*) 'Total of types <> total # of atoms'
          call exit(0)
        endif
      endif

c setup rectangular or cubic box

      volume = natoms * sigma*sigma*sigma / rhostar
      alattice = (volume/(nx*ny*nz))**(1.0/3.0)

      if (isetup.eq.0) then
        xboundlo = -(nx*alattice)/2.0
        xboundhi = -xboundlo
        yboundlo = -(ny*alattice)/2.0
        yboundhi = -yboundlo
        zboundlo = -(nz*alattice)/2.0
        zboundhi = -zboundlo
      else
        xboundlo = -(volume**(1.0/3.0))/2.0
        xboundhi = -xboundlo
        yboundlo = xboundlo
        yboundhi = xboundhi
        zboundlo = xboundlo
        zboundhi = xboundhi
      endif

c initial fcc lattice or random points

      if (isetup.eq.0) then

        m = 0
        do k = 1,nz*2
          do j = 1,ny*2
            do i = 1,nx*2
              if (mod(i+j+k,2).eq.1) then
                m = m + 1
                type(m) = 1
                x(1,m) = xboundlo + (i-1)*alattice/2.0
                x(2,m) = yboundlo + (j-1)*alattice/2.0
                x(3,m) = zboundlo + (k-1)*alattice/2.0
              endif
            enddo
          enddo
        enddo
        
      else

        do m = 1,natoms
          type(m) = 1
          x(1,m) = xboundlo + random(iseed)*(xboundhi-xboundlo)
          x(2,m) = yboundlo + random(iseed)*(yboundhi-yboundlo)
          x(3,m) = zboundlo + random(iseed)*(zboundhi-zboundlo)
        enddo

      endif

c convert some type 1 atoms to mixture types as needed

      do i = 2,ntypes
        do j = 1,target(i)
 10       k = natoms*random(iseed) + 1
          if (k.le.0) k = 1
          if (k.gt.natoms) k = natoms
          if (type(k).ne.1) goto 10
          type(k) = i
        enddo
      enddo

c write out file

      call fopen(polyfile(1:length(polyfile)))

      call fprints('LAMMPS Description\n')
      call fprintn

      call fprintsi('%d atoms\n',natoms)
      call fprintsi('%d bonds\n',0)
      call fprintsi('%d angles\n',0)
      call fprintsi('%d dihedrals\n',0)
      call fprintsi('%d impropers\n',0)
      call fprintn
      call fprintsi('%d atom types\n',ntypes)
      call fprintn
      call fprintf(xboundlo)
      call fprintf(xboundhi)
      call fprints('xlo xhi\n')
      call fprintf(yboundlo)
      call fprintf(yboundhi)
      call fprints('ylo yhi\n')
      call fprintf(zboundlo)
      call fprintf(zboundhi)
      call fprints('zlo zhi\n')

      call fprintn
      call fprints('Masses\n')
      call fprintn

      do i = 1,ntypes
        call fprinti(i)
        call fprintf(1.0)
        call fprintn
      enddo

      call fprintn
      call fprints('Atoms\n')
      call fprintn

      do i = 1,natoms
        call fprinti(i)
        call fprinti(0)
        call fprinti(type(i))
        call fprintf(0.0)
        call fprintf(x(1,i))
        call fprintf(x(2,i))
        call fprintf(x(3,i))
        call fprintn
      enddo

      call fclose

      end


c ************
c Subroutines
c ************

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


c returns the actual length of str
c backtracks from end of string skipping over spaces
c break if test into 2 to avoid compiler evaluation of str(0:0)

      integer function length(str)
      character*(*) str
      
      n = len(str)
 10   if (n.gt.0) then
        if (str(n:n).eq.' ') then
          n = n - 1
          goto 10
        endif
      endif
      length = n

      return
      end
