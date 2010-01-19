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
c string and number routines

c -----------------------------------------------------------------------
c returns TRUE if str1 matches 1st chars of str2, FALSE otherwise
c also returns m = loc of next char in str2
c could test if next char is a space or tab

      logical function match(str1,str2,m)
      implicit none
      character*(*) str1,str2
      integer m

      match = .FALSE.
      m = len(str1) + 1
      if (len(str1).gt.len(str2)) return
      if (str1.eq.str2(1:len(str1))) match = .TRUE.

      return
      end


c -----------------------------------------------------------------------
c return the integer that starts at loc m in str
c skip initial spaces and tabs
c scans to end of number marked by space or tab
c also return m = loc of char after end of number in str

      integer function intread(str,m)
      implicit none
      character*(*) str
      integer m

      integer n1,n2,nmax

      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.
     $     (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' '
     $     .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      read (str(n1:n2-1),*) intread
      m = n2

      return
      end
      
      
c -----------------------------------------------------------------------
c return the real*8 that starts at loc m in str
c skip initial spaces and tabs
c scans to end of number marked by space or tab
c also return m = loc of char after end of number in str

      real*8 function realread(str,m)
      implicit none
      character*(*) str
      integer m

      integer n1,n2,nmax
      
      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.
     $     (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' '
     $     .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      read (str(n1:n2-1),*) realread
      m = n2

      return
      end
      
      
c -----------------------------------------------------------------------
c return the substring that starts at loc m in str
c skip initial spaces and tabs
c substring is any chars up to next space or tab
c also return m = loc of char after end of substr in str

      subroutine strread(str,m,substr)
      implicit none
      character*(*) str,substr
      integer m

      integer n1,n2,nmax
      
      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.
     $     (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' '
     $     .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      substr = str(n1:n2-1)
      m = n2

      return
      end
      
      
c -----------------------------------------------------------------------
c returns the actual length of str
c backtracks from end of string skipping over spaces

      integer function length(str)
      implicit none
      character*(*) str

      integer n

      n = len(str)
      do while (n.gt.0.and.str(n:n).eq.' ')
        n = n - 1
      enddo
      length = n

      return
      end
