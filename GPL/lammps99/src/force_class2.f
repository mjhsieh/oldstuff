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
c Class II routines authored by Eric Simon, Cray Research
c       Tue Jan 21 1997
c	- testing code removed
c	- newton and ifactor added
c	- all virial code commented out
c	- Wilson-out-of-plane forces added
c

c 3- and 4-body forces for class 2 force fields

c many-body forces
c  force components stored for atoms I own
c  variable force = (-1/r * d(phi)/dr)
c  variable virial = (-r * d(phi)/dr)

c ---------------------------------------------------------------------

c class 2 angles including
c	- quartic angle
c	- bond/bond coupling
c	- bond/angle coupling
c written by Eric J. Simon, Cray Research Inc.  Tue Jul 16 11:49:44 CDT 1996
c modified by ejs 1/21/97: removed test code, added newton/ifactor, commented
c			   out virial calculation
c
c  assumes array anglecoeff contains:
c	1	theta_0 	(equilibrium angle value)
c	2 	2_K_theta 	(coefficient for quadratic term)
c	3	3_K_theta	(coefficient for cubic term)
c	4	4_K_theta	(coefficient for quartic term)
c  assumes array bondbondcoeff contains:
c	1	K_b1/b2		(coefficient for bond_AB/bond_BC coupling
c	2	r1_0 		(equilibrium bond 1 value)
c	3	r2_0 		(equilibrium bond 2 value)
c  assumes array bondanglecoeff contains:
c	1	K_b1/theta	(coefficient for coupling of b1/angle)
c	2	K_b2/theta	(coefficient for coupling of b2/angle)
c	3	r1_0 		(equilibrium bond 1 value)
c	4	r2_0 		(equilibrium bond 2 value)

      subroutine angle_class2(iflag)
      include "lammps.h"
      parameter (zero=0.0,two=2.0,three=3.0,four=4.0,small=0.00001)
      real*8 fi(3),fj(3),fk(3)

      e_angle = 0.0
      e_bondbond = 0.0
      e_bondangle = 0.0
      
      do m = 1,nanglelocal
        i1 = anglelist(1,m)
        i2 = anglelist(2,m)
        i3 = anglelist(3,m)
        itype = anglelist(4,m)

c 1st bond
        delx1 = x(1,i1) - x(1,i2)
        dely1 = x(2,i1) - x(2,i2)
        delz1 = x(3,i1) - x(3,i2)
        call minimg(delx1,dely1,delz1)
        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1
        r1 = sqrt(rsq1)

c 2nd bond
        delx2 = x(1,i3) - x(1,i2)  
        dely2 = x(2,i3) - x(2,i2)
        delz2 = x(3,i3) - x(3,i2)
        call minimg(delx2,dely2,delz2)
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2
        r2 = sqrt(rsq2)

c angle (cos and sin)
        c = delx1*delx2 + dely1*dely2 + delz1*delz2
        c = c / (r1*r2)
        
        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0
        
        s = sqrt(1.0 - c*c)
        if (s.lt.small) s = small
        s = 1.0/s
        
c energy
        dtheta = acos(c) - anglecoeff(1,itype)
	dtheta2 = dtheta * dtheta
	dtheta3 = dtheta2 * dtheta
	dtheta4 = dtheta3 * dtheta

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.maxown) ifactor = ifactor + 1
          if (i2.le.maxown) ifactor = ifactor + 1
          if (i3.le.maxown) ifactor = ifactor + 1
        else
          ifactor = 3
        endif

	e_angle = e_angle + (ifactor/3.0 * 
     1      (anglecoeff(2,itype) * dtheta2 +
     2       anglecoeff(3,itype) * dtheta3 +
     3       anglecoeff(4,itype) * dtheta4))

	de_angle = two   * anglecoeff(2,itype) * dtheta +
     1       three * anglecoeff(3,itype) * dtheta2 +
     2       four  * anglecoeff(4,itype) * dtheta3

c force on each of 3 atoms
        a = de_angle * s
        a11 = a*c / rsq1
        a12 = -a / (r1*r2)
        a22 = a*c / rsq2
        
        fi(1) = a11*delx1 + a12*delx2
        fk(1) = a22*delx2 + a12*delx1
        fi(2) = a11*dely1 + a12*dely2
        fk(2) = a22*dely2 + a12*dely1
        fi(3) = a11*delz1 + a12*delz2
        fk(3) = a22*delz2 + a12*delz1

c----------------- bond/bond coupling -------------------------------------
        
        r1_0 = bondbondcoeff(2,itype)
        r2_0 = bondbondcoeff(3,itype)
        dr1 = (r1 - r1_0)
        dr2 = (r2 - r2_0)
        tk1 = bondbondcoeff(1,itype) * dr1
        tk2 = bondbondcoeff(1,itype) * dr2

        e_bondbond = e_bondbond + (ifactor/3.0 * 
     1               (bondbondcoeff(1,itype)*dr1*dr2))
        
        fi(1) = fi(1) + delx1*tk2/r1
        fi(2) = fi(2) + dely1*tk2/r1
        fi(3) = fi(3) + delz1*tk2/r1

        fk(1) = fk(1) + delx2*tk1/r2
        fk(2) = fk(2) + dely2*tk1/r2
        fk(3) = fk(3) + delz2*tk1/r2

c----------------- bond/angle coupling -------------------------------------

	e_bondangle = e_bondangle + (ifactor/3.0 * 
     1       ( (bondanglecoeff(1,itype)*dr1*dtheta) +
     2       (bondanglecoeff(2,itype)*dr2*dtheta) ))

	aa1 = s * dr1 * bondanglecoeff(1,itype)
	aa2 = s * dr2 * bondanglecoeff(2,itype)

	aa11 = aa1 * c / rsq1
	aa12 = -aa1 / (r1 * r2)
	aa21 = aa2 * c / rsq1
	aa22 = -aa2 / (r1 * r2)

	vx11 = (aa11 * delx1) + (aa12 * delx2)
	vx12 = (aa21 * delx1) + (aa22 * delx2)
	vy11 = (aa11 * dely1) + (aa12 * dely2)
	vy12 = (aa21 * dely1) + (aa22 * dely2)
	vz11 = (aa11 * delz1) + (aa12 * delz2)
	vz12 = (aa21 * delz1) + (aa22 * delz2)

	aa11 = aa1 * c / rsq2
	aa21 = aa2 * c / rsq2

	vx21 = (aa11 * delx2) + (aa12 * delx1)
	vx22 = (aa21 * delx2) + (aa22 * delx1)
	vy21 = (aa11 * dely2) + (aa12 * dely1)
	vy22 = (aa21 * dely2) + (aa22 * dely1)
	vz21 = (aa11 * delz2) + (aa12 * delz1)
	vz22 = (aa21 * delz2) + (aa22 * delz1)

	b1 = bondanglecoeff(1,itype) * dtheta / r1
	b2 = bondanglecoeff(2,itype) * dtheta / r2

        fi(1) = fi(1) + (vx11 + (b1 * delx1) + vx12) 
        fi(2) = fi(2) + (vy11 + (b1 * dely1) + vy12)
        fi(3) = fi(3) + (vz11 + (b1 * delz1) + vz12)

        fk(1) = fk(1) + (vx21 + (b2 * delx2) + vx22)
        fk(2) = fk(2) + (vy21 + (b2 * dely2) + vy22)
        fk(3) = fk(3) + (vz21 + (b2 * delz2) + vz22)

C-------------------------------------------------------
C
C  Forces on J are found by translational symmetry from
C  the forces on I and K then all the forces are scattered
C  into the overall force array and the virial
C

        if (newton_bond.eq.1.or.i1.le.maxown) then
          f(1,i1) = f(1,i1) - fi(1)
          f(2,i1) = f(2,i1) - fi(2)
          f(3,i1) = f(3,i1) - fi(3)
	endif
        
        if (newton_bond.eq.1.or.i2.le.maxown) then
          f(1,i2) = f(1,i2) + fi(1) + fk(1)
          f(2,i2) = f(2,i2) + fi(2) + fk(2)
          f(3,i2) = f(3,i2) + fi(3) + fk(3)
	endif
        
        if (newton_bond.eq.1.or.i3.le.maxown) then
          f(1,i3) = f(1,i3) - fk(1)
          f(2,i3) = f(2,i3) - fk(2)
          f(3,i3) = f(3,i3) - fk(3)
	endif

        if (iflag.eq.1) then

          virial(1) = virial(1) - (ifactor/3.0)*(delx1*fi(1) + 
     1                                           delx2*fk(1))
          virial(2) = virial(2) - (ifactor/3.0)*(dely1*fi(2) + 
     1                                           dely2*fk(2))
          virial(3) = virial(3) - (ifactor/3.0)*(delz1*fi(3) + 
     1                                           delz2*fk(3))
          virial(4) = virial(4) - (ifactor/3.0)*(delx1*fi(2) + 
     1                                           delx2*fk(2))
          virial(5) = virial(5) - (ifactor/3.0)*(delx1*fi(3) + 
     1                                           delx2*fk(3))
          virial(6) = virial(6) - (ifactor/3.0)*(dely1*fi(3) + 
     1                                           dely2*fk(3))
        endif

      enddo

      return
      end

c class 2 dihedral torsion including
c	- quartic dihedral
c	- bond/torsion coupling (divided into midbond and endbond)
c	- angle/torsion coupling
c	- angle/angle/torsion coupling
c written by Eric J. Simon, Cray Research Inc.  Tue Jul 16 11:49:44 CDT 1996
c modified by ejs 1/21/97: removed test code, added newton/ifactor, commented
c			   out virial calculation
c
c dihedcoeff contains:
c	1	V1		force constant for the 1st term
c	2	Phi1_0		equilibrium dihedral value for 1st term
c	3	V2		force constant for the 2nd term
c	4	Phi2_0		equilibrium dihedral value for 2nd term
c	5	V3		force constant for the 3rd term
c	6	Phi3_0		equilibrium dihedral value for 3rd term
c
c midbondtorsioncoeff contains:
c	1	K1_BC		force constant 1 for mid bond BC
c	2	K2_BC		force constant 2 for mid bond BC
c	3	K3_BC		force constant 3 for mid bond BC
c	4	b0_BC		equilibrium bond length for mid bond BC
c
c endbondtorsioncoeff contains:
c	1	K1_AB		force constant 1 for end bond AB
c	2	K2_AB		force constant 2 for end bond AB
c	3	K3_AB		force constant 3 for end bond AB
c	4	b0_AB		equilibrium bond length for end bond AB
c	5	K1_CD		force constant 1 for end bond CD
c	6	K2_CD		force constant 2 for end bond CD
c	7	K3_CD		force constant 3 for end bond CD
c	8	b0_CD		equilibrium bond length for end bond CD
c
c angletorsioncoeff contains:
c	1	K1_ABC		force constant 1 for 1st angle ABC
c	2	K2_ABC		force constant 2 for 1st angle ABC
c	3	K3_ABC		force constant 3 for 1st angle ABC
c	4	Theta0_ABC	equilibrium angle value for 1st angle ABC
c	5	K1_BCD		force constant 1 for 2nd angle BCD
c	6	K2_BCD		force constant 2 for 2nd angle BCD
c	7	K3_BCD		force constant 3 for 2nd angle BCD
c	8	Theta0_BCD	equilibrium angle value for 2nd angle BCD
c
c angleangletorsioncoeff contains:
c	1	K		force constant
c	2	Theta1_0 	equilibrium angle value for 1st angle
c	3	Theta2_0 	equilibrium angle value for 2nd angle
c
c The real*8 array dphidr(i,j) holds the force derivative d(phi)/d(r)
c	for each coordinate i (X,Y,Z for i=1,2,3 respectively)
c	in each atom j (a,b,c,d for j=1,2,3,4 respectively)
c The real*8 array dbonddr(i,j,k) holds the force derivative d(bond)/d(r)
c	for each coordinate i (X,Y,Z) in each atom j (a,b,c,d) for each
c	bond (AB,BC,CD)
c The real*8 array dthetadr(i,j,k) holds the force derivatrive d(theta)/d(r)
c	for each coordinate i (X,Y,Z) in each atom j (a,b,c,d) for each
c	angle (ABC,BCD)

      subroutine dihedral_class2(iflag)
      include "lammps.h"
      parameter (zero=0.0,one=1.0,two=2.0,three=3.0,small=0.0000001)
      real*8 dcosphidr(3,4)
      real*8 dphidr(3,4)
      real*8 dbonddr(3,4,3)
      real*8 dthetadr(3,4,2)
      real*8 fabcd(3,4)
      integer iindex(4)
      
      tolerance = 0.05

      e_dihedral = 0.0
      e_midbondtorsion = 0.0
      e_endbondtorsion = 0.0
      e_angletorsion = 0.0
      e_angleangletorsion = 0.0
      e_bondbond13 = 0.0

      do m = 1,ndihedlocal
          i1 = dihedlist(1,m)
          i2 = dihedlist(2,m)
          i3 = dihedlist(3,m)
          i4 = dihedlist(4,m)
          iindex(1) = i1
          iindex(2) = i2
          iindex(3) = i3
          iindex(4) = i4
          itype = dihedlist(5,m)

          delx1 = x(1,i1) - x(1,i2)
          dely1 = x(2,i1) - x(2,i2)
          delz1 = x(3,i1) - x(3,i2)
          call minimg(delx1,dely1,delz1)

c 2nd bond in dihedral
          delx2 = x(1,i3) - x(1,i2)
          dely2 = x(2,i3) - x(2,i2)
          delz2 = x(3,i3) - x(3,i2)
          call minimg(delx2,dely2,delz2)
          delx2m = -1.0 * delx2
          dely2m = -1.0 * dely2
          delz2m = -1.0 * delz2
          call minimg(delx2m,dely2m,delz2m)

c 3rd bond in dihedral
          delx3 = x(1,i4) - x(1,i3)
          dely3 = x(2,i4) - x(2,i3)
          delz3 = x(3,i4) - x(3,i3)
          call minimg(delx3,dely3,delz3)
          
          r1mag2 = delx1*delx1 + dely1*dely1 + delz1*delz1
          r1 = sqrt(r1mag2)
          r2mag2 = delx2*delx2 + dely2*dely2 + delz2*delz2
          r2 = sqrt(r2mag2)
          r3mag2 = delx3*delx3 + dely3*dely3 + delz3*delz3
          r3 = sqrt(r3mag2)

          sb1 = one / r1mag2
	  rb1 = one / r1
          
          sb2 = one / r2mag2
	  rb2 = one / r2

          sb3 = one / r3mag2
	  rb3 = one / r3

          c0 = (delx1*delx3 + dely1*dely3 + delz1*delz3) * rb1*rb3

c angles

          r12c1 = rb1*rb2
          r12c2 = rb2*rb3
c          r12c1 = one / (r1 * r2)
c          r12c2 = one / (r2 * r3)
          costh12 = (delx1*delx2 + dely1*dely2 + delz1*delz2)*
     1         r12c1
          costh13 = c0
c          costh13 = (delx1*delx3 + dely1*dely3 + delz1*delz3)*
c     1         (rb1*rb3)
          costh23 = (delx2m*delx3 + dely2m*dely3 + delz2m*delz3)*
     1         r12c2
          
c cos and sin of 2 angles and final c
          sc1 = sqrt(one - costh12*costh12)
          if (sc1.lt.small) sc1 = small
          sc1 = one/sc1
          
          sc2 = sqrt(one - costh23*costh23)
          if (sc2.lt.small) sc2 = small
          sc2 = one/sc2
          
          s1 = sc1 * sc1
          s2 = sc2 * sc2
          s12 = sc1 * sc2
          c = (c0 + costh12*costh23) * s12
          
c error check
          if (c.gt.1.0+tolerance.or.c.lt.(-1.0-tolerance)) then
            write (6,*) 'Dihedral:',node,ntimestep,
     $           tag(i1),tag(i2),tag(i3),tag(i4)
            write (6,*) ' Error1:',node,x(1,i1),x(2,i1),x(3,i1)
            write (6,*) ' Error2:',node,x(1,i2),x(2,i2),x(3,i2)
            write (6,*) ' Error3:',node,x(1,i3),x(2,i3),x(3,i3)
            write (6,*) ' Error4:',node,x(1,i4),x(2,i4),x(3,i4)
          endif

          if (c.gt.1.0) c = 1.0
          if (c.lt.-1.0) c = -1.0
          cosphi = c
          phi = acos(c)
c	  if (phi .eq. 0.0) phi = small
          sinphi = sqrt(1.0 - c*c)
          sinphi = max(sinphi,small)
          

C bug check

          a11 = (-c*sb1*s1)
          a22 = sb2*(2.0*costh13*s12 - c*(s1+s2))
          a33 = (-c*sb3*s2)
          a12 = r12c1*(costh12*c*s1 + costh23*s12)
          a13 = rb1*rb3*s12
          a23 = r12c2*(-costh23*c*s2 - costh12*s12)
          
          sx1  = a11*delx1 + a12*delx2 + a13*delx3
          sx2  = a12*delx1 + a22*delx2 + a23*delx3
          sx12 = a13*delx1 + a23*delx2 + a33*delx3
          sy1  = a11*dely1 + a12*dely2 + a13*dely3
          sy2  = a12*dely1 + a22*dely2 + a23*dely3
          sy12 = a13*dely1 + a23*dely2 + a33*dely3
          sz1  = a11*delz1 + a12*delz2 + a13*delz3
          sz2  = a12*delz1 + a22*delz2 + a23*delz3
          sz12 = a13*delz1 + a23*delz2 + a33*delz3

c set up d(cos(phi))/d(r) and dphi/dr arrays

          dcosphidr(1,1) = -1*sx1
          dcosphidr(2,1) = -1*sy1
          dcosphidr(3,1) = -1*sz1
          dcosphidr(1,2) = (sx2 + sx1)
          dcosphidr(2,2) = (sy2 + sy1)
          dcosphidr(3,2) = (sz2 + sz1)
          dcosphidr(1,3) = (sx12 - sx2)
          dcosphidr(2,3) = (sy12 - sy2)
          dcosphidr(3,3) = (sz12 - sz2)
          dcosphidr(1,4) = -1*sx12
          dcosphidr(2,4) = -1*sy12
          dcosphidr(3,4) = -1*sz12

          do j=1,4
           do i=1,3
            dphidr(i,j) = -dcosphidr(i,j)/sinphi
           end do
          end do

c energy
          dphi1 = phi - dihedcoeff(2,itype)
          dphi2 = two*phi - dihedcoeff(4,itype)
          dphi3 = three*phi - dihedcoeff(6,itype)

          if (newton_bond.eq.0) then
            ifactor = 0
            if (i1.le.maxown) ifactor = ifactor + 1
            if (i2.le.maxown) ifactor = ifactor + 1
            if (i3.le.maxown) ifactor = ifactor + 1
            if (i4.le.maxown) ifactor = ifactor + 1
          else
            ifactor = 4
         endif

        e_dihedral = e_dihedral + (ifactor/4.0 * (
     1         dihedcoeff(1,itype)*(one - cos(dphi1)) +
     2         dihedcoeff(3,itype)*(one - cos(dphi2)) +
     3         dihedcoeff(5,itype)*(one - cos(dphi3))))
          
          de_dihedral = dihedcoeff(1,itype)*sin(dphi1) +
     1         two*dihedcoeff(3,itype)*sin(dphi2) +
     2         three*dihedcoeff(5,itype)*sin(dphi3)

c force on all 4 atoms
c torsion forces
          do j = 1, 4
              do i = 1,3
                fabcd(i,j) = de_dihedral*dphidr(i,j)
              enddo
          enddo

c--------------------------------------------------------------------------
c set up d(bond)/d(r) array, dbonddr(i,j,k) = coordinate i, atom j, bond k
          do i = 1, 3
            do j = 1,4
	      do k = 1,3
		dbonddr(i,j,k) = zero
	      enddo
            enddo
          enddo

c bond1:
          dbonddr(1,1,1) = delx1 / r1
          dbonddr(2,1,1) = dely1 / r1
          dbonddr(3,1,1) = delz1 / r1
          dbonddr(1,2,1) = -1.0 * delx1 / r1
          dbonddr(2,2,1) = -1.0 * dely1 / r1
          dbonddr(3,2,1) = -1.0 * delz1 / r1
c bond2:
          dbonddr(1,2,2) = delx2 / r2
          dbonddr(2,2,2) = dely2 / r2
          dbonddr(3,2,2) = delz2 / r2
          dbonddr(1,3,2) = -1.0*delx2 / r2
          dbonddr(2,3,2) = -1.0*dely2 / r2
          dbonddr(3,3,2) = -1.0*delz2 / r2
c bond3:
          dbonddr(1,3,3) = delx3 / r3
          dbonddr(2,3,3) = dely3 / r3
          dbonddr(3,3,3) = delz3 / r3
          dbonddr(1,4,3) = -1.0 * delx3 / r3
          dbonddr(2,4,3) = -1.0 * dely3 / r3
          dbonddr(3,4,3) = -1.0 * delz3 / r3

c--------------------------------------------------------------------------
c set up d(theta)/d(r) array, dthetadr(i,j,k) = coordinate i, atom j, angle k

          do i = 1, 3
            do j = 1,4
	      do k = 1,2
		dthetadr(i,j,k) = zero
	      enddo
            enddo
          enddo

          t1 = costh12 / r1mag2
          t2 = costh23 / r2mag2
          t3 = costh12 / r2mag2
          t4 = costh23 / r3mag2

c angle12:
          dthetadr(1,1,1) = sc1 * ((t1 * delx1) - (delx2 * r12c1)) 
          dthetadr(2,1,1) = sc1 * ((t1 * dely1) - (dely2 * r12c1)) 
          dthetadr(3,1,1) = sc1 * ((t1 * delz1) - (delz2 * r12c1)) 

          dthetadr(1,2,1) = sc1 *
     1         ((-1.0*t1 * delx1) + (delx2 * r12c1) +
     2         (-1.0*t3 * delx2) + (delx1 * r12c1))
          dthetadr(2,2,1) = sc1 *
     1         ((-1.0 * t1 * dely1) + (dely2 * r12c1) +
     2         (-1.0 * t3 * dely2) + (dely1 * r12c1))
          dthetadr(3,2,1) = sc1 *
     1         ((-1.0 * t1 * delz1) + (delz2 * r12c1) +
     2         (-1.0 * t3 * delz2) + (delz1 * r12c1))

          dthetadr(1,3,1) = sc1 * ((t3 * delx2) - (delx1 * r12c1)) 
          dthetadr(2,3,1) = sc1 * ((t3 * dely2) - (dely1 * r12c1)) 
          dthetadr(3,3,1) = sc1 * ((t3 * delz2) - (delz1 * r12c1)) 

c angle23:
          dthetadr(1,2,2) = sc2 * ((t2 * delx2) + (delx3 * r12c2)) 
          dthetadr(2,2,2) = sc2 * ((t2 * dely2) + (dely3 * r12c2)) 
          dthetadr(3,2,2) = sc2 * ((t2 * delz2) + (delz3 * r12c2)) 
          
          dthetadr(1,3,2) = sc2 *
     1         ((-1.0*t2 * delx2) - (delx3 * r12c2) +
     2         (t4 * delx3) + (delx2 * r12c2))
          dthetadr(2,3,2) = sc2 *
     1         ((-1.0 * t2 * dely2) - (dely3 * r12c2) +
     2         (t4 * dely3) + (dely2 * r12c2))
          dthetadr(3,3,2) = sc2 *
     1         ((-1.0 * t2 * delz2) - (delz3 * r12c2) +
     2         (t4 * delz3) + (delz2 * r12c2))

          dthetadr(1,4,2) = -1.0*sc2 * ((t4 * delx3) + (delx2 * r12c2))
          dthetadr(2,4,2) = -1.0*sc2 * ((t4 * dely3) + (dely2 * r12c2))
          dthetadr(3,4,2) = -1.0*sc2 * ((t4 * delz3) + (delz2 * r12c2))

c----------- bond / torsion coupling ------------------------------

c mid-bond/torsion coupling

c energy on bond2
          cos2phi = cos(two*phi)
          cos3phi = cos(three*phi)

          bt1 = midbondtorsioncoeff(1,itype) * cosphi
          bt2 = midbondtorsioncoeff(2,itype) * cos2phi
          bt3 = midbondtorsioncoeff(3,itype) * cos3phi
          sumbte = bt1 + bt2 + bt3
          db = r2 - midbondtorsioncoeff(4,itype)
          e_midbondtorsion = e_midbondtorsion + (ifactor/4.0 *
     1		(db * sumbte))

c force on bond2
          bt1 = -1.0 * midbondtorsioncoeff(1,itype) * sinphi
          bt2 = -2.0 * midbondtorsioncoeff(2,itype) * sin(two*phi)
          bt3 = -3.0 * midbondtorsioncoeff(3,itype) * sin(three*phi)
          sumbtf = bt1 + bt2 + bt3

          do j = 1, 4
              do i = 1, 3
                fabcd(i,j) = fabcd(i,j) +
     1             (db * sumbtf * dphidr(i,j) + 
     2                   sumbte * dbonddr(i,j,2)     )
              enddo
          enddo

c end-bond/torsion coupling

c energy on bond1
          bt1 = endbondtorsioncoeff(1,itype) * cosphi
          bt2 = endbondtorsioncoeff(2,itype) * cos2phi
          bt3 = endbondtorsioncoeff(3,itype) * cos3phi
          sumbte = bt1 + bt2 + bt3

          db = r1 - endbondtorsioncoeff(7,itype)
          e_endbondtorsion = e_endbondtorsion + (ifactor/4.0 *
     1				(db * (bt1+bt2+bt3)))

          bt1 = 1.0 * endbondtorsioncoeff(1,itype) * sinphi
          bt2 = 2.0 * endbondtorsioncoeff(2,itype) * sin(two*phi)
          bt3 = 3.0 * endbondtorsioncoeff(3,itype) * sin(three*phi)
          sumbtf = bt1 + bt2 + bt3

c force on bond1
          do j = 1, 4
              do i = 1, 3
	        fabcd(i,j) = fabcd(i,j) -
     1             (db * sumbtf * dphidr(i,j)   + 
     2                   sumbte * dbonddr(i,j,1) )
              enddo
          enddo

c energy on bond3

          bt1 = endbondtorsioncoeff(4,itype) * cosphi
          bt2 = endbondtorsioncoeff(5,itype) * cos2phi
          bt3 = endbondtorsioncoeff(6,itype) * cos3phi
          sumbte = bt1 + bt2 + bt3

          db = r3 - endbondtorsioncoeff(8,itype)
          e_endbondtorsion = e_endbondtorsion + (ifactor/4.0 *
     1                          (db * (bt1+bt2+bt3)))

          bt1 = -1.0 * endbondtorsioncoeff(4,itype) * sinphi
          bt2 = -2.0 * endbondtorsioncoeff(5,itype) * sin(two*phi)
          bt3 = -3.0 * endbondtorsioncoeff(6,itype) * sin(three*phi)
          sumbtf = bt1 + bt2 + bt3

c force on bond3
          do j = 1, 4
              do i = 1, 3
	        fabcd(i,j) = fabcd(i,j) +
     1             (db * sumbtf * dphidr(i,j)   +
     2                   sumbte * dbonddr(i,j,3) )
              enddo
          enddo

c
c----------- angle / torsion coupling --------------------------

c angle/torsion coupling

c energy on angle1
          at1 = angletorsioncoeff(1,itype) * cosphi
          at2 = angletorsioncoeff(2,itype) * cos2phi
          at3 = angletorsioncoeff(3,itype) * cos3phi
          sumbte = at1 + at2 + at3

          da = acos(costh12) - angletorsioncoeff(7,itype)
          e_angletorsion = e_angletorsion + (ifactor/4.0 *
     1                          (da * (at1+at2+at3)))

          bt1 = 1.0 * angletorsioncoeff(1,itype) * sinphi
          bt2 = 2.0 * angletorsioncoeff(2,itype) * sin(two*phi)
          bt3 = 3.0 * angletorsioncoeff(3,itype) * sin(three*phi)
          sumbtf = bt1 + bt2 + bt3

c force on angle1
          do j = 1, 4
              do i = 1, 3
                fabcd(i,j) = fabcd(i,j) - 
     1                       (da * sumbtf * dphidr(i,j) + 
     2                             sumbte * dthetadr(i,j,1) )
              enddo
          enddo

c energy on angle2
          at1 = angletorsioncoeff(4,itype) * cosphi
          at2 = angletorsioncoeff(5,itype) * cos2phi
          at3 = angletorsioncoeff(6,itype) * cos3phi
          sumbte = at1 + at2 + at3

          da = acos(costh23) - angletorsioncoeff(8,itype)
          e_angletorsion = e_angletorsion + (ifactor/4.0 *
     1                          (da * (at1+at2+at3)))

          bt1 = -1.0 * angletorsioncoeff(4,itype) * sinphi
          bt2 = -2.0 * angletorsioncoeff(5,itype) * sin(two*phi)
          bt3 = -3.0 * angletorsioncoeff(6,itype) * sin(three*phi)
          sumbtf = bt1 + bt2 + bt3

c force on angle2
          do j = 1, 4
              do i = 1, 3
                fabcd(i,j) = fabcd(i,j) + 
     1                       (da * sumbtf * dphidr(i,j) + 
     2                             sumbte * dthetadr(i,j,2) )
              enddo
          enddo


          
c----------- angle / angle / torsion coupling -----------------------

          da1 = acos(costh12) - angleangletorsioncoeff(2,itype)
          da2 = acos(costh23) - angleangletorsioncoeff(3,itype)
          
c     energy ejss
          e_angleangletorsion = e_angleangletorsion + 
     1          (ifactor/4.0 * 
     1      (angleangletorsioncoeff(1,itype)*da1*da2*cosphi))

c     force
          do j = 1, 4
              do i = 1, 3
                fabcd(i,j) = fabcd(i,j) - 
     1             angleangletorsioncoeff(1,itype) * (
     2             cosphi * (da2 * dthetadr(i,j,1) - 
     3                       da1 * dthetadr(i,j,2)) +
     4             sinphi * da1 * da2 * dphidr(i,j) )
              enddo
          enddo

c----------- Bond1 / Bond3 coupling -----------------------

        if (abs(bondbond13coeff(1,itype)).gt.small) then

        r1_0 = bondbond13coeff(2,itype)
        r3_0 = bondbond13coeff(3,itype)
        dr1 = (r1 - r1_0)
        dr2 = (r3 - r3_0)
        tk1 = -1.0 * bondbond13coeff(1,itype) * dr1 / r3
        tk2 = -1.0 * bondbond13coeff(1,itype) * dr2 / r1

        e_bondbond13 = e_bondbond13 + (ifactor/4.0 * 
     1               (bondbond13coeff(1,itype)*dr1*dr2))
        
        fabcd(1,1) = fabcd(1,1) + tk2 * delx1
        fabcd(2,1) = fabcd(2,1) + tk2 * dely1
        fabcd(3,1) = fabcd(3,1) + tk2 * delz1
      
        fabcd(1,2) = fabcd(1,2) - tk2 * delx1
        fabcd(2,2) = fabcd(2,2) - tk2 * dely1
        fabcd(3,2) = fabcd(3,2) - tk2 * delz1
        
        fabcd(1,3) = fabcd(1,3) - tk1 * delx3
        fabcd(2,3) = fabcd(2,3) - tk1 * dely3
        fabcd(3,3) = fabcd(3,3) - tk1 * delz3

        fabcd(1,4) = fabcd(1,4) + tk1 * delx3
        fabcd(2,4) = fabcd(2,4) + tk1 * dely3
        fabcd(3,4) = fabcd(3,4) + tk1 * delz3

        end if

c------------------------------------------------------------
c
c Scatter the forces for this torsion into the force array
c
          do j = 1, 4
            if (newton_bond.eq.1.or.iindex(j).le.maxown) then
              do i = 1, 3
	        f(i,iindex(j)) = f(i,iindex(j)) + fabcd(i,j)
              enddo
	    endif
          enddo
c
c Compute the virial contribution
c
        if (iflag.eq.1) then

          vx1 = fabcd(1,1)
          vx2 = fabcd(1,3) + fabcd(1,4)
          vx3 = fabcd(1,4)

          vy1 = fabcd(2,1)
          vy2 = fabcd(2,3) + fabcd(2,4)
          vy3 = fabcd(2,4)

          vz1 = fabcd(3,1)
          vz2 = fabcd(3,3) + fabcd(3,4)
          vz3 = fabcd(3,4)


          virial(1) = virial(1) + (ifactor/4.0)*
     1                    (delx1*vx1 + delx2*vx2 + delx3*vx3)
          virial(2) = virial(2) + (ifactor/4.0)*
     1                    (dely1*vy1 + dely2*vy2 + dely3*vy3)
          virial(3) = virial(3) + (ifactor/4.0)*
     1                    (delz1*vz1 + delz2*vz2 + delz3*vz3)
          virial(4) = virial(4) + (ifactor/4.0)*
     1                    (delx1*vy1 + delx2*vy2 + delx3*vy3)
          virial(5) = virial(5) + (ifactor/4.0)*
     1                    (delx1*vz1 + delx2*vz2 + delx3*vz3)
          virial(6) = virial(6) + (ifactor/4.0)*
     1                    (dely1*vz1 + dely2*vz2 + dely3*vz3)

        endif

      enddo

      return
      end

c  assumes array anglecoeff contains:
c class 2 out-of-plane improper 
c written by Eric J. Simon, Cray Research, Inc. Jul/31/1996-Jan/18/1997
c modified by ejs 1/21/97: removed test code, added newton/ifactor, commented
c			   out virial calculation, added forces
c
c for atoms ABCD, it is assumed that atom B is the central atom in the
c       trigonal form
c
c chi is calculated as the average of three chi values ABCD, CBAD, DBAC
c
c   sin_X(ABCD) =        r_CB cross r_DB                r_AB
c                -------------------------------  dot  -------
c                 sin_theta(CBD) |r_CB| |r_DB|          |r_AB|
c 
c The energy is calculated as:  
c    e_improper = K_chi * (chi chi_0)**2 
c
c impcoeff contains:
c       1       K               coefficient for chi 
c       2       chi_0           equilibrium value for chi 
c
      subroutine improper_class2(iflag)
      include "lammps.h"
      parameter (zero=0.0,one=1.0,two=2.0,three=3.0,
     1     four=4.0,negone=-1.0)

      integer iindex(4)
c matrix to hold the coordinates, where i=X/Y/Z, j=A/B/C/D
      real*8 r(3,4)
c difference vectors between coordinates, i=X/Y/Z, j=AB/CB/BD
      real*8 delr(3,3)
c magnitude of bond lengths (i=AB/CB/BD), inverse, and 2nd and 3rd power
      real*8 rmag(3), rinvmag(3), rmag2(3), rmag3(3)
c matrix to hold theta, cos(theta), sin(theta), i=ABC/CBD/ABD
      real*8 theta(3), costheta(3), sintheta(3), cossqtheta(3),
     1     sinsqtheta(3), invstheta(3)
c matrix that represtents the cross product vectors, where 1=X, 2=Y, 3=Z
      real*8 rABxrCB(3), rDBxrAB(3), rCBxrDB(3)
c derivatives of delr, where i=atomA/B/C/D, j=bond AB/CB/DB
      real*8 ddelr(4,3)
c derivaties of r and 1/r, where i=X/Y/Z, j=atomA/B/C/D, k=bond AB/CB/DB
      real*8 dr(3,4,3), dinvr(3,4,3)
c derivatives of theta, where i=X/Y/Z, j=atomA/B/C/D, k=thetaABC/CBD/ABD
      real*8 dthetadr(3,4,3)
c derivative of 1/sin(theta), where i=X/Y/Z, j=atomA/B/C/D, k=thetaABC/CBD/ABD
      real*8 dinvsth(3,4,3)
c derivatives of 1 / (|r_AB| * |r_CB| * |r_DB|)
      real*8 dinv3r(3,4), dinvs3r(3,4,3)
      real*8 drCBxrDB(3), rCBxdrDB(3), dot1, dot2, dd(3)
      real*8 drDBxrAB(3), rDBxdrAB(3)
      real*8 drABxrCB(3), rABxdrCB(3)
      real*8 fdot(3,4,3), fdot2(3,4,3)
      real*8 f1(3,4,3), f2(3,4,3)
      real*8 invs3r(3), inv3r
      real*8 tt1, tt3, sc1
      real*8 dotCBDBAB, dotDBABCB, dotABCBDB
      real*8 chi, deltachi, d2chi, cossin2
      real*8 drAB(3,4,3), drCB(3,4,3), drDB(3,4,3)
      real*8 dchi(3,4,3), dtotalchi(3,4)
      real*8 schiABCD, chiABCD
      real*8 schiCBDA, chiCBDA
      real*8 schiDBAC, chiDBAC
      real*8 fabcd(3,4)

      data r /12*0.0/
      data delr, ddelr /21*0.0/
      data rmag, rinvmag, rmag2, rmag3 /12*0.0/
      data theta, costheta, sintheta, cossqtheta, sinsqtheta,
     1     invstheta /18*0.0/

      data rABxrCB /3*0.0/
      data rDBxrAB /3*0.0/
      data rCBxrDB /3*0.0/
      data dr, dinvr /72*0.0/
      data dthetadr, dinvsth /72*0.0/
      data dinv3r, dinvs3r /48*0.0/
      data f1, f2 /72*0.0/
      data drCBxrDB, rCBxdrDB, dot1, dot2, dd /11*0.0/
      data drDBxrAB, rDBxdrAB /6*0.0/
      data drABxrCB, rABxdrCB /6*0.0/
      data fdot, fdot2 /72*0.0/
      data sc1, dotCBDBAB, dotDBABCB, dotABCBDB /4*0.0/
      data schiABCD, chiABCD /2*0.0/
      data schiCBDA, chiCBDA /2*0.0/
      data schiDBAC, chiDBAC /2*0.0/
      data cossin2, chi, deltachi, d2chi /4*0.0/
      data inv3r, invs3r /4*0.0/
      data drAB, drCB, drDB /108*0.0/
      data tt1, tt3 /2*0.0/
      data dchi, dtotalchi /48*0.0/

      e_improper = 0.0

      do m = 1,nimprolocal
        iindex(1) = improlist(1,m)
        iindex(2) = improlist(2,m)
        iindex(3) = improlist(3,m)
        iindex(4) = improlist(4,m)
        itype = improlist(5,m)

        if (improcoeff(1,itype).ne.0.0) then

c gather coordinates
        do i = 1,3
          do j = 1,4
            r(i,j) = x(i,iindex(j))
          enddo
        enddo

c calculate difference vectors
        do i = 1, 3
          delr(i,1) = r(i,1) - r(i,2)
          delr(i,2) = r(i,3) - r(i,2)
          delr(i,3) = r(i,4) - r(i,2)
        enddo
        do i = 1, 3
          call minimg(delr(1,i),delr(2,i),delr(3,i))
        enddo

c calculate bond lengths and associates values
        do i = 1, 3
          rmag2(i) = delr(1,i)*delr(1,i) + 
     1         delr(2,i)*delr(2,i) + delr(3,i)*delr(3,i)
          rmag(i) = sqrt(rmag2(i))
          rinvmag(i) = one / rmag(i)
          rmag3(i) = rmag2(i) * rmag(i)
        enddo

c angle ABC
        costheta(1) =(delr(1,1)*delr(1,2)+delr(2,1)*delr(2,2)+ 
     1       delr(3,1)*delr(3,2))/(rmag(1)*rmag(2))
c angle CBD
        costheta(2) =(delr(1,2)*delr(1,3)+delr(2,2)*delr(2,3)+ 
     2       delr(3,2)*delr(3,3))/(rmag(2)*rmag(3))
c angle ABD
        costheta(3) =(delr(1,1)*delr(1,3)+delr(2,1)*delr(2,3)+ 
     3       delr(3,1)*delr(3,3))/(rmag(1)*rmag(3))

c angle error check
        do i = 1, 3
          if (costheta(i) .eq. -1.0) then
            write (6,*) 'ERROR in improper: linear angle'
            write (6,*) 'Improper:',node,ntimestep,
     1            tag(iindex(1)),tag(iindex(2)),tag(iindex(3))
     2            ,tag(iindex(4))
            write (6,*) ' Error1:',node,x(1,iindex(1)),
     1           x(2,iindex(1)),
     2           x(3,iindex(1))
            write (6,*) ' Error2:',node,x(1,iindex(2)),x(2,iindex(2)),
     1           x(3,iindex(2))
            write (6,*) ' Error3:',node,x(1,iindex(3)),x(2,iindex(3)),
     1           x(3,iindex(3))
            write (6,*) ' Error4:',node,x(1,iindex(4)),x(2,iindex(4)),
     1           x(3,iindex(4))
          endif
        enddo

        do i = 1, 3
          if (costheta(i) .gt. one)  costheta(i) = one
          if (costheta(i) .lt. negone) costheta(i) = negone
          theta(i) = acos(costheta(i))
          cossqtheta(i) = costheta(i)*costheta(i)
          sintheta(i) = sin(theta(i))
          invstheta(i) = one / sintheta(i)
          sinsqtheta(i) = sintheta(i)*sintheta(i)
        enddo

c cross products
        call cross(delr(1,1),delr(2,1),delr(3,1),
     1       delr(1,2),delr(2,2),delr(3,2),
     2       rABxrCB(1),rABxrCB(2),rABxrCB(3))
        call cross(delr(1,3),delr(2,3),delr(3,3),
     1       delr(1,1),delr(2,1),delr(3,1),
     2       rDBxrAB(1),rDBxrAB(2),rDBxrAB(3))
        call cross(delr(1,2),delr(2,2),delr(3,2),
     1       delr(1,3),delr(2,3),delr(3,3),
     2       rCBxrDB(1),rCBxrDB(2),rCBxrDB(3))

c dot products
        call dot(rCBxrDB(1),rCBxrDB(2),rCBxrDB(3),
     1       delr(1,1),delr(2,1),delr(3,1),dotCBDBAB)
        call dot(rDBxrAB(1),rDBxrAB(2),rDBxrAB(3),
     1       delr(1,2),delr(2,2),delr(3,2),dotDBABCB)
        call dot(rABxrCB(1),rABxrCB(2),rABxrCB(3),
     1       delr(1,3),delr(2,3),delr(3,3),dotABCBDB)

        t = rmag(1) * rmag(2) * rmag(3)
        inv3r = one / t
        invs3r(1) =  invstheta(2) * inv3r
        invs3r(2) =  invstheta(3) * inv3r
        invs3r(3) =  invstheta(1) * inv3r

c chi ABCD
        schiABCD = dotCBDBAB * invs3r(1)
        chiABCD = asin(schiABCD)
c chi CBDA
        schiCBDA = dotDBABCB * invs3r(2)
        chiCBDA = asin(schiCBDA)
c chi DBAC
        schiDBAC = dotABCBDB * invs3r(3)
        chiDBAC = asin(schiDBAC)

c calculate final chi by averaging over three values
        chi = (chiABCD + chiCBDA + chiDBAC) / three
        deltachi = chi - improcoeff(2,itype)
        d2chi = deltachi * deltachi

        if (newton_bond.eq.0) then
          ifactor = 0
          if (iindex(1).le.maxown) ifactor = ifactor + 1
          if (iindex(2).le.maxown) ifactor = ifactor + 1
          if (iindex(3).le.maxown) ifactor = ifactor + 1
          if (iindex(4).le.maxown) ifactor = ifactor + 1
        else
          ifactor = 4
        endif
        e_improper = e_improper + (ifactor/four *
     1       (improcoeff(1,itype) * d2chi))

c forces

c define d(delr) for each value of delr, where i=atomA/B/C/D, j=bondAB/CB/DB
	do i = 1, 4
	  do j = 1, 3
	    ddelr(i,j) = zero
	  enddo
	enddo

        ddelr(1,1) = one
        ddelr(2,1) = negone
        ddelr(2,2) = negone
        ddelr(3,2) = one
        ddelr(2,3) = negone
        ddelr(4,3) = one

c calculate d(|r|)/dr and d(1/|r|)/dr for each direction, atom and bond
c define d(r) for each r, where i=X/Y/Z, j=atomA/B/C/D, k=bondAB/CB/DB
        do i = 1, 3
          do j = 1, 4
            do k = 1,3
                dr(i,j,k) = delr(i,k) * ddelr(j,k) / rmag(k)
                dinvr(i,j,k) = negone * dr(i,j,k) / rmag2(k)
            enddo
          enddo
        enddo

c calculate d(1 / (|r_AB| * |r_CB| * |r_DB|)/dr for each direction and atom
        do i = 1,3
          do j = 1,4
            dinv3r(i,j) = rinvmag(2) * (rinvmag(3) * dinvr(i,j,1) +
     1                                  rinvmag(1) * dinvr(i,j,3)) +
     2                    rinvmag(3) * rinvmag(1) * dinvr(i,j,2)
          enddo
        enddo

c calculate d(theta)/d(r)

c angleABC:
	tt1 = costheta(1) / rmag2(1)
	tt3 = costheta(1) / rmag2(2)
	sc1 = one / sqrt(one - cossqtheta(1))
        dthetadr(1,1,1) = sc1 * ((tt1 * delr(1,1)) - (delr(1,2) * 
     1                            rinvmag(1) * rinvmag(2)))
        dthetadr(2,1,1) = sc1 * ((tt1 * delr(2,1)) - (delr(2,2) * 
     1                            rinvmag(1) * rinvmag(2)))
        dthetadr(3,1,1) = sc1 * ((tt1 * delr(3,1)) - (delr(3,2) * 
     1                            rinvmag(1) * rinvmag(2)))
        dthetadr(1,2,1) = -1.0 * sc1 * ((tt1 * delr(1,1)) - (delr(1,2) * 
     1                            rinvmag(1) * rinvmag(2)) +
     2                           (tt3 * delr(1,2)) - (delr(1,1) * 
     3                            rinvmag(1) * rinvmag(2)))
        dthetadr(2,2,1) = -1.0 * sc1 * ((tt1 * delr(2,1)) - (delr(2,2) * 
     1                            rinvmag(1) * rinvmag(2)) +
     2                           (tt3 * delr(2,2)) - (delr(2,1) * 
     3                            rinvmag(1) * rinvmag(2)))
        dthetadr(3,2,1) = -1.0 * sc1 * ((tt1 * delr(3,1)) - (delr(3,2) * 
     1                            rinvmag(1) * rinvmag(2)) +
     2                           (tt3 * delr(3,2)) - (delr(3,1) * 
     3                            rinvmag(1) * rinvmag(2)))
        dthetadr(1,3,1) = sc1 * ((tt3 * delr(1,2)) - (delr(1,1) * 
     1                            rinvmag(1) * rinvmag(2)))
        dthetadr(2,3,1) = sc1 * ((tt3 * delr(2,2)) - (delr(2,1) * 
     1                            rinvmag(1) * rinvmag(2)))
        dthetadr(3,3,1) = sc1 * ((tt3 * delr(3,2)) - (delr(3,1) * 
     1                            rinvmag(1) * rinvmag(2)))

c angleCBD:
	tt1 = costheta(2) / rmag2(2)
	tt3 = costheta(2) / rmag2(3)
	sc1 = one / sqrt(one - cossqtheta(2))

        dthetadr(1,3,2) = sc1 * ((tt1 * delr(1,2)) - (delr(1,3) * 
     1                            rinvmag(2) * rinvmag(3)))
        dthetadr(2,3,2) = sc1 * ((tt1 * delr(2,2)) - (delr(2,3) * 
     1                            rinvmag(2) * rinvmag(3)))
        dthetadr(3,3,2) = sc1 * ((tt1 * delr(3,2)) - (delr(3,3) * 
     1                            rinvmag(2) * rinvmag(3)))
        dthetadr(1,2,2) = -1.0 * sc1 * ((tt1 * delr(1,2)) - (delr(1,3) * 
     1                            rinvmag(2) * rinvmag(3)) +
     2                           (tt3 * delr(1,3)) - (delr(1,2) * 
     3                            rinvmag(3) * rinvmag(2)))
        dthetadr(2,2,2) = -1.0 * sc1 * ((tt1 * delr(2,2)) - (delr(2,3) * 
     1                            rinvmag(2) * rinvmag(3)) +
     2                           (tt3 * delr(2,3)) - (delr(2,2) * 
     3                            rinvmag(3) * rinvmag(2)))
        dthetadr(3,2,2) = -1.0 * sc1 * ((tt1 * delr(3,2)) - (delr(3,3) * 
     1                            rinvmag(2) * rinvmag(3)) +
     2                           (tt3 * delr(3,3)) - (delr(3,2) * 
     3                            rinvmag(3) * rinvmag(2)))
        dthetadr(1,4,2) = sc1 * ((tt3 * delr(1,3)) - (delr(1,2) * 
     1                            rinvmag(3) * rinvmag(2)))
        dthetadr(2,4,2) = sc1 * ((tt3 * delr(2,3)) - (delr(2,2) * 
     1                            rinvmag(3) * rinvmag(2)))
        dthetadr(3,4,2) = sc1 * ((tt3 * delr(3,3)) - (delr(3,2) * 
     1                            rinvmag(3) * rinvmag(2)))

c angleABD:
	tt1 = costheta(3) / rmag2(1)
	tt3 = costheta(3) / rmag2(3)
	sc1 = one / sqrt(one - cossqtheta(3))
        dthetadr(1,1,3) = sc1 * ((tt1 * delr(1,1)) - (delr(1,3) * 
     1                            rinvmag(1) * rinvmag(3)))
        dthetadr(2,1,3) = sc1 * ((tt1 * delr(2,1)) - (delr(2,3) * 
     1                            rinvmag(1) * rinvmag(3)))
        dthetadr(3,1,3) = sc1 * ((tt1 * delr(3,1)) - (delr(3,3) * 
     1                            rinvmag(1) * rinvmag(3)))
        dthetadr(1,2,3) = -1.0 * sc1 * ((tt1 * delr(1,1)) - (delr(1,3) * 
     1                            rinvmag(1) * rinvmag(3)) +
     2                           (tt3 * delr(1,3)) - (delr(1,1) * 
     3                            rinvmag(3) * rinvmag(1)))
        dthetadr(2,2,3) = -1.0 * sc1 * ((tt1 * delr(2,1)) - (delr(2,3) * 
     1                            rinvmag(1) * rinvmag(3)) +
     2                           (tt3 * delr(2,3)) - (delr(2,1) * 
     3                            rinvmag(3) * rinvmag(1)))
        dthetadr(3,2,3) = -1.0 * sc1 * ((tt1 * delr(3,1)) - (delr(3,3) * 
     1                            rinvmag(1) * rinvmag(3)) +
     2                           (tt3 * delr(3,3)) - (delr(3,1) * 
     3                            rinvmag(3) * rinvmag(1)))
        dthetadr(1,4,3) = sc1 * ((tt3 * delr(1,3)) - (delr(1,1) * 
     1                            rinvmag(3) * rinvmag(1)))
        dthetadr(2,4,3) = sc1 * ((tt3 * delr(2,3)) - (delr(2,1) * 
     1                            rinvmag(3) * rinvmag(1)))
        dthetadr(3,4,3) = sc1 * ((tt3 * delr(3,3)) - (delr(3,1) * 
     1                            rinvmag(3) * rinvmag(1)))

c calculate d( 1 / sin(theta))/dr for each direction and atom and angle
	do k = 1,3
	  cossin2 = negone * costheta(k) / sinsqtheta(k)
	  do i = 1,3
	    do j = 1,4
		dinvsth(i,j,k) = cossin2 * dthetadr(i,j,k)
	    enddo
	  enddo
	enddo

c calculate d(1 / sin(theta) * |r_AB| * |r_CB| * |r_DB|)/dr for each direction,
c                                       atom and angle
	do i = 1,3
	  do j = 1,4
	    dinvs3r(i,j,1) = (invstheta(2) * dinv3r(i,j)) +
     1                       (inv3r * dinvsth(i,j,2))
	    dinvs3r(i,j,2) = (invstheta(3) * dinv3r(i,j)) +
     1                       (inv3r * dinvsth(i,j,3))
	    dinvs3r(i,j,3) = (invstheta(1) * dinv3r(i,j)) +
     1                       (inv3r * dinvsth(i,j,1))
	  enddo
	enddo

C------------- d( (r_CB x r_DB) dot r_AB) -----------------------------------

c drCB(i,j,k) where i=direction X/Y/Z, j=atom A/B/C/D, k=vector X'/Y'/Z'
	do i = 1, 3
		drCB(i,2,i) = negone
		drAB(i,2,i) = negone
		drDB(i,2,i) = negone
		drDB(i,4,i) = one
		drCB(i,3,i) = one
		drAB(i,1,i) = one
	enddo

	do i = 1,3
	  do j = 1,4
c r_CB x d(r_DB)
	      call cross (delr(1,2),delr(2,2),delr(3,2),
     1                  drDB(1,j,i), drDB(2,j,i), drDB(3,j,i),
     2			rCBxdrDB(1),rCBxdrDB(2),rCBxdrDB(3))
c d(r_CB) x r_DB
	      call cross (drCB(1,j,i), drCB(2,j,i), drCB(3,j,i),
     1                  delr(1,3),delr(2,3),delr(3,3),
     2			drCBxrDB(1),drCBxrDB(2),drCBxrDB(3))
		do k = 1, 3
c (r_CB x d(r_DB)) + (d(r_CB) x r_DB)
		  dd(k) = rCBxdrDB(k) + drCBxrDB(k)
		enddo
c (r_CB x d(r_DB)) + (d(r_CB) x r_DB) dot r_AB
	      call dot(dd(1),dd(2),dd(3),
     1                 delr(1,1),delr(2,1),delr(3,1),dot1)
c d(r_AB) dot (r_CB x r_DB)
	      call dot(rCBxrDB(1), rCBxrDB(2),rCBxrDB(3),
     1		       drAB(1,j,i),drAB(2,j,i),drAB(3,j,i),
     2           dot2)
	      fdot(i,j,1) = dot1 + dot2
	  enddo
	enddo

C ----------------------- d( (r_DB x r_AB) dot r_CB) -------------------------
	do i = 1,3
	  do j = 1,4
c r_DB x d(r_AB)
	      call cross (delr(1,3),delr(2,3),delr(3,3),
     1                  drAB(1,j,i), drAB(2,j,i), drAB(3,j,i),
     2			rDBxdrAB(1),rDBxdrAB(2),rDBxdrAB(3))
c d(r_DB) x r_AB
	      call cross (drDB(1,j,i), drDB(2,j,i), drDB(3,j,i),
     1                  delr(1,1),delr(2,1),delr(3,1),
     2			drDBxrAB(1),drDBxrAB(2),drDBxrAB(3))
		do k = 1, 3
c (r_DB x d(r_AB)) + (d(r_DB) x r_AB)
		  dd(k) = rDBxdrAB(k) + drDBxrAB(k)
		enddo
c (r_DB x d(r_AB)) + (d(r_DB) x r_AB) dot r_CB
	      call dot(dd(1),dd(2),dd(3),
     1                 delr(1,2),delr(2,2),delr(3,2),dot1)
c d(r_CB) dot (r_DB x r_AB)
	      call dot(rDBxrAB(1), rDBxrAB(2),rDBxrAB(3),
     1                 drCB(1,j,i),drCB(2,j,i),drCB(3,j,i),
     2           dot2)
	      fdot(i,j,2) = dot1 + dot2
	  enddo
	enddo

C----------------------------------------------------------------------------
c calculate d( (r_AB x r_CB) dot r_DB) for each direction and atom
	do i = 1,3
	  do j = 1,4
c r_AB x d(r_CB)
	      call cross (delr(1,1),delr(2,1),delr(3,1),
     1                  drCB(1,j,i), drCB(2,j,i), drCB(3,j,i),
     2			rABxdrCB(1),rABxdrCB(2),rABxdrCB(3))
c d(r_AB) x r_CB
	      call cross (drAB(1,j,i), drAB(2,j,i), drAB(3,j,i),
     1                  delr(1,2),delr(2,2),delr(3,2),
     2			drABxrCB(1),drABxrCB(2),drABxrCB(3))
		do k = 1, 3
c (r_AB x d(r_CB)) + (d(r_AB) x r_CB)
		  dd(k) = rABxdrCB(k) + drABxrCB(k)
		enddo
c (r_AB x d(r_CB)) + (d(r_AB) x r_CB) dot r_DB
	      call dot(dd(1),dd(2),dd(3),
     1                 delr(1,3),delr(2,3),delr(3,3),dot1)
c d(r_DB) dot (r_AB x r_CB)
	      call dot(rABxrCB(1), rABxrCB(2),rABxrCB(3),
     1                 drDB(1,j,i),drDB(2,j,i),drDB(3,j,i),
     2           dot2)
	      fdot(i,j,3) = dot1 + dot2
	  enddo
	enddo

	do i = 1,3
	  do j = 1,4

	    f2(i,j,1) = (fdot(i,j,1) * invs3r(1)) + 
     1                  (dinvs3r(i,j,1) * dotCBDBAB)
	    dchi(i,j,1) = f2(i,j,1) / cos(chiABCD)

	    f2(i,j,2) = (fdot(i,j,2) * invs3r(2)) + 
     1                  (dinvs3r(i,j,2) * dotDBABCB)
	    dchi(i,j,2) = f2(i,j,2) / cos(chiCBDA)

	    f2(i,j,3) = (fdot(i,j,3) * invs3r(3)) + 
     1                  (dinvs3r(i,j,3) * dotABCBDB)
	    dchi(i,j,3) = f2(i,j,3) / cos(chiDBAC)

	    dtotalchi(i,j) = (dchi(i,j,1)+dchi(i,j,2)+
     1                           dchi(i,j,3)) / three

	  enddo
	enddo

          do j = 1, 4
              do i = 1,3
                fabcd(i,j) = 
     1             -two*improcoeff(1,itype)*deltachi*dtotalchi(i,j)
              enddo
          enddo
c
c Scatter the forces into the overall force and virial arrays
c
          do j = 1, 4
            if (newton_bond.eq.1.or.iindex(j).le.maxown) then
              do i = 1,3
                f(i,iindex(j)) = f(i,iindex(j)) + fabcd(i,j)
              enddo
            endif
          enddo

          if (iflag.eq.1) then

           virial(1) = virial(1) + (ifactor/4.0)*(
     1                              delr(1,1)*fabcd(1,1) +
     2                              delr(1,2)*fabcd(1,3) +
     3                              delr(1,3)*fabcd(1,4)  )
           virial(2) = virial(2) + (ifactor/4.0)*(
     1                              delr(2,1)*fabcd(2,1) +
     2                              delr(2,2)*fabcd(2,3) +
     3                              delr(2,3)*fabcd(2,4)  )
           virial(3) = virial(3) + (ifactor/4.0)*(
     1                              delr(3,1)*fabcd(3,1) +
     2                              delr(3,2)*fabcd(3,3) +
     3                              delr(3,3)*fabcd(3,4)  )
           virial(4) = virial(4) + (ifactor/4.0)*(
     1                              delr(1,1)*fabcd(2,1) +
     2                              delr(1,2)*fabcd(2,3) +
     3                              delr(1,3)*fabcd(2,4)  )
           virial(5) = virial(5) + (ifactor/4.0)*(
     1                              delr(1,1)*fabcd(3,1) +
     2                              delr(1,2)*fabcd(3,3) +
     3                              delr(1,3)*fabcd(3,4)  )
           virial(6) = virial(6) + (ifactor/4.0)*(
     1                              delr(2,1)*fabcd(3,1) +
     2                              delr(2,2)*fabcd(3,3) +
     3                              delr(2,3)*fabcd(3,4)  )

          endif

        endif

      enddo

      return
      end
c-----------------------------------------------------------------------
c  subroutine to perform cross-product of 2 vectors (In1 and In2) and write
c       the result into a third vector (Out)
c       written by Eric J. Simon, Cray Research, Inc. Fri Aug  2 1996
c
c         i       j       k
c       In1x    In1y    In1z
c       In2x    In2y    In2z
c
c	In1XIn2 = In1Y*In2z-In1z*In2y, In1z*In2x-In1x*In2z, In1x*In2y-In1y*In2x
c
      subroutine cross(In1x, In1y, In1z, In2x, In2y, In2z, 
     1                 Outx, Outy, Outz)

      real*8 In1x, In1y, In1z, In2x, In2y, In2z, Outx
      real*8 Outy, Outz

      Outx = In1y*In2z - In2y*In1z
      Outy = In2x*In1z - In1x*In2z
      Outz = In1x*In2y - In2x*In1y

      return
      end
c------------------------------------------------------------------------
c  subroutine to perform dot-product of 2 vectors (In1 and In2) and write
c       the result into the scalar Out
c       written by Eric J. Simon, Cray Research, Inc. Fri Aug  2 1996
c
      subroutine dot(Inx1, Iny1, Inz1,Inx2, Iny2, Inz2, Out)
      real*8 Inx1, Iny1, Inz1, Inx2, Iny2, Inz2, Out
      Out = Inx1*Inx2 + Iny1*Iny2 + Inz1*Inz2
      return
      end

c class 2 angle/angle interaction
c written by Eric J. Simon, Cray Research Inc.  Tue Jul 16 11:49:44 CDT 1996
c modified by ejs 1/21/97: removed test code, added newton/ifactor, commented
c			   out virial calculation
c
      subroutine angleangle_class2(iflag)
      include "lammps.h"
      real*8 dthetadr(3,4,3)
      real*8 fabcd(3,4)
      integer iindex(4)
      parameter (zero=0.0,one=1.0,two=2.0,three=3.0)
      parameter (negone=-1.0,small=0.00001)

      e_angleangle = 0.0
      
      do m = 1,nimprolocal
C	do raj = 1,2
        i1 = improlist(1,m)
        i2 = improlist(2,m)
        i3 = improlist(3,m)
        i4 = improlist(4,m)
        iindex(1) = i1
        iindex(2) = i2
        iindex(3) = i3
        iindex(4) = i4
        itype = improlist(5,m)

c     bond AB
        delxAB = x(1,i1) - x(1,i2)
        delyAB = x(2,i1) - x(2,i2)
        delzAB = x(3,i1) - x(3,i2)
        call minimg(delxAB,delyAB,delzAB)
        rABmag2 = delxAB*delxAB + delyAB*delyAB + delzAB*delzAB
        rAB = sqrt(rABmag2)
        
c bond BC
        delxBC = x(1,i3) - x(1,i2)
        delyBC = x(2,i3) - x(2,i2)
        delzBC = x(3,i3) - x(3,i2)
        call minimg(delxBC,delyBC,delzBC)
        rBCmag2 = delxBC*delxBC + delyBC*delyBC + delzBC*delzBC
        rBC = sqrt(rBCmag2)

c bond BD
        delxBD = x(1,i4) - x(1,i2) 
        delyBD = x(2,i4) - x(2,i2)
        delzBD = x(3,i4) - x(3,i2)
        call minimg(delxBD,delyBD,delzBD)
        rBDmag2 = delxBD*delxBD + delyBD*delyBD + delzBD*delzBD
        rBD = sqrt(rBDmag2)
        
c angle ABC
	costhABC = (delxAB*delxBC + delyAB*delyBC +delzAB*delzBC) / 
     1       (rAB * rBC)
	if (costhABC .gt. one)  costhABC = one
	if (costhABC .lt. -1.0) costhABC = -1.0
	thetaABC = acos(costhABC)

c angle ABD
	costhABD = (delxAB*delxBD + delyAB*delyBD +delzAB*delzBD) / 
     1       (rAB * rBD)
	if (costhABD .gt. one)  costhABD = one
	if (costhABD .lt. -1.0) costhABD = -1.0
	thetaABD = acos(costhABD)

c angle CBD
	costhCBD = (delxBC*delxBD + delyBC*delyBD +delzBC*delzBD) /
     1       (rBC * rBD)
	if (costhCBD .gt. one)  costhCBD = one
	if (costhCBD .lt. -1.0) costhCBD = -1.0
	thetaCBD = acos(costhCBD)

	dthABC = thetaABC - angleanglecoeff(4,itype)
        dthCBD = thetaCBD - angleanglecoeff(5,itype)
        dthABD = thetaABD - angleanglecoeff(6,itype)

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.maxown) ifactor = ifactor + 1
          if (i2.le.maxown) ifactor = ifactor + 1
          if (i3.le.maxown) ifactor = ifactor + 1
          if (i4.le.maxown) ifactor = ifactor + 1
        else
          ifactor = 4
        endif

        e_angleangle = e_angleangle + (ifactor/4.0 * (
     1       (angleanglecoeff(2,itype) * dthABC * dthABD) +
     2       (angleanglecoeff(1,itype) * dthABC * dthCBD) +
     3       (angleanglecoeff(3,itype) * dthABD * dthCBD)))

c set up d(theta)/d(r) array, dthetadr(i,j,k) = coordinate i, atom j, angle k
        do i = 1, 3
          do j = 1, 4
            do k = 1, 3
              dthetadr(i,j,k) = zero
            enddo
          enddo
        enddo

c angle ABC:
	sc1 = sqrt(one/(one - costhABC*costhABC))
        t1 = costhABC / rABmag2
        t3 = costhABC / rBCmag2
	r12 = one / (rAB * rBC)

        dthetadr(1,1,1) = sc1 * ((t1 * delxAB) - (delxBC * r12))
        dthetadr(2,1,1) = sc1 * ((t1 * delyAB) - (delyBC * r12))
        dthetadr(3,1,1) = sc1 * ((t1 * delzAB) - (delzBC * r12))

        dthetadr(1,2,1) = sc1 * ((negone*t1 * delxAB) + 
     1       (delxBC * r12)+
     2       (negone*t3 * delxBC) + (delxAB * r12))
        dthetadr(2,2,1) = sc1 * ((negone * t1 * delyAB) + 
     1       (delyBC * r12)+
     2       (negone * t3 * delyBC) + (delyAB * r12))
        dthetadr(3,2,1) = sc1 * ((negone * t1 * delzAB) + 
     1       (delzBC * r12)+
     2       (negone * t3 * delzBC) + (delzAB * r12))

        dthetadr(1,3,1) = sc1 * ((t3 * delxBC) - (delxAB * r12))
        dthetadr(2,3,1) = sc1 * ((t3 * delyBC) - (delyAB * r12))
        dthetadr(3,3,1) = sc1 * ((t3 * delzBC) - (delzAB * r12))

c angle CBD:
	sc1 = sqrt(one/(one - costhCBD*costhCBD))
        t1 = costhCBD / rBCmag2
        t3 = costhCBD / rBDmag2
	r12 = one / (rBC * rBD)

        dthetadr(1,3,2) = sc1 * ((t1 * delxBC) - (delxBD * r12))
        dthetadr(2,3,2) = sc1 * ((t1 * delyBC) - (delyBD * r12))
        dthetadr(3,3,2) = sc1 * ((t1 * delzBC) - (delzBD * r12))

        dthetadr(1,2,2) = sc1 * ((negone*t1 * delxBC) + 
     1       (delxBD * r12)+(negone*t3 * delxBD) + 
     2       (delxBC * r12))
        dthetadr(2,2,2) = sc1 * ((negone * t1 * delyBC) + 
     1       (delyBD * r12)+
     2       (negone* t3 * delyBD) + (delyBC * r12))
        dthetadr(3,2,2) = sc1 * ((negone* t1 * delzBC) + 
     1       (delzBD * r12)+
     2       (negone* t3 * delzBD) + (delzBC * r12))

        dthetadr(1,4,2) = sc1 * ((t3 * delxBD) - (delxBC * r12))
        dthetadr(2,4,2) = sc1 * ((t3 * delyBD) - (delyBC * r12))
        dthetadr(3,4,2) = sc1 * ((t3 * delzBD) - (delzBC * r12))

c angle ABD:
	sc1 = sqrt(one/(one - costhABD*costhABD))
        t1 = costhABD / rABmag2
        t3 = costhABD / rBDmag2
	r12 = one / (rAB * rBD)

        dthetadr(1,1,3) = sc1 * ((t1 * delxAB) - (delxBD * r12))
        dthetadr(2,1,3) = sc1 * ((t1 * delyAB) - (delyBD * r12))
        dthetadr(3,1,3) = sc1 * ((t1 * delzAB) - (delzBD * r12))

        dthetadr(1,2,3) = sc1 * ((negone*t1 * delxAB) + 
     1       (delxBD * r12)+
     2       (negone*t3 * delxBD) + (delxAB * r12))
        dthetadr(2,2,3) = sc1 * ((negone * t1 * delyAB) + 
     1       (delyBD * r12)+
     2       (negone * t3 * delyBD) + (delyAB * r12))
        dthetadr(3,2,3) = sc1 * ((negone * t1 * delzAB) + 
     1       (delzBD * r12)+
     2       (negone * t3 * delzBD) + (delzAB * r12))

        dthetadr(1,4,3) = sc1 * ((t3 * delxBD) - (delxAB * r12))
        dthetadr(2,4,3) = sc1 * ((t3 * delyBD) - (delyAB * r12))
        dthetadr(3,4,3) = sc1 * ((t3 * delzBD) - (delzAB * r12))

c angleangle forces
	do j = 1, 4
	    do i = 1,3
              fabcd(i,j) = negone*((angleanglecoeff(1,itype)*
     2                                (dthABC * dthetadr(i,j,2)+ 
     3                                 dthCBD * dthetadr(i,j,1)))+ 
     4                             (angleanglecoeff(2,itype) *
     5                                (dthABC * dthetadr(i,j,3)+ 
     6                                 dthABD * dthetadr(i,j,1)))+ 
     7                             (angleanglecoeff(3,itype) *
     8                                (dthABD * dthetadr(i,j,2)+ 
     9                                 dthCBD * dthetadr(i,j,3))))
	    enddo
        enddo

c Scatter forces in the overall force array and into the virial
c
        do j = 1, 4
          if (newton_bond.eq.1.or.iindex(j).le.maxown) then
            do i = 1, 3
              f(i,iindex(j)) = f(i,iindex(j)) + fabcd(i,j)
            enddo
          endif
        enddo

        if (iflag.eq.1) then

          virial(1) = virial(1) + (ifactor/4.0)*
     1                            ( delxAB*fabcd(1,1) + 
     2                              delxBC*fabcd(1,3) + 
     3                              delxBD*fabcd(1,4)  )
          virial(2) = virial(2) + (ifactor/4.0)*
     1                            ( delyAB*fabcd(2,1) + 
     2                              delyBC*fabcd(2,3) + 
     3                              delyBD*fabcd(2,4)  )
          virial(3) = virial(3) + (ifactor/4.0)*
     1                            ( delzAB*fabcd(3,1) + 
     2                              delzBC*fabcd(3,3) + 
     3                              delzBD*fabcd(3,4)  )
          virial(4) = virial(4) + (ifactor/4.0)*
     1                            ( delxAB*fabcd(2,1) + 
     2                              delxBC*fabcd(2,3) + 
     3                              delxBD*fabcd(2,4)  )
          virial(5) = virial(5) + (ifactor/4.0)*
     1                            ( delxAB*fabcd(3,1) + 
     2                              delxBC*fabcd(3,3) + 
     3                              delxBD*fabcd(3,4)  )
          virial(6) = virial(6) + (ifactor/4.0)*
     1                            ( delyAB*fabcd(3,1) + 
     2                              delyBC*fabcd(3,3) + 
     3                              delyBD*fabcd(3,4)  )
        end if

      enddo

     
      return
      end
