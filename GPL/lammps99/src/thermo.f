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
c compute all thermodynamic info

      subroutine thermo(iflag)
      include "lammps.h"
      include "mpif.h"
      integer iflag

 900  format (' ',a,i9,a,f11.4,a)
 901  format (' ',a,f15.4,a,f15.4,a,f15.4)
 902  format(i8,6f10.6,f11.4)

      call temperature
      call pressure
      call energy

c sum previous computed energy components
c  except e_long is already totaled in Ewald or PPPM

      tmp = e_vdwl
      call mpi_allreduce(tmp,e_vdwl,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_coul
      call mpi_allreduce(tmp,e_coul,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      tmp = e_bond
      call mpi_allreduce(tmp,e_bond,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_angle
      call mpi_allreduce(tmp,e_angle,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_dihedral
      call mpi_allreduce(tmp,e_dihedral,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      tmp = e_improper
      call mpi_allreduce(tmp,e_improper,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      if (anglestyle.eq.2) then
        tmp = e_bondbond
        call mpi_allreduce(tmp,e_bondbond,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_bondangle
        call mpi_allreduce(tmp,e_bondangle,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif
      if (dihedstyle.eq.2) then
        tmp = e_midbondtorsion
        call mpi_allreduce(tmp,e_midbondtorsion,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        tmp = e_endbondtorsion
        call mpi_allreduce(tmp,e_endbondtorsion,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        tmp = e_angletorsion
        call mpi_allreduce(tmp,e_angletorsion,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        tmp = e_angleangletorsion
        call mpi_allreduce(tmp,e_angleangletorsion,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        tmp = e_bondbond13
        call mpi_allreduce(tmp,e_bondbond13,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      endif
      if (improstyle.eq.3) then
        tmp = e_angleangle
        call mpi_allreduce(tmp,e_angleangle,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      endif

c add in cross-terms if class II is computed but not displayed

      if (thermostyle.le.1) then
        if (anglestyle.eq.2)
     $       e_angle = e_angle + e_bondbond + e_bondangle
        if (dihedstyle.eq.2)
     $       e_dihedral = e_dihedral + e_midbondtorsion +
     $       e_endbondtorsion + e_angletorsion +
     $       e_angleangletorsion + e_bondbond13
        if (improstyle.eq.3)
     $       e_improper = e_improper + e_angleangle
      endif

c total potential and kinetic energies

      if (thermostyle.le.1) then
        e_potential = e_long + e_vdwl + e_coul +
     $       e_bond + e_angle + e_dihedral + e_improper
      else
        e_potential = e_long + e_vdwl + e_coul +
     $       e_bond + e_angle + e_dihedral + e_improper +
     $       e_bondbond + e_bondangle + e_midbondtorsion +
     $       e_endbondtorsion + e_angletorsion +
     $       e_angleangletorsion + e_bondbond13 + e_angleangle
      endif

      e_kinetic = idimension*boltz*t_current/2.0 * natoms
      e_total = e_kinetic + e_potential
      
c normalize for LJ reduced units

      if (units.eq.1) then
        rnorm = 1.0/natoms
        e_long = e_long*rnorm
        e_vdwl = e_vdwl*rnorm
        e_coul = e_coul*rnorm
        e_bond = e_bond*rnorm
        e_angle = e_angle*rnorm
        e_dihedral = e_dihedral*rnorm
        e_improper = e_improper*rnorm
        if (anglestyle.eq.2) then
          e_bondbond = e_bondbond*rnorm
          e_bondangle = e_bondangle*rnorm
        endif
        if (dihedstyle.eq.2) then
          e_midbondtorsion = e_midbondtorsion*rnorm
          e_endbondtorsion = e_endbondtorsion*rnorm
          e_angletorsion = e_angletorsion*rnorm
          e_angleangletorsion = e_angleangletorsion*rnorm
          e_bondbond13 = e_bondbond13*rnorm
        endif
        if (improstyle.eq.3) e_angleangle = e_angleangle*rnorm
        e_potential = e_potential*rnorm
        e_kinetic = e_kinetic*rnorm
        e_total = e_total*rnorm
      endif

c set time based on start, middle, or end

      if (iflag.eq.-1) then
        time = 0.0
      else if (iflag.eq.0) then
        time = mpi_wtime()
        time = time-time_loop
      else if (iflag.eq.1) then
        time = time_loop
      endif

c screen and file output

      if (thermostyle.eq.0) then

        if (node.eq.0) then
          write (6,900) '--------------- Step',ntimestep,
     $         ' ---- CPU = ',time,' (sec) ---------------'
          write (6,901) 'Total E =',e_total,
     $         ' Total KE=',e_kinetic,' Temp    =',t_current
          write (6,901) 'Total PE=',e_potential,
     $         ' E_bond  =',e_bond,' E_angle =',e_angle
          write (6,901) 'E_dihed =',e_dihedral,
     $         ' E_impr  =',e_improper,' E_vdwl  =',e_vdwl
          write (6,901) 'E_coul  =',e_coul,' E_long  =',e_long,
     $         ' Press   =',p_total*pfactor
          if (ensemble.eq.3)
     $         write (6,901) 'Volume  =',xprd*yprd*zprd
          write (6,*)

          write (1,900) '--------------- Step',ntimestep,
     $         ' ---- CPU = ',time,' (sec) ---------------'
          write (1,901) 'Total E =',e_total,
     $         ' Total KE=',e_kinetic,' Temp    =',t_current
          write (1,901) 'Total PE=',e_potential,
     $         ' E_bond  =',e_bond,' E_angle =',e_angle
          write (1,901) 'E_dihed =',e_dihedral,
     $         ' E_impr  =',e_improper,' E_vdwl  =',e_vdwl
          write (1,901) 'E_coul  =',e_coul,' E_long  =',e_long,
     $         ' Press   =',p_total*pfactor
          if (ensemble.eq.3)
     $         write (1,901) 'Volume  =',xprd*yprd*zprd
        endif

      else if (thermostyle.eq.1) then

        if (node.eq.0) then
          write (6,902) ntimestep,t_current,e_vdwl+e_coul,
     $         e_bond+e_angle+e_dihedral+e_improper,e_long,
     $         e_total,p_total*pfactor,xprd*yprd*zprd
          write (1,902) ntimestep,t_current,e_vdwl+e_coul,
     $         e_bond+e_angle+e_dihedral+e_improper,e_long,
     $         e_total,p_total*pfactor,xprd*yprd*zprd
        endif

      else if (thermostyle.eq.2) then

        if (node.eq.0) then
          write (6,900) '--------------- Step',ntimestep,
     $         ' ---- CPU = ',time,' (sec) ---------------'
          write (6,901) 'Total E =',e_total,
     $         ' Total KE=',e_kinetic,' Temp    =',t_current
          write (6,901) 'Total PE=',e_potential,
     $         ' E_bond  =',e_bond,' E_angle =',e_angle
          write (6,901) 'E_dihed =',e_dihedral,
     $         ' E_impr  =',e_improper,' E_vdwl  =',e_vdwl
          write (6,901) 'E_coul  =',e_coul,' E_long  =',e_long,
     $         ' Press   =',p_total*pfactor
          write (6,901) 'BondBond=',e_bondbond,
     $         ' BondAng =',e_bondangle,
     $         ' MidBond =',e_midbondtorsion
          write (6,901) 'EndBond =',e_endbondtorsion,
     $         ' AngTors =',e_angletorsion,
     $         ' AA-Tors =',e_angleangletorsion
          write (6,901) 'BnBn13  =',e_bondbond13,
     $         ' Ang-Ang =',e_angleangle
          if (ensemble.eq.3)
     $         write (6,901) 'Volume  =',xprd*yprd*zprd
          write (6,*)

          write (1,900) '--------------- Step',ntimestep,
     $         ' ---- CPU = ',time,' (sec) ---------------'
          write (1,901) 'Total E =',e_total,
     $         ' Total KE=',e_kinetic,' Temp    =',t_current
          write (1,901) 'Total PE=',e_potential,
     $         ' E_bond  =',e_bond,' E_angle =',e_angle
          write (1,901) 'E_dihed =',e_dihedral,
     $         ' E_impr  =',e_improper,' E_vdwl  =',e_vdwl
          write (1,901) 'E_coul  =',e_coul,' E_long  =',e_long,
     $         ' Press   =',p_total*pfactor
          write (1,901) 'BondBond=',e_bondbond,
     $         ' BondAng =',e_bondangle,
     $         ' MidBond =',e_midbondtorsion
          write (1,901) 'EndBond =',e_endbondtorsion,
     $         ' AngTors =',e_angletorsion,
     $         ' AA-Tors =',e_angleangletorsion
          write (1,901) 'BnBn13  =',e_bondbond13,
     $         ' Ang-Ang =',e_angleangle
          if (ensemble.eq.3)
     $         write (1,901) 'Volume  =',xprd*yprd*zprd
          write (1,*)
        endif

      endif

c check for lost atoms
c "lost" means an atom is outside the periodic box
c atom(s) can also be lost if sum of nlocal across all procs != global natoms

      call mpi_allreduce(nlocal,nsum,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)

      nlost = 0
      i = atompnt
      do ii = 1,nlocal
        jflag = 0
        if (perflagx.eq.0.and.
     $       (x(1,i).lt.xboundlo.or.x(1,i).ge.xboundhi)) jflag = 1
        if (perflagy.eq.0.and.
     $       (x(2,i).lt.yboundlo.or.x(2,i).ge.yboundhi)) jflag = 1
        if (perflagz.eq.0.and.
     $       (x(3,i).lt.zboundlo.or.x(3,i).ge.zboundhi)) jflag = 1
        if (jflag.eq.1) nlost = nlost + 1
        i = list(i)
      enddo

      itmp = nlost
      call mpi_allreduce(itmp,nlost,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)

      if (nsum.ne.natoms.or.nlost.gt.0) then
        if (node.eq.0) write (6,*) 'Atom counts =',nsum,nlost
        call error('Incorrect number of atoms')
      endif

      nthermo_next = min(ntimestep+nthermo,ntime_last)

      return
      end
      

c -------------------------------------------------------------------------
c compute current temperature
      
      subroutine temperature
      include "lammps.h"
      include "mpif.h"
      
      t_current = 0.0
      i = atompnt
      do ii = 1,nlocal
        t_current = t_current +
     $       (v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i)) *
     $       mass(type(i))
        i = list(i)
      enddo

      call mpi_allreduce(t_current,tmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      t_current = tmp / (boltz*idimension*natoms)
      
      return
      end
      

c -------------------------------------------------------------------------
c compute pressure from virial
c do not sum virial array in place since this routine may be called more
c  than once per timestep

      subroutine pressure
      include "lammps.h"
      include "mpif.h"

      dimension p_virial(3),p_virial_all(3)

      if (nrespa.eq.0) then
        p_virial(1) = virial(1) 
        p_virial(2) = virial(2)
        p_virial(3) = virial(3)
      else
        p_virial(1) = vir_stretch(1) + vir_intra(1) +
     $       vir_short(1) + vir_long(1) 
        p_virial(2) = vir_stretch(2) + vir_intra(2) +
     $       vir_short(2) + vir_long(2) 
        p_virial(3) = vir_stretch(3) + vir_intra(3) +
     $       vir_short(3) + vir_long(3)
      endif

      call mpi_allreduce(p_virial,p_virial_all,3,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

      tempk = natoms*boltz*t_current
      vol = xprd*yprd*zprd
      
      p_current(1) = (tempk + p_virial_all(1)) /vol
      p_current(2) = (tempk + p_virial_all(2)) /vol
      p_current(3) = (tempk + p_virial_all(3)) /vol

      p_total = (p_current(1) + p_current(2) + p_current(3))/3.0

      return
      end


c -------------------------------------------------------------------------
c compute nonbond energy for all LJ and Coulombic styles

      subroutine energy
      include "lammps.h"

      data pi /3.1415926/
      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/
      
      e_vdwl = 0.0
      e_coul = 0.0

      factor = 0.0
      if (coulstyle.ge.3) factor = 2.0

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich) - factor
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          call minimg(delx,dely,delz)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            ffactor = 1.0
            if (newton_nonbond.eq.0.and.j.gt.maxown) ffactor = 0.5
            r2inv = 1.0/rsq

            if (rsq.lt.cutcoulsq) then
              if (coulstyle.eq.1) then
                phicoul = coulpre*qtmp*q(j)*sqrt(r2inv)
                e_coul = e_coul + ffactor*spfactor*phicoul
              else if (coulstyle.eq.2) then
                if (rsq.lt.cutcoulintsq) then
                  phicoul = coulpre*qtmp*q(j)*sqrt(r2inv)
                else
                  r = sqrt(rsq)
                  phicoul = (coulpre*qtmp*q(j) / r) *
     $                 (cutcoul - r)**2 *
     $                 (cutcoul + (2.0D0*r) - (3.0D0*cutcoulint)) /
     $                 ( (cutcoul - cutcoulint)**3 )
                endif
                e_coul = e_coul + ffactor*spfactor*phicoul
              else
                r = sqrt(rsq)
                grij = gewsr*r
                expm2 = exp(-grij*grij)
                t = 1.0/(1.0+ewald_p*grij)
                erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $               (ewald_a4+t*ewald_a5))))*expm2
                prefactor = coulpre*qtmp*q(j)/r
                phicoul = prefactor*erfc
                if (spfactor.lt.1.0)
     $               phicoul = phicoul - (1.0-spfactor)*prefactor
                e_coul = e_coul + ffactor*phicoul
              endif
            endif

            if (rsq.lt.cutljsq(itype,jtype)) then
              if (nonstyle.eq.1) then
                r6inv = r2inv*r2inv*r2inv
                philj = r6inv *
     $               (lj3(itype,jtype)*r6inv - lj4(itype,jtype)) -
     $               offset(itype,jtype)
              else if (nonstyle.eq.2) then
                if (rsq.lt.cutljinnersq(itype,jtype)) then
                  r6inv = r2inv*r2inv*r2inv
                  philj = r6inv *
     $                 (lj3(itype,jtype)*r6inv - lj4(itype,jtype)) -
     $                 offset(itype,jtype)
                else
                  r = sqrt(rsq) 
                  t = r - cutljinner(itype,jtype)
                  tsq = t*t
                  philj = ljsw0(itype,jtype) -
     $                 ljsw1(itype,jtype)*t -
     $                 ljsw2(itype,jtype)*tsq/2.0 -
     $                 ljsw3(itype,jtype)*tsq*t/3.0 -
     $                 ljsw4(itype,jtype)*tsq*tsq/4.0 -
     $                 offset(itype,jtype)
                endif
              else if (nonstyle.eq.3) then
                r = sqrt(rsq)
                rshift = r - lj5(itype,jtype)
                rshiftsq = rshift*rshift
                r2inv = 1.0/rshiftsq
                r6inv = r2inv*r2inv*r2inv
                philj = r6inv *
     $               (lj3(itype,jtype)*r6inv - lj4(itype,jtype)) -
     $               offset(itype,jtype)
              else if (nonstyle.eq.4) then
                r = sqrt(rsq)
                arg = pi*r/lj3(itype,jtype)
                philj = lj4(itype,jtype)*(1.0+cos(arg))
              else if (nonstyle.eq.5) then
                rinv = sqrt(r2inv)
                r3inv = r2inv*rinv
                r6inv = r3inv*r3inv
                philj = r6inv *
     $               (lj3(itype,jtype)*r3inv - lj4(itype,jtype)) -
     $               offset(itype,jtype)
              endif
              e_vdwl = e_vdwl + ffactor*spfactor*philj
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end
