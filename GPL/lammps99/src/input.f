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

c ------------------------------------------------------------------------
c read all info from text input file one line at a time

      subroutine input(iopflag)
      include "lammps.h"
      include "mpif.h"
      integer iopflag

      external length
      logical match
      character*200 str
 900  format (a)
      
      do infinity = 1,100000000
        
        if (node.eq.0) read (5,900,end=999,err=999) str
        call mpi_bcast(str,200,mpi_character,0,mpi_comm_world,ierror)

        if (length(str).eq.0) then
          continue
        else if (match('#',str,m)) then
          if (node.eq.0) write (1,900) str(1:length(str))
        else if (match('units',str,m)) then
          call in_units(str,m)
        else if (match('dimension',str,m)) then
          call in_dimension(str,m)
        else if (match('processor grid',str,m)) then
          call in_processor_grid(str,m)
        else if (match('periodicity',str,m)) then
          call in_periodicity(str,m)
        else if (match('newton flag',str,m)) then
          call in_newton_flag(str,m)
        else if (match('timestep',str,m)) then
          call in_timestep(str,m)
        else if (match('respa',str,m)) then
          call in_respa(str,m)
        else if (match('neighbor',str,m)) then
          call in_neighbor(str,m)
        else if (match('special bonds',str,m)) then
          call in_special_bonds(str,m)
        else if (match('thermo flag',str,m)) then
          call in_thermo_flag(str,m)
        else if (match('thermo style',str,m)) then
          call in_thermo_style(str,m)
        else if (match('true flag',str,m)) then
          call in_true_flag(str,m)
        else if (match('dump atoms',str,m)) then
          call in_dump_atoms(str,m)
        else if (match('dump velocities',str,m)) then
          call in_dump_velocities(str,m)
        else if (match('dump forces',str,m)) then
          call in_dump_forces(str,m)
        else if (match('restart',str,m)) then
          call in_restart(str,m)
        else if (match('diagnostic',str,m)) then
          call in_diagnostic(str,m)
        else if (match('temp control',str,m)) then
          call in_temp_control(str,m)
        else if (match('press control',str,m)) then
          call in_press_control(str,m)
        else if (match('press_x control',str,m)) then
          call in_press_x_control(str,m)
        else if (match('press_y control',str,m)) then
          call in_press_y_control(str,m)
        else if (match('press_z control',str,m)) then
          call in_press_z_control(str,m)
        else if (match('nonbond style',str,m)) then
          call in_nonbond_style(str,m)
        else if (match('nonbond coeff',str,m)) then
          call in_nonbond_coeff(str,m)
        else if (match('mixing style',str,m)) then
          call in_mixing_style(str,m)
        else if (match('coulomb style',str,m)) then
          call in_coulomb_style(str,m)
        else if (match('pppm mesh',str,m)) then
          call in_pppm_mesh(str,m)
        else if (match('pppm order',str,m)) then
          call in_pppm_order(str,m)
        else if (match('dielectric',str,m)) then
          call in_dielectric(str,m)
        else if (match('bond style',str,m)) then
          call in_bond_style(str,m)
        else if (match('bond coeff',str,m)) then
          call in_bond_coeff(str,m)
        else if (match('angle style',str,m)) then
          call in_angle_style(str,m)
        else if (match('dihedral style',str,m)) then
          call in_dihedral_style(str,m)
        else if (match('improper style',str,m)) then
          call in_improper_style(str,m)
        else if (match('min style',str,m)) then
          call in_min_style(str,m)
        else if (match('min file',str,m)) then
          call in_min_file(str,m)
        else if (match('read data',str,m)) then
          iopflag = 1
          call in_read_data(str,m)
          return
        else if (match('create group',str,m)) then
          call in_create_group(str,m)
        else if (match('create temp',str,m)) then
          iopflag = 2
          call in_create_temp(str,m)
          return
        else if (match('read restart',str,m)) then
          iopflag = 3
          call in_read_restart(str,m)
          return
        else if (match('fix style',str,m)) then
          call in_fix_style(str,m)
        else if (match('assign fix',str,m)) then
          iopflag = 4
          call in_assign_fix(str,m)
          return
        else if (match('reset timestep',str,m)) then
          call in_reset_timestep(str,m)
        else if (match('run',str,m)) then
          iopflag = 5
          call in_run(str,m)
          return
        else if (match('minimize',str,m)) then
          iopflag = 6
          call in_minimize(str,m)
          return
        else
          call error('STRING: '//str(1:length(str)))
        endif
        
      enddo
      
 999  continue
      
      str = 'All Done'
      call mpi_bcast(str,200,mpi_character,0,mpi_comm_world,ierror)
      call error(str(1:8))
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_units(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      character*16 work
      logical match

      readflag = 0
      call strread(str,m,work)
      if (length(work).le.0) call error('Bad units param')
      if (match('real',work,n)) then
        units = 0
        boltz = 0.001987191
        dtfactor = 48.88821
        pfactor = 68589.796
        efactor = 332.054
      else if (match('lj',work,n)) then
        units = 1
        boltz = 1.0
        dtfactor = 1.0
        pfactor = 1.0
        efactor = 1.0
      else
        call error('Bad units param')
      endif
      if (node.eq.0) then
        if (units.eq.0) write (1,*) 'Units real'
        if (units.eq.1) write (1,*) 'Units lj'
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dimension(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      idimension = intread(str,m)
      readflag = 0
      if (idimension.lt.2.or.idimension.gt.3)
     $     call error('Bad dimension parameter')
      if (node.eq.0)
     $     write (1,*) 'Dimension',idimension

      return
      end

c ------------------------------------------------------------------------

      subroutine in_processor_grid(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      pgrid(1) = intread(str,m)
      pgrid(2) = intread(str,m)
      pgrid(3) = intread(str,m)
      if (pgrid(1).le.0.or.pgrid(2).le.0.or.pgrid(3).le.0)
     $     call error('Bad processor grid parameter')
      if (pgrid(1)*pgrid(2)*pgrid(3).ne.nprocs)
     $     call error('Specified processor grid not equal Nprocs')
      readflag = 0
      if (node.eq.0)
     $     write (1,*) 'Processor grid',pgrid(1),pgrid(2),pgrid(3)

      return
      end

c ------------------------------------------------------------------------

      subroutine in_periodicity(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      perflagx = intread(str,m)
      perflagy = intread(str,m)
      perflagz = intread(str,m)
      if (perflagx.lt.0.or.perflagx.gt.1.or.
     $     perflagy.lt.0.or.perflagy.gt.1.or.
     $     perflagz.lt.0.or.perflagz.gt.1)
     $     call error('Bad periodicity parameter')
      readflag = 0
      if (node.eq.0)
     $     write (1,*) 'Periodicity',perflagx,perflagy,perflagz

      return
      end

c ------------------------------------------------------------------------

      subroutine in_newton_flag(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      newton = intread(str,m)
      if (newton.lt.0.or.newton.gt.3)
     $     call error('Bad newton flag parameter')
      readflag = 0
      newton_bond = mod(newton,2)
      newton_nonbond = newton/2
      if (node.eq.0)
     $     write (1,*) 'Newton flag',newton

      return
      end

c ------------------------------------------------------------------------

      subroutine in_timestep(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      dt = realread(str,m)
      if (dt.le.0.0) call error('Bad timestep parameter')
      if (node.eq.0)
     $     write (1,*) 'Timestep',dt
      dt = dt/dtfactor

      return
      end

c ------------------------------------------------------------------------

      subroutine in_respa(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      nstretch = intread(str,m)
      nintra = intread(str,m)
      nshort = intread(str,m)
      if (nstretch.le.0.or.nintra.le.0.or.nshort.le.0)
     $     call error('Bad respa parameter')
      nrespa = 1
      if (nstretch.eq.1.and.nintra.eq.1.and.nshort.eq.1) nrespa = 0
      if (node.eq.0)
     $     write (1,*) 'Respa',nstretch,nintra,nshort
      if (nrespa.eq.0) then
        peratom_comm = 19
      else
        peratom_comm = 25
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_neighbor(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      skin = realread(str,m)
      neighstyle = intread(str,m)
      neighfreq = intread(str,m)
      neighdelay = intread(str,m)
      neightrigger = intread(str,m)
      if (skin.lt.0.0.or.
     $     (neighstyle.lt.0.or.neighstyle.gt.1).or.
     $     neighfreq.le.0.or.neighdelay.lt.0.or.
     $     (neightrigger.lt.0.or.neightrigger.gt.1))
     $     call error('Bad neighbor parameter')
      if (node.eq.0)
     $     write (1,*) 'Neighbor',skin,neighstyle,
     $     neighfreq,neighdelay,neightrigger

      return
      end

c ------------------------------------------------------------------------

      subroutine in_special_bonds(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      special(1) = realread(str,m)
      special(2) = realread(str,m)
      special(3) = realread(str,m)
      if (special(1).lt.0.0.or.special(1).gt.1.0.or.
     $     special(2).lt.0.0.or.special(2).gt.1.0.or.
     $     special(3).lt.0.0.or.special(3).gt.1.0)
     $     call error('Bad special bonds parameter')
      if (node.eq.0)
     $     write (1,*) 'Special bonds',
     $     special(1),special(2),special(3)

      return
      end

c ------------------------------------------------------------------------

      subroutine in_thermo_flag(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      nthermo = intread(str,m)
      if (nthermo.lt.0)
     $     call error('Bad thermo flag parameter')
      if (node.eq.0)
     $     write (1,*) 'Thermo flag',nthermo

      return
      end

c ------------------------------------------------------------------------

      subroutine in_thermo_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      thermostyle = intread(str,m)
      if (thermostyle.lt.0.or.thermostyle.gt.2)
     $     call error('Bad thermo style parameter')
      if (node.eq.0)
     $     write (1,*) 'Thermo style',thermostyle
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_true_flag(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length

      trueflag = intread(str,m)
      if (trueflag.lt.0.or.trueflag.gt.3)
     $     call error('Bad true flag parameter')
      if (node.eq.0)
     $     write (1,*) 'True flag',trueflag

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dump_atoms(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length

      ndumpatom = intread(str,m)
      if (ndumpatom.lt.0) call error('Bad dump atoms parameter')

      if (ndumpatom.eq.0) then
        if (dumpatomfileflag.gt.0) call dumpatomclose
        dumpatomfileflag = 0
        if (node.eq.0) write (1,*) 'Dump atoms',ndumpatom
      else
        call strread(str,m,dumpatomfile)
        if (length(dumpatomfile).eq.0)
     $       call error('Bad dump atoms param')
        if (dumpatomfileflag.gt.0) call dumpatomclose
        call dumpatomopen
        dumpatomfileflag = 1
        if (node.eq.0) write (1,*) 'Dump atoms',ndumpatom,' ',
     $       dumpatomfile(1:length(dumpatomfile))
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dump_velocities(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length

      ndumpvel = intread(str,m)
      if (ndumpvel.lt.0) call error('Bad dump velocities parameter')

      if (ndumpvel.eq.0) then
        if (dumpvelfileflag.gt.0) call dumpvelclose
        dumpvelfileflag = 0
        if (node.eq.0) write (1,*) 'Dump velocities',ndumpvel
      else
        call strread(str,m,dumpvelfile)
        if (length(dumpvelfile).eq.0)
     $       call error('Bad dump velocities param')
        if (dumpvelfileflag.gt.0) call dumpvelclose
        call dumpvelopen
        dumpvelfileflag = 1
        if (node.eq.0) write (1,*) 'Dump velocities',ndumpvel,' ',
     $       dumpvelfile(1:length(dumpvelfile))
      endif
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_dump_forces(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length

      ndumpforce = intread(str,m)
      if (ndumpforce.lt.0) call error('Bad dump forces parameter')

      if (ndumpforce.eq.0) then
        if (dumpforcefileflag.gt.0) call dumpforceclose
        dumpforcefileflag = 0
        if (node.eq.0) write (1,*) 'Dump forces',ndumpforce
      else
        call strread(str,m,dumpforcefile)
        if (length(dumpforcefile).eq.0)
     $       call error('Bad dump forces param')
        if (dumpforcefileflag.gt.0) call dumpforceclose
        call dumpforceopen
        dumpforcefileflag = 1
        if (node.eq.0) write (1,*) 'Dump forces',ndumpforce,' ',
     $       dumpforcefile(1:length(dumpforcefile))
      endif
      
      return
      end

c ------------------------------------------------------------------------


      subroutine in_restart(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length

      nrestart = intread(str,m)
      if (nrestart.lt.0) call error('Bad restart param')

      if (nrestart.eq.0) then
        if (node.eq.0) write (1,*) 'Restart',nrestart
      else
        call strread(str,m,restart_out1)
        call strread(str,m,restart_out2)
        if (length(restart_out1).eq.0.or.
     $       length(restart_out2).eq.0)
     $       call error('Bad restart param')
        if (node.eq.0) write (1,*) 'Restart',nrestart,
     $       restart_out1(1:length(restart_out1)),' ',
     $       restart_out2(1:length(restart_out2))
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_diagnostic(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad diagnostic parameter')

      call diag_lookup(work(1:length(work)),idiag)
      if (idiag.eq.0) call error('Non-existent diagnostic')

      ndiag(idiag) = intread(str,m)
      if (ndiag(idiag).lt.0) call error('Bad diagnostic param')

c close out this diagnostic

      if (ndiag(idiag).eq.0) then
        if (diagfileflag(idiag).gt.0) call diag_close(idiag)
        diagfileflag(idiag) = 0
        if (node.eq.0) write (1,*) 'Diagnostic ',
     $       work(1:length(work)),' ',ndiag(idiag)
      else

c initialize this diagnostic

        call strread(str,m,diagfile(idiag))
        if (length(diagfile(idiag)).eq.0)
     $       call error('Bad diagnostic param')
        if (diagfileflag(idiag).gt.0) call diag_close(idiag)
        if (diagfile(idiag)(1:length(diagfile(idiag))).ne.'none') then
          call diag_open(idiag)
          diagfileflag(idiag) = 1
        endif
        diagnparams(idiag) = intread(str,m)
        if (diagnparams(idiag).lt.0.or.diagnparams(idiag).gt.5)
     $       call error('Bad diagnostic param')
        do i = 1,5
          diagparam(i,idiag) = 0.0
          if (diagnparams(idiag).ge.i)
     $         diagparam(i,idiag) = realread(str,m)
        enddo
        if (node.eq.0) write (1,*) 'Diagnostic',
     $       work(1:length(work)),' ',ndiag(idiag),
     $       diagfile(idiag)(1:length(diagfile(idiag))),
     $       diagparam(1,idiag),diagparam(2,idiag),diagparam(3,idiag),
     $       diagparam(4,idiag),diagparam(5,idiag)
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_temp_control(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad temp control parameter')
      if (match('none',work,mm)) then
        tempflag = 0
      else if (match('rescale',work,mm)) then
        tempflag = 1
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_every = intread(str,m)
        t_window = realread(str,m)
        if (t_start.lt.0.0.or.t_stop.lt.0.0.or.
     $       t_every.le.0.or.t_window.lt.0.0)
     $       call error('Bad temp control parameter')
      else if (match('replace',work,mm)) then
        tempflag = 2
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_every = intread(str,m)
        iseed = realread(str,m)
        if (t_start.lt.0.0.or.t_stop.lt.0.0.or.t_every.le.0.or.
     $       iseed.le.0.or.iseed.ge.1000000000)
     $       call error('Bad temp control parameter')
        tmp = ranmars(iseed+node)
      else if (match('langevin',work,mm)) then
        tempflag = 3
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_freq = realread(str,m)
        iseed = intread(str,m)
        if (t_start.lt.0.0.or.t_stop.lt.0.0.or.
     $       t_freq.lt.0.0.or.iseed.le.0.or.iseed.ge.1000000000)
     $       call error('Bad temp control parameter')
        tmp = ranmars(iseed+node)
      else if (match('nose/hoover',work,mm)) then
        tempflag = 4
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_freq = realread(str,m)
        if (t_start.lt.0.0.or.t_stop.lt.0.0.or.t_freq.lt.0.0)
     $       call error('Bad temp control parameter')
      else
        call error('Bad temp control parameter')
      endif
      if (node.eq.0) then
        if (tempflag.eq.0) then
          write (1,*) 'Temp control none'
        else if (tempflag.eq.1) then
          write (1,*) 'Temp control rescale',
     $         t_start,t_stop,t_every,t_window
        else if (tempflag.eq.2) then
          write (1,*) 'Temp control replace',
     $         t_start,t_stop,t_every,iseed
        else if (tempflag.eq.3) then
          write (1,*) 'Temp control langevin',
     $         t_start,t_stop,t_freq,iseed
        else if (tempflag.eq.4) then
          write (1,*) 'Temp control nose/hoover',
     $         t_start,t_stop,t_freq
        endif
      endif
      if (tempflag.eq.3) t_freq = t_freq*dtfactor
      if (tempflag.eq.4) t_freq = t_freq*dtfactor

      return
      end

c ------------------------------------------------------------------------

      subroutine in_press_control(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad press control parameter')
      if (match('none',work,mm)) then
        pressflag = 0
      else if (match('nose/hoover',work,mm)) then
        pressflag = 1
        p_start(1) = realread(str,m)
        p_stop(1) = realread(str,m)
        p_freq = realread(str,m)
        p_start(2) = p_start(1)
        p_start(3) = p_start(1)
        p_stop(2) = p_stop(1)
        p_stop(3) = p_stop(1)
        p_freq2(1) = p_freq*p_freq
        p_freq2(2) = p_freq2(1)
        p_freq2(3) = p_freq2(1)
        if (p_freq.lt.0.0)
     $       call error('Bad press control parameter')
      else
        call error('Bad press control parameter')
      endif
      if (node.eq.0) then
        if (pressflag.eq.0) then
          write (1,*) 'Press control none'
        else if (pressflag.eq.1) then
          write (1,*) 'Press control nose/hoover',
     $         p_start(1),p_stop(1),p_freq
        endif
      endif
      if (pressflag.eq.1) then
        p_start(1) = p_start(1)/pfactor
        p_stop(1) = p_stop(1)/pfactor
        p_freq2(1) = p_freq*p_freq*dtfactor*dtfactor
        p_start(2) = p_start(2)/pfactor
        p_stop(2) = p_stop(2)/pfactor
        p_freq2(2) = p_freq*p_freq*dtfactor*dtfactor
        p_start(3) = p_start(3)/pfactor
        p_stop(3) = p_stop(3)/pfactor
        p_freq2(3) = p_freq*p_freq*dtfactor*dtfactor
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_press_x_control(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad press control parameter')
      if (pressflag.eq.1) 
     $     call error('Cannot use press control & press_x control') 
      if (match('none',work,mm)) then
        xpressflag = 0
      else if (match('nose/hoover',work,mm)) then
        xpressflag = 1
        p_start(1) = realread(str,m)
        p_stop(1) = realread(str,m)
        p_freq = realread(str,m)
        p_freq2(1) = p_freq*p_freq
        if (p_freq.lt.0.0)
     $       call error('Bad press control parameter')
      else
        call error('Bad press control parameter')
      endif
      if (node.eq.0) then
        if (xpressflag.eq.0) then
          write (1,*) 'x-pressure control none'
        else if (pressflag.eq.1) then
          write (1,*) 'x-pressure control nose/hoover',
     $         p_start(1),p_stop(1),p_freq
        endif
      endif
      if (xpressflag.eq.1) then
        p_start(1) = p_start(1)/pfactor
        p_stop(1) = p_stop(1)/pfactor
        p_freq2(1) = p_freq2(1)*dtfactor*dtfactor
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_press_y_control(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad press control parameter')
      if (pressflag.eq.1) 
     $     call error('Cannot use press control & press_y control') 
      if (match('none',work,mm)) then
        ypressflag = 0
      else if (match('nose/hoover',work,mm)) then
        ypressflag = 1
        p_start(2) = realread(str,m)
        p_stop(2) = realread(str,m)
        p_freq = realread(str,m)
        p_freq2(2) = p_freq*p_freq
        if (p_freq.lt.0.0)
     $       call error('Bad press control parameter')
      else
        call error('Bad press control parameter')
      endif
      if (node.eq.0) then
        if (ypressflag.eq.0) then
          write (1,*) 'y-pressure control none'
        else if (pressflag.eq.1) then
          write (1,*) 'y-pressure control nose/hoover',
     $         p_start(2),p_stop(2),p_freq
        endif
      endif
      if (ypressflag.eq.1) then
        p_start(2) = p_start(2)/pfactor
        p_stop(2) = p_stop(2)/pfactor
        p_freq2(2) = p_freq2(2)*dtfactor*dtfactor
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_press_z_control(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad press control parameter')
      if (pressflag.eq.1) 
     $     call error('Cannot use press control & press_z control') 
      if (match('none',work,mn)) then
        zpressflag = 0
      else if (match('nose/hoover',work,mm)) then
        zpressflag = 1
        p_start(3) = realread(str,m)
        p_stop(3) = realread(str,m)
        p_freq = realread(str,m)
        p_freq2(3) = p_freq*p_freq
        if (p_freq.lt.0.0)
     $       call error('Bad press control parameter')
      else
        call error('Bad press control parameter')
      endif
      if (node.eq.0) then
        if (zpressflag.eq.0) then
          write (1,*) 'z-pressure control none'
        else if (zpressflag.eq.1) then
          write (1,*) 'z-pressure control nose/hoover',
     $         p_start(3),p_stop(3),p_freq
        endif
      endif
      if (zpressflag.eq.1) then
        p_start(3) = p_start(3)/pfactor
        p_stop(3) = p_stop(3)/pfactor
        p_freq2(3) = p_freq2(3)*dtfactor*dtfactor
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_nonbond_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad nonbond style parameter')
      if (match('none',work,mm)) then
        nonstyle = 0
        do i = 1,maxtype
          do j = i,maxtype
            nontypeflag(i,j) = nonstyle
          enddo
        enddo
      else if (match('lj/cutoff',work,mm)) then
        nonstyle = 1
        cutlj = realread(str,m)
        offsetflag = intread(str,m)
        if (cutlj.le.0.0.or.offsetflag.lt.0.or.offsetflag.gt.1)
     $       call error('Bad nonbond style parameter')
        do i = 1,maxtype
          do j = i,maxtype
            noncoeff3(i,j) = cutlj
          enddo
        enddo
        if (mixflag.eq.0) mixstyle = 1
      else if (match('lj/smooth',work,mm)) then
        nonstyle = 2
        cutljinterior = realread(str,m)
        cutlj = realread(str,m)
        if (cutljinterior.le.0.0.or.cutljinterior.ge.cutlj)
     $       call error('Bad nonbond style parameter')
        do i = 1,maxtype
          do j = i,maxtype
            noncoeff3(i,j) = cutljinterior
            noncoeff4(i,j) = cutlj
          enddo
        enddo
        if (mixflag.eq.0) mixstyle = 1
      else if (match('lj/shift',work,mm)) then
        nonstyle = 3
        cutlj = realread(str,m)
        offsetflag = intread(str,m)
        if (cutlj.le.0.0.or.offsetflag.lt.0.or.offsetflag.gt.1)
     $       call error('Bad nonbond style parameter')
        do i = 1,maxtype
          do j = i,maxtype
            noncoeff4(i,j) = cutlj
          enddo
        enddo
        if (mixflag.eq.0) mixstyle = 1
      else if (match('soft',work,mm)) then
        nonstyle = 4
        cutlj = realread(str,m)
        if (cutlj.le.0.0) call error('Bad nonbond style parameter')
        do i = 1,maxtype
          do j = i,maxtype
            noncoeff3(i,j) = cutlj
          enddo
        enddo
        mixstyle = 1
      else if (match('class2/cutoff',work,mm)) then
        nonstyle = 5
        cutlj = realread(str,m)
        offsetflag = intread(str,m)
        if (cutlj.le.0.0.or.offsetflag.lt.0.or.offsetflag.gt.1)
     $       call error('Bad nonbond style parameter')
        do i = 1,maxtype
          do j = i,maxtype
            noncoeff3(i,j) = cutlj
          enddo
        enddo
        if (mixflag.eq.0) mixstyle = 3
      else
        call error('Bad nonbond style parameter')
      endif

      if (node.eq.0) then
        if (nonstyle.eq.0) then
          write (1,*) 'Nonbond style none'
        else if (nonstyle.eq.1) then
          write (1,*) 'Nonbond style lj/cutoff',
     $         cutlj,offsetflag
        else if (nonstyle.eq.2) then
          write (1,*) 'Nonbond style lj/smooth',
     $         cutljinterior,cutlj
        else if (nonstyle.eq.3) then
          write (1,*) 'Nonbond style lj/shift',
     $         cutlj,offsetflag
        else if (nonstyle.eq.4) then
          write (1,*) 'Nonbond style soft',
     $         cutlj
        else if (nonstyle.eq.5) then
          write (1,*) 'Nonbond style class2/cutoff',
     $         cutlj,offsetflag
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_nonbond_coeff(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      if (nonstyle.eq.0) call error('Cannot define nonbond '//
     $     'coeffs when nonbond style = none')
      i = intread(str,m)
      j = intread(str,m)
      if (i.le.0.or.i.gt.maxtype.or.j.lt.i.or.j.gt.maxtype)
     $     call error('Bad nonbond coeff parameter')
      nontypeflag(i,j) = nonstyle
      if (nonstyle.eq.1) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
      else if (nonstyle.eq.2) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
        noncoeff4(i,j) = realread(str,m)
      else if (nonstyle.eq.3) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
        noncoeff4(i,j) = realread(str,m)
      else if (nonstyle.eq.4) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
      else if (nonstyle.eq.5) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
      endif
      if (nonstyle.eq.1) then
        if (noncoeff1(i,j).lt.0.0.or.noncoeff2(i,j).lt.0.0.or.
     $       noncoeff3(i,j).lt.0.0)
     $       call error('Bad nonbond coeff parameter')
      else if (nonstyle.eq.2) then
        if (noncoeff1(i,j).lt.0.0.or.noncoeff2(i,j).lt.0.0.or.
     $       noncoeff3(i,j).lt.0.0.or.noncoeff4(i,j).lt.0.0)
     $       call error('Bad nonbond coeff parameter')
        if (noncoeff3(i,j).eq.0.0.and.noncoeff4(i,j).eq.0.0) then
          continue
        else
          if (noncoeff3(i,j).eq.0.0.or.
     $         noncoeff3(i,j).ge.noncoeff4(i,j))
     $         call error('Bad nonbond coeff parameter')
        endif
      else if (nonstyle.eq.3) then
        if (noncoeff1(i,j).lt.0.0.or.noncoeff2(i,j).lt.0.0.or.
     $       noncoeff4(i,j).lt.0.0)
     $       call error('Bad nonbond coeff parameter')
      else if (nonstyle.eq.4) then
        if (noncoeff1(i,j).lt.0.0.or.noncoeff2(i,j).lt.0.0.or.
     $       noncoeff3(i,j).lt.0.0)
     $       call error('Bad nonbond coeff parameter')
      else if (nonstyle.eq.5) then
        if (noncoeff1(i,j).lt.0.0.or.noncoeff2(i,j).lt.0.0.or.
     $       noncoeff3(i,j).lt.0.0)
     $       call error('Bad nonbond coeff parameter')
      endif
      if (node.eq.0) then
        if (nonstyle.eq.1) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
        else if (nonstyle.eq.2) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),
     $         noncoeff3(i,j),noncoeff4(i,j)
        else if (nonstyle.eq.3) then
          write (1,*) 'Nonbond coeff',
     $         noncoeff1(i,j),noncoeff2(i,j),
     $         noncoeff3(i,j),noncoeff4(i,j)
        else if (nonstyle.eq.4) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
        else if (nonstyle.eq.5) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_mixing_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad mixing style parameter')
      mixflag = 1
      if (match('geometric',work,mm)) then
        mixstyle = 1
      else if (match('arithmetic',work,mm)) then
        mixstyle = 2
      else if (match('sixthpower',work,mm)) then
        mixstyle = 3
      else
        call error('Bad mixing style parameter')
      endif
      if (node.eq.0) then
        if (mixstyle.eq.1) then
          write (1,*) 'Mixing style geometric'
        else if (mixstyle.eq.2) then
          write (1,*) 'Mixing style arithmetic'
        else if (mixstyle.eq.3) then
          write (1,*) 'Mixing style sixthpower'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_coulomb_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad coulomb style parameter')
      if (match('none',work,mm)) then
        coulstyle = 0
        cutcoul = 0.0
      else if (match('cutoff',work,mm)) then
        coulstyle = 1
        cutcoul = realread(str,m)
        if (cutcoul.le.0.0)
     $       call error('Bad coulomb style parameter')
      else if (match('smooth',work,mm)) then
        coulstyle = 2
        cutcoulint = realread(str,m)
        cutcoul = realread(str,m)
        if (cutcoulint.le.0.0.or.cutcoulint.ge.cutcoul)
     $       call error('Bad coulomb style parameter')
      else if (match('ewald',work,mm)) then
        coulstyle = 3
        cutcoul = realread(str,m)
        long_prec = realread(str,m)
        if (cutcoul.le.0.0.or.
     $       long_prec.le.0.0.or.long_prec.gt.1.0)
     $       call error('Bad coulomb style parameter')
      else if (match('pppm',work,mm)) then
        coulstyle = 4
        cutcoul = realread(str,m)
        long_prec = realread(str,m)
        if (cutcoul.le.0.0.or.
     $       long_prec.le.0.0.or.long_prec.gt.1.0)
     $       call error('Bad coulomb style parameter')
      else
        call error('Bad coulomb style parameter')
      endif
      if (node.eq.0) then
        if (coulstyle.eq.0) then
          write (1,*) 'Coulomb style none'
        else if (coulstyle.eq.1) then
          write (1,*) 'Coulomb style cutoff',cutcoul
        else if (coulstyle.eq.2) then
          write (1,*) 'Coulomb style smooth',cutcoulint,cutcoul
        else if (coulstyle.eq.3) then
          write (1,*) 'Coulomb style ewald',cutcoul,long_prec
        else if (coulstyle.eq.4) then
          write (1,*) 'Coulomb style pppm',cutcoul,long_prec
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_pppm_mesh(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      nx_pppm_input = intread(str,m)
      if (nx_pppm_input.eq.0) then
        meshflag = 0
        if (node.eq.0) write (1,*) 'Pppm mesh 0'
        return
      endif
        
      meshflag = 1
      ny_pppm_input = intread(str,m)
      nz_pppm_input = intread(str,m)
      if (nx_pppm_input.le.0.or.ny_pppm_input.le.0.or.
     $     nz_pppm_input.le.0)
     $     call error('Bad pppm mesh parameter')
      if (node.eq.0) write (1,*) 'Pppm mesh',
     $     nx_pppm_input,ny_pppm_input,nz_pppm_input

      return
      end

c ------------------------------------------------------------------------

      subroutine in_pppm_order(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      orderflag = intread(str,m)
      if (orderflag.le.0) call error('Bad pppm order parameter')
      if (node.eq.0)
     $     write (1,*) 'Pppm order',orderflag

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dielectric(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      dielectric = realread(str,m)
      if (dielectric.le.0.0)
     $     call error('Bad dielectric parameter')
      if (node.eq.0)
     $     write (1,*) 'Dielectric',dielectric

      return
      end

c ------------------------------------------------------------------------

      subroutine in_bond_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad bond style parameter')
      if (match('none',work,mm)) then
        bondstyle = 0
        do i = 1,maxbondtype
          bondtypeflag(i) = 0
        enddo
      else if (match('harmonic',work,mm)) then
        bondstyle = 1
      else if (match('fene/standard',work,mm)) then
        bondstyle = 2
      else if (match('fene/shift',work,mm)) then
        bondstyle = 3
      else if (match('nonlinear',work,mm)) then
        bondstyle = 4
      else if (match('class2',work,mm)) then
        bondstyle = 5
      else
        call error('Bad bond style parameter')
      endif
      if (node.eq.0) then
        if (bondstyle.eq.0) then
          write (1,*) 'Bond style none'
        else if (bondstyle.eq.1) then
          write (1,*) 'Bond style harmonic'
        else if (bondstyle.eq.2) then
          write (1,*) 'Bond style fene/standard'
        else if (bondstyle.eq.3) then
          write (1,*) 'Bond style fene/shift'
        else if (bondstyle.eq.4) then
          write (1,*) 'Bond style nonlinear'
        else if (bondstyle.eq.5) then
          write (1,*) 'Bond style class2'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_bond_coeff(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      if (bondstyle.eq.0) call error(
     $     'Cannot define bond coeffs when bond style = none')
      if (bondstyle.eq.5) call error(
     $     'Cannot define bond coeffs when bond style = class2')
      i = intread(str,m)
      if (i.le.0.or.i.gt.maxbondtype)
     $     call error('Bad bond coeff parameter')
      bondtypeflag(i) = bondstyle
      if (bondstyle.eq.1) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
      else if (bondstyle.eq.2) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
        bondcoeff(3,i) = realread(str,m)
        bondcoeff(4,i) = realread(str,m)
      else if (bondstyle.eq.3) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
        bondcoeff(3,i) = realread(str,m)
        bondcoeff(4,i) = realread(str,m)
        bondcoeff(5,i) = realread(str,m)
      else if (bondstyle.eq.4) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
        bondcoeff(3,i) = realread(str,m)
      endif
      if (bondstyle.eq.1) then
        if (bondcoeff(1,i).lt.0.0.or.
     $       bondcoeff(2,i).lt.0.0)
     $       call error('Bad bond coeff parameter')
      else if (bondstyle.eq.2) then
        if (bondcoeff(1,i).lt.0.0.or.
     $       bondcoeff(2,i).le.0.0.or.
     $       bondcoeff(3,i).lt.0.0.or.
     $       bondcoeff(4,i).le.0.0)
     $       call error('Bad bond coeff parameter')
      else if (bondstyle.eq.3) then
        if (bondcoeff(1,i).lt.0.0.or.
     $       bondcoeff(2,i).le.0.0.or.
     $       bondcoeff(3,i).lt.0.0.or.
     $       bondcoeff(4,i).le.0.0)
     $       call error('Bad bond coeff parameter')
      else if (bondstyle.eq.4) then
        if (bondcoeff(1,i).lt.0.0.or.
     $       bondcoeff(2,i).le.0.0.or.
     $       bondcoeff(3,i).le.0.0)
     $       call error('Bad bond coeff parameter')
      endif
      if (node.eq.0) then
        if (bondstyle.eq.1) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,2)
        else if (bondstyle.eq.2) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,4)
        else if (bondstyle.eq.3) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,5)
        else if (bondstyle.eq.4) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,3)
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_angle_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad angle style parameter')
      if (match('none',work,mm)) then
        anglestyle = 0
        do i = 1,maxangletype
          angletypeflag(i) = 0
        enddo
      else if (match('harmonic',work,mm)) then
        anglestyle = 1
      else if (match('class2',work,mm)) then
        anglestyle = 2
      else
        call error('Bad angle style parameter')
      endif
      if (node.eq.0) then
        if (anglestyle.eq.0) then
          write (1,*) 'Angle style none'
        else if (anglestyle.eq.1) then
          write (1,*) 'Angle style harmonic'
        else if (anglestyle.eq.2) then
          write (1,*) 'Angle style class2'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dihedral_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad dihedral style parameter')
      if (match('none',work,mm)) then
        dihedstyle = 0
        do i = 1,maxdihedtype
          dihedtypeflag(i) = 0
        enddo
      else if (match('harmonic',work,mm)) then
        dihedstyle = 1
      else if (match('class2',work,mm)) then
        dihedstyle = 2
      else
        call error('Bad dihedral style parameter')
      endif
      if (node.eq.0) then
        if (dihedstyle.eq.0) then
          write (1,*) 'Dihedral style none'
        else if (dihedstyle.eq.1) then
          write (1,*) 'Dihedral style harmonic'
        else if (dihedstyle.eq.2) then
          write (1,*) 'Dihedral style class2'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_improper_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad improper style parameter')
      if (match('none',work,mm)) then
        improstyle = 0
        do i = 1,maximprotype
          improtypeflag(i) = 0
        enddo
      else if (match('harmonic',work,mm)) then
        improstyle = 1
      else if (match('cvff',work,mm)) then
        improstyle = 2
      else if (match('class2',work,mm)) then
        improstyle = 3
      else
        call error('Bad improper style parameter')
      endif
      if (node.eq.0) then
        if (improstyle.eq.0) then
          write (1,*) 'Improper style none'
        else if (improstyle.eq.1) then
          write (1,*) 'Improper style harmonic'
        else if (improstyle.eq.2) then
          write (1,*) 'Improper style cvff'
        else if (improstyle.eq.3) then
          write (1,*) 'Improper style class2'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_min_style(str,m)
      include "lammps.h"
      character*(*)  str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (match('hftn',work,mm)) then
        optstyle = 1
      else
        call error('Bad min style parameter')
      endif
      if (node.eq.0) then
        if (optstyle.eq.1) write (1,*) 'Min style htfn'
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_min_file(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      call strread(str,m,opt_outfile)
      if (length(opt_outfile).le.0)
     $     call error('Bad min file parameter')
      optfileflag = 1
      if (node.eq.0)
     $     write (1,*) 'Min file ',
     $     opt_outfile(1:length(opt_outfile))

      return
      end

c ------------------------------------------------------------------------

      subroutine in_read_data(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      call strread(str,m,datafile)
      if (length(datafile).le.0)
     $     call error('Bad read data parameter')
      if (node.eq.0)
     $     write (1,*) 'Read data ',datafile(1:length(datafile))

      return
      end

c ------------------------------------------------------------------------

      subroutine in_create_group(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad create group parameter')
      if (match('types',work,mm)) then
        creategroup = 1
        createtypelo = intread(str,m)
        createtypehi = intread(str,m)
        if (createtypelo.le.0.or.createtypehi.lt.createtypelo)
     $       call error('Bad create style parameter')
      else if (match('region',work,mm)) then
        creategroup = 2
        toggle = -1.0
        do i = 1,6
          mm = m
          call strread(str,m,work)
          if (work(1:3).eq.'INF') then
            createregion(i) = toggle*1.0E20
          else
            m = mm
            createregion(i) = realread(str,m)
          endif
          toggle = -toggle
        enddo
        if (createregion(1).gt.createregion(2).or.
     $       createregion(3).gt.createregion(4).or.
     $       createregion(5).gt.createregion(6))
     $       call error('Bad create group parameter')
      else if (match('remainder',work,mm)) then
        creategroup = 3
      else
        call error('Bad create group parameter')
      endif
      if (node.eq.0) then
        if (creategroup.eq.1) then
          write (1,*) 'Create group types',
     $         createtypelo,createtypehi
        else if (creategroup.eq.2) then
          write (1,*) 'Create group region',
     $         createregion(1),createregion(2),createregion(3),
     $         createregion(4),createregion(5),createregion(6)
        else if (creategroup.eq.3) then
          write (1,*) 'Create group remainder'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_create_temp(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad create temp parameter')
      if (match('uniform',work,mm)) then
        createstyle = 1
        t_create = realread(str,m)
        iseed = intread(str,m)
        if (t_create.lt.0.0.or.iseed.le.0.or.iseed.ge.1000000000)
     $       call error('Bad create temp parameter')
      else if (match('gaussian',work,mm)) then
        createstyle = 2
        t_create = realread(str,m)
        iseed = intread(str,m)
        if (t_create.lt.0.0.or.iseed.le.0.or.iseed.ge.1000000000)
     $       call error('Bad create temp parameter')
      else if (match('velocity',work,mm)) then
        createstyle = 3
        createvec(1) = realread(str,m)
        createvec(2) = realread(str,m)
        createvec(3) = realread(str,m)
      else
        call error('Bad create temp parameter')
      endif
      if (node.eq.0) then
        if (createstyle.eq.1) then
          write (1,*) 'Create temp uniform ',t_create,iseed
        else if (createstyle.eq.2) then
          write (1,*) 'Create temp gaussian ',t_create,iseed
        else if (createstyle.eq.2) then
          write (1,*) 'Create temp velocity ',
     $         createvec(1),createvec(2),createvec(3)
        endif
      endif
      if (createstyle.eq.3) then
        createvec(1) = createvec(1)*dtfactor
        createvec(2) = createvec(2)*dtfactor
        createvec(3) = createvec(3)*dtfactor
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_read_restart(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      call strread(str,m,restart_in)
      if (length(restart_in).le.0)
     $     call error('Bad read restart parameter')
      if (node.eq.0)
     $     write (1,*) 'Read restart ',
     $     restart_in(1:length(restart_in))

      return
      end

c ------------------------------------------------------------------------

      subroutine in_fix_style(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      mm = m
      call strread(str,mm,work)
      if (length(work).le.0)
     $     call error('Bad fix style parameter')
      if (match('none',work,mm)) then
        nfixes = 0
        do i = 1,maxfix
          fixstyle(i) = 0
        enddo
        if (readflag.eq.1) then
          i = atompnt
          do ii = 1,nlocal
            fix(i) = 0
            i = list(i)
          enddo
        endif
      else
        iwhich = intread(str,m)
        if (iwhich.le.0) call error('Bad fix style parameter')
        if (iwhich.ge.maxfix)
     $       call error('Too many fixes - boost maxfix')
        nfixes = max(nfixes,iwhich)
        call strread(str,m,work)
        if (match('setforce',work,mm)) then
          fixstyle(iwhich) = 1
          do i = 1,3
            mm = m
            call strread(str,m,work)
            if (work(1:4).eq.'NULL') then
              fixflag(i,iwhich) = 1
              fixcoeff(i,iwhich) = 0.0
            else
              fixflag(i,iwhich) = 0
              m = mm
              fixcoeff(i,iwhich) = realread(str,m)
            endif
          enddo
        else if (match('addforce',work,mm)) then
          fixstyle(iwhich) = 2
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)
        else if (match('aveforce',work,mm)) then
          fixstyle(iwhich) = 3
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)
        else if (match('rescale',work,mm)) then
          fixstyle(iwhich) = 4
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = intread(str,m)
          fixcoeff(4,iwhich) = realread(str,m)
          if (fixcoeff(1,iwhich).lt.0.0.or.
     $         fixcoeff(2,iwhich).lt.0.0.or.
     $         nint(fixcoeff(3,iwhich)).le.0.or.
     $         fixcoeff(4,iwhich).lt.0.0)
     $         call error('Bad fix style parameter')
        else if (match('langevin',work,mm)) then
          fixstyle(iwhich) = 5
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)
          iseed = intread(str,m)
          fixflag(1,iwhich) = intread(str,m)
          fixflag(2,iwhich) = intread(str,m)
          fixflag(3,iwhich) = intread(str,m)
          if (fixcoeff(1,iwhich).lt.0.0.or.
     $         fixcoeff(2,iwhich).lt.0.0.or.
     $         fixcoeff(3,iwhich).lt.0.0.or.
     $         iseed.le.0.or.iseed.ge.1000000000)
     $         call error('Bad fix style parameter')
          if (fixflag(1,iwhich).lt.0.or.fixflag(1,iwhich).gt.1.or.
     $         fixflag(2,iwhich).lt.0.or.fixflag(2,iwhich).gt.1.or.
     $         fixflag(3,iwhich).lt.0.or.fixflag(3,iwhich).gt.1)
     $         call error('Bad fix style parameter')
          tmp = ranmars(iseed+node)
        else if (match('nose/hoover',work,mm)) then
          fixstyle(iwhich) = 6
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)
          if (fixcoeff(1,iwhich).lt.0.0.or.
     $         fixcoeff(2,iwhich).lt.0.0.or.
     $         fixcoeff(3,iwhich).lt.0.0)
     $         call error('Bad fix style parameter')
        else if (match('springforce',work,mm)) then
          fixstyle(iwhich) = 7
          do i = 1,3
            mm = m
            call strread(str,m,work)
            if (work(1:4).eq.'NULL') then
              fixflag(i,iwhich) = 1
              fixcoeff(i,iwhich) = 0.0
            else
              fixflag(i,iwhich) = 0
              m = mm
              fixcoeff(i,iwhich) = realread(str,m)
            endif
          enddo
          fixcoeff(4,iwhich) = realread(str,m)
          if (fixcoeff(4,iwhich).lt.0.0)
     $         call error('Bad fix style parameter')
        else if (match('dragforce',work,mm)) then
          fixstyle(iwhich) = 8
          do i = 1,3
            mm = m
            call strread(str,m,work)
            if (work(1:4).eq.'NULL') then
              fixflag(i,iwhich) = 1
              fixcoeff(i,iwhich) = 0.0
            else
              fixflag(i,iwhich) = 0
              m = mm
              fixcoeff(i,iwhich) = realread(str,m)
            endif
          enddo
          fixcoeff(4,iwhich) = realread(str,m)
          fixcoeff(5,iwhich) = realread(str,m)
          if (fixcoeff(4,iwhich).lt.0.0.or.fixcoeff(5,iwhich).lt.0.0)
     $         call error('Bad fix style parameter')
        else
          call error('Bad fix style parameter')
        endif
      endif
      if (node.eq.0) then
        if (nfixes.eq.0) then
          write (1,*) 'Fix style none'
        else if (fixstyle(iwhich).eq.1) then
          write (1,*) 'Fix style',iwhich,' setforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich).eq.2) then
          write (1,*) 'Fix style',iwhich,' addforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich).eq.3) then
          write (1,*) 'Fix style',iwhich,' aveforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich).eq.4) then
          write (1,*) 'Fix style',iwhich,' rescale',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         nint(fixcoeff(3,iwhich)),fixcoeff(4,iwhich)
        else if (fixstyle(iwhich).eq.5) then
          write (1,*) 'Fix style',iwhich,' langevin',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich),iseed,
     $         fixflag(1,iwhich),fixflag(2,iwhich),fixflag(3,iwhich)
        else if (fixstyle(iwhich).eq.6) then
          write (1,*) 'Fix style',iwhich,' nose/hoover',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich).eq.7) then
          write (1,*) 'Fix style',iwhich,' springforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich),fixcoeff(4,iwhich)
        else if (fixstyle(iwhich).eq.8) then
          write (1,*) 'Fix style',iwhich,' dragforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich),fixcoeff(4,iwhich),
     $         fixcoeff(5,iwhich)
        endif
      endif
      if (nfixes.gt.0) then
        if (fixstyle(iwhich).eq.1.or.fixstyle(iwhich).eq.2.or.
     $       fixstyle(iwhich).eq.3) then
          fixcoeff(1,iwhich) = fixcoeff(1,iwhich)*dtfactor*dtfactor
          fixcoeff(2,iwhich) = fixcoeff(2,iwhich)*dtfactor*dtfactor
          fixcoeff(3,iwhich) = fixcoeff(3,iwhich)*dtfactor*dtfactor
        endif
        if (fixstyle(iwhich).eq.5)
     $       fixcoeff(3,iwhich) = fixcoeff(3,iwhich)*dtfactor
        if (fixstyle(iwhich).eq.6)
     $       fixcoeff(3,iwhich) = fixcoeff(3,iwhich)*dtfactor
        if (fixstyle(iwhich).eq.7)
     $       fixcoeff(4,iwhich) = fixcoeff(4,iwhich)*dtfactor*dtfactor
        if (fixstyle(iwhich).eq.8)
     $     fixcoeff(4,iwhich) = fixcoeff(4,iwhich)*dtfactor*dtfactor
      endif
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_assign_fix(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      fixwhich = intread(str,m)
      if (fixwhich.le.0.or.fixwhich.ge.maxfix)
     $     call error('Bad assign fix parameter')
      call strread(str,m,work)
      if (length(work).le.0)
     $     call error('Bad assign fix parameter')
      if (match('atom',work,mm)) then
        fixgroup = 1
        fixatom = intread(str,m)
        if (fixatom.le.0.or.fixatom.gt.natoms) 
     $       call error('Bad assign fix parameter')
      else if (match('molecule',work,mm)) then
        fixgroup = 2
        fixatom = intread(str,m)
        if (fixatom.le.0) 
     $       call error('Bad assign fix parameter')
      else if (match('type',work,mm)) then
        fixgroup = 3
        fixtype = intread(str,m)
        if (fixtype.le.0) call error('Bad assign fix parameter')
      else if (match('region',work,mm)) then
        fixgroup = 4
        toggle = -1.0
        do i = 1,6
          mm = m
          call strread(str,m,work)
          if (work(1:3).eq.'INF') then
            fixregion(i) = toggle*1.0E20
          else
            m = mm
            fixregion(i) = realread(str,m)
          endif
          toggle = -toggle
        enddo
        if (fixregion(1).gt.fixregion(2).or.
     $       fixregion(3).gt.fixregion(4).or.
     $       fixregion(5).gt.fixregion(6))
     $       call error('Bad assign fix parameter')
      else if (match('remainder',work,mm)) then
        fixgroup = 5
      else
        call error('Bad assign fix parameter')
      endif
      if (node.eq.0) then
        if (fixgroup.eq.1) then
          write (1,*) 'Assign fix',fixwhich,' atom',fixatom
        else if (fixgroup.eq.2) then
          write (1,*) 'Assign fix',fixwhich,' molecule',fixatom
        else if (fixgroup.eq.3) then
          write (1,*) 'Assign fix',fixwhich,' type',fixtype
        else if (fixgroup.eq.4) then
          write (1,*) 'Assign fix',fixwhich,' region',
     $         (fixregion(i),i=1,6)
        else if (fixgroup.eq.5) then
          write (1,*) 'Assign fix',fixwhich,' remainder'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_reset_timestep(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      ntimestep = intread(str,m)
      if (node.eq.0)
     $     write (1,*) 'Reset timestep',ntimestep

      return
      end

c ------------------------------------------------------------------------

      subroutine in_run(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match

      nsteps = intread(str,m)
      if (nsteps.lt.0) call error('Bad run parameter')
      if (node.eq.0)
     $     write (1,*) 'Run',nsteps

      return
      end

c ------------------------------------------------------------------------

      subroutine in_minimize(str,m)
      include "lammps.h"
      character*(*) str
      integer m

      external length
      logical match
      character*16 work

      opt_stop_tol = realread(str,m)
      opt_max_iters = intread(str,m)
      opt_max_fns = intread(str,m)
      nsteps = 1
      if (opt_stop_tol.le.0.0.or.
     $     opt_max_iters.lt.0.or.opt_max_fns.lt.0)
     $     call error('Bad minimize parameter')
      if (node.eq.0)
     $     write (1,*) 'Minimize',
     $     opt_stop_tol,opt_max_iters,opt_max_fns

      return
      end
