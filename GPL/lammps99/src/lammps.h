c LAMMPS 99 - Molecular Dynamics Simulator, Copyright 1998, 1999
c Authored by Steve Plimpton
c   (505) 845-7873, sjplimp@cs.sandia.gov
c   Dept 9221, MS 1111, Sandia National Labs, Albuquerque, NM 87185-1111
c See the README file for more information

      implicit real*8 (a-h,o-z)

      include 'param.h'

      parameter (maxperatom=25)
      parameter (maxatom=maxown+maxother)
      parameter (maxatomp1=maxatom+1)
      parameter (maxexchtot=maxexch*(maxperatom+maxsneigh+
     $     3*maxbondper+4*maxangleper+5*maxdihedper+5*maximproper))
      parameter (maxrestot=maxown*(maxperatom-3+3*maxbondper+
     $     4*maxangleper+5*maxdihedper+5*maximproper)+1)
      parameter (maxspec=2*maxsneigh*maxown)

      real*8 buf1(maxexchtot),buf2(2*maxexchtot)
      real*8 buf3(3*maxsendone),buf4(8*maxown),buf5(maxrestot)
      real*8 buf6(maxsendone)

      integer ibuf1(maxsendone),ibuf2(maxsendone)
      integer ibuf3(maxspec),ibuf4(maxspec)

      real*8 anglecoeff(4,maxangletype)
      real*8 angleanglecoeff(6,maximprotype)
      real*8 angleangletorsioncoeff(3,maxdihedtype)
      real*8 angletorsioncoeff(8,maxdihedtype)
      real*8 bondanglecoeff(4,maxangletype)
      real*8 bondbondcoeff(3,maxangletype)
      real*8 bondbond13coeff(3,maxdihedtype)
      real*8 bondcoeff(5,maxbondtype)
      real*8 boundhi(maxswap),boundlo(maxswap)
      real*8 cutforcesq(maxtype,maxtype),cutljsq(maxtype,maxtype)
      real*8 cutljinner(maxtype,maxtype)
      real*8 cutljinnersq(maxtype,maxtype)
      real*8 cutneighsq(maxtype,maxtype)
      real*8 dihedcoeff(6,maxdihedtype)
      real*8 endbondtorsioncoeff(8,maxdihedtype)
      real*8 f(3,maxatom),fixcoeff(8,maxfix),fixstore(5*maxfix)
      real*8 improcoeff(3,maximprotype)
      real*8 lj1(maxtype,maxtype),lj2(maxtype,maxtype)
      real*8 lj3(maxtype,maxtype),lj4(maxtype,maxtype)
      real*8 lj5(maxtype,maxtype)
      real*8 ljsw0(maxtype,maxtype)
      real*8 ljsw1(maxtype,maxtype),ljsw2(maxtype,maxtype)
      real*8 ljsw3(maxtype,maxtype),ljsw4(maxtype,maxtype)
      real*8 mass(maxtype)
      real*8 midbondtorsioncoeff(4,maxdihedtype)
      real*8 noncoeff1(maxtype,maxtype),noncoeff2(maxtype,maxtype)
      real*8 noncoeff3(maxtype,maxtype),noncoeff4(maxtype,maxtype)
      real*8 offset(maxtype,maxtype)
      real*8 q(maxatom)
      real*8 v(3,maxown),x(3,maxatom),xhold(3,maxown)

      integer angleatom1(maxangleper,maxown)
      integer angleatom2(maxangleper,maxown)
      integer angleatom3(maxangleper,maxown)
      integer anglelist(4,maxangle)
      integer angletype(maxangleper,maxown)
      integer angletypeflag(maxangletype)
      integer bin(maxatom),binpnt(maxbin)
      integer bondatom1(maxbondper,maxown)
      integer bondatom2(maxbondper,maxown)
      integer bondlist(3,maxbond)
      integer bondtype(maxbondper,maxown)
      integer bondtypeflag(maxbondtype)
      integer dihedatom1(maxdihedper,maxown)
      integer dihedatom2(maxdihedper,maxown)
      integer dihedatom3(maxdihedper,maxown)
      integer dihedatom4(maxdihedper,maxown)
      integer dihedlist(5,maxdihed)
      integer dihedtype(maxdihedper,maxown)
      integer dihedtypeflag(maxdihedtype)
      integer fix(maxown),fixflag(3,maxfix)
      integer fixptr(maxfix),fixstyle(maxfix)
      integer improatom1(maximproper,maxown)
      integer improatom2(maximproper,maxown)
      integer improatom3(maximproper,maxown)
      integer improatom4(maximproper,maxown)
      integer improlist(5,maximpro)
      integer improtype(maximproper,maxown)
      integer improtypeflag(maximprotype)
      integer list(maxown),localptr(maxtotal),molecule(maxown)
      integer nlist(maxown*maxneigh+maxneigh)
      integer nliststart(maxown),nliststop(maxown)
      integer nontypeflag(maxtype,maxtype),nrlist(maxswap+1)
      integer nslist(maxswap+1),numangle(maxown),numbond(maxown)
      integer num1bond(maxown),num2bond(maxown),num3bond(maxown)
      integer numdihed(maxown),numimpro(maxown)
      integer rpart(maxswap),slist(maxsend)
      integer spart(maxswap),specbond(maxsneigh,maxown)
      integer tag(maxatom),true(maxown),type(maxatom)
      integer velflag(maxown)

      double precision time_angle,time_bond,time_comm,time_current
      double precision time_dihedral,time_exch,time_fcomm
      double precision time_improper
      double precision time_io
      double precision time_loop,time_neigh1,time_neigh2
      double precision time_nonbond,time_other,time_total

      real*8 binsizex,binsizey,binsizez,boltz,border(2,3)
      real*8 coulpre,createregion(6),createvec(3),cutcoul
      real*8 cutcoulint,cutcoulintsq,cutcoulsq,cutforce
      real*8 cutlj,cutljinterior,cutneigh
      real*8 dielectric,dt,dtfactor,dthalf
      real*8 efactor,e_angle,e_bond,e_coul,e_dihedral
      real*8 e_improper,e_total,e_vdwl,fixregion(6)
      real*8 skin,special(3),triggersq,two_1_3
      real*8 t_create,t_current
      real*8 t_nph,t_start,t_stop,t_window
      real*8 xmc,ymc,zmc,xms,yms,zms,xboundlo,yboundlo,zboundlo
      real*8 xboundhi,yboundhi,zboundhi,xprd,yprd,zprd
      real*8 e_bondbond,e_bondbond13,e_bondangle,e_endbondtorsion
      real*8 e_midbondtorsion,e_angletorsion
      real*8 e_angleangletorsion,e_angleangle

      integer anglestyle,atompnt
      integer bondstyle,boxflag,coulstyle,creategroup,createstyle
      integer createtypehi,createtypelo
      integer dihedstyle,fixatom,fixgroup
      integer dumpatomfileflag,dumpvelfileflag,dumpforcefileflag
      integer fixnum,fixtype,fixwhich,freepnt
      integer idimension,improstyle,iseed,itime,iversion
      integer max_angle,max_angleper,max_bond,max_bondper
      integer max_dihed,max_dihedper,max_exch,max_impro
      integer max_improper,max_nlocal
      integer max_neigh,max_nother,max_slist,max_swap
      integer mbinx,mbiny,mbinz,mbinxlo,mbinylo,mbinzlo
      integer me(3),mixflag,mixstyle,mpart(2,3)
      integer nanglelocal,nangles,nangletypes,natoms
      integer nbinx,nbiny,nbinz,nbondlocal,nbonds,nbondtypes
      integer ndanger
      integer ndihedlocal,ndihedrals,ndihedtypes
      integer need(3),neighago,neighdelay,neighfreq
      integer neighstyle,neightop,neightrigger
      integer newton,newton_bond,newton_nonbond
      integer nfixes,nimprolocal
      integer nimpropers,nimprotypes,nlocal,nother,node,nonstyle
      integer nprocs,nsteps,nswap,ntimestep
      integer nthermo,nthermo_next,thermostyle
      integer noutput_next,ntime_last
      integer ndumpatom,ndumpatom_prev,ndumpatom_next
      integer ndumpvel,ndumpvel_prev,ndumpvel_next
      integer ndumpforce,ndumpforce_prev,ndumpforce_next
      integer nrestart,nrestart_next
      integer ntypes,numneigh,offsetflag
      integer peratom_comm,peratom_io,perflagx,perflagy,perflagz
      integer pgrid(3)
      integer readflag,restartfileflag
      integer t_every,tempflag
      integer trueflag,units

      character*80 datafile,dumpatomfile,dumpvelfile,dumpforcefile
      character*80 restart_in,restart_out1,restart_out2

      common /bk00/ buf1
      common /bk01/ buf2
      common /bk02/ buf3
      common /bk03/ buf4
      common /bk04/ buf5
      common /bk05/ buf6
      common /bk07/ ibuf1
      common /bk08/ ibuf2
      common /bk09/ ibuf3
      common /bk10/ ibuf4

      common /bk20/ anglecoeff,angleanglecoeff,
     $     angleangletorsioncoeff,angletorsioncoeff,
     $     bondanglecoeff,bondbondcoeff,bondbond13coeff,bondcoeff
      common /bk21/ boundhi,boundlo
      common /bk22/ cutforcesq,cutljsq,cutljinner,cutljinnersq
      common /bk23/ cutneighsq
      common /bk24/ dihedcoeff,endbondtorsioncoeff
      common /bk25/ f,fixcoeff,fixstore,improcoeff
      common /bk26/ lj1,lj2,lj3,lj4,lj5
      common /bk27/ ljsw0,ljsw1,ljsw2,ljsw3,ljsw4
      common /bk28/ mass,midbondtorsioncoeff,
     $     noncoeff1,noncoeff2,noncoeff3,noncoeff4
      common /bk29/ offset,q
      common /bk30/ v
      common /bk31/ x,xhold

      common /bk40/ angleatom1,angleatom2,angleatom3
      common /bk41/ anglelist,angletype,angletypeflag
      common /bk42/ bin,binpnt
      common /bk43/ bondatom1,bondatom2,bondlist
      common /bk44/ bondtype,bondtypeflag
      common /bk46/ dihedatom1,dihedatom2,dihedatom3,dihedatom4
      common /bk47/ dihedlist,dihedtype,dihedtypeflag
      common /bk48/ fix,fixflag,fixptr,fixstyle
      common /bk49/ improatom1,improatom2,improatom3,improatom4
      common /bk50/ improlist,improtype,improtypeflag
      common /bk51/ list,localptr,molecule
      common /bk52/ nlist,nliststart,nliststop,nontypeflag
      common /bk53/ nrlist,nslist
      common /bk54/ numangle,numbond,num1bond,num2bond,num3bond
      common /bk55/ numdihed,numimpro,rpart,slist,spart
      common /bk56/ specbond,tag,true,type
      common /bk57/ velflag

      common /bk90/ time_angle,time_bond,time_comm,time_current,
     $     time_dihedral,time_exch,time_fcomm,
     $     time_improper,time_io,time_loop,
     $     time_neigh1,time_neigh2,time_nonbond,
     $     time_other,time_total

      common /bk91/ binsizex,binsizey,binsizez,boltz,border,coulpre,
     $     createregion,createvec,cutcoul,cutcoulint,cutcoulintsq,
     $     cutcoulsq,cutforce,cutlj,cutljinterior,
     $     cutneigh,dielectric,dt,dtfactor,dthalf,efactor,
     $     e_angle,e_bond,e_coul,e_dihedral,e_improper,e_total,
     $     e_vdwl,e_bondbond,e_bondbond13,
     $     e_bondangle,e_endbondtorsion,
     $     e_midbondtorsion,e_angletorsion,
     $     e_angleangletorsion,e_angleangle,
     $     fixregion,
     $     skin,special,triggersq,two_1_3,
     $     t_create,t_current,
     $     t_nph,t_start,t_stop,t_window,
     $     xmc,ymc,zmc,xms,yms,zms,xboundlo,yboundlo,zboundlo,
     $     xboundhi,yboundhi,zboundhi,xprd,yprd,zprd

      common /bk92/ anglestyle,atompnt,
     $     bondstyle,boxflag,coulstyle,creategroup,createstyle,
     $     createtypehi,createtypelo,
     $     dihedstyle,fixatom,fixgroup,
     $     dumpatomfileflag,dumpvelfileflag,dumpforcefileflag,
     $     fixnum,fixtype,fixwhich,freepnt,
     $     idimension,improstyle,iseed,itime,iversion,
     $     max_angle,max_angleper,max_bond,max_bondper,
     $     max_dihed,max_dihedper,max_exch,max_impro,
     $     max_improper,max_nlocal,max_neigh,max_nother,
     $     max_slist,max_swap,mbinx,mbiny,mbinz,mbinxlo,
     $     mbinylo,mbinzlo,me,mixflag,mixstyle,mpart,
     $     nanglelocal,nangles,nangletypes
      common /bk93/ natoms,nbinx,nbiny,nbinz,nbondlocal,nbonds,
     $     nbondtypes,
     $     ndanger,ndihedlocal,ndihedrals,
     $     ndihedtypes,need,neighago,neighdelay,
     $     neighfreq,neighstyle,neightop,neightrigger,
     $     newton,newton_bond,newton_nonbond,nfixes,
     $     nimprolocal,nimpropers,nimprotypes,nlocal,nother,
     $     node,nonstyle,
     $     nprocs,nsteps,nswap,ntimestep,
     $     nthermo,nthermo_next,thermostyle,
     $     noutput_next,ntime_last,
     $     ndumpatom,ndumpatom_prev,ndumpatom_next,
     $     ndumpvel,ndumpvel_prev,ndumpvel_next,
     $     ndumpforce,ndumpforce_prev,ndumpforce_next,
     $     nrestart,nrestart_next,
     $     ntypes,numneigh,
     $     offsetflag,peratom_comm,peratom_io,
     $     perflagx,perflagy,perflagz,pgrid,
     $     readflag,restartfileflag,
     $     t_every,tempflag,
     $     trueflag,units

      common /bk94/ datafile,dumpatomfile,dumpvelfile,
     $     dumpforcefile,restart_in,restart_out1,restart_out2

c Ewald data structures

      integer kxvecs(kmaxdim),kyvecs(kmaxdim),kzvecs(kmaxdim)
      real*8 ug(kmaxdim),eg(3,kmaxdim),ek(3,maxown)
      real*8 sfacrl(kmaxdim),sfacim(kmaxdim)
      real*8 sfacrl_all(kmaxdim),sfacim_all(kmaxdim)
      real*8 cs(-nkmax:nkmax,3,maxown),sn(-nkmax:nkmax,3,maxown)
      real*8 time_long,time_rho,time_poiss,time_field,long_prec

      common /bk100/ kxvecs,kyvecs,kzvecs,kcount,kmax
      common /bk101/ ug,eg,ek,sfacrl,sfacim,sfacrl_all,sfacim_all
      common /bk102/ cs,sn
      common /bk103/ e_long,ewald_g,gewsr,long_prec
      common /bk104/ qsum,qsqsum
      common /bk105/ time_long,time_rho,time_poiss,time_field

c PPPM data structures

      integer orderflag,nfft,nlower,nupper
      integer nx_pppm,ny_pppm,nz_pppm
      integer nx_pppm_input,ny_pppm_input,nz_pppm_input
      integer nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in
      integer nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out
      integer nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,
     $     nzlo_ghost,nzhi_ghost
      integer nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft
      real*8 plan1_fft,plan2_fft,plan_remap

      common /bk110/ meshflag,orderflag,nfft,nlower,nupper
      common /bk111/ nx_pppm,ny_pppm,nz_pppm,
     $     nx_pppm_input,ny_pppm_input,nz_pppm_input
      common /bk112/ nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
     $     nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
     $     nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,
     $     nzlo_ghost,nzhi_ghost,
     $     nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft
      common /bk113/ plan1_fft,plan2_fft,plan_remap

      integer partgrid(3,maxown)
      real*8 particle(4,maxown)
      real*8 density_fft(maxfft),density_brick(maxgrid)
      real*8 vdx_brick(maxgrid),vdy_brick(maxgrid),vdz_brick(maxgrid)
      real*8 fkvecs(3,0:8192-1),greensfn(maxfft)
      complex*16 workvec1(maxfft),workvec2(maxfft)

      common /bk120/ partgrid
      common /bk121/ particle
      common /bk122/ density_fft,density_brick
      common /bk123/ vdx_brick,vdy_brick,vdz_brick
      common /bk124/ fkvecs
      common /bk125/ greensfn
      common /bk126/ workvec1,workvec2

c Ensemble data structures

      integer ensemble,pressflag,xpressflag,ypressflag,zpressflag
      real*8 masssum

      common /bk130/ ensemble,pressflag,
     $     xpressflag,ypressflag,zpressflag
      common /bk131/ t_freq,t_target,p_freq2(3),p_total,
     $     p_freq,p_current(3),p_start(3),p_stop(3),p_target(3)
      common /bk133/ eta_dot,omega(3),omega_dot(3),
     $     masssum,pfactor,virial(6)

c RESPA data structures

      integer nrespa,nstretch,nintra,nshort
      real*8 f_stretch(3,maxatom),f_intra(3,maxatom)
      real*8 f_short(3,maxatom),f_long(3,maxatom)
      real*8 vir_stretch(6),vir_intra(6)
      real*8 vir_short(6),vir_long(6)
      real*8 dthalf_intra,dthalf_short,dthalf_long

      common /bk140/ nrespa,nstretch,nintra,nshort
      common /bk141/ f_stretch,f_intra,f_short,f_long
      common /bk142/ vir_stretch,vir_intra,vir_short,vir_long
      common /bk143/ dthalf_intra,dthalf_short,dthalf_long

c Optimization data structures

      character*80 opt_outfile
      integer optstyle,optfileflag,opt_max_iters,opt_max_fns
      real*8 opt_stop_tol,opt_time1,opt_time2,opt_time3

      common /bk150/ opt_outfile
      common /bk151/ optstyle,optfileflag,opt_max_iters,opt_max_fns
      common /bk152/ opt_stop_tol,opt_time1,opt_time2,opt_time3

c Diagnostic data structures

      character*80 diagfile(maxdiag)
      character*16 diagnames(maxdiag)
      integer ndiag(maxdiag),diagfileflag(maxdiag)
      integer diagnparams(maxdiag),diagcall(maxdiag),
     $     diagnext(maxdiag),diagprev(maxdiag)
      real*8 diagparam(5,maxdiag)
      integer numdiag,ndiag_next

      common /bk160/ diagfile,diagnames
      common /bk161/ ndiag,diagfileflag,
     $     diagnparams,diagcall,diagnext,diagprev
      common /bk162/ diagparam
      common /bk163/ numdiag,ndiag_next
