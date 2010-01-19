*=====================================================================
*     MDynaMix v.5.0
*     Part 2
*
*     Input data
C
C     This file contains subroutines responsible for input of data:
C
C     1. INPUT  - input of simulation parameters
C     2. READMOL - reads molecular structure and force field parameters
C     3. TAKESTR (char*128) - read line as character line
C     4. Other units are supplementary functions:
C
C        ANGLE   - calculate angles
C        PARSES  - perform character input in free format  
*
*    1. Input of simulations parameters
*    ----------------------------------
C================ INPUT ===============================================
	subroutine INPUT   
	include "prcm.h"
C   maximal length of input lines (the same numbers in the next 2 lines) 
	parameter (LSTR=128)   
	character*128 STR,TAKESTR,AUX,FPOT,AUX2
	character*32 KEYW
	integer ISPT(MAXPOT),JSPT(MAXPOT)
	IE=0                 ! error index
	IO=4                 ! input channel (fortran standard input)
	open(unit=4,file='md.input',status='old')
*  1.1 Defaults
*  -----------
	IPRINT = 5
	IFNAME=0
	PATHDB='.'
	LRST=.true.
	LDMP=.true.
	LFORMI=.false.
	LFORMO=.false.
	LINPUT=.false.
	L0AVS=.false.
	L0CPU=.false.
	L0VS=.false.
	LNVT=.false.
	LSCLT=.false.
	LSCTD=.false.
	LNPT=.false.
	LSEP=.false.
	NTYPES=0
	LHEX=.false.
	LOCT=.false.
	ICELL=0
	RHO=0.8
	IDST=0
	IDTS=0
	ISHK=0
	IBOXS=0
	BOXL=0.
	BOYL=0.
	BOZL=0.
	DT=2.
	NSTEPS=100000
	NFREQ=1
	IAVER=1
	NAVER=10000
	LSHEJK=.false.
	RCUT=12.
	SHORT=5.
	ICHNB=10
	ALPHA=2.8
	FEXP=9.
	TRTEMP=298.
	TRPRES=1.
	do I=1,NTPS
	  ISHEJK(I)=1
	  LMOVE(I)=.true.
	  IINIT(I)=0
	  LIST(I)=1
	  L14NB(I)=.true.
	  L15NB(I)=.true.
	  C14EL(I)=1.
	  C14LJ(I)=1.
	  IPOT(I)=0
	end do
	ICHECK=0
	LRR=.false.
	LCHT=.false.
	LCHP=.false.
	LSCFT=.false.
	LXMOL=.false.
	ISTAV=0
	IEXT=0
	LBONDL=.false.
	LCOMD=.false.
	LKONG = .false.
	LGEOM=.false.
	LGATHER=.false.
	LVISUAL=.false.
	LGR=.false.
	LGRST=.false.
	LGDMP=.false.
	ICMOM=0
	LCMOM=.true.          
*  1.2 Proceedind input file
*  ------------------------
*  1.2.1 Start parsing input file
C  this subroutine reads a line from the input file 
C   ignorint line starting with "#"
 10	STR=TAKESTR(IO,IE) 
	read(STR,*,end=10)KEYW
*  1.2.2 Basic simulation parameters
	if(KEYW(1:12).eq.'Main_filenam')then
	  read(STR,*,err=89)KEYW,AUX
	  do I=1,128
	    if(AUX(I:I).ne.' ')go to 1
	  end do
 1	  LD=LENG(AUX,128)
	  FNAME=AUX(I:LD)    
	  IFNAME=1
	  LDF=LD-I+1
	else if(KEYW(1:13).eq.'Verbose_level')then
	  read(STR,*,err=89,end=89)KEYW,IPRINT 
	else if(KEYW(1:7).eq.'Path_DB')then
          do J=9,112
	    if(STR(J:J).ne.' ')then
	      read(STR(J:128),'(a)',err=89)PATHDB
	      go to 7
	    end if
	  end do
 7	  LPD=LENG(PATHDB,112)
	else if(KEYW(1:6).eq.'Visual')then
          read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')LVISUAL=.true.
	else if(KEYW(1:9).eq.'Read_rest')then  
	  read(STR,*,err=8,end=8)KEYW,AUX,AUX2
	  if(AUX2(1:5).eq.'ASCII')LFORMI=.true.
 8	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'n')LRST=.false.
	else if(KEYW(1:9).eq.'Dump_rest')then  
	  read(STR,*,err=9,end=9)KEYW,NDUMP,AUX2
	  if(AUX2(1:5).eq.'ASCII')LFORMO=.true.
 9	  read(STR,*,err=89,end=89)KEYW,NDUMP
	  if(NDUMP.le.0.and.TASKID.eq.MAST)then
	    write(*,*)'!Warning: Number of steps for restart is ',NDUMP
	    write(*,*)' Restart file will not be written!'
	  end if
	else if(KEYW(1:10).eq.'Check_only')then  
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')LINPUT=.true.
	else if(KEYW(1:7).eq.'Zero_av')then  
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')L0AVS=.true.
	else if(KEYW(1:8).eq.'Zero_CPU')then  
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')L0CPU=.true.
	else if(KEYW(1:8).eq.'Zero_vel')then  
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')L0VS=.true.
	else if(KEYW(1:13).eq.'Molecule_type')then
	  read(STR,*,err=89,end=89)KEYW,NTYPES
	  do ITYP=1,NTYPES
 	    STR=TAKESTR(IO,IE) 
	    read(STR,*,end=15,err=15)NAMOL(ITYP),NSPEC(ITYP),AUX
	    if(AUX(1:5).eq.'fixed')LMOVE(ITYP)=.false.
 15	    read(STR,*,err=89)NAMOL(ITYP),NSPEC(ITYP)
	  end do
	else if(KEYW(1:3).eq.'PBC')then
	  read(STR,*,err=89)KEYW,AUX
	  if(AUX(1:3).eq.'hex')then 
            LHEX=.true.
	    ICELL=2
	  else if(AUX(1:3).eq.'oct')then
	    LOCT=.true.
	    ICELL=1
	  else if(AUX(1:3).eq.'rec')then
	    LHEX=.false.
	    LOCT=.false.
	    ICELL=0
	  else
	    if(TASKID.eq.MAST)then
	write(*,*)' !!!Warning: Unknown periodic boundary conditions.',
     +' Rectangular PBC implied'
	    end if
	  end if
	else if(KEYW(1:3).eq.'Box')then
	  read(STR,*)KEYW,BOXL,BOYL,BOZL
	  IBOXS=1
	else if(KEYW(1:7).eq.'Density')then
	  read(STR,*,err=89,end=89)KEYW,RHO
	  IDST=1
	else if(KEYW(1:9).eq.'Nose_ther')then
	  read(STR,*,err=89,end=89)KEYW,TRTEMP,QT
	  LNVT=.true.
	else if(KEYW(1:12).eq.'Velocity_sca')then
	  read(STR,*,err=89,end=89)KEYW,TRTEMP,TDELT
	  LSCLT=.true.
	else if(KEYW(1:14).eq.'Separate_therm')then  
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')LSCTD=.true.
	else if(KEYW(1:9).eq.'COM_check')then  
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')then
	    LCMOM=.true.
	    read(STR,*,err=89,end=89)KEYW,AUX,ICMOM
	  else
	    LCMOM=.false.
	  end if
	else if(KEYW(1:12).eq.'Barostate_NH')then
	  read(STR,*,err=89,end=89)KEYW,TRPRES,QPR
	  LNPT=.true.
	else if(KEYW(1:13).eq.'Barostate_ani')then  
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')LSEP=.true.
	else if(KEYW(1:9).eq.'Time_step')then
	  read(STR,*,err=89,end=89)KEYW,DT
	else if(KEYW(1:11).eq.'Number_step')then
	  read(STR,*,err=89,end=89)KEYW,NSTEPS
	else if(KEYW(1:10).eq.'Double_tim')then
	  read(STR,*,err=89,end=89)KEYW,NFREQ
	  LSHEJK=.false.
	  IDTS=1
	else if(KEYW(1:6).eq.'Constr')then
	  read(STR,*,err=89,end=89)KEYW,TOLER
	  LSHEJK=.true.
	  ISHK=1
	  read(STR,*,err=18,end=18)KEYW,TOLER,(ISHEJK(I),I=1,NTYPES)
 18	  continue
	else if(KEYW(1:6).eq.'Output')then
	  read(STR,*,err=89,end=89)KEYW,IAVER
	else if(KEYW(1:8).eq.'Serie_av')then
	  read(STR,*,err=89,end=89)KEYW,NAVER
	else if(KEYW(1:7).eq.'R_cutof')then
	  read(STR,*,err=89,end=89)KEYW,RCUT
	else if(KEYW(1:7).eq.'R_short')then
	  read(STR,*,err=89,end=89)KEYW,SHORT
	else if(KEYW(1:14).eq.'Neighbour_list')then
	  read(STR,*,err=89,end=89)KEYW,ICHNB
	else if(KEYW(1:9).eq.'Electrost')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:3).eq.'Cut')then
	    ALPHA=0.d0
	    FEXP=0.
	  else if(AUX(1:2).eq.'Ew')then
	    read(STR,*,err=89,end=89)KEYW,AUX,ALPHA,FEXP
	  else if(AUX(1:2).eq.'RF')then
	    read(STR,*,err=89,end=89)KEYW,AUX,ALPH,FEXP
	    ALPHA=-abs(ALPH)
	  else 
	    write(*,*)'!!!Warning: unknown electrostatics :',AUX
	  end if
	else if(KEYW(1:7).eq.'Startup')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:4).eq.'xmol')then
	    ICHECK=0
	  else if(AUX(1:7).eq.'Mol_COM')then
	    ICHECK=-1
	  else if(AUX(1:3).eq.'xyz')then
	    ICHECK=-2
	  else if(AUX(1:3).eq.'FCC')then
	    ICHECK=1
	  else if(AUX(1:8).eq.'Cyl_hole')then
	    ICHECK=2
	    read(STR,*,err=89,end=89)KEYW,AUX,RDNA
	  else if(AUX(1:8).eq.'Sph_hole')then
	    ICHECK=3
	    read(STR,*,err=89,end=89)KEYW,AUX,RDNA
	  else if(AUX(1:5).eq.'Cubic')then
	    ICHECK=4
	  else
	    write(*,*)'!!! Warning: Sturtup parameter unknown ',
     +  '   XMOL input file assumed'
	    ICHECK=0
	  end if
	else if(KEYW(1:9).eq.'Start_rot')then
	  read(STR,*,err=89,end=89)KEYW,(IINIT(I),I=1,NTYPES)
	else if(KEYW(1:9).eq.'Bind_atom')then
	  read(STR,*,err=89,end=89)KEYW,FILREF,RLR
	  LRR=.true.
	else if(KEYW(1:8).eq.'Change_T')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')LCHT=.true.
	else if(KEYW(1:8).eq.'Change_V')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')LCHP=.true.
	else if(KEYW(1:9).eq.'Cut_force')then
	  read(STR,*,err=89,end=89)KEYW,FTSCL
	  LSCFT=.true.
	else if(KEYW(1:12).eq.'Average_from')then
	  read(STR,*,err=89,end=89)KEYW,IHIST
	else if(KEYW(1:11).eq.'Average_int')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')ISTAV=1
	else if(KEYW(1:9).eq.'Dump_XMOL')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:3).eq.'com')LCOMD=.true.
	  if(AUX(1:1).ne.'n')LXMOL=.true.
	else if(KEYW(1:8).eq.'Trajecto')then
	  read(STR,*,err=89,end=89)KEYW,AUX,TTREK,LTREK,AUX2
	  if(AUX(1:6).eq.'bincrd')then
	    ITREK=1
	  else if(AUX(1:6).eq.'binvel')then
	    ITREK=3
	  else if(AUX(1:6).eq.'asccrd')then
	    ITREK=2
	  else if(AUX(1:6).eq.'ascvel')then
	    ITREK=4
	  else
	    write(*,*)'!!! Warning: unknown type of trajectory ',
     +' Binary coordinate-only is implied'
	  end if
	  if(AUX2(1:3).eq.'all')then
	    do I=1,NTYPES
	      LIST(I)=1
	    end do
	  else
	    read(STR,*,err=89,end=89)
     +      KEYW,AUX,TTREK,LTREK,(LIST(I),I=1,NTYPES)
	  end if
	else if(KEYW(1:8).eq.'El_field')then
	  read(STR,*,err=89,end=89)KEYW,EXAMPL,EXFREQ
	  IEXT=1
	else if(KEYW(1:9).eq.'Bond_list')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:1).eq.'y')LBONDL=.true.
	else if(KEYW(1:11).eq.'Combination')then
	  read(STR,*,err=89,end=89)KEYW,AUX
	  if(AUX(1:4).eq.'Kong')then
            LKONG=.true.
	  else if(AUX(1:5).eq.'Sigma')then
            LGEOM=.true.
	  else if(AUX(1:4).eq.'Lore')then
	    continue
	  else
	    write(*,*)'Unclear combination rule. Use Lorentz-Berthelot'
	  end if
	else if(KEYW(1:6).eq.'Gather')then
	  LGATHER=.true.
	else if(KEYW(1:8).eq.'RDF_calc')then
	  read(STR,*,err=89,end=89)KEYW,AUX,RDFCUT,MAX
	  LGR=.true.
	  if(AUX(1:2).eq.'RW')then
	    LGRST=.true.
	    LGDMP=.true.
	  else if(AUX(1:2).eq.'RO')then
	    LGRST=.true.
	    LGDMP=.false.
	  else if(AUX(1:2).eq.'WO')then
	    LGRST=.false.
	    LGDMP=.true.
	  else if(AUX(1:3).eq.'NOR')then
	    LGRST=.false.
	    LGDMP=.false.
	  else  
	    write(*,*)'!!! Warning: unknown option in RDFcalc',
     + ' RW is implied'
	    LGRST=.true.
	    LGDMP=.true.
	  end if
*  Reading of RDF information
          do IS=1,NSITES
	    NGRI(IS)=0
	    do J=1,MXGRS
	      IGRI(IS,J)=0
	      MGRI(IS,J)=0
	    end do
          end do
          N          = 0    
	  STR=TAKESTR(IO,IE) 
          READ(STR,*,err=89,end=89) MAXRDF
	  if(MAXRDF.le.0)then
	    write(*,*)'!!! Warning: No RDFs in the list!'
	  else
            DO K       = 1,MAXRDF
              N          = N+1  
	      IE=-1
	      STR=TAKESTR(IO,IE)
	      if(IE.lt.0)go to 19
	      if(str(1:1).eq.'&')then
                READ(str(2:24),*)IREP
	        do IR=1,IREP
	          STR=TAKESTR(IO,IE)
	          read(STR,*,err=89,end=89)I,J  
	          NGRI(I)=NGRI(I)+1 
	          if(J.ne.I)NGRI(J)=NGRI(J)+1
	          if(NGRI(I).gt.MXGRS)then
 	          write(*,*)' List of RDF exceeded for site ',I,' ',NM(I)
	          write(*,*)' increase MXGRS in prcm.h'
		  write(*,*)' RDF-s input interrupted'          
	          go to 19
	          end if 
	          if(NGRI(J).gt.MXGRS)then
 	          write(*,*)' List of RDF exceeded for site ',J,' ',NM(J)
	          write(*,*)' increase MXGRS in prcm.h'
		  write(*,*)' RDF-s input interrupted'          
	          go to 19
	          end if 
*  NGRI - number of RDF-pairs for this site
*  IGRI - 2-nd site for this RDF
*  MGRI - RDF number
	          IGRI(I,NGRI(I))=J
	          MGRI(I,NGRI(I))=N
	          IGRI(J,NGRI(J))=I
	          MGRI(J,NGRI(J))=N
	        end do
	      else      ! single pair 
	        IPER=1
	        read(str,*,err=89,end=89)I,J
	        NGRI(I)=NGRI(I)+1 
	        if(J.ne.I)NGRI(J)=NGRI(J)+1
	        if(NGRI(I).gt.MXGRS)then
 	        write(*,*)' List of RDF exceeded for site ',I,' ',NM(I)
	        write(*,*)' increase MXGRS in prcm.h'
		write(*,*)' RDF-s input interrupted'          
	        go to 19
	        end if 
	        if(NGRI(J).gt.MXGRS)then
 	        write(*,*)' List of RDF exceeded for site ',J,' ',NM(J)
	        write(*,*)' increase MXGRS in prcm.h'
	        write(*,*)' RDF-s input interrupted'          
	        go to 19
	        end if 
	        IGRI(I,NGRI(I))=J
	        MGRI(I,NGRI(I))=N
	        IGRI(J,NGRI(J))=I
	        MGRI(J,NGRI(J))=N
	      end if
            END DO
	  end if     ! if(MAXRDF      
 19	  MAXRDF     = N
*  TCF input
	else if(KEYW(1:8).eq.'TCF_calc')then
	  read(STR,*,err=89,end=89)KEYW,AUX,NSTEG,JUMP
	  LCF=.true.
	  if(AUX(1:2).eq.'RW')then
	    LCFRST=.true.
	    LCFDMP=.true.
	  else if(AUX(1:2).eq.'RO')then
	    LCFRST=.true.
	    LCFDMP=.false.
	  else if(AUX(1:2).eq.'WO')then
	    LCFRST=.false.
	    LCFDMP=.true.
	  else if(AUX(1:3).eq.'NOR')then
	    LCFRST=.false.
	    LCFDMP=.false.
	  else  
	    write(*,*)'!!! Warning: unknown option in TCFcalc',
     + ' RW is implied'
	    LCFRST=.true.
	    LCFDMP=.true.
	  end if
	  STR=TAKESTR(IO,IE)
	  read(STR,*,err=89,end=89)(ITCF(I),I=1,NTCF)
	  if(ITCF(5)*ITCF(6).ne.0)then
	    STR=TAKESTR(IO,IE)
            read(STR,*,err=89,end=89)(N1(I),I=1,NTYPES)
	    STR=TAKESTR(IO,IE)
            read(STR,*,err=89,end=89)(N2(I),I=1,NTYPES)
	  end if
	else if(KEYW(1:3).eq.'End'.or.KEYW(1:3).eq.'end')then
	  if(IPRINT.ge.6)write(*,*)' Ending reading input file'
	  go to 25
	else
	  write(*,*)'!!! Warning: Wrong keyword: ',KEYW
	end if
	if(IPRINT.ge.8)write(*,'(a)')STR
	go to 10
*
*   1.3. Control output
*   -------------------
 25	continue
	if(IPRINT.ge.7.and.TASKID.eq.MAST)then
	   write(*,*)' The following input parameters are read input: '
	   write(*,*)' Basis file name ',FNAME
	   write(*,*)' Molecular database path ',PATHDB(1:LPD)
	   if(LVISUAL)write(*,*)' Interface to visual shell is on'
	   if(LRST)then
	     write(*,*)' Read from restart file ',FNAME(1:LDF)//'.dmp'
	   else
	     write(*,*)' Do not read from restart file'
	   end if
 	   if(LDMP)then
	     write(*,*)' Dump restart file each ',NDUMP,' steps'
	   else
	     write(*,*)' Do not dump restart file'
	   end if
	   if(LINPUT)then
	     write(*,*)' Only check input / results. Do not run MD'
	   else
	     write(*,*)' Run MD simulations'
	   end if
	   if(L0AVS)then
	     write(*,*)' Begin new averaging'
	   else
	     write(*,*)' Continue collecting averages'
	   end if
	   if(L0CPU)then
	     write(*,*)' Begin new counting of CPU time'
	   else
	     write(*,*)' Continue counting of CPU time'
	   end if
	   if(L0VS)then
	     write(*,*)' Set all velocities to zero'
	   end if
	   write(*,*)' Molecule types:',NTYPES
	   do I=1,NTYPES
	     if(LMOVE(I))then 
	       write(*,*)NAMOL(I),NSPEC(I)
	     else
	       write(*,*)NAMOL(I),NSPEC(I),' FIXED'
	     end if
	   end do
	   if(LHEX)write(*,*)' Hexagonal periodioc boundary conditions'
	   if(LOCT)write(*,*)' Truncated octahedron PBC'
	   if(LKONG)write(*,*)' Kong combination rules'
	   if(LGEOM)write(*,*)' Geometrical sigma combination rules'
	   if(IDST.ne.0)write(*,*)' Density ',RHO
	   if(IBOXS.ne.0)write(*,*)' Box sizes:',BOXL,BOYL,BOZL
           if(IEXT.eq.1)write(*,*)
     +'External filed ',EXAMPL,' Frequency (hz) ',EXFREQ
	   if(LNVT)then
	     write(*,*)' Nose thermostate:   T=',TRTEMP,' Qt=',QT
	   end if
	   if(LSCLT)then
	     write(*,*)' Velocity scaling:   T=',TRTEMP,' DeltaT=',TDELT
	   end if
	   if(LSCTD)then
	     write(*,*)
     +' Temperature control for each molecular species separately'
	   else
	     write(*,*)' Temperature control for the whole system'
	   end if
	   if(LCMOM)then
	     write(*,*)' Excess COM momenta removed each ',ICMOM,' steps'
	   else
	     write(*,*)' Excess COM momenta are kept'
	   end if
	   if(LNPT)then
	     write(*,*)' Nose barostate:   P=',TRPRES,' Qpr=',QPR
	   end if
	   if(LSEP)then
	     write(*,*)
     +' Pressure control for each direction separately'
	   else
	     write(*,*)' Isotropic pressure controll'
	   end if
	   write(*,*)' Time step:',DT,' fs'
	   write(*,*)' Number of steps to run:',NSTEPS
	   if(LSHEJK)then
	     write(*,*)' Constrain bonds with tolerance ',
     +TOLER,' for molecules:'
	     do I=1,NTYPES
		if(ISHEJK(I).ne.0)write(*,*)I,NAMOL(I)
	     end do
	   end if
	   write(*,*)' Output step and energies each ',IAVER,' steps'
	   write(*,*)' Intermediate averaging in series of ',
     * NAVER,' steps'
	   write(*,*)' Rcutoff ',RCUT
	   write(*,*)' Rshort ',SHORT
	   if(ALPHA.eq.0.d0)then
	      write(*,*)' No treatment of long-range electrostatics'
	   else if(ALPHA.lt.0.)then
	      write(*,*)' Reaction field method with ',-ALPHA,FEXP
	   else
	      write(*,*)' Ewald method with ',ALPHA,FEXP
	   end if
	   if(ICHECK.eq.-2)then
      write(*,*)'Start atom coordinates from .inp file in xyz format'
	     else if(ICHECK.eq.-1)then
      write(*,*)'Start molecular COM coordinates from .inp file'
	     else if(ICHECK.eq.0)then
      write(*,*)'Start atom coordinates from .inp file in xmol format'
	     else if(ICHECK.eq.1)then
      write(*,*)'Start atom coordinates from FCC lattice'
	     else if(ICHECK.eq.2)then
      write(*,*)'Distribute molecules around cylindrical hole'
      write(*,*)'Radius ',RDNA	
	     else if(ICHECK.eq.3)then
      write(*,*)'Distribute molecules around spherical hole'
      write(*,*)'Radius ',RDNA	
	     else if(ICHECK.eq.4)then
      write(*,*)'Start atom coordinates from cubic lattice'
	   end if
	   write(*,*)' Initial orientation parameter',
     +    (IINIT(I),I=1,NTYPES)
	   if(LRR)then
	     write(*,*)' Standard structure taken from file ',filref
	     write(*,*)' Deviation from "standard" structure',RLR
	   end if
	   if(LCHT)write(*,*)' Change temperature after restart'
	   if(LCHP)write(*,*)' Change volume after restart'
	   if(LSCFT)write(*,*)' Cut large forces at level ',FTSCL
	 if(LGATHER)write(*,*)' Gather each molecule at one place'
	 write(*,*)' final averaging from ',IHIST,' serie '
	 if(LXMOL)then
	   if(LCOMD)then
	     write(*,*)' Write down molecular COM in final configuration'
	   else
	     write(*,*)' Write down final configuration'
	   end if
	 end if
	 if(ITREK.gt.0)write(*,*)
     +' Dump configurations each ',TTREK,' fs ',LTREK,' conf. per file'
	 if(ITREK.gt.0)write(*,*)
     +' Include molecules: ',(LIST(I),I=1,NTYPES)
	 if(LFORMI)write(*,*)' ASCII format of input restart file'
	 if(LFORMO)write(*,*)' ASCII format of output restart file'
	 if(ISTAV.eq.1)write(*,*)' Calculate average bond lenghts,',
     + ' angles and their energies'
	if(LGR)then
	   write(*,*)' calculate RDF'
	   if(LGRST)write(*,*)' restart RDF file'
	   if(LGDMP)write(*,*)' dump RDF restart file'
	   write(*,*)' RDF cut off ',RDFCUT 
	   write(*,*)' RDF resolution ',RDFCUT/MAX
	end if
	if(LCF)then
          write(*,*)' calculate TCF'
 	  if(LCFRST)write(*,*)' restart TCF file'
	  if(LCFDMP)write(*,*)' dump TCF file'
	  if(LCF)write(*,*)' TCF for ',NSTEG,
     + ' steps with jump = ',JUMP
	  write(*,*)' which tcf calculate:',
     +(ITCF(I),I=1,NTCF)
	  if(ITCF(5)*ITCF(6).ne.0)then
            write(*,*)' unit vector for reorientation tcf:'
	    write(*,*) (N1(I),I=1,NTYPES)
	    write(*,*) (N2(I),I=1,NTYPES)
	  end if
	end if
      end if   !  if(IPRINT
*
*  1.4. Some controls
*  ------------
      IF(NTYPES.GT.NTPS) STOP '!!! INCREASE NTPS !'
      if(NTYPES.eq.0)then
	 write(*,*)'!!! No molecular types have been given!'
	 stop
      end if
      if(IFNAME.eq.0)then
	write(*,*)'!!! No Main_filename has been given !'
	stop
      end if
      if(IBOXS.eq.1.and.IDST.eq.1)then
	write(*,*)'!!! Error: both box sizes and density are specified'
	call FINAL
      end if
      if(ISHK.eq.1.and.IDTS.eq.1)then
	write(*,*)'!!! Error: both double time step and SHAKE ',
     +' algorithms are specified'
	call FINAL
      end if
      if(LHEX)BOYL=cos30*BOXL
      if(IPRINT.ge.7)then
	  if(BOXL*BOYL*BOZL.le.0)then
	    write(*,*)' density ',RHO
	  else
	    write(*,*)' box sizes: ',BOXL,BOYL,BOZL
	  end if
      end if
      return
*   1.5  Diagnos               
 89   IE=99
      STR=TAKESTR(IO,IE)
      write(*,*)' Are you running Windows ?  '
      return
      end
*
C=====================================================================
*
*    2. Read molecular structure and force field
*    -------------------------------------------
*
C    This subroutine read the complete information about the molecule 
C    of type NTYP from the .mol file located in molecular database. 
C    It containes initial configuration of the molecule, atom types, 
C    names, charges, LJ parameters, list of covalent bonds, angles, 
C    dihedral angles with corresponding force field parameters.
C    Only in a few cases these data may be overwritten by the data
C    from the input file (see subroutine INPUT)
C
C    Information how to write .mol files are given in README file
C
C================ READMOL =============================================
C
      subroutine READMOL(NTYP)
      include "prcm.h"
      DIMENSION CM(3)
      character namn*64,STR*128,FIL*128,TAKESTR*128,NAMH*4,NMF*5
      data NX/0/
*
*  2.1 Open the .mol file
*  ----------------------
      FIL=' '
C  molecule name                    
      NAMN=NAMOL(NTYP)
      IE=0
      NX0=NX
      LD=LENG(PATHDB,112)
C  path
      FIL(1:LD)=PATHDB(1:LD)
      FIL(LD+1:LD+1)='/'
      LN=LENG(NAMN,64)
C  add file name
      FIL(LD+2:LD+LN+1)=NAMN(1:LN)
C  add extension .mol
      FIL(LD+LN+2:LD+LN+6)='.mmol'
      open(unit=12,file=fil,status='old',err=999)
      if(IPRINT.ge.5.and.TASKID.eq.MAST)
     +write(*,*)'*** MOLECULE ',NAMN,'   Type  No.',NTYP
*
*  2.2 Read information about the atoms
*  ------------------------------------
        IE=0
	STR=TAKESTR(12,IE)
C  Number of atoms
	read(STR,*,err=99,end=99)NSTS
	if(IPRINT.ge.6)PRINT * ,'*** NR OF SITES:                ',NSTS
	NSITS(NTYP)=NSTS
	NSBEG=NX+1 
C  check unit vector from the input file             
	if(LCF)then 
	  if(N1(NTYP).gt.NSTS)then  
	    N1(NTYP)=1
	    write(*,*)
     +    '!!! Too large site number for unit vector. Set to 1'
	  end if  
	  if(N2(NTYP).gt.NSTS)then  
	    N2(NTYP)=NSTS
	    write(*,*)
     +    '!!! Too large site number for unit vector. Set to ',NSTS
	  end if  
      end if
C   Cycle over atoms
      do I=1,NSTS
	  STR=TAKESTR(12,IE)
	  IX=NX+I
	  NM(IX)=STR(1:4)
	  read(STR(5:80),*,err=99,end=99)(R(IX,J),J=1,3),
     +    MASS(IX),CHARGE(IX),SIGMA(IX),EPSIL(IX)
          if(IPRINT.ge.6)
     +      write(*,'(a5,a4,a4,3f7.3,a3,f8.4,a3,f6.3,a5,f8.4,a4,f9.4)')
     +      'atom ',NM(IX),'  R=',(R(IX,J),J=1,3),
     +      ' M=',MASS(IX),' Q=',CHARGE(IX),' si=',
     +      SIGMA(IX),' eps=',EPSIL(IX)
      end do
      NSEND=NX+NSTS
*
*  2.3 Make local coordinate system COM at zero
*  -------------------------------------------- 
      SUMM       = 0.D0
      DO IS      = NSBEG,NSEND
        SUMM       = SUMM+MASS(IS)
      END DO! OF IS
C  This is total mass of the molecule in atomic units
      SUMMAS(NTYP) = SUMM 
      CM(1)      = 0.D0
      CM(2)      = 0.D0
      CM(3)      = 0.D0
      DO IS      = NSBEG,NSEND
        CM(1)      = CM(1)+R(IS,1)*MASS(IS)
        CM(2)      = CM(2)+R(IS,2)*MASS(IS)
        CM(3)      = CM(3)+R(IS,3)*MASS(IS)
      END DO! OF IS
      CM(1)      = CM(1)/SUMM
      CM(2)      = CM(2)/SUMM
      CM(3)      = CM(3)/SUMM
C  Set local atom coordinate so that molecular center of mass was zero.
C  Exceptions: molecules with IINIT=1 (see p. 1.18 this file) which are
C  put in the cell as they are
      if(IINIT(NTYP).ne.1)then
        if(IPRINT.ge.6)write(*,*)' Setting local COM to 0.  Type  ',
     &    NAMOL(NTYP)
        DO K       = 1,3
          DO IS      = NSBEG,NSEND
            R(IS,K)    = R(IS,K)-CM(K)
          END DO! OF IS
        END DO! OF I
      end if
C  Print extra information
      if(IPRINT.ge.8)then
        write(*,*)'*** Molecular GEOMETRY (C.O.M. IN ORIGIN): '
        DO IS      = NSBEG,NSEND
          write(*,'(I4,3x,A4,3X,3(1X,F7.3))')
     +    IS,NM(IS),R(IS,1),R(IS,2),R(IS,3)
        END DO! OF IS
      end if
      if(IPRINT.ge.10)then
        write(*,*)
        write(*,*)'*** INTER ATOMIC DISTANCES '
        do IS=NSBEG,NSEND
	  do JS=IS+1,NSEND
	    RR2=0.
	    do K=1,3
	      RR2=RR2+(R(IS,K)-R(JS,K))**2
	    end do
 	    RR=sqrt(RR2)
            write(*,'(3a4,2x,f9.4)')NM(IS),' -> ',NM(JS),RR
	  end do
	end do
      end if
*
*  2.4 Output of the reference
*  ---------------------------
      STR=TAKESTR(12,IE)
      read(STR,*,err=99,end=99)NSR
      do I=1,NSR
	STR=TAKESTR(12,IE)
	if(IPRINT.ge.6)write(*,'(a80)')STR
      end do
*
*  2.5 Input bonds information
*  --------------------------- 
      STR=TAKESTR(12,IE)
C  Number of bonds
      read(STR,*,err=99,end=99)NBD
      NRB(NTYP)=NBD
      do K=NRBB+1,NRBB+NBD
	STR=TAKESTR(12,IE)
	read(STR,*,err=36,end=36)ID(K),IB(K),JB(K),RB(K),FB(K),
     + DB(K),RO(K)
	go to 37
 36	read(STR,*,err=36,end=36)ID(K),IB(K),JB(K),RB(K),FB(K)
	DB(K)=0.
	RO(K)=0.
	ID(K)=0
 37	continue
C  If bond length is zero in the input file, it is set from the
C  atom coordinates
	if(RB(K).le.0.d0)RB(K)=sqrt((R(IB(K),1)-R(JB(K),1))**2+
     +  (R(IB(K),2)-R(JB(K),2))**2+(R(IB(K),3)-R(JB(K),3))**2)
      end do
      if(IPRINT.ge.8.and.NBD.gt.0)then
        write(*,*) 
        PRINT "('*** PARAMETERS FOR ',I3,'  BONDS')",NBD
        PRINT *,'*** HARMONIC POTENTIALS:'
        IDTOT      = 0
        DO M       = NRBB+1,NRBB+NBD
          if(ID(M).ne.0)IDTOT=IDTOT+1
           write(*,'(a5,a4,a1,a4,2I6,a,f8.4,a,f8.2)')
     +    'Bond ',NM(IB(M)+NX),'-',NM(JB(M)+NX),IB(M),JB(M),' Req.=',
     +    RB(M),'   Force=',FB(M)
        END DO! OF M
        IF(IDTOT.GT.0) THEN
          write(*,*)
          write(*,*)'*** MORSE POTENTIALS:'
          DO M       = NRBB+1,NRBB+NRB(NTYP)
            if(ID(M).ne.0)
     +      PRINT "('*** BOND: ',I3,'  - ',I3,'  RO: ',F7.3,'  DISS  ',
     X      'ENERGY  :',F9.3)",IB(M),JB(M),RO(M),DB(M)
          END DO! OF M
        END IF
      end if
*
*  2.6 Input angles information
*  ----------------------------
      STR=TAKESTR(12,IE)
      read(STR,*,err=99,end=99)NAN
      NRA(NTYP)=NAN   
      do K=NRAA+1,NRAA+NAN
	STR=TAKESTR(12,IE)
	read(STR,*,err=99,end=99)IA(K),JA(K),KA(K),RA(K),FA(K)
	if(RA(K).le.0)RA(K)=ANGLG(R,IA(K),JA(K),KA(K),NS)*TODGR
      end do
      if(IPRINT.ge.8.and.NAN.gt.0)then
	PRINT *
        PRINT *,'*** COVALENT ANGLES '
	do K=NRAA+1,NRAA+NAN
          PRINT "('***',3i6,3(1X,A4),'  ->',F9.3,'  Force c.',f9.3)",
     +       IA(K) ,   JA(K) ,   KA(K) ,
     +    NM(IA(K)+NX),NM(JA(K)+NX),NM(KA(K)+NX),RA(K),FA(K)        
	end do
      end if
*
*  2.7 Input torsional angles information
*  -------------------------------------- 
C  (potential of standard 1+cos(al+n*f) type )
      STR=TAKESTR(12,IE)
      read(STR,*,err=99,end=99)NTT       
      NRT(NTYP)=NTT
      if(NTT.gt.0)then
	do K=NRTT+1,NRTT+NTT
	  ITORS(K)=0
	  STR=TAKESTR(12,IE)
C  Input AMBER torsion angle parameters
	  read(STR,*,err=99,end=99)
     +	  IT(K),JT(K),KT(K),LT(K),RT(K),FT(K),NMUL(K)
C  Set parameters to other types of torsions to zero.
	  FT2(K)=0.d0
	  FT3(K)=0.d0
	  FT4(K)=0.d0
	  FT5(K)=0.d0
	end do
	if(IPRINT.ge.8)then
	  PRINT *
          PRINT *,'*** DIHEDRALS '
	  write(*,*)'  1+cos(nf) type:'
	  do K=NRTT+1,NRTT+NTT
            PRINT"(4I6,4(1X,A4),' ->',I2,F8.2,'  F.c.',f9.3)"
     +      ,   IT(K),JT(K),KT(K),LT(K)
     +         ,NM(IT(K)+NX),NM(JT(K)+NX),NM(KT(K)+NX),NM(LT(K)+NX),
     +        NMUL(K),RT(K),FT(K)        
	  end do
	end if
      end if                 
      NRBB         = NRBB+NRB(NTYP)
      NRAA         = NRAA+NRA(NTYP)
      NRTT         = NRTT+NRT(NTYP)
*
*  2.8 Optional force field terms
*  ------------------------------
C  Optional force field terms may be specified after the list of torsions.
C  The corresponding fields begin from a key word, followed by the number
C  of the external parameters of this type and their list.
C  The number of types of external parameters may be increased in 
C  future releases   
 100  IE=-1
      STR=TAKESTR(12,IE)  
      if(IE.lt.0)go to 20
      do I=1,72
	if(STR(I:I).ne.' ')go to 105
      end do
      go to 100
 105  NMF=STR(I:I+4)
*  2.8.1  MM3 type torsion angle  (tors1)
C         Sum(Kn*cos(n*f)) , n=1,2,3
	if(NMF.eq.'Tors1'.or.NMF.eq.'tors1'.or.NMF.eq.'TORS1')then	
        STR=TAKESTR(12,IE)
	  read(STR,*,err=99,end=99)NTT       
	  NRT(NTYP)=NRT(NTYP)+NTT
	  if(NTT.gt.0)then
	    do K=NRTT+1,NRTT+NTT  
	      ITORS(K)=1
	      STR=TAKESTR(12,IE)
	      read(STR,*,err=99,end=99)
     +	    IT(K),JT(K),KT(K),LT(K),FT1(K),FT2(K),FT3(K)
	      NMUL(K)=0   
	      FT(K)=0.d0 
	      RT(K)=0.d0 
	      RT(K)=0.d0 
	      FT4(K)=0.d0
	      FT5(K)=0.d0
	    end do
	    if(IPRINT.ge.8)then
	      PRINT *
            PRINT *,'*** DIHEDRALS '
	      write(*,*)'  Sum(1+Cos(nf)) type'
	      do K=NRTT+1,NRTT+NTT
              PRINT"(4I6,4(1X,A4),' -> ',3f9.3)",
     +        IT(K),JT(K),KT(K),LT(K),
     +        NM(IT(K)+NX),NM(JT(K)+NX),NM(KT(K)+NX),NM(LT(K)+NX),
     +        FT1(K),FT2(K),FT3(K)        
	        FT1(K)=0.5*FT1(K)
	        FT2(K)=0.5*FT2(K)
	        FT3(K)=0.5*FT3(K)
	      end do
	    end if
          NRTT         = NRTT+NTT
	  end if
*  2.8.2   Dihedral of Sum(cos(f)**n) type   (tors5)
	else if(NMF.eq.'Tors5'.or.NMF.eq.'tors5'.or.NMF.eq.'TORS5')then	
        STR=TAKESTR(12,IE)
	  read(STR,*,err=99,end=99)NTT       
	  NRT(NTYP)=NRT(NTYP)+NTT
	  if(NTT.gt.0)then
	    do K=NRTT+1,NRTT+NTT  
	      ITORS(K)=5
	      STR=TAKESTR(12,IE)
	      read(STR,*,err=99,end=99)
     +      IT(K),JT(K),KT(K),LT(K),FT1(K),FT2(K),FT3(K),FT4(K),FT5(K)
	      NMUL(K)=0   
	      FT(K)=-(FT1(K)+FT2(K)+FT3(K)+FT4(K)+FT5(K)) 
	      RT(K)=180.d0 
	    end do
	    if(IPRINT.ge.8)then
	      PRINT *
            PRINT *,'*** DIHEDRALS '
	      write(*,*)'  Sum(Cosf)**n type'
	      do K=NRTT+1,NRTT+NTT
              PRINT"(4I6,4(1X,A4),' ->',5f8.2)",
     +        IT(K),JT(K),KT(K),LT(K),
     +        NM(IT(K)+NX),NM(JT(K)+NX),NM(KT(K)+NX),NM(LT(K)+NX),
     +        FT1(K),FT2(K),FT3(K),FT4(K),FT5(K)        
	      end do
	    end if
          NRTT         = NRTT+NTT
	end if
      else if(NMF.eq.'Impro'.or.NMF.eq.'impro'.or.NMF.eq.'IMPRO')then	
*   2.8.3 Improper dihedrals  (impro)
        STR=TAKESTR(12,IE)
	read(STR,*,err=99,end=99)NTT       
	NRT(NTYP)=NRT(NTYP)+NTT
	if(NTT.gt.0)then
	  do K=NRTT+1,NRTT+NTT  
	    ITORS(K)=-1
	    STR=TAKESTR(12,IE)
	    read(STR,*,err=99,end=99)
     +      IT(K),JT(K),KT(K),LT(K),RT(K),FT(K)
	    NMUL(K)=0   
	  end do
	  if(IPRINT.ge.8)then
	    PRINT *
            PRINT *,'*** Improper DIHEDRALS '
	    do K=NRTT+1,NRTT+NTT
              PRINT"(4I6,4(1X,A4),' ->',5f8.2)",
     +        IT(K),JT(K),KT(K),LT(K),
     +        NM(IT(K)+NX),NM(JT(K)+NX),NM(KT(K)+NX),NM(LT(K)+NX),
     +        RT(K),FT(K)        
	    end do
	  end if
          NRTT         = NRTT+NTT
        end if
*  2.8.5.    1-4 interactions
*  2.8.5.1   1-4 exclusion
      else if(NMF.eq.'no_14'.or.NMF.eq.'NO_14')then 
	L14NB(NTYP)=.false.
*  2.8.5.2   All intramolecular LJ and electrostatic excluded
      else if(NMF.eq.'noi15'.or.NMF.eq.'NOI15')then 
	L15NB(NTYP)=.false.
*  2.6.5.3   1-4 electrostatic scaling factors
      else if(NMF.eq.'sel14'.or.NMF.eq.'SEL14')then
	read(STR,*,err=99)NMF,C14EL(NTYP)
	if(TASKID.eq.MAST.and.IPRINT.ge.5)write(*,*)
     &' Scale 14 electrostatics by ',C14EL(NTYP)
*  2.6.5.4   1-4 LJ scaling factors
      else if(NMF.eq.'slj14'.or.NMF.eq.'SLJ14')then
	read(STR,*,err=99)NMF,C14LJ(NTYP)
	if(TASKID.eq.MAST.and.IPRINT.ge.5)write(*,*)
     &' Scale 14 LJ interactions by ',C14LJ(NTYP)
      else if(NMF.eq.'Speci'.or.NMF.eq.'speci'.or.NMF.eq.'SPECI')then
*  2.8.5.5 Specific 1-4 interactions
C     These special LJ parameters are used only if specified
C     below atoms are "1-4" bound 
 	STR=TAKESTR(12,IE)
	read(STR,*,err=99,end=99)IADLJ 
	if(TASKID.eq.MAST.and.IPRINT.ge.6)write(*,*)
     &' Overwrite ',IADLJ,' LJ parameters'
	do I=NNAD+1,NNAD+IADLJ
	  STR=TAKESTR(12,IE)
	  read(STR,*,err=99,end=99)ILJ(I),SIGAD(I),EPAD(I) 
	  ILJ(I)=ILJ(I)+NX
	  if(IPRINT.ge.8)write(*,'(i4,2x,a4,2x,2f12.3)')
     +    ILJ(I),NM(ILJ(I)),SIGAD(I),EPAD(I)
	end do   
	NNAD=NNAD+IADLJ
	if(NNAD.gt.NNADM)stop ' increase NNADM'
*   2.8.6  Signals flexible SPC water
      else if(NMF(1:4).eq.'fSPC')then
	if(NSTS.eq.3)then
	  IPOT(NTYP)=2
	  if(IPRINT.ge.5)write(*,*)' Flexible SPC water'
	else
	  write(*,*)'!!! Error: flexible SPC water in .mmol file ',
     +   ' with ',NSTS,' molecules'
          write(*,*)'File ',FIL
	  call FINAL
	end if

*   2.8.7  Other cases
      else if(NMF(1:3).eq.'End'.or.NMF(1:3).eq.'end'.or.
     +  NMF(1:3).eq.'END')then
	  go to 20	
      else
	  LD=LENG(STR,73)
	  write(*,*)' !!! Unknown additional keyword ',STR(I:LD)
	  write(*,*)' !!! for Molecule ',NAMOL(NTYP)
      end if
      go to 100
*
*  2.9   Addresses for bonds, angles and dihedral angles
*  -----------------------------------------------------
 20   IADB(NTYP)   = NRBB-NRB(NTYP)+1
      IADA(NTYP)   = NRAA-NRA(NTYP)+1
      IADT(NTYP)   = NRTT-NRT(NTYP)+1
      NX=NX+NSTS
      close (12)
      return                      
 99	IE=99
	STR=TAKESTR(12,IE)
      write(*,*)' This message should never appear on the screen ...'
 999  write(*,*)'!!! Molecule ',NAMN,' not found in the data base'
      write(*,*)'    File not found: ',FIL
      stop
      end
*
*================================================================
*
*    3. Read and check a line from input file
*    ----------------------------------------
*
C    This subroutine reads a line from Fortran input file 
C    defined by "KAN" (i.e. opened with parameter unit=KAN)
C    It returns the line if it is not started with # (commentaries); 
C    otherwise it reads and check the next line.
C    Parameter IE:
C    If input value of IE is NOT 99, the subroutine performs normal action
C    (see above) and return IE=0 (normal exit) or IE=-1 (read error)
C    If input value of IE=99, the subroutine is used for diagnostics:
C    it reports the file name and the line number where input error occur.
C    
C
C============== TAKESTR =========================================
C
*   3.1 Definitions
      character*128 function TAKESTR(KAN,IE)
      character*128 AUX, fn*32 
      integer LCOUNT(256)   
      data LCOUNT /256*0/
      save AUX,LCOUNT
*   3.2 Normal action           
      if(IE.eq.99)go to 10
C  LCOUNT variable count line number (including commentaries) 
 1	LCOUNT(KAN)=LCOUNT(KAN)+1
 	read(KAN,'(a128)',err=10,end=20)AUX
	if(AUX(1:1).eq.'#')go to 1
	TAKESTR=AUX 
	IE=0
	return
*   3.3 Error diagnostic
 10	write(*,*)'!!! error in input file '   
	if(KAN.eq.5)then
	  write(*,*)' Standard input ' 
	else
	  inquire(unit=KAN,name=fn) 
	  write(*,*)' File ',fn
	end if
	write(*,*)' in line ',LCOUNT(KAN)
	write(*,*)AUX
	stop
*    3.4 End of file case (is not signaled if input IE < 0 )
 20	if(IE.ge.0)then
 	  write(*,*)'!!! end of file reached' 
	  if(KAN.eq.5)then
	    write(*,*)' Standard input ' 
	  else
	    inquire(unit=KAN,name=fn) 
	    write(*,*)' File ',fn
	  end if
	  stop
	end if  
	TAKESTR='  ' 
	IE=-1
	return
	end 
*
*   4. A few specific procedures
*   ----------------------------
*   4.1 Calculate covalent angles from configuration in .mol file
C    
C=============== ANGLG =============================================
C
      FUNCTION ANGLG(R,I,J,K,N)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION R(N,3)
      RIJ=0.
      RJK=0.
      SCL=0.
      DO M=1,3
	RIJ=RIJ+(R(I,M)-R(J,M))**2
        RJK=RJK+(R(J,M)-R(K,M))**2
        SCL =SCL +(R(I,M)-R(J,M))*(R(K,M)-R(J,M))
      end do
      ANGLG = DACOS(SCL/sqrt(RIJ*RJK))
      RETURN
      END
*
