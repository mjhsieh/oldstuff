*   MDynaMix v 4.3
*
*   file dimpar.h
*
C     This file defines sizes of working arrays 
C     Edit parameters of this part to adjust to you SYSTEMS:
C     MOLECULAR and COMPUTER
C     Of course, the first one prefers as large as possible values and
C     the second one would like as low as possible values, so you should
C     compromise...       
C
C     Parameters which affect the memory most:
C     NPART, NTOT, MAXCF, NRQS, LHIST, NBLMX
C     (NBLMX scales as 1/NPROCS for parallel jobs)  
C
C     You must recompile the sources each time you change this file.
C
C
C  These parameters set maximum number of:              
C       NPART -  molecules
      PARAMETER (NPART  = 5000                             )
C       NTPS - different types of molecules
      PARAMETER (NTPS   =   6                              )
C       NS - sites (i.e., sum of atoms if one take one molecule of each type)
C       NTOT - total atoms
      PARAMETER (NS     = 1200  , NTOT   = 20000             ) 
C       NBSMX - max neighbours per atom inside cutoff of fast forces (SHORT)
C       (not essential if SHAKE algorithm used)
C       NBLMX - max neighbours per atom inside cutoff Rcut
C       These parameters may be decreased for parallel jobs as 1/NPROCS
      PARAMETER (NBSMX=64, NBLMX=800)
      PARAMETER (NBSMAX = NBSMX*NTOT,NBLMAX=NBLMX*NTOT          )  
C       NB     different covalent bonds 
C             (goes over ALL bonds on ONE molecule of EACH type)  
C       NBO    total covalent bonds per processor
      PARAMETER (NB     = 2000  , NBO     = 50000             )
C       NA,NAO  the same for covalent angles
      PARAMETER (NA    = 3000 , NAO     = 50000           )    
C       NT,NTO  the same for dihedral angles
      PARAMETER (NT     = 5000 , NTO   = 50000       )  
C       NHYD  - different types of hydrogen bonds
C       NNADM - number of specific Lennard-Jones parameters for 1-4
C               interactions 
C       NPLM  - max number of "linked" atoms 
C               (attached to some points by harmonic forces) 
      PARAMETER (NHYD=20, NNADM=1000, NPLM=100 )                 
C              RDF resolution    num of RDFs     Max RDF for 1 site 
      PARAMETER (MAXS   = 400  , MAXGR  = 30 ,   MXGRS=80)        
C                    points for tcf calculations  ; total tcf
C   set MAXCF=1 if you do not want calculate tcf and want to save memory
      PARAMETER (MAXCF  = 1    , NTCF   = 12               )
      PARAMETER (MAXTCF = (MAXCF+2)*NPART*3                ) 
C                   different averages
      PARAMETER (NRQS   = 12000 , NTPP =(NTPS*(NTPS+1))/2+1 ) 
C   max number of "bound" neighbours 
      PARAMETER (NBDMAX    = 40 )           
C       intermediate averaging points
      PARAMETER (LHIST  = 500)                        
C       maximum number of nodes   
      parameter (nodes=96)
C       maximum number of extra potentials
      parameter (MAXPOT=10)
*

