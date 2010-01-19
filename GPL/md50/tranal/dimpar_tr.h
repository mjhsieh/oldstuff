* This file sets up array boundaries for TRANAL programm
      implicit real*8 (A-H,O-Z)
* max number of molecules  and atoms
      PARAMETER (NPART  = 5000, NTOT=20000)
* max number of  mol.types          bonds    
      PARAMETER (NTPS   =   8  , NBM=2000 )
*                sites            torsions 
      PARAMETER (NS     = 1000 , NT=2000 )
*       NMX,NMY - max sizes of 2D plots for 3-body cor.func.
*       NOMAX - max number of "to" sites for SDF calculations    
*       NDUPM - max number of "from" sites for SDF calculations    
      parameter (NMX=200,NMY=200,NOMAX=100,NDUPM=30,NOX2=2*NOMAX)
*   max sizes of 3D distributions (SDF)
      parameter (NOMXM=240,NOMYM=240,NOMZM=240)
*   MAXRDF - max number of different RDFs
*   MAXGR - max number of site pairs in one RDF
*   MAXCF - max number of time points in TCF calculations 
*   NAAM - grid size for RDF calculations
      parameter (MAXRDF=20,MAXGR=50,NAAM=400,MAXCF=500)
* number of pixels in 2D plots
      parameter (NTY=200,NTX=200)
* max number of "from" and "to" atoms in residence time calculations 
      parameter (NRT1M=260,NRT2M=1100)
* MAXMOL - max number of molecules in a single diffusion calculation
* MTRM - max number of time points for diffusion calculation
*     (Obs: arrays MAXMOL*MTRM are used!)
* NITTM - max number oftime points for residence time calculations
      parameter (MAXMOL=1000,MTRM=1000,NTCD=MTRM,NITTM=1000)
* NTOR  - number of different torsions distributions
* NTA - angular grid for torsions distributions
      parameter (NTOR=50,NTA=360)
* NEQHM - max number of bonds for order parameter calculations
* IODM - resolution for distribution of orientational order parameter
      parameter(NEQHM=100,IODM=1000)

