/* replicate

   Read LAMMPS 99 data file as input including class II support
   Creates new file as output (with new molecule tags as needed)

   New file can be translated, dilated, remapped to periodic box
     and/or replicated to create a larger system

   Syntax:  replicate [options] < infile > outfile

     Options:

     -2d
          LAMMPS data files are for 2-d problems, not 3-d
          this option must be specified before -replicate or 
	    -inbox or -outbox

     -repeat nx ny nz

          use the input data as a unit cell and replicate the unit
	    cell by (nx,ny,nz) in the +x, +y, +z directions
	  this creates new atoms, bonds, angles, etc
	  if not specified, nx = ny = nz = 1
	  do not specify nz if -2d option is set

     -inbox xlo xhi ylo yhi zlo zhi

          use these values for the bounding box of the input system
	    and for the unit cell used for replicating purposes
	  if not specified use the box bounds read in from infile
	  do not specify zlo and zhi if -2d option is set

     -outbox xlo xhi ylo yhi zlo zhi

          use these values for the bounding box of the output system
	  if not specified the output bounding box is computed from the input
	    bounding box by holding the "lo" values constant and creating 
	    new "hi" values by multiplying by the replicating factors
	  do not specify zlo and zhi if -2d option is set

      -remap

          remap all final atom positions into the output bounding box
	    in a periodic sense
	  works no matter how far away from the box the atoms are

      -truein

          true flags are specified in infile (after atom coords)


      -trueout

          true flags will be written to outfile (after atom coords)
	  if remapping is done, true flags will be modified accordingly

      Notes:
      
      Can do a system dilation/contraction by using "-repeat 1 1 1"
        and specifying a slightly larger/smaller output box than input
        box.

      Can add/delete true flags by using -truein and -trueout.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct box {
  double xlo,xhi,xsize;
  double ylo,yhi,ysize;
  double zlo,zhi,zsize;
};

main(int argc, char *argv[])

{
  char line[1000];            /* strings for reading/parsing input file */
  char linetoken[1000];
  char *token;             
                             /* numbers of various quantities as read in */

  int natoms,nbonds,nangles,ndihedrals,nimpropers;
  int ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes;

  int nmolecules;            /* total # of mols specified by mol tags */

                             /* vectors of various quantities as read in */

  int *molecule,*type,*btype,*atype,*dtype,*itype;
  int *truex,*truey,*truez;
  double *q,*x,*y,*z;
  int *bond1,*bond2,*angle1,*angle2,*angle3;
  int *dihed1,*dihed2,*dihed3,*dihed4,*impro1,*impro2,*impro3,*impro4;

  int three_d;               /* 0 for 2-d problem, 1 for 3-d */
  int nx,ny,nz;              /* replication factors */
  int inbox,outbox;          /* flags for whether in/out boxes were explicitly
				specified */
  int remap;                 /* 0 = no remapping, 1 = remap to output box */
  int truein;                /* 0/1 = no/yes input of true flags */
  int trueout;               /* 0/1 = no/yes output of true flags */

  struct box in,out;         /* input and output bounding boxes */
  int ntotal;                /* total replication factor = nx*ny*nz */
  double newx,newy,newz;     /* replicated atom position */
  int tx,ty,tz;              /* temporary true flags */

  int i,j,k,m,n,tmp,del,offset;   /* temp variables */

  char *gettoken(char *);    /* function defs */
  void copyline(int); 
  void scale(struct box, int, int, int, 
	     struct box, double *, double *, double *, int);
  void pbc(struct box, double *, double *, double *, int *, int *, int *, int);

/* default input values */

  three_d = 1;
  nx = ny = nz = 1;
  inbox = outbox = 0;
  remap = 0;
  truein = 0;
  trueout = 0;
  
/* read input options from command line
   should do more error checking for missing args */

  for (i = 1; i < argc; i++) {

    if (!strcmp(argv[i],"-2d"))
      three_d = 0;
    else if (!strcmp(argv[i],"-repeat")) {
      sscanf(argv[i+1],"%d",&nx);
      sscanf(argv[i+2],"%d",&ny);
      i += 2;
      if (three_d) {
	sscanf(argv[i+1],"%d",&nz);
	i++;
      }
      else
	nz = 1;
    }
    else if (!strcmp(argv[i],"-inbox")) {
      inbox = 1;
      sscanf(argv[i+1],"%lg",&in.xlo);
      sscanf(argv[i+2],"%lg",&in.xhi);
      sscanf(argv[i+3],"%lg",&in.ylo);
      sscanf(argv[i+4],"%lg",&in.yhi);
      i += 4;
      if (three_d) {
	sscanf(argv[i+1],"%lg",&in.zlo);
	sscanf(argv[i+2],"%lg",&in.zhi);
	i += 2;
      }
    }
    else if (!strcmp(argv[i],"-outbox")) {
      outbox = 1;
      sscanf(argv[i+1],"%lg",&out.xlo);
      sscanf(argv[i+2],"%lg",&out.xhi);
      sscanf(argv[i+3],"%lg",&out.ylo);
      sscanf(argv[i+4],"%lg",&out.yhi);
      i += 4;
      if (three_d) {
	sscanf(argv[i+1],"%lg",&out.zlo);
	sscanf(argv[i+2],"%lg",&out.zhi);
	i += 2;
      }
    }
    else if (!strcmp(argv[i],"-remap"))
      remap = 1;
    else if (!strcmp(argv[i],"-truein"))
      truein = 1;
    else if (!strcmp(argv[i],"-trueout"))
      trueout = 1;
    else {
      fprintf(stderr,"Syntax error: replicate [options] < infile > outfile\n");
      exit(1);
    }
  }

/* ntotal = total replication factor */

  ntotal = nx*ny*nz;

/* read/write header */

  copyline(2);

  gets(line); sscanf(line,"%d",&natoms);
  gets(line); sscanf(line,"%d",&nbonds);
  gets(line); sscanf(line,"%d",&nangles);
  gets(line); sscanf(line,"%d",&ndihedrals);
  gets(line); sscanf(line,"%d",&nimpropers);

  printf("%d atoms\n",natoms*ntotal);
  printf("%d bonds\n",nbonds*ntotal);
  printf("%d angles\n",nangles*ntotal);
  printf("%d dihedrals\n",ndihedrals*ntotal);
  printf("%d impropers\n",nimpropers*ntotal);

  copyline(1);

  gets(line); sscanf(line,"%d",&ntypes);
  if (nbonds) {gets(line); sscanf(line,"%d",&nbondtypes);}
  if (nangles) {gets(line); sscanf(line,"%d",&nangletypes);}
  if (ndihedrals) {gets(line); sscanf(line,"%d",&ndihedtypes);}
  if (nimpropers) {gets(line); sscanf(line,"%d",&nimprotypes);}

  printf("%d atom types\n",ntypes);
  if (nbonds) printf("%d bond types\n",nbondtypes);
  if (nangles) printf("%d angle types\n",nangletypes);
  if (ndihedrals) printf("%d dihedral types\n",ndihedtypes);
  if (nimpropers) printf("%d improper types\n",nimprotypes);

/* read input box bounds, setup input box */

  copyline(1);

  if (inbox) {
    gets(line);
    gets(line);
    if (three_d) gets(line);
  }
  else {
    gets(line);
    sscanf(line,"%lg %lg",&in.xlo,&in.xhi);
    gets(line);
    sscanf(line,"%lg %lg",&in.ylo,&in.yhi);
    if (three_d) {
      gets(line);
      sscanf(line,"%lg %lg",&in.zlo,&in.zhi);
    }
  }
  
  in.xsize = in.xhi - in.xlo;
  in.ysize = in.yhi - in.ylo;
  if (three_d) 
    in.zsize = in.zhi - in.zlo;
  else
    in.zsize = 0.0;

/* derive output box bounds, write output box */

  if (!outbox) {
    out.xlo = in.xlo;
    if (nx > 1)
      out.xhi = out.xlo + (in.xhi-in.xlo)*nx;
    else
      out.xhi = in.xhi;
    out.ylo = in.ylo;
    if (ny > 1)
      out.yhi = out.ylo + (in.yhi-in.ylo)*ny;
    else
      out.yhi = in.yhi;
    if (three_d) {
      out.zlo = in.zlo;
      if (nz > 1) 
	out.zhi = out.zlo + (in.zhi-in.zlo)*nz;
      else
	out.zhi = in.zhi;
    }
  }

  out.xsize = out.xhi - out.xlo;
  out.ysize = out.yhi - out.ylo;
  if (three_d) out.zsize = out.zhi - out.zlo;
  
  printf("%g %g xlo xhi\n",out.xlo,out.xhi);
  printf("%g %g ylo yhi\n",out.ylo,out.yhi);
  if (three_d) printf("%g %g zlo zhi\n",out.zlo,out.zhi);

/* malloc space for input atoms and topology */

  molecule = (int *) malloc(natoms*sizeof(int));
  type = (int *) malloc(natoms*sizeof(int));
  q = (double *) malloc(natoms*sizeof(double));
  x = (double *) malloc(natoms*sizeof(double));
  y = (double *) malloc(natoms*sizeof(double));
  z = (double *) malloc(natoms*sizeof(double));
  truex = (int *) malloc(natoms*sizeof(int));
  truey = (int *) malloc(natoms*sizeof(int));
  truez = (int *) malloc(natoms*sizeof(int));

  if (molecule == NULL || type == NULL || q == NULL || 
      x == NULL || y == NULL || z == NULL ||
      truex == NULL || truey == NULL || truez == NULL) {
    fprintf(stderr,"Error in atom malloc - no space\n");
    exit(1);
  }

  btype = (int *) malloc(nbonds*sizeof(int));
  bond1 = (int *) malloc(nbonds*sizeof(int));
  bond2 = (int *) malloc(nbonds*sizeof(int));

  if (btype == NULL || bond1 == NULL || bond2 == NULL) {
    fprintf(stderr,"Error in bond malloc - no space\n");
    exit(1);
  }

  atype = (int *) malloc(nangles*sizeof(int));
  angle1 = (int *) malloc(nangles*sizeof(int));
  angle2 = (int *) malloc(nangles*sizeof(int));
  angle3 = (int *) malloc(nangles*sizeof(int));

  if (atype == NULL || angle1 == NULL || angle2 == NULL || angle3 == NULL) {
    fprintf(stderr,"Error in angle malloc - no space\n");
    exit(1);
  }

  dtype = (int *) malloc(ndihedrals*sizeof(int));
  dihed1 = (int *) malloc(ndihedrals*sizeof(int));
  dihed2 = (int *) malloc(ndihedrals*sizeof(int));
  dihed3 = (int *) malloc(ndihedrals*sizeof(int));
  dihed4 = (int *) malloc(ndihedrals*sizeof(int));

  if (dtype == NULL || dihed1 == NULL || dihed2 == NULL || 
      dihed3 == NULL || dihed4 == NULL) {
    fprintf(stderr,"Error in dihedral malloc - no space\n");
    exit(1);
  }

  itype = (int *) malloc(nimpropers*sizeof(int));
  impro1 = (int *) malloc(nimpropers*sizeof(int));
  impro2 = (int *) malloc(nimpropers*sizeof(int));
  impro3 = (int *) malloc(nimpropers*sizeof(int));
  impro4 = (int *) malloc(nimpropers*sizeof(int));

  if (itype == NULL || impro1 == NULL || impro2 == NULL || 
      impro3 == NULL || impro4 == NULL) {
    fprintf(stderr,"Error in improper malloc - no space\n");
    exit(1);
  }

/* read identifier strings one by one in free-form part of data file */

  while (token = gettoken(linetoken)) {

/* read atoms */

    if (!strcmp(token,"Atoms")) {

      for (i = 0; i < natoms; i++) {
	gets(line);
	if (truein)
	  sscanf(line,"%d %d %d %lg %lg %lg %lg %d %d %d",
		 &tmp,&molecule[i],&type[i],&q[i],&x[i],&y[i],&z[i],
		 &truex[i],&truey[i],&truez[i]);
	else {
	  sscanf(line,"%d %d %d %lg %lg %lg %lg",
		 &tmp,&molecule[i],&type[i],&q[i],&x[i],&y[i],&z[i]);
	  truez[i] = truey[i] = truex[i] = 500;
	}
      }
      
/* find number of molecules */

      nmolecules = 0;
      for (i=0; i < natoms; i++)
        if (molecule[i] > nmolecules) nmolecules = molecule[i];

/* replicate set of N atoms as many times as requested
   generate new coords and molecule number and true flags as needed */

      n = 0;
      for (k = 0; k < nz; k++) {
	for (j = 0; j < ny; j++) {
	  for (i = 0; i < nx; i++) {
	    for (m = 0; m < natoms; m++) {
	      n++;
	      newx = x[m] + i*in.xsize;
	      newy = y[m] + j*in.ysize;
	      newz = z[m] + k*in.zsize;
	      tx = truex[m];
	      ty = truey[m];
	      tz = truez[m];
	      if (outbox) scale(in,nx,ny,nz,out,&newx,&newy,&newz,three_d);
	      if (remap) pbc(out,&newx,&newy,&newz,&tx,&ty,&tz,three_d);
              offset = k*ny*nx*nmolecules + j*nx*nmolecules + i*nmolecules;
	      if (trueout) 
		printf("%d %d %d %g %g %g %g %d %d %d\n",
		       n,molecule[m]+offset,type[m],q[m],
		       newx,newy,newz,tx,ty,tz);
	      else
		printf("%d %d %d %g %g %g %g\n",
		       n,molecule[m]+offset,type[m],q[m],newx,newy,newz);
	    }
	  }
	}
      }
    }

/* read bonds and replicate */

  else if (!strcmp(token,"Bonds")) {

    for (i = 0; i < nbonds; i++) {
      gets(line);
      sscanf(line,"%d %d %d %d",&n,&btype[i],&bond1[i],&bond2[i]);
    }
    
    n = 0;
    for (k = 0; k < nz; k++)
      for (j = 0; j < ny; j++)
	for (i = 0; i < nx; i++)
	  for (m = 0; m < nbonds; m++) {
	    n++;
	    del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	    printf("%d %d %d %d\n",n,btype[m],bond1[m]+del,bond2[m]+del);
	  }
  }

/* read angles and replicate */

  else if (!strcmp(token,"Angles")) {

    for (i = 0; i < nangles; i++) {
      gets(line);
      sscanf(line,"%d %d %d %d %d",&n,&atype[i],
	     &angle1[i],&angle2[i],&angle3[i]);
    }
    
    n = 0;
    for (k = 0; k < nz; k++)
      for (j = 0; j < ny; j++)
	for (i = 0; i < nx; i++)
	  for (m = 0; m < nangles; m++) {
	    n++;
	    del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	    printf("%d %d %d %d %d\n",n,atype[m],
		   angle1[m]+del,angle2[m]+del,angle3[m]+del);
	  }
  }

/* read dihedrals and replicate */

  else if (!strcmp(token,"Dihedrals")) {

    for (i = 0; i < ndihedrals; i++) {
      gets(line);
      sscanf(line,"%d %d %d %d %d %d",&n,&dtype[i],&dihed1[i],&dihed2[i],
	     &dihed3[i],&dihed4[i]);
    }

    n = 0;
    for (k = 0; k < nz; k++)
      for (j = 0; j < ny; j++)
	for (i = 0; i < nx; i++)
	  for (m = 0; m < ndihedrals; m++) {
	    n++;
	    del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	    printf("%d %d %d %d %d %d\n",n,dtype[m],
		   dihed1[m]+del,dihed2[m]+del,dihed3[m]+del,dihed4[m]+del);
	  }
  }

/* read impropers and replicate */

  else if (!strcmp(token,"Impropers")) {

    for (i = 0; i < nimpropers; i++) {
      gets(line);
      sscanf(line,"%d %d %d %d %d %d",&n,&itype[i],&impro1[i],&impro2[i],
	     &impro3[i],&impro4[i]);
    }

    n = 0;
    for (k = 0; k < nz; k++)
      for (j = 0; j < ny; j++)
	for (i = 0; i < nx; i++)
	  for (m = 0; m < nimpropers; m++) {
	    n++;
	    del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	    printf("%d %d %d %d %d %d\n",n,itype[m],
		   impro1[m]+del,impro2[m]+del,impro3[m]+del,impro4[m]+del);
	  }
  }

/* non-replicated sections - just copy proper number of lines from in to out */

  else if (!strcmp(token,"Masses"))
    copyline(ntypes);
  else if (!strcmp(token,"Nonbond Coeffs"))
    copyline(ntypes);   /* version 4.0 */
    /* copyline(ntypes*(ntypes+1)/2); */  /* version 3.0 */
  else if (!strcmp(token,"Bond Coeffs"))
    copyline(nbondtypes);
  else if (!strcmp(token,"Angle Coeffs"))
    copyline(nangletypes);
  else if (!strcmp(token,"Dihedral Coeffs"))
    copyline(ndihedtypes);
  else if (!strcmp(token,"Improper Coeffs"))
    copyline(nimprotypes);
  else if (!strcmp(token,"BondBond Coeffs"))
    copyline(nangletypes);
  else if (!strcmp(token,"BondAngle Coeffs"))
    copyline(nangletypes);
  else if (!strcmp(token,"MiddleBondTorsion Coeffs"))
    copyline(ndihedtypes);
  else if (!strcmp(token,"EndBondTorsion Coeffs"))
    copyline(ndihedtypes);
  else if (!strcmp(token,"AngleTorsion Coeffs"))
    copyline(ndihedtypes);
  else if (!strcmp(token,"AngleAngleTorsion Coeffs"))
    copyline(ndihedtypes);
  else if (!strcmp(token,"BondBond13 Coeffs"))
    copyline(ndihedtypes);
  else if (!strcmp(token,"AngleAngle Coeffs"))
    copyline(nimprotypes);
  else {
    fprintf(stderr,"Error in input data file - unknown identifier %s\n",token);
    exit(1);
  }

  }

}


/* ------------------------------------------------------------------- */

/* return a LAMMPS keyword from data file - keyword is entire 2nd line
   echo the 3 lines to stdout */

char *gettoken(char *line)
{
  char *a, *b, *c;
  char dummy[1000];

  a = gets(dummy);
  b = gets(line);
  c = gets(dummy);
  if (b) printf("\n%s\n\n",b);
  return b;
}

/* copy n lines from stdin to stdout */

void copyline(int n)
{
  char line[1000];

  while (n) {
    gets(line);
    puts(line);
    n--;
  }
}

/* point (x,y,z) is somewhere in (or near) the replicated box "in"
   rescale (x,y,z) so it is in the same relative location in box "out"
   do this by 1st converting to 0-1 coord system in "in" space */
   
void scale(struct box in, int nx, int ny, int nz,
	   struct box out, double *x, double *y, double *z, int three_d)
{
  double rel;

  rel = (*x-in.xlo) / (nx*in.xsize);
  *x = out.xlo + rel*out.xsize;
  rel = (*y-in.ylo) / (ny*in.ysize);
  *y = out.ylo + rel*out.ysize;
  if (three_d) {
    rel = (*z-in.zlo) / (nz*in.zsize);
    *z = out.zlo + rel*out.zsize;
  }
}

/* remap the point (x,y,z) so it is guaranteed to be inside the box 
   using standard periodic imaging, update true flags */

void pbc(struct box box, double *x, double *y, double *z, 
	 int *truex, int *truey, int *truez, int three_d)
{
  while (*x < box.xlo) {
    *x += box.xsize;
    *truex -= 1;
  }
  while (*x >= box.xhi) {
    *x -= box.xsize;
    *truex += 1;
  }
  while (*y < box.ylo) {
    *y += box.ysize;
    *truey -= 1;
  }
  while (*y >= box.yhi) {
    *y -= box.ysize;
    *truey += 1;
  }
  if (three_d) {
    while (*z < box.zlo) {
      *z += box.zsize;
      *truez -= 1;
    }
    while (*z >= box.zhi) {
      *z -= box.zsize;
      *truez += 1;
    }
  }
}

