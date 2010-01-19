#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define ns (100)
#define maxstr (80)
#define npHmax (1000)

void error (char * s)
{
  fputs (s, stderr);
  exit(1);
}

main(int argc, char* argv[])
{
  int nsite;		/* Number of titratable sites			*/
  float pKint[ns];	/* Intrinsic pK of sites			*/
  int ida[ns];		/* ida=1 for cation 2 for anion			*/
  float g[ns*ns];	/* site - site interaction matrix 		*/
  int nsr;		/* Reduced number of sites			*/
  float pK[ns];	/* Effective pK from titration curves		*/
  float pH;
  float beta;
  int i, j, nsitetest;
  float T;			/* Temperature		*/
  float cs[ns];
  FILE *gf;
  FILE *pKf;
  float tot;
  float tcur[ns][npHmax];
  float pHcur[npHmax];
  float pHmax = 30.1;
  float pHmin = -20.0;
  float pHstep = 0.2;
  int npH;
  char sname[ns][maxstr];
  char molname[maxstr];
  char filename[maxstr];
  char stype[ns];
  char line[ns];
  int dry = 0;
  
  char *ptr;
  int iarg;
  char ustring[] = "Usage: redti [-cutoff f] [-pHrange f f] molname\n";
  float cutoff = 0.05;  /* Default value of pfix cutoff */
  
  
  T = 300;
  beta = 1.67112e5 / T;	/* for g in electrons**2 / Angs, T in kelvin */
  
  for (iarg = 1; iarg<argc; ++iarg) {
    if (iarg+1<argc) {
      if (!strcmp(argv[iarg], "-cutoff")) {
	cutoff = (float) strtod(argv[++iarg],&ptr);
	if (ptr<=argv[iarg]) error (ustring);
	if (cutoff > 1.0 || cutoff < 0.0) {
	  fprintf (stderr, "%f ", cutoff);
	  error ("is absurd cutoff value.\n");
	}
      }
      else if (!strcmp(argv[iarg], "-pHrange")) {
	pHmin = (float) strtod(argv[++iarg],&ptr);
	if (ptr<=argv[iarg]) error (ustring);
	pHmax = (float) strtod(argv[++iarg],&ptr);
	if (ptr<=argv[iarg]) error (ustring);
	if (pHmin >= pHmax) {
	  fprintf (stderr, "%f to %f ", pHmin, pHmax);
	  error ("is absurd pH range.\n");
	}
      }
      else if (!strcmp(argv[iarg], "-dry"))
	dry = 1;
      else
	error(ustring);
    }
    else {
      if (strlen (argv[iarg]) > maxstr - 7)
	error ("molname is too long\n");
      strcpy (molname, argv[iarg]);
    }
  }
  
  strcpy (filename, molname);
  strcat (filename, ".pkint");
  if ((gf = fopen(filename, "r")) == NULL) {
    printf ("Could not open file %s.\n", filename);
    error ("\n");
  }
  
  for (i=0; fgets(line,80,gf) != NULL; ++i) {
    if (i >= ns) {
      printf ("WARNING: more than %d sites in %s\n", ns, filename);
      break;
    }
    if (sscanf(line, "%f %1s %15[^\n]", pKint+i, stype+i, sname[i]) != 3)
      printf ("WARNING: error reading line %d of %s\n", i+1, filename);
    nsite = i+1;
  }
  fclose (gf);
  
  strcpy (filename, molname);
  strcat (filename, ".g");
  getg (filename, &nsitetest, g);
  if (nsitetest != nsite) printf ("WARNING getg gets different nsite\n");
  
  /* Call to aveg probably not needed.  Leave it in anyway. */
  aveg (nsite, g);
  
  
  for (i=0; i<nsite; ++i)
    switch (stype[i])
      {
      case 'C' : ida[i] = 1; break;
      case 'A' : ida[i] = 2; break;
	default : 
	  printf ("WARNING stype[%d] = %c\n", i, stype[i]);
	ida[i] = 200;
	break;
      }
  
  printf (" i   pKint   ida\n");
  for (i=0; i<nsite; ++i)
    {
      printf (" %d  %f  %d  \n", i, pKint[i], ida[i]);
    }
  
  for (i=0; i<nsite; ++i)
    pK[i] = pHmin-1.0;
  
  npH = 0;
  if ((pHmax-pHmin)/pHstep >= npHmax) 
    error ("ERROR number of pH steps > npHmax\n");
  
  for (pH=pHmin; pH<pHmax; pH += pHstep)
    {
      if (npH > npHmax) {
	printf ("ERROR npH > npHmax\n");
	exit (1);
      }
      if (dry) {
	nsr = rfixprosee (nsite, pH, beta, pKint, ida, g, cutoff, cs, &tot);
	printf ("pH = %f, nsr = %d\n", pH, nsr);
	continue;
      }

      pHcur[npH] = pH;
      nsr = rfixpro (nsite, pH, beta, pKint, ida, g, cutoff, cs, &tot);
      printf ("pH = %f, nsr = %d tot = %f\n", pH, nsr, tot);
      
      for (i=0; i<nsite; ++i)
	{
	  tcur[i][npH] = cs[i];
	  if (npH>0 && tcur[i][npH] < 0.5 && tcur[i][npH-1] >= 0.5)
	    pK[i] = (0.5 - tcur[i][npH])/(tcur[i][npH] - tcur[i][npH-1])
	      * (pHcur[npH] - pHcur[npH-1])  + pHcur[npH];
	}
      /*
	fprintf (tctot,"%f %f\n", pH, tot);
	*/
      ++npH;
    }
  
  if (dry)
    exit (0);

  for(i=0; i<nsite; ++i)
    if (tcur[i][npH-1] > 0.5) pK[i] = pHmax+1.0;
  
  
  strcpy (filename, molname);
  strcat (filename, ".pkout");
  if ((pKf = fopen(filename, "w")) == NULL) {
    printf ("Could not open file %s.\n", filename);
    error ("\n");
  }
  
  fprintf (pKf, "             pKint          pK         pKexp\n");
  for (i=0; i<nsite; ++i)
    {
      fprintf (pKf, "%7s %12.3f ", sname[i], pKint[i]);
      if (pK[i] < pHmin)
	fprintf (pKf, " <%f ", pHmin);
      else if (pK[i] > pHmax)
	fprintf (pKf, " >%f", pHmax);
      else
	fprintf (pKf, "%12.3f ", pK[i]);
      fprintf (pKf, "\n");
    }
  
  /*  A quick hack to write out some individual titration curves 
   */
  
  for (i=0; i<nsite; ++i) {
    fprintf (pKf , "Site %s\n", sname[i]);
    for (j=0; j<npH; ++j)
      if (tcur[i][j] > 0.05 && tcur[i][j] < 0.95) 
	fprintf (pKf, "%f %f\n", pHcur[j], tcur[i][j]);
  }
  
  /* A bit of code to print out a plottable full titration curve on the tail
     of the pK output file
     */
  fprintf (pKf, "Whole protein\n");
  for (j=0; j<npH; ++j) {
    for (i=0, tot=0.0; i<nsite; ++i)
      tot += tcur[i][j];
    if (tot > 0.1 && tot < ((float) nsite) - 0.1)
      fprintf (pKf, "%f %f  \n", pHcur[j], tot);
  }
  
  /*
    {
    FILE *tcdat;
    FILE *tctot;
    char fname[20];
    
    printf ("npH = %d, writing out curves for each site\n", npH);
    for (i=0; i<nsite; ++i)
    {
    sprintf (fname, "C1S%d.DAT", i);
    if (NULL == (tcdat = fopen (fname, "w")))
    printf ("WARNING error opening file %s\n", fname);
    for (j=0; j<npH; ++j)
    fprintf (tcdat, "%f %f\n", pHcur[j], tcur[i][j]);
    fclose (tcdat);
    }
    }
    */
}


/* Read in g matrix from file produced by matmake with file name fname.
   if nsite == 0 on input nsite is determined from the file, otherwise
   nsite is taken as the correct value and program complains if file
   doesn't conform.	*/

getg (fname, nsitep, g)
     char *fname;		/* input */
     int *nsitep;		/* input/output */
     float *g;		/* output */
{
  int i, j, teof, ir, jr, nsite;
  FILE *gf;
  if ((gf = fopen(fname, "r")) == NULL) {
    printf ("Could not open file %s.\n", fname);
    error ("\n");
  }
  nsite = *nsitep;
  i = j = 0;
  while (1)
    {
      teof = fscanf (gf, " %d %d", &ir, &jr);
      if (teof != 2) break;
      if (j>0 && jr == 1) 
	{
	  if (nsite == 0)
	    nsite = j;
	  else if (j != nsite)
	    printf ("WARNING row i ends too soon\n");
	  j = 0;
	  ++i;
	}
      if (ir != i+1)
	printf("WARNING ir = %d read when %d expected\n", ir, i+1);
      if (jr != j+1)
	printf("WARNING jr = %d read when %d expected\n", jr, j+1);
      if (1 != fscanf(gf, " %f ", g + i*nsite + j)) 
	printf ("ERROR reading gfile\n");
      ++j;
    }
  if (i != nsite-1 || j != nsite)
    printf ("WARNING error reading gfile or unexpected end of file\n");
  fclose (gf);
  *nsitep = nsite;
}


/* Average i,j and j,i components, zero the diagonal */
aveg (nsite, g)
     int nsite;
     float *g;
{
  int i, j;
  float gave;
  
  for (i=0; i<nsite; ++i)
    {
      g[i*nsite + i] = 0.0;
      for (j=0; j<i; ++j)
	{
	  gave = (g[i*nsite + j] + g[j*nsite + i]) / 2.0;
	  g[i*nsite + j] = gave;
	  g[j*nsite + i] = gave;
	}
    }
}

