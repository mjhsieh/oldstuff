#include <stdio.h>
#include <stdlib.h>

void * emalloc (n)
size_t n;
{
    void *p;

    p = malloc(n);
    if (p == NULL)
	printf ("out of memory\n");
    return p;
}

/*
rfixpro
This routine calculates the protonation state of a multi site
molecule at a given pH by Reducing the number of sites by FIXing
some of them prior to the PROtonation calculation.
It uses pfix to set up the fixed sites and tc to do the protonation
calculation.
Return value is the number of "titratable" sites included in the reduced set.
*/

int rfixpro (ns, pH, beta, pKint, ida, g, cutoff, cs, tot)
     int ns;		/* Number of titratable sites			*/
     float pH;
     float beta;        /* 1 / kT */
     float pKint[];	/* Intrinsic pK of sites (input)		*/
     int ida[];         /* ida=1 for cation 2 for anion	 (input)	*/
     float g[];   	/* site - site interaction matrix (input)	*/
     float cutoff;       /* Cutoff value for fixing sites (input) */
     float cs[];        /* degree of protonation of sites (output) */
     float *tot;        /* total protonation (output) */
{
  int nsr;	/* Number of sites in reduced set 	*/
  int *fprot;	/* = 1 or 0 if fixed prot or deprot = -1 if unfixed  */
  int *rmap;	/* rmap[i] = full set index of reduced set site i  */
  int *fmap;	/* fmap[i] = reduced set index of full set site i  */
   		/* fmap[i] = -1 if i not in reduced set	*/
  float *pKintr;	/* Intr. pKs for reduced set 	*/
  int *idar;	/* Reduced form of ida 	*/
  float *gr;	/* Reduced form of g 	*/
  float *csr;   /* Reduced form of cs */
  int i;

  fprot = (int *) emalloc ((unsigned) ns * sizeof (int));
  rmap = (int *) emalloc ((unsigned) ns * sizeof (int));
  fmap = (int *) emalloc ((unsigned) ns * sizeof (int));
  pKintr = (float *) emalloc ((unsigned) ns * sizeof (float));
  idar = (int *) emalloc ((unsigned) ns * sizeof (int));
  gr = (float *) emalloc ((unsigned) ns*ns * sizeof (float));
  csr = (float *) emalloc ((unsigned) ns * sizeof (float));

  pfix (ns, pH, beta, pKint, ida, g, cutoff,
	&nsr, fprot, rmap, fmap, pKintr, idar, gr);
  tc (nsr, pKintr, idar, gr, pH, beta, csr);
  *tot = 0.0;
  for (i=0; i<ns; ++i)
    {
      if (fmap[i] != -1)
	cs[i] = csr[fmap[i]];
      else 
	cs[i] = (float) fprot[i];
      *tot += cs[i];
    }

  free ( (char *) fprot);
  free ( (char *) rmap);
  free ( (char *) fmap);
  free ( (char *) pKintr);
  free ( (char *) idar);
  free ( (char *) gr);
  free ( (char *) csr);

  return (nsr);
}

