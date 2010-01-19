#include <stdio.h>
#include <math.h>

pfix (nsite, pH, beta, pKint, ida, g, cutoff,
   	nsrp, fprot, rmap, fmap, pKintr, idar, gr)
int nsite;	/* Full number of sites (input)		*/
float pH;	/* (input) */
float beta;	/* 1/kT in units of 1/g's units (input)		*/
float pKint[];	/* Intrinsic pKs for full set (input)	*/
int ida[];	/* =1 for cationic site, =2 anionic (input)	*/
float g[];	/* Full g matrix in unit compatable with beta (input)	*/
float cutoff;       /* Cutoff value for fixing sites (input) */
int *nsrp;	/* Number of sites in reduced set (output)	*/
int fprot[];	/* = 1 or 0 if fixed prot or deprot = -1 if unfixed (output) */
int rmap[];	/* rmap[i] = full set index of reduced set site i (output) */
int fmap[];	/* fmap[i] = reduced set index of full set site i (output) */
   		/* fmap[i] = -1 if i not in reduced set	*/
float pKintr[];	/* Intr. pKs for reduced set (output)	*/
int idar[];	/* Reduced form of ida (output)	*/
float gr[];	/* Rduced form of g (output)	*/

{
   int i,j, nsr;
   float emax, emin, prmax, prmin, wt, sum;

/* Choose the sites to be fixed by calculating the maximum and
   minimum possible protonation for each site			*/

   nsr = 0;
   for (i=0; i<nsite; ++i)
      {
      float arg;
      emax = 0.0;
      emin = 0.0;
      for (j=0; j<nsite; ++j)
	 {
	 if (ida[j] == 1) emax += g[i*nsite + j];
	 if (ida[j] == 2) emin -= g[i*nsite + j];
	 }

      /* Protect against overflow of exp() */
      arg = 2.3025851 * (pKint[i] - pH) - beta * emax;
      if (arg > 20.0) {
	prmin = 1;
      }
      else if (arg < -20) {
	prmin = 0;
      }
      else {
	wt = exp(arg);
	prmin = wt / (1 + wt);
      }

      arg = 2.3025851 * (pKint[i] - pH) - beta * emin;
      if (arg > 20.0) {
	prmax = 1;
      }
      else if (arg < -20) {
	prmax = 0;
      }
      else {
	wt = exp(arg);
	prmax = wt / (1 + wt);
      }

      if (prmax < prmin)
	 printf ("WARNING prmin = %f, prmax = %f\n", prmin, prmax);
      if (prmax < cutoff) fprot[i] = 0;
      else if (prmin > 1.0-cutoff) fprot[i] = 1;
      else 
	 {
	 fprot[i] = -1;
	 ++nsr;
	 }
      }

   for (i=0, j=0; i<nsite; ++i)
      {
      if (fprot[i] == -1)
	 {
	 pKintr[j] = pKint[i];
	 idar[j] = ida[i];
	 rmap[j] = i;
	 fmap[i] = j;
	 ++j;
	 }
      else fmap[i] = -1;
      }

/* Calculate pKint shifts due to fixed sites */
   for(j=0; j<nsr; ++j)
      {
      sum = 0.0;
      for (i=0; i<nsite; ++i)
	 if (fprot[i] != -1)
	    {
	    if (ida[i] == 1 && fprot[i] == 1) /* Site fixed positive */
	       sum += g[nsite*rmap[j] + i];
	    if (ida[i] == 2 && fprot[i] == 0) /* Site fixed negative */
	       sum -= g[nsite*rmap[j] + i];
	    }
      pKintr[j] += -beta * sum / 2.3025851;
      }

/* Make reduced g matrix */
   for (i=0; i<nsr; ++i)
      for(j=0; j<nsr; ++j)
	 {
	 gr[i*nsr + j] = g [nsite*rmap[i] + rmap[j]];
	 }
   *nsrp = nsr;
}
