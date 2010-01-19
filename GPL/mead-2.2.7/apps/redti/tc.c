/*
  Given the number of sites, the intrinsic pK of each site, the site-site
  interaction matrix, and the pH, caluclate the Boltzmann averaged
  protonation of each site.
*/
#include <stdio.h>
#include <math.h>
tc (nsite, pKint, ida, g, pH, beta, cs)
int nsite;		/* Number of titratable sites (input)		*/
float pKint[];		/* Intrinsic pK of sites (pK units) (input)	*/
int ida[];		/* ida=1 for cation 2 for anion	(input)		*/
float *g;		/* site - site interaction matrix 
			   such that g[i*nsite+j] is the potential
			   at site j due to a charge at site i.
			   g must be symmetric and positive definite
			   and have diagonal elements = 0.
			   (usually units of charge^2 / length) (input)	*/
float pH;
float beta;		/* 1/kT (must have inverse units of g) (input)	*/
float cs[];		/* protonation fraction of site	(output)	*/
/* pKint, ida, and cs must have dimension nsite. g must have nsite*nsite */

{
/*
BASIC IDEAS:
    We must do a Boltzmann average over the 2^nsite possible protonation
    configurations of the molecule.  If "state" is a binary number
    representing a protonation state then the ratio of concentrations
    [state]/[0] ("0" is the unprotonated state) is:
    exp (-beta * free_energy + nprot*ln([H]))
    where free_energy = 2.3/beta * pKint * p  +  qs*G*qs / 2 - q0*Gq0
    where pKint, p, qs and q0 are dimension nsite vectors of intrinsic pK,
    protonation and charges of "state" and "0", respectively;
    and G is a matrix of site-site interactions.
    First we calculate each the fractional population of each state by the
    formula above, then the contribution of each state to the protonation
    of each site is summed to give the resulting fractional protonation
    vector, cs.

NOTES FOR IMPROVEMENTS:
    It would be more sensible to have ida=0 for cation and 1 (or -1)
    for anion.
    It would speed and slightly simplify the main loop if sint included
    ln pH in the way that z does.

*/

   int i, j;
   float c;         /* Un-Normalized population of "state" */
   float ctot;      /* Normalazation factor for cs's (sum of all c's) */
   float lnc;       /* ln(c) */
   float lncmax, lnctruemax, lncmin, wsum;
   float z, zmax, zmin;
   float cnone;     /* c of state "0" */
   float lncnone;   /* ln(cnone) */
   float lncH;      /* ln([H]) */
   unsigned int top, state, x;
   unsigned int *mask;
                    /* top, state, x, and mask are used in a binary scheme
		       for keeping track of protonation states. */
   float esum;      /* beta * free_energy */
   float nprot;     /* Number of protons in a given state */
   float *sint, *w;
                    /* Used in scheme to simplify free_energy formula.
		       sint*e + e*w*e = beta times free energy	
		       relative to deprotonated state where	
		       e=protonation vector			
		     */
   extern void *calloc();

/*
WARNING: It is assumed here that sizeof * 8 gives the number of bits
         in an unsigned integer (i.e. char are 8 bits).
	 This routine cannot handle more sites than there are bits in an
	 unsigned integer (It would be too slow even if it could.)
*/
   if (nsite > 8 * (sizeof (unsigned int)))
      printf ("ERROR nsite > size of unsigned integer = %d\n", 
   		8 * (sizeof (unsigned int)));


/*  SET UP SOME PRELIMINARIES:
    Calculation of free_energy can be simplified to a formula of the form
    beta * free_energy = sint * p + p*W*p
    where p is the protonation vector, sint = -2.3*pKint + beta*g*q0
    (that is, sint[i] = beta * free-en of adding proton i to deprot state)
    and W = beta*g/2.
    This way sums can be done over protation vectors rather than charge
    vectors.
    PROTONATION STATES
    are represented by the variable "state".  Each bit of state represents
    a site which is protonated/deprotonated if the bit is one/zero.
    To go through all states, state takes values from 0 to top = 2^nsite.
    (Actually state=0 is the reference and isn't done explicitly.
    "Mask" is used to pull the protonation information out of state.
    mask[i]'s i'th bit is one and all other bits are zero so that
    mask[i]&state == 1 if site i is protonated in state.
    OVERFLOW/UNDERFLOW PREVENTION:
    If the logarithms get too high or low there will be trouble.  To
    prevent this, precalculate the lowest and highest logs that can
    occur and subtract lncmax from all lnc's calculated so the lnc
    never exceeds zero.  Still could be trouble if lnc is too low,
    but that can be handled easily by setting exp(lnc) = 0. 
*/

   if (nsite == 0) return (1);
   if (nsite < 0) {
     printf ("WARNING: tc called with nsite < 0\n");
     return (0);
   }

   lncH = -2.3025851*pH;
   mask = (unsigned int *) calloc ((unsigned) nsite, sizeof (unsigned int));
   sint = (float *) calloc ((unsigned) nsite, sizeof (float));
   w = (float *) calloc ((unsigned) nsite*nsite, sizeof (float));

   x=1;
   top = 1;
   lncmax = lncmin = 0.0;
   wsum = 0.0;
   for (i=0; i<nsite; ++i)
      {
      cs[i] = 0;
      mask[i] = x << i;
      top *= 2;
      sint[i] = -2.3025851 * pKint[i];
      for (j=0; j<nsite; ++j)
	 {
	 if (ida[j] == 2) 	/* if site is anionic	*/
	    sint[i] += -beta * g[i*nsite + j];
	 w[i*nsite+j] = beta * g[i*nsite + j] / 2.0;
	 if (w[i*nsite+j] < 0.0) 
	   printf ("ERROR w[%d,%d] = %f < 0. Bad g matrix?\n",
		   i, j, w[i*nsite+j]);
	 else wsum += w[i*nsite+j];
	 }
      z = lncH - sint[i];
      if (z > 0.0) lncmax += z;
      else lncmin += z;
      if (i==0) zmax = zmin = z;
      else {
	if (z > zmax) zmax = z;
	if (z < zmin) zmin = z;
      }
   }

/*  If lncmax/min is still zero it means all z's were negative/positive,
    so set lncmax/min to highest/lowest value since at least one of the
    z's will occur in all the summations
*/
   if (lncmax == 0) lncmax = zmax;
   if (lncmin == 0) lncmin = zmin;

/* W term can only decrease lnc, so subtract wsum from lncmin */
   lncmin -= wsum;

/* Some little safety margins: */
   lncmax += 0.0001;
   lncmin -= 0.0001;

/* This loop terminates when lncmax estimate is ok */
 while(1)
 {

/*  THE MAIN LOOP OVER ALL PROTONATION STATES
*/
   ctot = 0;
   lnctruemax = 0.0;
   for (state=1; state<top; ++state)
      {
/*  Loop over sites gathering free_energy and number of proton sums.
*/
      esum = 0;
      nprot = 0;
      for (i=0; i<nsite; ++i)
	 if (state & mask[i])
	    {
	    ++nprot;
	    esum += sint[i];
	    for (j=0; j<nsite; ++j)
	       if (state & mask[j])
	          {
	          esum += w[i*nsite + j];
		  }
	    }
/* Calculate and rescale the ln of concentration of this state.
   Keep track of the maximium log that actually occurs so re-rescaling
   can be done if needed.
*/
      lnc = nprot*lncH - esum;
      if (lnc > lnctruemax) lnctruemax = lnc;
      lnc -= lncmax;

/* Check whether lnc has exceeded precalculated bounds */
      if (lnc > 0.0) {
	printf ("\nWARNING lnc > 0 (lnc = %f, lncmax = %f)\n", lnc, lncmax);
	printf ("pH = %f, state = %d\n", pH, (int) state);
	for (i=0; i<nsite; ++i) {
	  printf ("  pKint[%d] = %f ", i, pKint[i]);
	  if (state & mask[i])
	    printf (" (protonated)\n");
	  else
	    printf (" (unprotonated)\n");
	}
      }
      if (lnc + lncmax < lncmin) {
	printf ("\nWARNING lnc+lncmax < lncmin (lnc+lncmax=%f, lncmin=%f)\n",
		lnc+lncmax, lncmin);
	printf ("pH = %f, state = %d\n", pH, (int) state);
	for (i=0; i<nsite; ++i) {
	  printf ("  pKint[%d] = %f ", i, pKint[i]);
	  if (state & mask[i])
	    printf (" (protonated)\n");
	  else
	    printf (" (unprotonated)\n");
	}
      }

/* A "safe" exp function for calculating concentration of this state */
      if (lnc > -60.0) c = exp (lnc);
      else c = 0;

/* This state's contributions to the fractional protonations of each site */
      for (i=0; i<nsite; ++i)
	 if (state & mask[i])
	    {
	    cs[i] += c;
	    }

/* Keep track of total concentration for normalization later. */
      ctot += c;
      }
/* Do the fully unprotonated state.  */
   lncnone = -lncmax;
   if (lncnone > -60.0) cnone = exp(lncnone);
   else cnone = 0;
   ctot += cnone;

/* See if re-rescaling is needed.  If so, reset lncmax and repeat. */
   if (lncmax - lnctruemax > 30)
     {
       lncmax = lnctruemax + 0.001; /* safety margin */
     }
   else break;
 }

   free (mask);
   free (sint);
   free (w);

   if (ctot <= 0)
      {
      printf ("ERROR tc fails due to ctot <= 0\n");
      return (0);
      }
   for (i=0; i<nsite; ++i)
      {
      cs[i] = cs[i] / ctot;
      }
   return (1);
}
