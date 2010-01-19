/* ------------------------------------------------------------------------
 * MD Display Version 3.0
 * Copyright (C) 2002 Chris Moth and Terry Lybrand
 * Vanderbilt University
 *
 *    Derived from MD Display Version 2.0
 *
 *  Copyright (c) 1990  Timothy J. Callahan, Eric Swanson, and Terry Lybrand
 *
 *  Please cite:
 *
 *    Callahan, T.J., Swanson, E., Lybrand, T.P.:
 *       MD Display: An interactive graphics program for visualization
 *       of molecular dynamics trajectories,
 *    Journal of Molecular Graphics, 1996, 14:39-41
 *
 * Do not distribute this code without explicit permission.
 * Do not incorporate into other codes without explicit permission.


 preproc.c reads a variety of MD formats, and makes a binary movie
 file which can be displayed by mddisplay.  My contribution to this module
 has been quite thin - mainly converting it to ANSI C, making it Amber 7
 compatible, and fixing a few small bugs.

   6/25/90 -C option to read in CONECT cards added
      idea and code from  Mike Siani

   7/15/90 Added "#mol" extension; also from MAS

   7/22/90 Added "nojump()" to undo periodic boundaries
        (useful if you want to look at filtered trajectories)

   8/10/90 Added feature to display H-bonds.

   10/3/91 (E.S.)  Version 1.1
           Parse -f n1:n2 to select range of frames.
           Added header record to coord file.

    4/9/92 (E.S.)  Version 1.1b
           Fixed testfloat108 to detect end of file.
           Fixed testfloat108 and amberparm_readfloat108 to read 10F8.3 format
           correctly (fixed bug when field had -123.456).
           Added "-j" option; nojump is not called if -j specified.
           Fixed nojump to keep residues intact.
           Fixed bug when reading >100 frames.  Now only reads first
           100 frames unless more are requested (-f option).

   10/6/92 Modifications from the University of South Alabama
           Support for Charmm22 trajectories
           -m  <formatted PSF File>
           -c  <binary DCD File>
           No support for frame control
           No support for adding hydrogen bonds
           Error Checking is minimal.

   10/28/93  (E.S.)
           Handle conventional pdb names for H and D.
    1/26/94  (E.S.)
           Changed bonding array to int (from short) to handle
           large structures.
      6/94   (E.S.)
           Added support for GROMOS and Biosym trajectories.

  Jan/2001 Brought code to ANSI C compliance, improved error reporting
           removed some register declarations that might harm performance
           (Chris Moth)

           Updated Version # to 2.2

   BUGS: scanf, many shorts  Why not floats?
         stack usage.
         why not count (look at input file size) and set runtime
        maximums that way?

    BUG: Hbond parameters - how do they work?  What are excluded list?

 BUG: Do you use the exit() values somewhere?  Howbout errors on stderr?
      CWM took liberty of reducing all exit()s to
         exit(EXIT_SUCCESS) and exit(EXIT_FAILURE)....  may need to be undone

  BUG: Changed fread() to fread_or_exit_if_error()

  BUG: some dependence on sizeof(int)=32.  What about 64bit machines???

  BUG: upgraded filename lengths to FILENAME_MAX to allow for complex
       path names or shell expansions of path names to preprocessor.

  BUG: Desparately need format documentation.

  BUG: We use lots of floats - which causes Microsoft Visual C++ to complain about
      the need to typecast....  Are floats really more efficient in current times?

      Added -d Debug flag

*/



#include <stdio.h>
#include <stdlib.h>

#include <fcntl.h>
#include <math.h>

#include <errno.h>
/* #include <io.h>  THIS IS A BUG - NON ANSI C for open() call - see TERRY */


#include <ctype.h>
#include <string.h>

#include <assert.h>

#include "displib.h"

#define VERSION "3.0"     /* program version */
#define COORDVERSION 101  /* coord file version number */


typedef float float3[3]; /* A triplet in float format, often an x,y,z */

typedef int   int3[3];   /* A triplet in int format, often an x,y,z */

typedef short short3[3]; /* A triplet in short format, often an x,y,z */


static char bondFileName[FILENAME_MAX],
            attribFileName[FILENAME_MAX],
            binFileName[FILENAME_MAX],
            surfaceFileName[FILENAME_MAX];
/*          name[FILENAME_MAX];  Replaced by defaultFileName in main() */

#define MAXCOORDFILES 30
int nCoordFiles=0;

static char *pdbFileName=0,*coordFileNames[MAXCOORDFILES],parameterFileName[FILENAME_MAX]="parm.fmt",*insurf=0;
static char *psfFileName=0, *bshFileName=0, *gromos_mtFileName;
/* BUG - what are these variables */

static int frames=9999,fr1=0,fr2=0,frinc=1;
static int conpress=0,oneIfInBox=0,jumpok=0,gxscale;
static float fbox[3],hbcut=(-1.0);
static short nmol=0, /* # of Molecules.  Incremented differently in each format */

             mol[MAXATM];
static int conectflag=0;

typedef char rlabtype[5];
static rlabtype* rlab;

void allocate_global_rlab(int nres,const char* calling_function_name)
{
   char error_context[200];

   assert(rlab == 0); // Make sure this is called only once!

   if (nres <= 0)
      {
      printf("\npreproc:%s:  %d residues in source data.  At least one is needed\n",calling_function_name,nres);

      exit(EXIT_FAILURE);
      }
   sprintf(error_context,"preproc:allocate_global_rlab - %d residues function: %s",nres,calling_function_name);
   rlab = malloc_or_exit_if_error(nres*sizeof(rlabtype),error_context);
   memset(rlab,0,nres*sizeof(rlabtype));
}

static atomattrib* attribs;

void allocate_global_atomattrib(int natom,int nres,const char* calling_function_name)
{
   char error_context[300];
   assert(attribs == 0); // Make sure this is called only once!
   sprintf(error_context,"preproc:allocate_global_rlab %d atoms %d residues - function: %s",natom,nres,calling_function_name);
   attribs = malloc_or_exit_if_error((nres+natom+4)*sizeof(atomattrib),error_context);
   memset(attribs,0,(nres+natom+4)*sizeof(atomattrib));
}


/* Default atom names & radii.  Go ahead and twiddle with them. */
enum {nAtomTypes=16};

struct
   {
   const char* symbol;
   float radius;
   } atom[nAtomTypes] =
  {{"C", 0.85f},  /* Carbon */

   {"N", 0.85f},  /* Nitrogen */

   {"O", 0.75f},  /* Oxygen */

   {"S", 1.05f},  /* Sulfur */

   {"P", 0.96f},  /* Phosphorous */

   {"F", 0.60f},  /* Flourine */

   {"CL",1.00f},  /* Chlorine */

   {"BR",1.15f},  /* Bromine */

   {"I", 1.35f},  /* Iodine */

   {"H", 0.40f},  /* Hydrogen */

   {"A", 0.85f},  /* Ambiguity */

   {"M", 1.32f},
   {"ZN",0.01f},  /* Zinc */

   {"Ca",0.01f},  /* Calcium */

   {"FE",0.01f},  /* Iron */

   {"D", 0.40f}}; /* Deuterium */


#ifdef OLD_STUFF
static char atype[][5]={"C","N","O","S","P","F","CL","BR","I","H","A","M",
                    "ZN","Ca","FE","D"};

static float atrad[16]={0.85, 0.85, 0.75, 1.05, 0.96, 0.60, 1.00,
                       1.15, 1.35, 0.40, 0.85, 1.32, 0.01, 0.01, 0.01, 0.40};
#endif






/* Function prototypes */



/*------------------------------------------------------------*/

/* Read in a .pdb format file. */

/* The data format as of this update to ANSI C 1/8/2000 was found at: */

/* http://www.hhmi.umbc.edu/toolkit/structure/procheck/manappa.html */

/* For work in 'C', offsets start at 0 of course. */


/* A better source of file format documentation is at */

/* http://www.rcsb.org/pdb */


/*
Appendix A - Brookhaven file format

The table below shows the Brookhaven file format for the coordinate records (ie ATOM and HETATM) of your PDB file.
Each record holds the coordinates and other details of a single atom.

---------------------------------------------------------------------------
Field |    Column    | FORTRAN |
  No. |     range    | format  | Description
---------------------------------------------------------------------------
   1. |    1 -  6    |   A6    | Record ID (eg ATOM, HETATM)
   2. |    7 - 11    |   I5    | Atom serial number
   -  |   12 - 12    |   1X    | Blank
   3. |   13 - 16    |   A4    | Atom name (eg " CA " , " ND1")
   The following rules are used in assigning atom names.
      * Greek letter remoteness codes are transliterated as follows:
        alpha = A, beta = B, gamma = G, delta = D, epsilon = E, zeta = Z,
        eta = H, etc.

      * Atoms for which some ambiguity exists in the crystallographic results
        are designated A. This usually applies only to the terminal
        atoms of asparagine and glutamine and to the ring atoms of
        histidine.

      * The extra oxygen atom of the carboxy terminal amino acid
        is designated OXT.

      * Six characters (columns) are reserved for atom names, assigned as follows.

         COLUMN     VALUE
         -----------------------------------------------------------------------
         13 - 14    Chemical symbol - right justified, except for hydrogen atoms
         15         Remoteness indicator (alphabetic)
         16         Branch designator (numeric)

   4. |   17 - 17    |   A1    | Alternative location code (if any)
   5. |   18 - 20    |   A3    | Standard 3-letter amino acid code for residue
   -  |   21 - 21    |   1X    | Blank
   6. |   22 - 22    |   A1    | Chain identifier code
   7. |   23 - 26    |   I4    | Residue sequence number
   8. |   27 - 27    |   A1    | Insertion code (if any)
   -  |   28 - 30    |   3X    | Blank
   9. |   31 - 38    |  F8.3   | Atom's x-coordinate
  10. |   39 - 46    |  F8.3   | Atom's y-coordinate
  11. |   47 - 54    |  F8.3   | Atom's z-coordinate
  12. |   55 - 60    |  F6.2   | Occupancy value for atom
  13. |   61 - 66    |  F6.2   | B-value (thermal factor)
   -  |   67 - 67    |   1X    | Blank
  14. |   68 - 68    |   I3    | Footnote number
---------------------------------------------------------------------------

Each molecule is separated by a "TER" record as described here:

TER

Overview

The TER record indicates the end of a list of ATOM/HETATM records for a chain.

Record Format

COLUMNS         DATA TYPE         FIELD        DEFINITION
-------------------------------------------------------------------------
 1 -  6         Record name       "TER   "
 7 - 11         Integer           serial       Serial number.
18 - 20         Residue name      resName      Residue name.
22              Character         chainID      Chain identifier.
23 - 26         Integer           resSeq       Residue sequence number.
27              AChar             iCode        Insertion code.
*/

/*------------------------------------------------------------*/
/* This routine is called by qsort to compare two integers */
/* And returns which ever is lower */
/* See, from sort() below, that ap abd bp actual point to */
/* an array of two integers, the second of which is the key */


typedef struct {
   int index;    /* The associated index entry from original array */

   int keyvalue; /* The value that will be compared */

   } sortKeyType;

static int sortKeyCompareFunction(const void* ap,const void* bp)
{
   return ((sortKeyType*)ap)->keyvalue - ((sortKeyType*)bp)->keyvalue;

#ifdef OLDSTUFF
   int ans;

   if (ap[1]>bp[1]) ans=1; else ans=(-1);

return(ans);
#endif
}

/*------------------------------------------------------------*/

void sort(int* key,int* indx,int n)
{
   register int i;
/* register short t; */


   sortKeyType* sortKeys = (sortKeyType*)malloc_or_exit_if_error(n*sizeof(sortKeyType),"sort - sortkeys");

   sortKeyType* skPtr = sortKeys;

   /* Copy indices and keys to temp array */

   for (i=0;i<n;i++)
      {
      skPtr->index = i+1;
      skPtr->keyvalue = key[i+1];
      skPtr++;
      }

   qsort((void *)(sortKeys), n, sizeof(sortKeyType), sortKeyCompareFunction);

   skPtr = sortKeys;
   for (i=0;i<n;i++)
      {
      indx[i+1]=skPtr->index;
      skPtr++;
      }
   free(sortKeys);
}


/* Caveat: buf[bufsize] must be a legal memory reference */


int ParseInteger(char* buf,int bufsize)
{
   int retval;
   char temp = buf[bufsize]; /* Save buf[bufsize] */

   buf[bufsize]=0;           /* zero the end */

   retval =atoi(buf);           /* convert to integer */

   buf[bufsize]= temp;       /* Restore buf[bufsize] */

   return retval;            /* Return integer to caller */

}


float ParseFloat(char* buf,int bufsize)
{
   float retval;
   char temp = buf[bufsize]; /* Save buf[bufsize] */

   buf[bufsize]=0;           /* zero the end */

   retval = (float)atof(buf);         /* convert to integer */

   buf[bufsize]= temp;       /* Restore buf[bufsize] */

   return retval;            /* Return integer to caller */

}


static void ParseStringSupport(char* dest, int dest_size,
                             const char* src,int src_size,
                          int pad_with_spaces)
{
   /* Skip leading blanks */

   while (src_size && (*src == ' '))
      {
      src++;
      src_size--;
      }

   /* Skip trailing blanks */

   while (src_size)
      {
      if (src[src_size-1] != ' ')
         break;
      src_size--;
      }

   if (src_size == 0)
      *dest = 0;
   else
      {
      if (src_size >= dest_size)
         src_size = dest_size - 1;

      memcpy(dest,src,src_size);

      dest[src_size] = 0;
      }

   if (pad_with_spaces)
      {
      /* Remember destination size includes a 0 byte */

      int remaining = (dest_size-1) - src_size;
      if (remaining > 0)
         memset(dest+src_size,' ',remaining);
      }

   dest[dest_size-1] = 0;
}


#define macroParseStringAndPad(dest,src,src_size) ParseStringSupport(dest,sizeof(dest),src,src_size,1);


#define macroParseString(dest,src,src_size) ParseStringSupport(dest,sizeof(dest),src,src_size,0);


/*------------------------------------------------------------*/

void WriteAttributeFile(int nAttributes)
{
   /* BUG - changed to ANSI compatible fopen call.. - what about permissions??? */

   FILE* f=fopen_or_exit_if_error(attribFileName,"wb");

   printf("\nWriting attribute file: %s\n",attribFileName);
   if (debugFlag)
      printf("\tnAttributes=%d\n\n",nAttributes);

   fwrite_or_exit_if_error(attribs,sizeof(atomattrib),nAttributes,f);
   fclose_or_exit_if_error(attribFileName,f);
}

void WriteBondFile(int nBonds,int ib[][3])
{
   /* BUG - changed to ANSI compatible fopen - what about permissions */

   FILE* f = fopen_or_exit_if_error(bondFileName,"wb");

   printf("Writing bond file: %s\n",bondFileName);
   if (debugFlag)
      printf("\tnBonds=%d\n",nBonds);

   /* BUG - need error checking on write */

   fwrite_or_exit_if_error(ib,3*sizeof(int),nBonds,f);
   fclose_or_exit_if_error(bondFileName,f);
}

void WriteBinaryCoordinateFileWithAllFrames(
                            short* xx,
                            int natom,
                            int bytesperframe,
                            int frames)
{
   int header[5];
   FILE* f = fopen_or_exit_if_error(binFileName,"wb");

   printf("Writing Binary Coord File with all frames: %s\n",binFileName);
   if (debugFlag)
      printf("\tnatom=%d  bytesperframe=%d  frames=%d\n",
                        natom,bytesperframe,frames);

   header[0]=COORDVERSION;
   header[1]=natom;
   header[2]=bytesperframe;
   header[3]=frames;
   header[4]= -COORDVERSION;

   assert(sizeof(header) == 20);
   fwrite_or_exit_if_error(header,sizeof(header),1,f);

   fwrite_or_exit_if_error(xx,bytesperframe,frames,f);
   fclose_or_exit_if_error(binFileName,f);
}

void WriteSurfaceFile(int n,const short* surface)
{
   FILE* f = fopen_or_exit_if_error(surfaceFileName,"wb");

  fwrite_or_exit_if_error((char *)surface,sizeof(short),n,f);

   fclose_or_exit_if_error(surfaceFileName,f);
}



/*------------------------------------------------------------*/

int stcp(const char* a,const char* b)
{
   int i=0;

   /* While we have both a and b to compres, and they are the same... */

   /* BUG::: We need to document what this is all about. */

   /* What if a and b are eactly the same???? */


   while (i<6 && a[i]==b[i] && a[i])
      i++;

   /* Now that we're past the matching part... quick analysis... */

   if (!b[i] || b[i]==32 || b[i]=='*' || a[i]=='*')
      return 0;   /* match! */
   else
      return 1;         /* no match */
}


/*------------------------------------------------------------*/
/* convert residue label to number to speed comparisons */

int resnum(const char* res)
{
   register int ans; /* This requires a 32bit int, of course. */


   if (*res=='*')
      ans=0;
   else
      {
      ans = (int)res[0];
      ans <<= 8;
      ans += (int)res[1];
      ans <<= 8;
      ans += (int)res[2];
      }

return(ans);
}


/*------------------------------------------------------------*/
/* DECIDES IF ATOM i IS A POSSIBLE H-BOND PARTICIPANT         */
/*   USES BRUTE FORCE METHOD */

int righttype(int i,int* res,char atm[50][5],int n)
{
   int ans,resn,nn;

   resn = resnum(rlab[attribs[i].res]);
   ans = 0;
   nn = 0;

   while (!ans && nn<n)
      {
      if (resn == res[nn] || !res[nn])
         {    /* RIGHT RESIDUE TYPE? */
         if (! stcp(attribs[i].label,atm[nn]))
            ans = 1;  /* RIGHT ATOM TYPE? */
         }

      nn++;
      }

   return(ans);
}

/*------------------------------------------------------------*/

int read_hb_parms(int  donres[50],
                  char donatm[50][5],
                  int* ndonp,
                  int  accres[50],
                  char accatm[50][5],
                  int* naccp)

{

   FILE *fsp;
   char *not_eof,type;
   char res[5],atm[5],instr[80];
   int ndon,nacc;

   fsp=fopen("hb.parm","r");  /* FIRST TRY LOCAL DIRECTORY */

   if (fsp==0)
      {
      /* NOW TRY SYSTEM DIRECTORY (CHANGE IF YOU WANT) */
      const char* hbSystemDirFileName = "/usr/local/lib/hb.parm";
      fsp=fopen(hbSystemDirFileName,"r");

      if (fsp==0)
         {
         printf("Can't open Hydrogen Bond Parameter File:%s %s\n",hbSystemDirFileName,strerror(errno));
         return 1;
         }
      }

   ndon=0;
   nacc=0;

   do
      {
      not_eof=fgets(instr,79,fsp);
      if (not_eof && instr[0]!='#')
         {
         sscanf(instr,"%s %s %c",res,atm,&type);
         if (type=='D')
            {
            donres[ndon] = resnum(res);
            sprintf(donatm[ndon],"%-4s",atm);
            ndon++;
            }
         if (type=='A')
            {
            accres[nacc] = resnum(res);
            sprintf(accatm[nacc],"%-4s",atm);
            nacc++;
            }
         }
      }
   while (not_eof);

   *ndonp = ndon;
   *naccp = nacc;

   return 0;
}


/*------------------------------------------------------------*/
/* Given partially processed data, locate hydrogen bonds */


void findhb(int nat, /* Number of atoms */

            short* xx, /* Coordinates of atoms */

            int ib[MAXDRAW][3]) /* Bonded list */

{
   int i,j,ii,jj,iii,jjj,k,l,na=0,nd=0;
   int start,end,check;
   int curratom, possatom, ptr, iter;
   int accar[MAXATM],donar[MAXATM];
   int ihbcut2,dist2;
   int dx,dy,dz;
   int nhb=0,nbnd;
   int nabors[MAXATM][MAXBPA],nnabors[MAXATM];
   short excl[MAXATM*20],ex[300];
   int ipex[MAXATM];
   int ndon,nacc;
   int donres[50],accres[50];
   char donatm[50][5],accatm[50][5];
   int flag;

   printf("Searching for possible h-bonds.\n");
   nbnd = ib[0][2]+1; /* Number of bonds in the molecules */

   /* A distance within which a H-bond is possible */

   ihbcut2 = (int)(KAPPA * KAPPA * hbcut * hbcut);

   flag = read_hb_parms(donres,donatm,&ndon,accres,accatm,&nacc);
   if (flag)
      return; /* No parameter file was found - so we are done. */


   for (i=1;i<=nat;i++)
      nnabors[i]=0;

   for (i=1;i<nbnd;i++)
      {
      ii=(int)(ib[i][0]/3+1); /* Convert to atom #ii */

      jj=(int)(ib[i][1]/3+1); /* Atom #jj */

      nabors[ii][nnabors[ii]++]=jj;
      nabors[jj][nnabors[jj]++]=ii;
      }

   ptr=0;

   /* GENERATE EXCLUDED LIST */

   for (i=1;i<=nat;i++)
      {
      /*
      printf("EXCLUDE LIST FOR ATOM %3s-%3d: \n",attribs[i].label,attribs[i].res);
      */

      ipex[i]=ptr;
      start=0;
      end=1;
      ex[0]=(short)i; /* POtential BUG in type conversion */


      for (iter=0;iter<3;iter++)
         {
         int new_position=end;

         /* FOR EACH ATOM ALREADY IN EXCLUDED LIST */
         for (j=start;j<end;j++)
            {
            curratom=ex[j];

             for (k=0;k<nnabors[curratom];k++)
               {
               possatom = nabors[curratom][k];

               /* CHECK POSSATOM TO MAKE SURE IT IS NOT ALREADY IN LIST */
               for (l=0; l<new_position && ex[l]!=possatom; l++)
                  {
                  }

               if (l==new_position)
                  {
                  /* ADD TO END OF LIST */
                  ex[new_position++]=(short)possatom; /* Potential BUG in type conversion */

                  }
               }
            }
         start=end;
         end=new_position;
         }

   /* TRANSFER TO MAIN LIST */
   for (j=0;j<end;j++)
      {
      /* Possible BUG - could overrun array */

      excl[ptr++] = ex[j];
      }
   }

   ipex[i]=ptr;

   /* MAKE LISTS OF DONORS & ACCEPTORS */
   for (i=1;i<=nat;i++)
      {
      if (righttype(i,donres,donatm,ndon)) donar[nd++]=i*3-3;
      if (righttype(i,accres,accatm,nacc)) accar[na++]=i*3-3;
      }

   printf("%d donors, %d acceptors found.\n",nd,na);

   for (i=0;i<nd;i++)
      {
      ii=donar[i];
      iii=ii/3+1;
      for (j=0;j<na;j++)
         {
         jj = accar[j];
         jjj=jj/3+1;
         dx = (int)(xx[ii  ]-xx[jj  ]);
         dy = (int)(xx[ii+1]-xx[jj+1]);
         dz = (int)(xx[ii+2]-xx[jj+2]);
         dist2 = dx*dx+dy*dy+dz*dz;
         if (dist2 < ihbcut2)
            {
       /* MAKE SURE IT ISN'T IN EXCLUDED LIST */
             for (check=ipex[iii];
                    check<ipex[iii+1] &&
                    excl[check]!=jjj;
                       check++)
                  {
                  }

             if (check==ipex[iii+1])
               {
               ib[nhb+nbnd][0]=ii;
               ib[nhb+nbnd][1]=jj;
               ib[nhb+nbnd][2]=(-1);
               nhb++;
               if (nhb+nbnd==MAXDRAW)
                  {
                  i=nd;
                  j=na;
                  printf("Too many H-bonds found...\n");
                  }
               printf("HB: (%3d %3s) (%3d %3s)\n",attribs[iii].res,
                      attribs[iii].label, attribs[jjj].res, attribs[jjj].label);
               }
            }
         }
      }

   ib[0][2] += nhb;
   printf("Found %d H-bonds.\n",nhb);
}


/* Read in pdb file, generate connectivities, return number of atoms */

int gencon(float ave[3],
           int ib[MAXDRAW][3])
{
   int ires=1; /* # of Residues read so far... */

   int iprevres=0; /* # of residue on previous iteration */

   int iatom=1;    /* Index of "atom #" */

/* int header[5]; */

   int nres,rs; /* atomSerialNumber; */

   int i,j;
   int bonds = 1;


/* int typ,i,ii,jj,j,m, */

   int nat=0; /* # of atoms */

/*    found */

   short xx[MAXATM*3];
/* int jb[MAXDRAW][2],nb[MAXATM]; */


   float rad[MAXATM]; /* Radius of each atom */

/* register int temp,sdist;  BUG - what is the type of sdist?????  See below */

/* float dist,xf[3]; */

   FILE *pdbFILE;
   char *not_eof,instr[100]; /* ,aname[5]; */

   int key[MAXATM*2],lookseq[MAXATM*2];

   int nb[MAXATM]; /* Number of bonds (connectivity) from each atom */


   printf("Reading in PDB file.\n");

   pdbFILE=fopen_or_exit_if_error(pdbFileName,"r");

   do
      {
      int reclen;

      not_eof=fgets(instr,sizeof(instr),pdbFILE);

      if (not_eof)
         reclen = strlen(instr);

      /* If we got a line of data, and the first two characters */

      /* are "AT" then the record ID is an "Atom" */

      if (not_eof &&
          reclen >= 54 && instr[0] == 'A' && instr[1] == 'T')
         {
         /* atomSerialNumber = ParseInteger(instr+6,5); */

         /* BUG - in original code this was "atom" and not used */

         /* anywhere else below */


         /* Parse fields from the pdb record in instr */


         macroParseStringAndPad(attribs[iatom].label,instr+12,4);

         rs = ParseInteger(instr+22,4);

         macroParseString(rlab[rs],instr+17,3);

         attribs[iatom].res=(unsigned short)rs;
         attribs[iatom].color=1;       /* default white */
         attribs[iatom].showlabl=0;    /* default no labels */

         /* Get X,Y, Z coordinates in float format */

         /* and move to 2 byte "shorts" for speed? */

         /* Iterate over x,y, and z positions */

         /* Move instr+30, instr+38, instr+46 */

         /* into xx[iatom*3-3], xx[iatom*3-2], xx[iatom*3-1] */

         /* after multiplying my our scaling factor KAPPA */

         for (i=0;i<3;i++)
            xx[iatom*3-3+i] = (short)(ParseFloat(instr+30+i*8,8)*KAPPA);


         if (iatom>nat)
            nat=iatom;

         /* If residueSequenceNumber is not same as previous, then */

         /* increment cound */

         if(rs!=iprevres)
            ires++;

         nres=ires-1;
         allocate_global_rlab(nres,"gencon");

         /* Set previous residue # for next record read */

         iprevres=rs;

         /* Given an atom, we need to be able to answer the question */

         /* What molecule is that atom a part of? */

         mol[iatom]=nmol; /* Assign the atom to a molecule. */

         attribs[iatom].moln=nmol;


         iatom++;

         if (iatom>MAXATM) {
            printf ("Exceeded MAXATM.\n");
            exit(2);
         }
      }

      /* If we see "TER", then that is the end of this molecule. */

      /* There could easily be other molecules in the file. */

      /* So, we increment our molecule count and press on... */


      if (not_eof && instr[0]=='T' && instr[1]=='E' && instr[2]=='R')
         nmol++;
      }
   while (not_eof); /* End "do" loop */


/*  add info about residue labels at end of attrib array  */

   for (i=1;i<=nres;i++)
      {
      for (j=0;j<3;j++)
         attribs[nat+i].label[j]=rlab[i][j];
      attribs[nat+i].label[3]=0;
      attribs[nat+i].color=(-1);
      attribs[nat+i].res=(unsigned short)i;
      }


/*  move center of mass to origin */

   {
   int xc=0,yc=0,zc=0;
   for (i=0;i<nat;i++)
      {
      xc+=xx[i*3];
      yc+=xx[i*3+1];
      zc+=xx[i*3+2];
      }
   xc=xc/nat; yc=yc/nat; zc=zc/nat;

   /* BUG: Lots of type conversion here... Are they OK? */

   ave[0]=(float)(((float)xc)/KAPPA);
   ave[1]=(float)(((float)yc)/KAPPA);
   ave[2]=(float)(((float)zc)/KAPPA);
   for (i=0;i<nat;i++)
      {
      xx[i*3]-=   (short)xc;
      xx[i*3+1]-= (short)yc;
      xx[i*3+2]-= (short)zc;
      }
   }

/*
dorescolor(nat,rlab);
*/

/*    MAS 2-18-90
   option to read in CONECT files from the pdbfile. */

/*-----------------01/09/01 12:07pm-----------------
 From http://www.rcsb.org/pdb/ we have the following documentation 

 on CONECT records

CONECT

   Overview

   The CONECT records specify connectivity between atoms for which
   coordinates are supplied. The connectivity is described
   using the atom serial number as found in the entry.
   CONECT records are mandatory for HET groups (excluding water) and for
   other bonds not specified in the standard residue connectivity table
   which involve atoms in standard residues (see Appendix
   4 for the list of standard residues).
   These records are generated by the PDB.

   Record Format

COLUMNS         DATA TYPE        FIELD           DEFINITION
---------------------------------------------------------------------------------
 1 -  6         Record name      "CONECT"

 7 - 11         Integer          serial          Atom serial number

12 - 16         Integer          serial          Serial number of bonded atom

17 - 21         Integer          serial          Serial number of bonded atom

22 - 26         Integer          serial          Serial number of bonded atom

27 - 31         Integer          serial          Serial number of bonded atom

32 - 36         Integer          serial          Serial number of hydrogen bonded
                                                 atom
37 - 41         Integer          serial          Serial number of hydrogen bonded
                                                 atom
37 - 41         Integer          serial          Serial number of hydrogen bonded
                                                 atom
42 - 46         Integer          serial          Serial number of salt bridged
                                                 atom
47 - 51         Integer          serial          Serial number of hydrogen bonded
                                                 atom
52 - 56         Integer          serial          Serial number of hydrogen bonded
                                                 atom
57 - 61         Integer          serial          Serial number of salt bridged
--------------------------------------------------*/


   if (conectflag)
      {
      int imas,atom1,atom2;

/*    char *bufptr,linetype[10]; */


      printf("reading in CONECT records\n");

         fclose_or_exit_if_error(pdbFileName,pdbFILE);

      pdbFILE=fopen_or_exit_if_error(pdbFileName,"r");

      imas=0;
      while (imas<MAXATM)
         {
         int instr_len;
         char* bufptr;

         instr[0] = 0;
         if(NULL==fgets(instr,sizeof(instr),pdbFILE))
            break; /*EOF*/

         instr_len = strlen(instr);
         /* If this input line is not a CONECT record read next input */

         if ( memcmp(instr,"CONECT",6) || (instr_len < 11))
            {
            imas++;
            continue; /* Goes back to top of while loop */

            }

         /* a CONECT record, extract the bond information */
         atom1 = ParseInteger(instr+6,5);

         /* Make sure that instr is all zeroes to right of input data */

         memset(instr+instr_len,0,sizeof(instr)-instr_len);

         bufptr = instr+11;
         while (*bufptr)
            { /* get successive atom2s */
            atom2 = ParseInteger(bufptr,5);

            /* following is necessary because display does
               not expect the redundancy present in the
               CONECT records; that is, if 44 is bonded to
               45, it is only stated once here, whereas,
               the following CONECT records are expected:
               CONECT 44 45
               CONECT 45 44
            */
            if(atom1<atom2)
               {
               ib[bonds][0] = atom1; /* We have bond from atom1 to atom2 */

               ib[bonds][1] = atom2;
               nb[atom1]++;          /* Increment # of bonds out of atom1 */

               nb[atom2]++;          /* and out of atom2 */

               bonds++;
               }

            if (nb[atom1]>MAXBPA)
                printf("WARNING: Atom %d has %d bonds.\n",atom1,nb[atom1]);

            if (nb[atom2]>MAXBPA)
                printf("WARNING: Atom %d has %d bonds.\n",atom2,nb[atom2]);

            /* Point at next bonded atom */

            bufptr += 5;
            }
         imas++;
         }

      printf("read in %d CONECT bonds\n",bonds-1);
      } /* End if (conectflag) */

   else
      {
      int i;
      printf("Generating connectivity.  Please wait.\n");
      for (i=1;i<=nat;i++)
         {
         const char* atomNamePtr = attribs[i].label;
         nb[i]=0;

         /* Not sure about this????? (BUG??) */

         /* Apparently, if we have a Hydrogen or Dueterium, it could */

         /* be preceeded by a spurious??? digit */

         if (isdigit(atomNamePtr[0]) &&
               (atomNamePtr[1] == 'H' || atomNamePtr[1] == 'D'))
            {
            atomNamePtr++;
            }

#ifdef OLDSTUFF
      sprintf(aname,"%s",attribs[i].label);
      if (isdigit(aname[0]) && (aname[1]=='H' || aname[1]=='D'))
         sprintf(aname,"%s",&attribs[i].label[1]);
#endif

         /* The labels that come in from the .pdb file will have the */

         /* element symbol, then followed by, often, A,B,C */

         /* (alpha, beta, gamma) etc. */


         /* The code below effectively compares the element symbol */

         /* and ignores the ABC business..... */


         /* BUG??: */

         /* Care should be excercised in using this code */

         /* It could easily fail if new elements are added */

         /* which are ambiguous with the Alpha, Betta, etc lettering */

         {
         int found = (-1);
         int typ = 0;

         while (typ<nAtomTypes && found<0)
            {
            int j=0;
            const char* typeCandidate = atom[typ].symbol;  /* See table above */


            while (j<4 && (atomNamePtr[j] == typeCandidate[j]))
               j++;

            if (! typeCandidate[j])
               found=typ;

            typ++;
            }

         if (found<0)
            {
            printf("Unknown atom %s (#%d).\n",atomNamePtr,i);
            found=0;
            }

         rad[i] = atom[found].radius;
         }
      /*rad[i]=1.00;*/
         }
      } /* End "else" part of if (conectflag) */


   for (i=1;i<=nat;i++)
      key[i]=xx[i*3-3];

   /* BUG??: Why don't we sort by inter-atom distance?  That is quite doable. */

   printf("    -sorting by x coordinate.\n");
   sort(key,lookseq,nat);
   printf("    -checking inter-atomic distances.\n");

   /* Look at each atom, to numberOfAtoms-1 */

   for (i=1;i<nat;i++)
      {
      int ii=lookseq[i]; /* Now, look at "adjacent" atoms by x coordinate. */

      int jj;
      float dist;

      /* Check distance to jth atom.  Only consider combinations */

      /* not yet considered.... */

      /* so start at i+1... and interate down to the last atom */

      for (j=i+1;j<=nat;j++)
         {
         int m;

         jj=lookseq[j];
#ifdef SUSPECT_BUG
         /*float*/ short short_dist=0; /* It was declared int above - but inited as a float??? */

#endif

         /* Compute Euclidean distance (Sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2) */

         /* between atoms */

         /* BUG - I am very skeptical that we can square integers safely like this */

         /* hmmm.....  I'm changing this to just use "dist" instead of short_dist */

         dist = 0.0; /* Euclidean distance between atoms */


         for (m=0;m<3;m++)
            {
            float temp=(float)(xx[ii*3+m-3]-xx[jj*3+m-3]);
            temp *= temp;
            dist+=temp;

#ifdef OLD_STUFF
            int temp=xx[ii*3+m-3]-xx[jj*3+m-3];
            temp*=temp;
            short_dist+=temp;
#endif
            }

#ifdef OLD_STUFF
         dist=sqrt((float)sdist)/KAPPA;
#endif
         dist=(float)(sqrt(dist)/KAPPA);

         if (dist<0.01)
            {
            printf("WARNING: Overlapping atoms,");
            printf("probably caused by gap in sequence.\n");
            }

         /* If the distance between atoms is less than the */

         /* VDW radius??? (BUG - comment) */

         /* AND the atoms are part of the same molecule */

         /* THEN, they are indeed bonded! */

         if (dist<(rad[ii]+rad[jj]) && (mol[ii]==mol[jj]))
            {
            /* */

            ib[bonds][0]=(ii<jj ? ii : jj);
            ib[bonds][1]=(ii<jj ? jj : ii);
            nb[ii]++; /* Increment ## of bonds outof atom ii */

            nb[jj]++; /* Increment # of bonds out of atom jj */

            bonds++;
            if (nb[ii]>MAXBPA)
               printf("WARNING: Atom %d has %d bonds.\n",ii,nb[ii]);
            if (nb[jj]>MAXBPA)
               printf("WARNING: Atom %d has %d bonds.\n",jj,nb[jj]);
            }

         /* If our X coordinates are now too far apart, we can */

         /* forget about other atoms in the "lookseq" */

         /* I.e. by considering X coordinate, we have a minimum */

         /* distance, even if Y1=Y2 and Z1=Z2 */

         /* BUG issue:  Why don't we skip ahead as we approach */

         /* good atoms???  Would reduce runtime. */

         if ((float)(xx[jj*3-3]-xx[ii*3-3]) > MAXCHK*KAPPA)
            break;
         }

         /* Every 100th atom, output some info on the screen */

      if (!(i%100))
         printf("      >atom %d\n",i);
      }

   /* Same idea as above - but this time output the last atom. */

   /* (The last atom, of course, is not considered above) */

   printf("      >atom %d\n",nat);


   printf("\n\n%d bonds found.\n",bonds-1);
   printf("\nSorting bonds.\n");

   /* Preparation for sorting... set keys to interatom bonds */

   for (i=1;i<bonds;i++)
      key[i]=ib[i][1]; /* key[i] = # of second atom in this bonded pair */


   sort(key,lookseq,bonds-1);

   {
   int jb[MAXDRAW][2];
   /* Create a temp copy of ib in jb.  Only, jb has been sorted */

   /* by the # of the second atom in each bonded pair */

   for (i=1;i<bonds;i++)
      {
      jb[i][0]=ib[lookseq[i]][0];
      jb[i][1]=ib[lookseq[i]][1];
      }

   /* Now, we'll sort jb by the # of the first atom in each bonded pair */

   for (i=1;i<bonds;i++)
      key[i]=jb[i][0];

   sort(key,lookseq,bonds-1);

   /* Now, we recrate a sorted ib.  However, note that we index */

   /* from 0, and multiply each entry by 3.  Hmmm (BUG? comment?) */

   /* Perhaps this is something we need to do to satisfy output format?? */

   for (i=1;i<bonds;i++)
      {
      ib[i][0]=(jb[lookseq[i]][0]-1)*3;
      ib[i][1]=(jb[lookseq[i]][1]-1)*3;
      ib[i][2]=1; /* BUG comment -what is a 1 */

      }
   }

   /* Finally, lets store some "header" information in the array */

   ib[0][0]=nat;     /* Number of atoms */

   ib[0][1]=nres;    /* Number of residues */

   ib[0][2]=bonds-1; /* Number of bonds */


   WriteAttributeFile(nat+nres+1);

   /* FIND & ADD H-BONDS */
   if ((coordFileNames[0] == 0) && hbcut>0.0)
      {
      findhb(nat,xx,ib);
      bonds=ib[0][2]+1;
      }

   /* WRITE BONDFILE HERE IF NO COORDFILE GIVEN */
   if (coordFileNames[0] == 0)
      {
      int bytesperframe;

      WriteBondFile(bonds+1,ib);

/* And, if coordfile not specified, write 1 frame from the coords in the PDB file */

      bytesperframe=3*nat*sizeof(short);
      WriteBinaryCoordinateFileWithAllFrames(
               xx,
               nat,
               bytesperframe,
               1);
      }


   return nat;    /* return number of atoms */
}





/*------------------------------------------------------------*/

/* BUG what is this? */

void surfdude(int nat,float ave[3])
{
   short surface[MAXSRF*3],lo,ho;
   int atm=0,indx;
   char *not_eof,dum[82],instr[82];
   FILE *fsp;
   short x,y,z;
   float xf,yf,zf;
   char  card[6],alab[6],reslab[6];
   int   res; /* ,bytes; */


   printf("Doing surface file.\n");
   fsp=fopen_or_exit_if_error(insurf,"r");

   indx=sizeof(short)*nat+4;  /* position index to beginning of coord section*/
                           /*   of the array                              */

   do
      {
      not_eof=fgets(instr,80,fsp);
      if (not_eof)
         {
         if (instr[0]==32)
            {        /* read in another surface point */
            sscanf(instr," %f %f %f",&xf,&yf,&zf);
            x=(short)((xf-ave[0])*KAPPA);
            y=(short)((yf-ave[1])*KAPPA);
            z=(short)((zf-ave[2])*KAPPA);
            surface[indx++]=x;
            surface[indx++]=y;
            surface[indx++]=z;
            }
         else
            {        /* next atom */
            not_eof=fgets(dum,80,fsp);
            sscanf(instr,"%s %s %s %d",card,alab,reslab,&res);
            while(atm<nat &&
               (res!=attribs[atm].res || stcp(attribs[atm].label,alab))
               )
               {
               lo=(short)(indx%10000);     /* low order part of address */
               ho=(short)(indx/10000);     /* high order part of address */
               surface[2*atm]  =lo;
               surface[2*atm+1]=ho;
               atm++;
               }

            if (atm==nat)
               {
               printf("ERROR: Atom %s in residue %d found ",alab,res);
               printf("in surface file but not in PDB file!\n");
               exit(EXIT_FAILURE);
               }

            lo=(short)(indx%10000);     /* low order part of address */
            ho=(short)(indx/10000);     /* high order part of address */
            surface[2*atm]  =lo;
            surface[2*atm+1]=ho;
            atm++;
            }
         } /* End if (not_eof) */

      }
   while (not_eof);

/* OLD CODE:  bytes=sizeof(short)*indx+2; */


  WriteSurfaceFile(indx+1,surface);
}

int testfloat108(int n,
                 float arr[],
                 FILE* fsp)
{
/* float flt; */

   register int i;
   char *ptr,instr[82];

   /* BUG???  Looks like n could be passed in as 0 */

   /* In this case, i would increment below until eof... */

   /* Hmmm... What _is_ n? */

   if (! n)
      fgets(instr,sizeof(instr),fsp);

/* Reads up to 10 reals per line in every 8th column (i.e. 10F8.3) */
   for (i=0;i<n;i++)
      {
      if (i%10==0) /* Read input when we need a new line */

         {
         memset(instr,0,sizeof(instr));
         ptr=fgets(instr,82,fsp);
         if (!ptr)
            {
            printf("Insufficient data in coord file\n");
            return (-1);
            }
         }

      arr[i] = ParseFloat(ptr,8);

/*      sscanf(ptr,"%10f",&flt);  Bug- why was this a 10 before??? */

/*      Seems like we should be working with 8 cols at a time! */

/*        BUG??? */


/*      arr[i]=flt;  Replaced with above */


      ptr+=8;
      }

   /* BUG - need comment, example here... */

   /* LOOK FOR CONSTANT PRESSURE BOX DIMENSIONS */
   ptr=fgets(instr,82,fsp);
   if(ptr)
      {
      i=0;
      while (i<81 && instr[i])
         i++;

      if (i<27 && n>3 && ! conpress)
         {
         conpress=1;
         printf("Constant pressure box detected...\n");
         }
      }

   return 0;
}

/*------------------------------------------------------------*/


// Let's read a line of text from amber parameter file
// Skip over any comments in amber7 format which
// are lines starting with a '%'
static int amberparms_has_FORMAT;
static int amberparms_fieldcount;
static char amberparms_formatchar;
static int amberparms_fieldwidth;
static int amberparms_fieldright;

static char* amberparm_fgets(char instr[82],FILE* fsp,const char* callingFunction)
{
   char* ptr;
   do
      {
      instr[0] = 0;

      ptr = fgets(instr,82,fsp);

      if (debugFlag)
         {
         printf("%-15.15s",callingFunction);
         if (ptr)
            {
            printf("|%-50.50s",ptr);
            if (ptr[0] == '%')
               printf("*NotParsed*");
            }
         printf("\n");
         fflush(stdout);
         }
      if ((ptr) && 
          (! memcmp(ptr,"%FORMAT(",8)))
         {
         char* fmtstart = instr+8;
         char* fmtend = strchr(fmtstart,')');
         if (fmtend)
            {
            int nargs;
            *fmtend = 0;   
            amberparms_has_FORMAT = 1;
            nargs = sscanf(fmtstart,"%d%c%d.%d",
               &amberparms_fieldcount,
               &amberparms_formatchar,
               &amberparms_fieldwidth,
               &amberparms_fieldright);


            if (nargs == 3)
               amberparms_fieldright = 0;
            else if (nargs != 4)
               {
               printf("%%FORMAT error on input line %s(%s)\n",instr,fmtstart);
               exit(EXIT_FAILURE);
               }

            if (debugFlag)
               {
               printf("  New Format: %d%c%d.%d\n",
                        amberparms_fieldcount,
                        amberparms_formatchar,
                        amberparms_fieldwidth,
                        amberparms_fieldright);
               fflush(stdout);
               }

            }
         }
      }
   while ((ptr != 0) && (instr[0] == '%'));

   return ptr;
   
}


/* Returns number of floats read... or 0 if error */

int amberparm_readfloat108(int n,short arr[],float ave[3],FILE *fsp)
{
   float flt;
   register int i;
   char *ptr,instr[82];

   if (!n)
      amberparm_fgets(instr,fsp,"readfloat");

   for (i=0;i<n;i++)
      {
      if (i%10==0)
         {
         ptr = amberparm_fgets(instr,fsp,"readfloat");
         if (!ptr)
            return 0;
         }

      /* BUG: This was a scanf with a %10f....  Needs research. */

      flt = ParseFloat(ptr,8);
      arr[i]=(short)((flt-ave[i%3])*KAPPA);
      ptr+=8;
      }
   return i;
}

/*------------------------------------------------------------*/

void nojump(int nat,short* xx,short cpbox[MAXFRAMES][3])
{
/*  Remove discontinuities--if a residue "jumps" from one side  */
/*  of the box to the opposite side, put it back (it will then  */
/*  be outside the box).  If constant pressure box the box size */
/*  varies with each frame; if constant volume box the box size */
/*  is fixed (=size at frame 0) */

   int fr,co,ico,ncord,ires,prevres,xoff[3];

   ncord = 3*(nat+oneIfInBox);

   printf("Making trajectories continuous.\n");

   for (fr=1;fr<frames;fr++)
      {
      prevres = -1; /* no previous as yet. */

      for (co=0;co<ncord-oneIfInBox*3;co+=3)
         {
         ires = attribs[co/3+1].res;
         if (ires != prevres) /* Are we at a new residue? */

            {
            for (ico=co;ico<co+3;ico++)
               {
               /* BUG?? - is this comment OK? */

               int xx_this_frame = xx[fr*ncord+ico];
               int xx_0_frame = xx[ico];
               int distance = xx_this_frame - xx_0_frame;

               /* If we are outside of the constant pressure box */

               /* Then fixup xoff (BUG - is this right?) */

               if (distance > cpbox[fr*conpress][ico%3]/2)
                  xoff[ico%3] = -cpbox[fr*conpress][ico%3];
               else if (distance < -cpbox[fr*conpress][ico%3]/2)
                  xoff[ico%3] =  cpbox[fr*conpress][ico%3];
               else
                  xoff[ico%3] = 0;
               }
            prevres = ires;
            }

         /* Now apply xoff to our xx array of coordinates */

         /* BUG - comments and (short) below... help! */

         for (ico=co;ico<co+3;ico++)
            xx[fr*ncord+ico] += (short)xoff[ico%3];
         }
      }

}




/*------------------------------------------------------------*/

void dotraj(int nat,float ave[3],int ib[MAXDRAW][3])
{
   short* xx;
   /* short cpbox[MAXFRAMES][3],*xx; */

   /* float test[MAXATM*3]; */

   int nread,i,ii,j,last; /* ,val; */
   int nCoordFile;

   int bytesperframe,bonds,nbytes;
   // float x,y,z; // ,minx,maxx,miny,maxy,minz,maxz;
   char titl[82];
   FILE *fsp;
/* extern char *coord; */

   float zero[3];

   short3* cpbox = (short3*)malloc_or_exit_if_error(MAXFRAMES*3*sizeof(short),"dotraj - cpbox");
   float* test = (float*)malloc_or_exit_if_error(MAXATM*3*sizeof(float),"dotraj - test");
    
   for (i=0;i<3;i++)
      zero[i]=0.0;

   fsp=fopen_or_exit_if_error(coordFileNames[0],"r");

   // Skip over first blank line in AMBER crd file
   fgets(titl,82,fsp);

   i = testfloat108(nat*3,test,fsp);
   fclose_or_exit_if_error(coordFileNames[0],fsp);

   if(i<0)
      exit(EXIT_FAILURE);

/* This block finds geometric center */
/*   of the system and centers it at */
/*   the origin.                     */
/*     (using first frame only)      */

/*  move center of mass to origin */

   {
   float xc=0.0,yc=0.0,zc=0.0;
   for (i=0;i<nat;i++)
      {
      xc+=test[i*3];
      yc+=test[i*3+1];
      zc+=test[i*3+2];
      }
   xc=xc/nat; yc=yc/nat; zc=zc/nat;

   /* BUG: Lots of type conversion here... Are they OK? */

   ave[0]=(float)(((float)xc));
   ave[1]=(float)(((float)yc));
   ave[2]=(float)(((float)zc));
   }
#ifdef LATER

   minx=maxx=test[0];
   miny=maxy=test[1];
   minz=maxz=test[2];

   for (i=0;i<nat;i++)
      {
      x=test[i*3];
      y=test[i*3+1];
      z=test[i*3+2];

      if (x>maxx) maxx=x;
      if (x<minx) minx=x;
      if (y>maxy) maxy=y;
      if (y<miny) miny=y;
      if (z>maxz) maxz=z;
      if (z<minz) minz=z;
      }
      
   printf("\nDimension ranges in first frame x(%f,%f) y(%f,%f)z(%f,%f)\n",minx,maxx,miny,maxy,minz,maxz);
   printf("\nDimension multiplied by KAPPA x(%d,%d) y(%d,%d)z(%d,%d)\n",(int)(minx*KAPPA),(int)(maxx*KAPPA),(int)(miny*KAPPA),(int)(maxy*KAPPA),(int)(minz*KAPPA),(int)(maxz*KAPPA));

   ave[0]=(float)((minx+maxx)/2.0f);
   ave[1]=(float)((miny+maxy)/2.0f);
   ave[2]=(float)((minz+maxz)/2.0f);
#endif

   if (oneIfInBox && !conpress) // If period box and NOT constant pressure
      {
      for (i=0;i<3;i++)
         cpbox[0][i]=(short)(fbox[i]*KAPPA);

      printf(" CV Box dimensions: %d  %d  %d\n",
                                 cpbox[0][0],
                                       cpbox[0][1],
                                          cpbox[0][2]);
      }

   if (frames==9999)
     printf("Loading all frames...\n");
   else if(frinc>1)
            printf("Loading frames %d to %d, interval %d frames...\n",
                                 fr1,fr1+frames-1,  frinc);
         else
            printf("Loading frames %d to %d...\n",fr1,fr1+frames-1);


   /* Allocate memory for coords */
   if(frames==9999)
     nbytes=2*3*(nat+oneIfInBox)*101;
   else
     nbytes=2*3*(frames/frinc+1)*(nat+oneIfInBox);

   printf("Allocating %d bytes of memory\n",nbytes);
   xx  = (short*)malloc_or_exit_if_error(nbytes,"dotraj - xx");
   if(xx==0)
      {
      fprintf(stderr,"\nERROR:  Unable to allocate sufficient dynamic memory.\n");
      exit(EXIT_FAILURE);
      }

   if (fr1>1)
      printf("Skipping %d frames.\n",fr1-1);

   nCoordFile = 0;
   /* reset to top of file */
   printf("\nOpening coordinate file %s\n",coordFileNames[nCoordFile]);

   fsp=fopen_or_exit_if_error(coordFileNames[nCoordFile],"r");
   fgets(titl,82,fsp);                     /* discard title line   */

   for (i=0;i<fr1;i++)
      {               /* skip fr1 frames    */
      nread=amberparm_readfloat108(nat*3,xx,ave,fsp);

      if (nread<nat*3)
         {
         if (nCoordFile < (nCoordFiles-1))
            {
            fclose_or_exit_if_error(coordFileNames[nCoordFile],fsp);
            nCoordFile++;
            printf("\nOpening coordinate file %s\n",coordFileNames[nCoordFile]);
            fsp=fopen_or_exit_if_error(coordFileNames[nCoordFile],"r");
            // Skip over first blank line in AMBER crd file
            fgets(titl,82,fsp);
            nread=amberparm_readfloat108(nat*3,xx,ave,fsp);
            }
   
         if (nread < nat*3)
            {
            fprintf(stderr,"End of file in %s encountered after %d frames.\n\n",coordFileNames[nCoordFile],i);
            exit(EXIT_FAILURE);
            }
         }

      /* BUG::::  Check this line.  It was 3,cpbox,zero,fsp */

      /* before - but I think it should be cpbox[subscript] */

      /* Actually, since this is a skip, it probably does not matter */

      /* So, I'll just make it cpbox[0] */

      if (oneIfInBox && conpress)
          amberparm_readfloat108(3,cpbox[0],zero,fsp);
      }

   for (i=0;i<frames;i++)
      {               /* read in frames       */
      ii=(i+frinc-1)/frinc;
      if (i%frinc==0)
         {
         last=ii;
         }

      nread=amberparm_readfloat108(nat*3,&xx[ii*(nat+oneIfInBox)*3],ave,fsp);
#ifdef LAYI_TESTING
      {
      int i = ii*(nat+oneIfInBox)*3;
      printf("First atom: x=%f y=%f z=%f\n",ave[0] + (float)xx[i]/KAPPA,ave[1] + (float)xx[i+1]/KAPPA,ave[2] + (float)xx[i+2]/KAPPA);
      }
#endif
      if (nread<nat*3)
         {
         if (nCoordFile < (nCoordFiles-1))
            {
            fclose_or_exit_if_error(coordFileNames[nCoordFile],fsp);
            nCoordFile++;
            printf("\nOpening coordinate file %s\n",coordFileNames[nCoordFile]);
            fsp=fopen_or_exit_if_error(coordFileNames[nCoordFile],"r");
            // Skip over first blank line in AMBER crd file
            fgets(titl,82,fsp);
            nread=amberparm_readfloat108(nat*3,&xx[ii*(nat+oneIfInBox)*3],ave,fsp);
            }

         if (nread < nat*3)
            {
            frames=i;
            goto zz;
            }
         }

      // Let the user see we are making progress.
      printf(" %d\r",i+fr1);
      fflush(stdout);

      if (oneIfInBox && conpress)
         {                /* dimensions of CP box */
         amberparm_readfloat108(3,cpbox[ii],zero,fsp);   /* get from coord file  */
         for (j=0;j<3;j++)
            xx[(ii+1)*(nat+oneIfInBox)*3-3+j]=cpbox[ii][j];
         }

      /* BUG - what is a CV box... */

      if (oneIfInBox && !conpress)                 /* dimensions of CV box */
      for (j=0;j<3;j++)
         xx[(ii+1)*(nat+oneIfInBox)*3-3+j]=cpbox[0][j];

      if (ii>=100 && frames==9999)
         {
         printf("Reading terminated after frame %d; additional frames may exist.\n",i);
         printf("Use -f option to read all frames.\n\n");
         frames=100;
         break;
         }
      } // End for loop over frames
   frames = last+1;
zz:
   printf(" %d frames read.  \n\n",frames);
   
   fclose_or_exit_if_error(coordFileNames[nCoordFile],fsp);
   

   /* LET SOLVENT FLOAT OUT OF BOX (MAKE TRAJECTORIES CONTINUOUS) */
   if (oneIfInBox && !jumpok)
      nojump(nat,xx,cpbox);

   /* USE FIRST FRAME TO GET POSSIBLE H-BONDS */
   if (hbcut>0.0) findhb(nat,xx,ib);

   /* NOW WRITE OUT BONDS */
   bonds=ib[0][2]+1;

   WriteBondFile(bonds+1,ib);

#ifdef OLD_CODE_REPLACED_BY_LINE_ABOVE
   fdes=open(bondfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
   checkwrite(fdes,bondfile);
   printf("Writing bond file.\n");
   bytes=(bonds*3+3)*sizeof(int);
   write(fdes,(char *)ib,(unsigned)bytes);
   close(fdes);
#endif

   bytesperframe=3*(nat+oneIfInBox)*sizeof(short);

   WriteBinaryCoordinateFileWithAllFrames(
      xx,
      nat,
      bytesperframe,
      frames);

   free(test);
   free(cpbox);
}


/*------------------------------------------------------------*/

void amberparm_readlabel(int n,char arr[][5],FILE* fsp)
{
   int i,j;
   char instr[82];

   /* BUG - what happens on EOF? */


   if (!n)
         amberparm_fgets(instr,fsp,"readlabel");

   for (i=0;i<n;i++)
      {
      if (i%20==0)
         {
         do
            {
            amberparm_fgets(instr,fsp,"readlabel");
            }
         while (instr[0]==13 || instr[0]==10 || instr[0]==0);
         }

      arr[i][4]=0;
      for (j=0;j<4;j++)
         {
         arr[i][j]=instr[(i*4)%80+j];
         }
      }
}

/*------------------------------------------------------------*/

void amberparm_readint(int n,int arr[],FILE* fsp)
{
   int i;
   char *ptr,instr[82];

   /* Not at all sure what this is about???? */

   if (!n)
      amberparm_fgets(instr,fsp,"readint");

   if (! amberparms_has_FORMAT) // Amber 6 and earlier
      {
      amberparms_fieldcount = 12;
      amberparms_fieldwidth = 6;
      }

   for (i=0;i<n;i++)
      {
      if (i%amberparms_fieldcount==0)  /* Time for another line? */

         {
         do
            {
            /* BUG - what about EOF. */
            ptr = amberparm_fgets(instr,fsp,"readint");
            }
         while ((ptr != 0) &&
                (instr[0]==13 || instr[0]==10 || instr[0]==0));

         if ( ptr == 0 )
            {
            fprintf(stderr,"EOF reading integers\n");
            exit(EXIT_FAILURE);
            }

#ifdef OLD_CODE
         ptr=instr;
#endif
         }

      arr[i] = ParseInteger(ptr+1,amberparms_fieldwidth-1);
      ptr+=amberparms_fieldwidth;

#ifdef OLD_CODE
      *(ptr+6)=0;
      sscanf(ptr+1,"%d",arr+i);
#endif
      }
}

/*------------------------------------------------------------*/
void amberparm_readshort(int n,short arr[],FILE* fsp)
{
   int i;
   char *ptr,instr[82];

   /* Not at all sure what this is about???? */

   if (!n)
      amberparm_fgets(instr,fsp,"readshort");

   if (! amberparms_has_FORMAT) // Amber 6 and earlier
      {
      amberparms_fieldwidth = 6;
      amberparms_fieldcount = 12;
      }


   for (i=0;i<n;i++)
      {
      if (i%amberparms_fieldcount==0)  /* Time for another line? */

         {
         do
            {
            /* BUG - what about EOF. */

            ptr = amberparm_fgets(instr,fsp,"readshort");
            }
         while ((ptr != 0) &&
                (instr[0]==13 || instr[0]==10 || instr[0]==0));

         if ( ptr == 0 )
            {
            fprintf(stderr,"EOF reading integers\n");
            exit(EXIT_FAILURE);
            }

#ifdef OLD_CODE
         ptr=instr;
#endif
         }

      arr[i] = (short)ParseInteger(ptr+1,amberparms_fieldwidth-1);
      ptr+=amberparms_fieldwidth;

#ifdef OLD_CODE
      *(ptr+6)=0;
      sscanf(ptr+1,"%d",arr+i);
#endif
      }
}
/*------------------------------------------------------------*/

void amberparm_readfloat(int n,float arr[],FILE* fsp)
{
   int i;
   char *ptr,instr[82];

   if (!n) amberparm_fgets(instr,fsp,"readfloat");

   if (! amberparms_has_FORMAT) // Amber 6 and earlier
      {
      amberparms_fieldcount = 5;
      amberparms_fieldwidth = 16;
      }

   for (i=0;i<n;i++)
      {
      if (i%amberparms_fieldcount==0)
         {
         do
            {
            ptr = amberparm_fgets(instr,fsp,"readfloat");
            }
         while ((ptr != 0) &&
                (instr[0]==13 || instr[0]==10 || instr[0]==0));

#ifdef OLD_CODE
         ptr=instr;
#endif
         }
      arr[i] = ParseFloat(ptr+1,amberparms_fieldwidth-1);
      ptr+=amberparms_fieldwidth;
      }
}



/* CWM changed from ib[MAXDRAW*3] to ib[MAXDRAW][3] for consistency */

/* with other uses of ib */

int readparm(int ib[][3],int forceBox) /* ,int *bondp) BUG - not used??? */

{
   int parms[31];
   char *not_eof,instr[90];
   float beta;
   int delbond=0;
   int nat,i,j;
   int bonds; /* BUG - changed from short to int by CWM */

   int ifbox;
   int nres,nbona;
   int nbonh;
   int numbnd,numang,nptra,natyp,ntype,ntypes,nttyp;
   int ntheth,ntheta,nphih,nphia,nnb,nspm,nphb;
/* int fdes,bytes; */

   FILE *fsp;

   typedef char char5[5];

   short* inpt = (short*)malloc_or_exit_if_error(3*MAXDRAW * sizeof(short),"readparm - inpt");
   int* idum = (int*)malloc_or_exit_if_error(MAXATM*sizeof(int),"readparm - idum");
   char5* cdum = (char5*)malloc_or_exit_if_error(MAXATM*5*sizeof(char),"readparm - cdum");
   float* fdum = (float*)malloc_or_exit_if_error(2*MAXATM*sizeof(float),"readparm - fdum");
   memset(parms,0,sizeof(parms));
   memset(idum,0,MAXATM*sizeof(int));

   fsp=fopen_or_exit_if_error(parameterFileName,"r");

   amberparm_readlabel(3,cdum,fsp);
   amberparm_readint(31,parms,fsp);
   nat=parms[0];

   ntypes=parms[1];
   if (ntypes<=0 || ntypes>MAXATP)
      {
      printf("\n Number of atom types = %d; <=0 or exceeds max. \n\n",ntypes);
      exit(EXIT_FAILURE);
      }

   ntype=ntypes*ntypes;
   nres=parms[11];


   allocate_global_rlab(nres,"readparm");
   allocate_global_atomattrib(nres,nat,"readparm"); 


   /*
   printf("\n\nEnter type of bond you want deleted from display (0 for none):");
   scanf("%s",instr);
   sscanf(instr,"%d",&delbond);
   */
   if (delbond)
      {
      /* BUG: This is uncool.... Is the idea to set all bits to oneor what? */

      /*      -1 would be a better choice for signed char on all architectures */

      /*    255 is not a legal "char" value */

      attribs[0].label[0]=(char)0xff;              /* flag for bond removal */
      attribs[0].res=(unsigned short)delbond;
      }
   else
      {
      attribs[0].label[0]=0;
      }


#ifdef OLD_IDEA
   if (nres <=0 || nres > MAXRES)
      {
      printf("\n Number of residues = %d; <=0 or exceeds MAXRES (%d).\n\n",
         nres, MAXRES);
      exit(EXIT_FAILURE);
      }
#endif

   nttyp=ntypes*(ntypes+1)/2;
   nbonh=parms[2];
   nbona=parms[3];
   ntheth=parms[4];
   ntheta=parms[5];
   nphih=parms[6];
   nphia=parms[7];
   nnb=parms[10];
   numbnd=parms[15];
   numang=parms[16];
   nptra=parms[17];
   natyp=parms[18];
   nphb=parms[19];
   ifbox=parms[27]; // Amber 0=no box, 1=hexahedron (parallelpiped), 2=truncated octahedron
   bonds=(int)(nbonh+nbona); /* BUG: Always < 32K???  Hmmm */

   // if (oneIfInBox && getenv("KRISTINA"))
      // oneIfInBox = 0;


   amberparm_readlabel(nat,cdum,fsp);   /*igraph*/

   for (i=1;i<=nat;i++)
      {
      for (j=0;j<4;j++)
         {
         attribs[i].label[j]=cdum[i-1][j];
         }
      attribs[i].label[4]=0;
      attribs[i].moln=0;
      }

   oneIfInBox = 0;
   if (forceBox) // For Layi in Norway
      ifbox = 1;
   if (ifbox == 1) // A hexahedron (original "box")
      {
      sprintf(attribs[nat+1].label,"BOX");   /* DUMMY ATOM */
      attribs[nat+1].moln=999;
      oneIfInBox=1;
      }
   else if (ifbox == 2) // truncated octahedron
      {
      sprintf(attribs[nat+1].label,"OCT");   /* DUMMY ATOM */
      attribs[nat+1].moln=1000;
      oneIfInBox=1;
      }
   else if (ifbox)
      {
      printf("Sorry: Unrecognized ifbox parameter: %d in AMBER topology file\n",ifbox);
      exit(EXIT_FAILURE);
      }

   amberparm_readfloat(nat,fdum,fsp);   /*cg*/
   amberparm_readfloat(nat,fdum,fsp);   /*amass*/
   amberparm_readint(nat,idum,fsp);     /*iac*/
   amberparm_readint(nat,idum,fsp);     /*iblo*/
   amberparm_readint(ntype,idum,fsp);   /*ico*/
   amberparm_readlabel(nres,cdum,fsp);  /*lbres*/

   for (i=1;i<=nres;i++)
      {
      // Copy residue names to attribs array
      for (j=0;j<3;j++)
         {
         attribs[i+nat+oneIfInBox].label[j]=cdum[i-1][j];
         }

      attribs[i+nat+oneIfInBox].label[3]=0;
      attribs[i+nat+oneIfInBox].color=(-1); /* neg. color signifies this is residue */
      attribs[i+nat+oneIfInBox].res= (unsigned short)i; /* BUG: sizeof(i) > sizeof(.res) */

      }

   amberparm_readint(nres,idum,fsp);   /*ipres*/
   idum[nres]=nat+1;
   for (i=0,j=1;i<nres;i++)
      {
      do
         {
         attribs[j++].res=(unsigned short)(i+1); /* Bug... sizoef(i) > sizeof(res) */

         if (strcmp(attribs[i+1+nat+oneIfInBox].label,"WAT")==0)
            {
            attribs[j-1].moln=1;
            }
         attribs[j-1].color=1;
         attribs[j-1].showlabl=0;
         }
      while (j<idum[i+1]);
      }

   amberparm_readfloat(numbnd,fdum,fsp);   /*rk*/
   amberparm_readfloat(numbnd,fdum,fsp);   /*req*/
   amberparm_readfloat(numang,fdum,fsp);   /*tk*/
   amberparm_readfloat(numang,fdum,fsp);   /*teq*/
   amberparm_readfloat(nptra,fdum,fsp);    /*pk*/
   amberparm_readfloat(nptra,fdum,fsp);    /*pn*/
   amberparm_readfloat(nptra,fdum,fsp);    /*phase*/
   amberparm_readfloat(natyp,fdum,fsp);    /*solty*/
   amberparm_readfloat(nttyp,fdum,fsp);    /*cn1*/
   amberparm_readfloat(nttyp,fdum,fsp);    /*cn2*/

   ib[0][0]=nat+oneIfInBox;
   ib[0][1]=nres;
   ib[0][2]=bonds;

   /* BUG - CWM changed subscripting to make format like */

   /* other uses of ib[] elsewhere... */

   amberparm_readint(nbonh*3,ib[1],fsp);      /*ibh,jbh,icbh -- this is what we're after*/
   amberparm_readint(nbona*3,ib[1+nbonh],fsp);  /*iba,jba,icba*/


   /* WRITE OUT THE BOND FILE NOW IF COORDFILE NOT GIVEN
   if (!coordfile) {
      fdes=open(bondfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
      checkwrite(fdes,bondfile);
      bytes=(bonds*3+3)*sizeof(int);
      write(fdes,(char *)ib,(unsigned)bytes);
      close(fdes);
   }
   */

   /* GET CV BOX DIMENSIONS IF NECESSARY */
   if (oneIfInBox && !conpress)
      {
      amberparm_readshort(ntheth*4,inpt,fsp);         /* angles */
      amberparm_readshort(ntheta*4,inpt,fsp);         /* angles */
      amberparm_readshort(nphih*5,inpt,fsp);          /* dihedrals */
      amberparm_readshort(nphia*5,inpt,fsp);          /* dihedrals */
      amberparm_readshort(nnb,inpt,fsp);              /* nonbonds */
      amberparm_readfloat(nphb,fdum,fsp);             /* h-bond parms */
      amberparm_readfloat(nphb,fdum,fsp);             /* h-bond parms */
      amberparm_readfloat(nphb,fdum,fsp);             /* h-bond parms */
      amberparm_readlabel(nat,cdum,fsp);              /* isymbl */
      amberparm_readlabel(nat,cdum,fsp);              /* itree */
      amberparm_readshort(nat,inpt,fsp);              /* join */
      amberparm_readshort(nat,inpt,fsp);              /* irotat */

      /*FINALLY THE BOUNDARY STUFF */
      amberparm_readshort(3,inpt,fsp);                /* iptres,nspm,nspsol*/
      nspm=inpt[1];
      amberparm_readshort(nspm,inpt,fsp);
      not_eof = amberparm_fgets(instr,fsp,"readparm");

      if (not_eof)
         sscanf(instr,"%e %e %e %e",&beta,fbox+0,fbox+1,fbox+2);
      else
         oneIfInBox=0;

      if (fbox[0]<1.0 || fbox[1]<1.0 || fbox[2]<1.0) oneIfInBox=0;

      /* BUG comment what kind of box is this? */

      if (oneIfInBox)
         printf(" %s dimensions: %7.3f %7.3f %7.3f\n",
                                 attribs[nat+1].label,
                                 fbox[0],fbox[1],fbox[2]);
      else
         attribs[nat+1].moln=998; /* BUG what is 998? */

      oneIfInBox=1;
      }

   WriteAttributeFile(nat+nres+1+oneIfInBox);

   /* free blocks in reverse order */

   free(fdum);
   free(cdum);
   free(idum);
   free(inpt);
   return nat;

#ifdef OLD_REPLACED_BY_LINE_ABOVE
   /* WRITE OUT THE ATTRIBUTE FILE */
   fdes=open(attribfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
   checkwrite(fdes,attribfile);
   bytes=(nat+nres+1+oneIfInBox)*sizeof(atomattrib);
   write(fdes,(char *)attribs,(unsigned)bytes);
   close(fdes);
#endif

}





/*------------------------------------------------------------*/
/* University of South Alabama                                */
/* Jeffry Madura and Rodney Hamilton                          */
/* Department of Chemistry                                    */
/* Chem Bldg Room 223                                         */
/* Mobile, AL  36688-0002                                     */
/* Routine: To read Charmm22 formatted PSF files              */

int readpsf(/* float ave[3], */ int ib[MAXDRAW][3])
{
   FILE *psf;
   char *not_eof,instr[82]; /* ,aname[5]; */

/* char *bufptr;  ,linetype[10]; */

   int atom1,atom2;
   int ires=1,iprevres=0,iatom=1;

   int nres,rs,i,j,nat=0,bonds=1;
/* int nres,rs,atom,typ,i,j,nat=0,bonds=1; */

   /* ,ii,jj,j,m,nat=0,found,bonds=1; */

/* int xc=0,yc=0,zc=0,fdes;  ,bytes,bytesperframe; */

   int nlines,natoms,nbonds;
/* short xx[MAXATM*3], */

   short jb[MAXDRAW][2],nb[MAXATM];
/* float rad[MAXATM],dist,xf[3]; */

   int key[MAXATM*2],lookseq[MAXATM*2];

   printf("Reading in PSF file.\n");
   psf=fopen_or_exit_if_error(psfFileName,"r");

/* if (!psf) { printf("%s does not exist.\n",psffile); exit(2); } */


   /* read psf header */
   fgets(instr,82,psf); /* Read Header PSF */
   fgets(instr,82,psf); /* Read Blank Line */
   not_eof=fgets(instr,82,psf); /* Read Number of Title Lines */

   if ( ! not_eof )
      {
      printf("Terminating due to lack of title lines in %s",psfFileName);
      exit(EXIT_FAILURE);
      }


   sscanf(instr,"%d",&nlines);

   /* read the title lines */
   printf("Reading the title data.......\n");
   for (i=0;i<nlines;i++)
      {
      fgets(instr,82,psf);
      }

   fgets(instr,82,psf); /* Read Blank Line */

   /*  get the number of atoms */
   not_eof=fgets(instr,82,psf);
   sscanf(instr,"%d",&natoms);

   if ( ! not_eof )
      {
      printf("Terminating due to lack of input in %s",psfFileName);
      exit(EXIT_FAILURE);
      }

   /* read atom data */
   printf("Reading the data for %d atoms........\n",natoms);

   for (j=0;j<natoms;j++)
      {
      not_eof=fgets(instr,82,psf);
      if ( ! not_eof )
         {
         printf("Terminating due to lack of atoms input in %s",psfFileName);
         exit(EXIT_FAILURE);
         }

#ifdef BUG_IN_OLD_CODE
/* OK - the "atom" variable in the old code was a local integer */

/* But, the old code did nothing with it - so I've pulled it */

/* out here completely */

      sscanf(instr,"%d",&atom);
#endif
      sscanf(instr+24,"%s",attribs[iatom].label);
      i=0;    /* PAD WITH SPACES */
      while (attribs[iatom].label[i] && i<4)
         i++;
      while (i<4)
         attribs[iatom].label[i++]=32;

      attribs[iatom].label[4]=0;
      sscanf(instr+14,"%d",&rs);
      sscanf(instr+19,"%s",rlab[rs]);

      attribs[iatom].res=(unsigned short)rs;
      attribs[iatom].color=1;       /* default white */
      attribs[iatom].showlabl=0;    /* default no labels */

      if(iatom>nat)
          nat=iatom;

      if(rs!=iprevres)
         ires++;

      nres=ires-1;
      allocate_global_rlab(nres,"readpsf");

      if (iprevres>rs)
         { nmol++;}  /* Molecule Changed */

      mol[iatom]=nmol;
      iprevres=rs;
      attribs[iatom].moln=nmol;
      iatom++;
      }

   /*  add info about residue labels at end of attrib array  */
   printf("Reading residue labels......\n");
   for (i=1;i<=nres;i++)
      {
      for (j=0;j<3;j++)
         attribs[nat+i].label[j]=rlab[i][j];

      attribs[nat+i].label[3]=0;
      attribs[nat+i].color=(-1);
      attribs[nat+i].res=(unsigned short)i; /* BUG: nres always short? */

      }

   /*    MAS 2-18-90
      Routine to read Bond Information from PSF file. */
   printf("Reading bond information.....\n");
   fgets(instr,82,psf);
   not_eof=fgets(instr,82,psf);

   if (! not_eof)
         {
         printf("EOF reading bonds header records in %s\n",psfFileName);
         exit(EXIT_FAILURE);
         }

   sscanf(instr,"%d",&nbonds);

   for (i=0;i<nbonds;i++)
      {
      int nFieldsRead = fscanf(psf,"%d %d",&atom1,&atom2);

      if ( nFieldsRead != 2 )
         {
         printf("EOF or input error reading bonds in %s\n",psfFileName);
         exit(EXIT_FAILURE);
         }
      ib[bonds][0] = atom1;
      ib[bonds][1] = atom2;
      nb[atom1]++; nb[atom2]++;
      bonds++;
      }

   printf("%d bonds found.\n",bonds-1);
   printf("Sorting bonds.\n");
   for (i=1;i<bonds;i++)
      key[i]=ib[i][1];

   sort(key,lookseq,bonds-1);

   for (i=1;i<bonds;i++)
      {
      jb[i][0]=(short)ib[lookseq[i]][0]; /* BUG: (short) OK? */

      jb[i][1]=(short)ib[lookseq[i]][1]; /* BUG: (short) OK? */

      }

   for (i=1;i<bonds;i++)
      key[i]=jb[i][0];

   sort(key,lookseq,bonds-1);

   for (i=1;i<bonds;i++)
      {
      ib[i][0]=(jb[lookseq[i]][0]-1)*3;
      ib[i][1]=(jb[lookseq[i]][1]-1)*3;
      ib[i][2]=1;
      }

   ib[0][0]=nat;
   ib[0][1]=nres;
   allocate_global_rlab(nres,"readpsf");
   ib[0][2]=bonds-1;

/* Write Attribute File */
   WriteAttributeFile(nat+nres+1);

#ifdef OLD_REPLACED_BY_LINE_ABOVE
   fdes=open(attribfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
   checkwrite(fdes,attribfile);
   printf("Writing attribute file.\n");
   bytes=(nat+nres+1)*sizeof(atomattrib);
   write(fdes,(char *)attribs,(unsigned)bytes);
   close(fdes);
#endif

   /* FIND & ADD H-BONDS
   if (!coordfile && hbcut>0.0) {
      findhb(nat,xx,ib);
      bonds=ib[0][2]+1;
   }*/

/* Write Bond File */
   WriteBondFile(bonds+1,ib);

#ifdef OLD_CODE_REPLACED_BY_LINE_ABOVE
   fdes=open(bondfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
   checkwrite(fdes,bondfile);
   printf("Writing bond file.\n");
   bytes=(bonds*3+3)*sizeof(int);
   write(fdes,(char *)ib,(unsigned)bytes);
   close(fdes);
#endif

   return nat;
}


/*------------------------------------------------------------*/
/* University of South Alabama                                */
/* Jeffry Madura and Rodney Hamilton                          */
/* Department of Chemistry                                    */
/* Chem Bldg Room 223                                         */
/* Mobile, AL  36688-0002                                     */
/* Routine: To read Charmm22 binary DCD files                 */

void dcdtraj(int natom)
{
   FILE *fsp;
   char title[82];
   short *xx;
   double filehdr;
   float ave[3];
   float *x,*y,*z,*tax,*tay,*taz;
   float xa,ya,za,minx,maxx,miny,maxy,minz,maxz;
   int i,icnt,ii,j;
   int nfreat,bytesperframe;
   int data,nlines,info[20];
   int *freelist;

   /* BUG: The "rb" parameter on fopen may not be */

   /* platform independent.  On DOS/Windows, you need this for binary */

   /* file access. */

   fsp = fopen_or_exit_if_error(coordFileNames[0],"rb");

   /* changed fread to fread_or_exit_if_error() */

   fread_or_exit_if_error(&filehdr,sizeof(double),1,fsp); /* Read Header */
   fread_or_exit_if_error(info,sizeof(int),20,fsp);       /* Read File Info */
   fread_or_exit_if_error(&data,2*sizeof(int),1,fsp);

   fread_or_exit_if_error(&nlines,sizeof(int),1,fsp);  /* Read Number of Title Lines */
   for(i=0;i<nlines;i++)
      {
      fread_or_exit_if_error(title,10*sizeof(double),1,fsp); /* Read Title Line */

      title[79] = '\0';
      }

   fread_or_exit_if_error(&data,2*sizeof(int),1,fsp); /* Read Junk */

   fread_or_exit_if_error(&natom,sizeof(int),1,fsp);  /* Read Number of Atoms */
   fread_or_exit_if_error(&data,sizeof(int),1,fsp);
   nfreat = natom;

   if(info[8]>0)
      {
      nfreat=natom-info[8];

      freelist=(int*)malloc_or_exit_if_error(nfreat*sizeof(int),"dcdtraj - freelist");

      fread_or_exit_if_error(&data,sizeof(int),1,fsp);
      fread_or_exit_if_error(freelist,sizeof(int),nfreat,fsp);
      fread_or_exit_if_error(&data,sizeof(int),1,fsp);
      }

   printf("Number of atoms: %d Number of Free: %d\n",natom,nfreat);

   if (frames==9999)
      { frames=info[0]; }

   printf("Reading %d frames.......\n",frames);
   if(info[8]>0)
      {
      tax = (float*)malloc_or_exit_if_error(nfreat*sizeof(float),"dcdtraj - tax");
      tay = (float*)malloc_or_exit_if_error(nfreat*sizeof(float),"dcdtraj - tay");
      taz = (float*)malloc_or_exit_if_error(nfreat*sizeof(float),"dcdtraj - taz");
      }

   x = (float*)malloc_or_exit_if_error(natom*sizeof(float),"dcdtraj - x");
   y = (float*)malloc_or_exit_if_error(natom*sizeof(float),"dcdtraj - y");
   z = (float*)malloc_or_exit_if_error(natom*sizeof(float),"dcdtraj - z");
   xx = (short*)malloc_or_exit_if_error(natom*3*sizeof(short)*frames,"dcdtraj - xx");

   j=0;                      /* Position in array for the xx pointer */
   for(icnt=0;icnt<frames;icnt++)
      {
      if(info[8]>0 && icnt>0)
         { /* There are fixed atoms. */
         fread_or_exit_if_error(&data,sizeof(int),1,fsp);
         fread_or_exit_if_error(tax,sizeof(float),nfreat,fsp);
         fread_or_exit_if_error(&data,2*sizeof(int),1,fsp);
         fread_or_exit_if_error(tay,sizeof(float),nfreat,fsp);
         fread_or_exit_if_error(&data,2*sizeof(int),1,fsp);
         fread_or_exit_if_error(taz,sizeof(float),nfreat,fsp);
         fread_or_exit_if_error(&data,sizeof(int),1,fsp);

/* Put the Free Atom values into the x,y,z arrays */
         for(ii=0;ii<nfreat;ii++)
            {
            x[freelist[ii]-1] = tax[ii];
            y[freelist[ii]-1] = tay[ii];
            z[freelist[ii]-1] = taz[ii];
            }
         }
      else
         { /* Read in First Frame or Trajectories with no FIXED */
         fread_or_exit_if_error(&data,sizeof(int),1,fsp);
         fread_or_exit_if_error(x,sizeof(float),natom,fsp);
         fread_or_exit_if_error(&data,2*sizeof(int),1,fsp);
         fread_or_exit_if_error(y,sizeof(float),natom,fsp);
         fread_or_exit_if_error(&data,2*sizeof(int),1,fsp);
         fread_or_exit_if_error(z,sizeof(float),natom,fsp);
         fread_or_exit_if_error(&data,sizeof(int),1,fsp);

         if (icnt==0)
            {  /* First Frame Only */
               /* This block finds geometric center */
               /*   of the system and centers it at */
               /*   the origin.                     */
               /*     (using first frame only)      */

            minx= maxx = x[0];
            miny= maxy = y[0];
            minz= maxz = z[0];

            for (ii=0;ii<natom;ii++)
               {
               xa=x[ii];
               ya=y[ii];
               za=z[ii];

               if (xa>maxx) maxx=xa;
               if (xa<minx) minx=xa;
               if (ya>maxy) maxy=ya;
               if (ya<miny) miny=ya;
               if (za>maxz) maxz=za;
               if (za<minz) minz=za;
               }

            /* BUG: More typecasting */

            ave[0]=(float)((minx+maxx)/2.0f);
            ave[1]=(float)((miny+maxy)/2.0f);
            ave[2]=(float)((minz+maxz)/2.0f);
            }

         } /* End of Else */

   /* This routine assigns the x,y,z coordiantes into the xx pointer */
      /* BUG: more typecasting */

      for(ii=0;ii<natom;ii++)
         {
         xx[j++]=(short)((x[ii]-ave[0])*KAPPA);
         xx[j++]=(short)((y[ii]-ave[1])*KAPPA);
         xx[j++]=(short)((z[ii]-ave[2])*KAPPA);
         }
      } /* End of For */

   fclose_or_exit_if_error(coordFileNames[0],fsp);

   /* WRITE BINARY COORDINATE FILE WITH ALL FRAMES */
   bytesperframe=3*(natom)*sizeof(short);
   WriteBinaryCoordinateFileWithAllFrames(
      xx,
      natom,
      bytesperframe,
      frames);
}

/*------------------------------------------------------------*/
/* The following routines support reading Fortran unformatted */
/* data files.                                                */
/*------------------------------------------------------------*/

/* CWM did a fair amount of "clean up" here to deal with */

/* buffer overruns, and make code more C-like.  It now uses */

/* ansi c fread() call instead of read() */


/* Read next Fortran unformatted record.  The byte count is  */
/* at the beginning and end of the record.                   */

/* BUG: Does FORTRAN always output this way???? Or, do we need */

/* a LSB/MSB parameter at runtime??? */


int readuf_core(FILE* fp, char* buff,size_t buff_size)
{
   int i, nb; /* ,i1,i2; */


   long i1,i2; /* CWM changed to long to guaranee 32 bits */


   /* CWM removed "static" here to save static data space... */

   /* static char buff[12*MAXATM+1000]; */

   /* BUG - what about big/little endian idssues here??? */

   nb = fread(&i1,4,1,fp);

   if ((nb <= 0) || (i1 <= 0))
      return 1;

   if ((size_t)i1 < buff_size) /* 12*MAXATM+1000) */

      {
      nb = fread(buff,i1,1,fp);
      if (nb != i1)
         return 2;
      }
   else
      {
      /* BUG - this code makes no sense to me. */

      /* If we read in chunks, do we have a reasonable return value??? */

      /* Why read in lots of 1024??? */

      /* We could dynamically allocate a buffer of needed size....? */


      for (i=0;i<i1/1024;i++)
         fread(buff,1024,1,fp);

      fread(buff,i1%1024,1,fp);

      nb = fread(&i2,4,1,fp);

      /* There is a place in code where we _do_ look for return value of 4! */

      if ((i1 == i2) || (nb < 1))
         return 4; /* File integrity is OK - but too large for buffer */

      else
         return 3; /* File has integrity problems. */

      }

  /* Normal read should end with record size repeated at end of file. */

  nb = fread(&i2,4,1,fp);

  if ((i1 != i2) || (nb < 1))  /* File has fundamental integrity problems */

      return 3;

/*  *ibpp = &buff[0]; */


  return 0; /* Successful read into buffer */

}

static int record_count = 0;

void readuf_error(int error_value)
{
   printf("Error %d reading history file at record %d\n",error_value,record_count);
   exit(EXIT_FAILURE);
}


void readuf(FILE* fp, char* buff,size_t buff_size)
{
   int i;
   record_count++;

   i = readuf_core(fp,buff,buff_size);
   if (i)
      readuf_error(i);
}

/* Get character data from buffer */
/* Dramatically simplified by CWM */

char* getufchar(char *bufptr, char *c, int nc)
{
   memcpy(bufptr,c,nc);
   bufptr+= nc;
   return bufptr;
}

/* Get integer data from buffer */
char* getuflong(char *bufptr, long *ia, int na)
{
   /* Streamlined by CWM to reduce runtime arithmetic */

   /* From looking at old code, clearly, FORTRAN stores data */

   /* with MSB first */


   while ( na > 0 )
      {
      *ia = (((int)(bufptr[0])) << 24) |
            (((int)(bufptr[1])) << 16) |
            (((int)(bufptr[2])) << 8)  |
            (((int)(bufptr[3])));

      ia++;
      na--;
      bufptr += 4;
      }

   return bufptr;

#ifdef OLD_CODE_REPLACED_ABOVE
int i,ii,k;
char *ibp;

  ibp = *ibpp;
  for (i=0; i<na; i++) {
    ia[i]=0;
    for (k=0;k<4;k++) {
      ii = *ibp++;
      ia[i] |= ii<<(3-k)*8;
    }
  }
  *ibpp = ibp;
  return 0;
#endif
}

/* Get float (real*4) data from buffer */
char* getuffloat(char *bufptr, float *a, int na)
{
   assert(sizeof(float) == 4);

   /* CWM notes that the old code below can replaced by */

   /* a call to getuflong() - assuming that sizeof(int) == 4 */

   /* and sizeof(float)==4.  Just cast away and run... */


   return getuflong(bufptr,(long*)a,na);


#ifdef OLD_CODE_REPLACED_ABOVE
int i,ii,k;
char *ibp;
union {
  int ia;   /* BUG: Again, ints are 64 soon! */

  float x;
} u;

  ibp = *ibpp;
  for (i=0; i<na; i++) {
    u.ia=0;
    for (k=0;k<4;k++) {
      ii = *ibp++;
      u.ia |= ii<<(3-k)*8;
    }
    a[i] = u.x;
  }
  *ibpp = ibp;
  return 0;
#endif
}

/* Get double (real*8) data from buffer */
char* getufdouble(char* bufptr, double *d, int na)
{
/* int i,ii,j,k; */

/* char *ibp; */

   union
      {
      int id[2];
      double d;
      } u;

   assert(sizeof(u.d) == 8);
   assert(sizeof(u.id) == 8);

   while ( na > 0 )
      {
      u.id[0] = (((int)(bufptr[0])) << 24) |
                (((int)(bufptr[1])) << 16) |
                (((int)(bufptr[2])) << 8)  |
                (((int)(bufptr[3])));

      u.id[1] = (((int)(bufptr[4])) << 24) |
                (((int)(bufptr[5])) << 16) |
                (((int)(bufptr[6])) << 8)  |
                (((int)(bufptr[7])));

      *d++ = u.d;
      na--;
      bufptr += 8;
      }

   return bufptr;

#ifdef OLD_CODE_REPLACED_ABOVE
  ibp = *ibpp;
  for (i=0; i<na; i++) {
    for (j=0;j<2;j++) {
      u.id[j]=0;
      for (k=0;k<4;k++) {
        ii = *ibp++;
        u.id[j] |= ii<<(3-k)*8;
      }
    }
    d[i] = u.d;
  }
  *ibpp = ibp;
  return 0;
#endif
}

/*------------------------------------------------------------*/
int GetYesOrNo(void)
{
   char ans[8];
   int retval = -1;

   while ( retval == -1 )
      {
      char testChar;
      scanf("%4s",ans);
      testChar = ans[0];

      if ( testChar == 'y' || testChar == 'Y')
         {
         retval = 1;
         break;
         }

      if ( testChar == 'n' || testChar == 'N')
         {
         retval = 0;
         break;
         }

      printf("\nA response of Y or N is required -->");
      }

   return retval;
}


/* Read GROMOS molecular topology file */
int read_gromos_mt(int* nsolv,int* isolv,int ib[MAXDRAW][3])
{
   int i,j,k,n;
   int nsol,nat,ires;
   FILE *filep;
   char line[256];
   long nline=0,ntun;
   long nratt; /* Number of atom types */

   long nraa2; /* Number of residues */

   long nrp,nbty;
   long nbond,nbonh, ibond[MAXDRAW][2];
   typedef char resnametype[4];
   char atpname[MAXATP][4]; /* Text description of atom types */

   resnametype* resname;


   long longArray[MAXDRAW][3];
   float dummy;

   struct atmstuff {
      float winv;
      char label[8];
      long mres;
      long iac;
      } atmstf[MAXATM];

/* char *ibp; */


/* First decide if this is a formatted or unformatted file */

   filep = fopen_or_exit_if_error(gromos_mtFileName,"r");
   if (fgets(line,120,filep)==NULL)
      goto bin;

   nline++;

   if (fgets(line,120,filep)==NULL)
      goto bin;

   nline++;
   if (sscanf(line,"%5d",&ntun)!=1)
      goto bin;

   if (ntun<1 || ntun>3)  goto bin;

   if (ntun>1)
      gxscale=10;
   else
      gxscale=1;

   if (fgets(line,120,filep)==NULL)
      goto bin;

   nline++;
   if (sscanf(line,"%5d",&nratt)!=1) /* # of atom types */

      goto bin;

   if (nratt<1 || nratt>MAXATP)
      goto err1;

/* OK, read remainder of formatted topology file */
/* I.e. - we don't worry about this being binary anymore. */


   for (k=0;k<nratt;k++)
      {
      if (k%16==0)
         nline++;

      /* Text description of atom types */

      if (fscanf(filep,"%5s",atpname[k])!=1)
         goto err1;
      }

   n = nratt*(nratt+1);

   for (k=0;k<n;k++)
      {
      if (k%8==0)
         nline++;

      if (fscanf(filep,"%10f",&dummy)!=1)
         goto err1;
      }

   for (k=0;k<n;k++)
      {
      if (k%8==0)
         nline++;

      if (fscanf(filep,"%10f",&dummy)!=1)
         goto err1;
      }

   if (fgets(line,120,filep)==NULL)
      goto err1;

   nline++;

   if (fgets(line,120,filep)==NULL)
      goto err1;

   if (sscanf(line,"%5d",&nraa2)!=1) /* Number of residues */

      goto err1;

   if (nraa2<1)
      goto err1;

   resname = malloc_or_exit_if_error(sizeof(resnametype)*nraa2,"read_gromos_mt");

   for (k=0;k<nraa2;k++)
      {
      if (k%16==0)
         nline++;
      if (fscanf(filep,"%5s",resname[k])!=1) /* Description of Residues */

         goto err1;
      }

   if (fgets(line,120,filep)==NULL)
      goto err1;
   nline++;

   if (fgets(line,120,filep)==NULL)
      goto err1;
   if (sscanf(line,"%5d",&nrp)!=1)
      goto err1;

   if (nrp<1 || nrp>MAXATM)
      goto err1;

   /* Read the "atom stuff" for each atom. */

   for (k=0;k<nrp;k++)
      {
      nline++;
      if (fgets(line,120,filep)==NULL)  goto err1;

      /* BUG: CWM changed sscanf input to "ld" from "d" */

      /* beause .mres and .iac are longs - not ints! */


      if (sscanf(line,"%5ld%5s%5ld%10f",
                       &atmstf[k].mres,
                           atmstf[k].label,
                               &atmstf[k].iac,
                                   &atmstf[k].winv)!=4)
         {
         goto err1;
         }
      nline++;

      if (fgets(line,120,filep)==NULL)
         goto err1;
      }

   /* Now take the atmstf from this specific format */

   /* and populate _our_ attribs array */

   for (i=0;i<nrp;i++)
      {
      attribs[i].res = (unsigned short)atmstf[i].mres; /* BUG: typecast OK? */

      attribs[i].moln = 0;
      strncpy(attribs[i].label,atmstf[i].label,4);
      attribs[i].label[4]=0;
      attribs[i].color=1;
      attribs[i].showlabl=0;
      }

   nline++;
   if (fgets(line,120,filep)==NULL)
      goto err1;

   /* Skip over the nbty - BUG: Comment what is nbty??? */


   if (sscanf(line,"%5d",&nbty)!=1)
      goto err1;

   n = (nbty+3)/4;

   for (k=0;k<n;k++)
      {
      nline++;
      if (fgets(line,82,filep)==NULL)
         goto err1;
      }

   /* Bonds */
   nline++;
   if (fgets(line,120,filep)==NULL)
      goto err1;
   /* Read about bonds involving hydrogen */

   if (sscanf(line,"%5d",&nbonh)!=1)
      goto err1;

   for (k=0;k<nbonh;k++)
      {
      nline++;

      if (fgets(line,120,filep)==NULL)
         goto err1;

      if (sscanf(line,"%5d%5d",&ibond[k][0],&ibond[k][1])!=2)
         goto err1;
      }

   nline++;

   if (fgets(line,120,filep)==NULL)
      goto err1;

   /* Now handle bonds which don't involve H */

   /* BUG: Notice that in binary version, bin: below, we */

   /* seem to throw around bonds involving H */

   if (sscanf(line,"%5d",&nbond)!=1)
      goto err1;

   {
   /* Is this an OK var name????  CWM does not want */

   /* to recompute on each iteration... */

   int totalBondCount = nbonh+nbond;
   for (k=nbonh;k<totalBondCount;k++)
      {
      nline++;
      if (fgets(line,120,filep)==NULL)
         goto err1;

      if (sscanf(line,"%5d%5d",&ibond[k][0],&ibond[k][1])!=2)
         goto err1;
      }

   /* Convert the Gromos bonds to our format */

   for (k=0;k<totalBondCount;k++)
      {
      ib[k+1][0] = (ibond[k][0]-1)*3;
      ib[k+1][1] = (ibond[k][1]-1)*3;
      ib[k+1][2] =1;
      }
   }

   /* Close the input FILE. */

   fclose_or_exit_if_error(gromos_mtFileName,filep); /* CWM moved this here to free a descriptor */


   filep = fopen_or_exit_if_error("bond.dat","w");
   fprintf(filep,"%5d %5d\n",nbonh,nbond);
   for (k=0;k<nbonh+nbond+1;k++)
      fprintf(filep,"%5d %5d %5d\n",ib[k][0],ib[k][1],ib[k][2]);
   fclose_or_exit_if_error("bond.dat",filep);

   /* Ask about solvent, and add if requested */
   printf(" Number of solvent (water) molecules: ");
   scanf("%d",&nsol);
   printf("\n");
   *nsolv = nsol;

   if (nsol>0)
      {
      printf(" Keep solvent? ");
      /* BUG: CWM changed to '4's to avoid buffer overrun of ans */

      if ( GetYesOrNo() )
         {
         *isolv = 1;
         ires = nraa2;

         /* Add solvent molecules at the end of the array */

         /* after the atoms */

         for (i=nrp;i<nrp+3*nsol;i+=3)
            {
            ires++;
            strncpy(attribs[i].label,"HW1 ",4);
            strncpy(attribs[i+1].label,"OW  ",4);
            strncpy(attribs[i+2].label,"HW2 ",4);
            for (k=0;k<3;k++)
               {
               attribs[i+k].label[4]=0;
               attribs[i+k].res = (unsigned short)ires; /* BUG: typecast OK? */

               attribs[i+k].moln = 1;
               attribs[i+k].color=3;
               attribs[i+k].showlabl=0;
               }
            }

         k = nrp; /* Set k to first solvent molecule after atoms */


         /* Now add these solvent molecules to the bond list */

         /* BUG: Need comment on what is going on here.... */

         /* Seems like we are linking every 3rd solvent molecule?? */

         for (i=nbond+nbonh+1;i<nbond+nbonh+2*nsol;i+=2)
            {
            ib[i][0] = 3*k;
            ib[i][1] = 3*(k+1);
            ib[i][2] = 1;
            ib[i+1][0] = 3*k;
            ib[i+1][1] = 3*(k+2);
            ib[i+1][2] = 1;
            k+=3;
            }
         }
      else
         {
         *isolv = 0;
         nsol = 0;
         }
      }

   nat = nrp+3*nsol;

/*  Add residue info at end of attrib array  */
   {
   j = nat+1;

   /* CWM optimzed to avoid some repeated subscript calculations */

   for (i=0;i<nraa2;i++)
      {
      strncpy(attribs[j].label,resname[j],3); /* Residue descriptions */

      attribs[j].label[3]=0;
      attribs[j].color=(-1);
      attribs[j].res=(unsigned short)(i+1); /* BUG - typecast OK? */

      j++;
      }

/* old code: for (j=nraa2;j<nraa2+nsol;j++) */


   j = nraa2+nat+1; /* <-- Yes, this is redundant -but readaboe */

   for (i=0;i<nsol;i++)
      {
      strncpy(attribs[j].label,"SOL",3);
      attribs[j].label[3]=0;
      attribs[j].color=(-1);
      attribs[j].res=(unsigned short)(i + nraa2 + 1); /* BUG - typecast OK? */

      j++;
      }
   }

  ib[0][0]=nrp+3*nsol;
  ib[0][1]=nraa2+nsol;
  ib[0][2]=nbond+nbonh+2*nsol;

  goto write;

err1:
  printf("\n Error reading gromos formatted molecular topology file at line %d\n\n",nline);
  fclose(filep);


/* Path for reading unformatted (binary) topology file */

bin:
   /* BUG: Need to verify count of readuf lines against original.c */


   {
   char buff[12*MAXATM+1000]; /* Was static in readuf earlier */


   /* BUG: Filename was missing in original version! */

   printf("Attempt to read %s as a binary topology file?\n",
      gromos_mtFileName);

   if ( ! GetYesOrNo() )
     exit(EXIT_FAILURE);

    /* CWM changed to ANSI C fopen here... */

   filep = fopen_or_exit_if_error(gromos_mtFileName,"rb");

   record_count=0;
   readuf(filep,buff,sizeof(buff));

   getuflong(buff,&nratt,1);        /* number of atom types */

   readuf(filep,buff,sizeof(buff));

   /* BUG: This assumes binary data is 8*nratt */

   /* However, atpname[x] is only size 4.  What gives? */

   getufchar(buff,atpname[0],8*nratt);    /* atom types text */

   readuf(filep,buff,sizeof(buff));

   getuflong(buff,&nraa2,1);        /* number of residues */

   readuf(filep,buff,sizeof(buff));
   /* BUG: This assumes binary data is 8*nraa2 */

   /* However, resname[x] is only size 4.  What gives? */

   getufchar(buff,resname[0],8*nraa2);    /* residue names */

   readuf(filep,buff,sizeof(buff));
   getuflong(buff,&nrp,1);          /* number of solute atoms */

/* *natp = nrp;  Old code */


   readuf(filep,buff,sizeof(buff));

   /* This is a major BUG... atmstf is not made of longs! */

   /* What gives here????? */

#ifdef BUG_NEEDS_FIX
   i = getuflong(buff,atmstf,5*nrp);    /* atom info */
#endif
   for (i=0;i<nrp;i++)
      {
      strncpy(attribs[i].label,atmstf[i].label,4);
      attribs[i].label[4]=0;
      attribs[i].color=0;
      attribs[i].showlabl=0;
      }


   readuf(filep,buff,sizeof(buff));

   readuf(filep,buff,sizeof(buff));

   readuf(filep,buff,sizeof(buff));
   getuflong(buff,&nbonh,1);        /* number of bonds with H */

   readuf(filep,buff,sizeof(buff));

   readuf(filep,buff,sizeof(buff));
   getuflong(buff,&nbond,1);        /* number of bonds - not H */

   readuf(filep,buff,sizeof(buff));

   /* This code was getuflong(&ipb,buff,3*nbond) */

   /* BUG - Need documentation of this longArray business... */

   getuflong(buff,longArray[0],3*nbond);

   /* Update our bond data with gromos binary data */

   for (i=0;i<nbond;i++)
     {
     ib[i+1][0] = (longArray[i][0]-1)*3;
     ib[i+1][1] = (longArray[i][1]-1)*3;
     ib[i+1][2] =1;
     }

   /* BUG::::  Why do we discard bonds involving hydrogen here */

   /* while keeping them in the text version of teh GROMOS fil */


   ib[0][0]=nrp;
   ib[0][1]=nraa2;
   ib[0][2]=nbond;
   ib[nbond+1][0]=0;
   ib[nbond+2][0]=0;
   ib[nbond+3][0]=0;
   nat = nrp;
   nsol = 0;

   fclose_or_exit_if_error(gromos_mtFileName,filep);
   } /* End bin: */


write:

  for (k=0;k<10;k++)
    printf("%s %d %d %d %d\n",attribs[k].label,attribs[k].moln,
          attribs[k].res,attribs[k].color,attribs[k].showlabl);

/* Write Attribute File */
   WriteAttributeFile(nat+nraa2+1);

#ifdef OLD_REPLACED_BY_LINE_ABOVE
  fdes=open(attribfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
  checkwrite(fdes,attribfile);
  printf("Writing attribute file.\n");
  nbytes=(nat+nraa2+1)*sizeof(atomattrib);
  write(fdes,(char *)attribs,(unsigned)nbytes);
  close(fdes);
#endif

  for (k=0;k<11;k++)
    printf("%5d %5d %5d\n",ib[k][0],ib[k][1],ib[k][2]);

/* Write Bond File */
   WriteBondFile((nbonh+nbond+2*nsol)+1,ib);

#ifdef OLD_CODE_REPLACED_BY_LINE_ABOVE
   fdes=open(bondfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
   checkwrite(fdes,bondfile);
   printf("Writing bond file.\n");
   nbytes=((nbonh+nbond+2*nsol)*3+3)*sizeof(int);
   write(fdes,(char *)ib,(unsigned)nbytes);
   close(fdes);
#endif

  return(nrp); /* I'm frankly not too convinced... */

}



void read_gromos_traj(int natom,int nsolv,int isolv)
{
   int i,k,nbytes,ifr,ifr2;
   int natot;
   float xyz[3*MAXATM], ave[3], dum;
   char line[128];
   short *xx;
   FILE *filep;
/* BUG: Need to verify count of readuf lines against original.c */


   if (isolv)
      natot = natom +3*nsolv;
   else
      natot = natom;

/* Allocate dynamic memory for coords */
   if (frames==9999)
      {
      nbytes=2*3*(natot)*101;
      frames = 100;
      }
   else
      nbytes=2*3*(frames/frinc+1)*(natot);

   printf("Allocating %d bytes of memory\n",nbytes);
   xx  = (short *)malloc_or_exit_if_error(nbytes,"read_gromos_traj - xx");

/* First try to read a formatted trajectory file, and if that fails */
/* try to read a binary file */

   printf("Attempting to read formatted trajectory file\n");
   filep = fopen_or_exit_if_error(coordFileNames[0],"r");

  /* BUG: It is not clear why fgets returning NULL means */

  /* this is a binary file. */


   if (fgets(line,120,filep)==NULL)
      goto bin;

/* ifr1 = 1; */

   ifr2 = MAXFRAMES;
   ifr=0;
   record_count=0;

   while (ifr<ifr2)
      {
      printf("%5d\n",ifr);
      for (i=0;i<natot;i++)
         for (k=0;k<3;k++)
            if (fscanf(filep,"%f",&xyz[i*3+k])!=1)
               {
               fclose_or_exit_if_error(coordFileNames[0],filep);
               if (ifr==0)
                  goto bin;
               else
                  {
                  ifr--;
                  goto write;
                  }
               }

      if (isolv==0)
         {
         for (i=0;i<nsolv;i++)
            for (k=0;k<9;k++)
               if (fscanf(filep,"%f",&dum)!=1)
                  {
                  fclose_or_exit_if_error(coordFileNames[0],filep);
                  goto write;
                  }
         }

      /* Appears to be read of dummy record... or finish of line */

      fgets(line,120,filep);

      for (k=0;k<3;k++)
         ave[k]=0;
      for (i=0;i<natot;i++)
         for (k=0;k<3;k++)
            ave[k] += xyz[3*i+k];

      for (k=0;k<3;k++)
         ave[k] /= natot;

      for (i=0;i<natot;i++)
         for (k=0;k<3;k++)
            xx[3*(ifr*natot+i)+k] = (short)(KAPPA*(xyz[3*i+k]-ave[k])*gxscale);
      ifr++;
      }

   fclose_or_exit_if_error(coordFileNames[0],filep);
   goto write;

bin:
   {
   char buff[12*MAXATM+1000]; /* Was static in readuf earlier */


   printf("Attempting to read binary trajectory file\n");

   /* BUG: IS rb OK for UNIX? */

   filep = fopen_or_exit_if_error(coordFileNames[0],"rb");

/* ifr1 = 1; */

   ifr2 = MAXFRAMES;
   ifr=0;
   record_count=0;

   while (ifr<ifr2)
      {
      readuf(filep,buff,sizeof(buff));
      getuffloat(buff,xyz,3*natom);

      for (k=0;k<3;k++)
         ave[k]=0;

      for (i=0;i<natom;i++)
         for (k=0;k<3;k++)
            ave[k] += xyz[3*i+k];

      for (k=0;k<3;k++)
         ave[k] /= natom;

      for (i=0;i<natom;i++)
         for (k=0;k<3;k++) /* BUG: Nasty typecast here: */

            xx[3*(ifr*natom+i)+k] = (short)(KAPPA*(xyz[3*i+k]-ave[k]));
      }
   fclose_or_exit_if_error(coordFileNames[0],filep);
   return;
   } /* end of bin: BUG: Why no binary coordiante file write here? */


write:
   /* WRITE BINARY COORDINATE FILE WITH ALL FRAMES */
   WriteBinaryCoordinateFileWithAllFrames(
            xx,     /* The coordinates */

            natot,  /* # of Atoms */

            3*(natot+oneIfInBox)*sizeof(short), /* bytesperframe */

            ifr+1); /* # of Frames */


   return;
}


/*------------------------------------------------------------*/


/* Read Biosym .his file */

#define MAXBND MAXATM

void read_bs_his(float ave[3],int ib[MAXDRAW][3])
{
   FILE *filep;
   int i,k,jres,mol;
   int ifr,ifr1,ifr2,nbytes,bytesperframe;
   short *xx;
   long icontrol,nattyp,nnres,nmol,natom,natmov,nres,idum;
   long iatx[MAXATM], imove[MAXATM];
   typedef long long2[2];
   int* istmol;
   long* natmol;
   long* nresmol;
   long* ixresname;
   long2* iatres;
   long nbond, ibond[MAXBND][2];
   float xyz[3*MAXATM];
   char version[80],title1[80],title2[80];
   typedef char char4[4];
   char atpname[MAXATP][4], atname[MAXATM][4];
   char4* resname;
   double dvers, atpmass[MAXATP];
   char buff[12*MAXATM+1000]; /* Was static in readuf earlier */

   char* bufptr;

   /* BUG: Need to verify count of readuf lines against original.c */



   filep = fopen_or_exit_if_error(bshFileName,"rb");

/* Header info */
   record_count = 0;
   readuf(filep,buff,sizeof(buff));
   getuflong(buff,&icontrol,1);

   readuf(filep,buff,sizeof(buff));
   bufptr = getufchar(buff,version,80);
   /* BUG: In original, this was not &dvers - but dvers */

   /*      Thank heaven for prototypes! */

   getufdouble(bufptr,&dvers,1);

   readuf(filep,buff,sizeof(buff));
   getufchar(buff,title1,80);

   readuf(filep,buff,sizeof(buff));
   getufchar(buff,title2,80);

/* Atom types and masses */
   readuf(filep,buff,sizeof(buff));
   bufptr = getuflong(buff,&nattyp,1);
   if (nattyp>MAXATP)
      {
      printf("\nError: (# of atom types) NATTYP=%d exceeds (Max # atom types) MAXATP=%d\n",nattyp,MAXATP);
      exit(EXIT_FAILURE);
      }

   /* BUG: Originl was just atpname....  Is this kosher??? */

   bufptr = getufchar(bufptr,atpname[0],4*nattyp);
   bufptr += 4; /* BUG: Skip over a long value... Need to document! */

   getufdouble(bufptr,atpmass,nattyp);

/* Residue names */
   readuf(filep,buff,sizeof(buff));
   bufptr = getuflong(buff,&nnres,1);
   allocate_global_rlab(nnres,"read_bs_his");
#ifdef OLD_IDEA
   if (nnres>MAXRES)
      {
      printf("\nError:  (# of residues) NNRES=%d exceeds (Max # of residues) MAXRES=%d\n\n",nnres,MAXRES);
      exit(EXIT_FAILURE);
      }
#endif

   /* BUG - originally just resname,...  Is this OK??  Need docuentation on format */

   resname = malloc_or_exit_if_error(sizeof(char4)*nnres,"read_bs_his:resname");
   getufchar(bufptr,resname[0],4*nnres);

/* Atoms:  index to atom types, and atom names */
   readuf(filep,buff,sizeof(buff));
   bufptr = getuflong(buff,&natom,1);

#ifdef OLD_IDEA
   if (natom>MAXATM)
      {
      printf("\nError:  NATOM=%d exceeds MAXATM=%d\n\n",natom,MAXATM);
      exit(EXIT_FAILURE);
      }
#endif

   bufptr = getuflong(bufptr,iatx,natom);
   /* BUG: Is it OK to do this with atname... Need doc on format */

   getufchar(bufptr,atname[0],4*natom);

/* Atoms which move */
   readuf(filep,buff,sizeof(buff));

   bufptr = getuflong(buff,&idum,1);
   bufptr = getuflong(bufptr,&natmov,1);
   if (natmov>natom)
      {
      printf("NATMOV=%d > NATOM=%d??\n\n",natmov,natom);
      exit(EXIT_FAILURE);
      }

   getuflong(bufptr,imove,natmov);

/* Breakdown of atoms, residues, mols */
   readuf(filep,buff,sizeof(buff));
   bufptr = getuflong(buff,&nmol,1);

#ifdef OLD_IDEA
   if (nmol>MAXMOL)
      {
      printf("\nError: (# of molecules) NMOL=%d exceeds MAXMOL=%d\n\n",nmol,MAXMOL);
      exit(EXIT_FAILURE);
      }
#endif

   natmol = malloc_or_exit_if_error(sizeof(long)*(nmol+1),"read_bs_his:natmol");
   bufptr = getuflong(bufptr,natmol,nmol);

   nresmol = malloc_or_exit_if_error(sizeof(long)*(nmol+1),"read_bs_his:nresmol");
   getuflong(bufptr,nresmol,nmol);

   istmol = malloc_or_exit_if_error(sizeof(int)*(nmol+1),"read_bs_his:istmol");

   {
   int na=0;
   for (i=0;i<=nmol;i++);
      {
      istmol[i] = na;
      na += natmol[i];
      }
   }

/* Residues */
   readuf(filep,buff,sizeof(buff));
   bufptr = getuflong(buff,&nres,1);
   allocate_global_rlab(nres,"read_bs_his");
      allocate_global_atomattrib(nres,natom,"read_bs_his"); 

   /* BUG - need to verify iatres is in expected format */

   iatres = malloc_or_exit_if_error(nres * sizeof(long2),"read_bs_his: ixiatres");
   bufptr = getuflong(bufptr,iatres[0],2*nres);

   ixresname = malloc_or_exit_if_error(nres * sizeof(long),"read_bs_his: ixresname");
   getuflong(bufptr,ixresname,nres);

/* Bonds */
   readuf(filep,buff,sizeof(buff));
   getuflong(buff,&nbond,1);
   if (nbond>MAXBND)
      {
      printf("\nError:  NBOND=%d exceeds MAXBND=%d\n\n",nbond,MAXBND);
      exit(EXIT_FAILURE);
      }

   readuf(filep,buff,sizeof(buff));
   /* BUG - need to check ibond format in record */

   getuflong(buff,ibond[0],2*nbond);
   for (i=0;i<nbond;i++)
      {
      ib[i+1][0] = (ibond[i][0]-1)*3;
      ib[i+1][1] = (ibond[i][1]-1)*3;
      ib[i+1][2] =1;
      }
   ib[0][0]=natom;
   ib[0][1]=nres;
   ib[0][2]=nbond;
   ib[nbond+1][0]=0;
   ib[nbond+2][0]=0;
   ib[nbond+3][0]=0;

/* Store atom attributes */
   for (i=0;i<natom;i++)
      {
      strncpy(attribs[i+1].label,atname[i],4);
      attribs[i+1].label[4]=0;
      attribs[i+1].color=1;
      attribs[i+1].showlabl=0;
      }

   for (mol=0;mol<nmol;mol++)
      {
      for (i=istmol[mol];i<istmol[mol+1];i++)
         attribs[i].moln=(short)mol; /* BUG: Is type conversion OK? */

      }

   for (jres=0;jres<nres;jres++)
      {
      for (i=iatres[jres][0];i<=iatres[jres][1];i++)
         attribs[i].res=(unsigned short)(jres+1); /* BUG: Is type conversion OK? */

      }

/*  Add other residue info at end of attrib array  */
   for (jres=0;jres<nres;jres++)
      {
      strncpy(attribs[natom+jres+1].label,resname[ixresname[jres]-1],3);
      attribs[natom+jres+1].label[3]=0;
      attribs[natom+jres+1].color=(-1);
      attribs[natom+jres+1].res=(unsigned short)(jres+1); /* BUG: Type conversion */

      }

/* Allocate dynamic memory for coords */
   if(frames==9999)
      nbytes=sizeof(short)*3*(natom)*102;
   else
      nbytes=sizeof(short)*3*(frames/frinc+2)*(natom);

   printf("Allocating %d bytes of memory\n",nbytes);
   xx  = (short*)malloc_or_exit_if_error(nbytes,"read_bs_his - xx");

/* Write Attribute File */
   WriteAttributeFile(natom+nres+1);

#ifdef OLD_REPLACED_BY_LINE_ABOVE
  fdes=open(attribfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
  checkwrite(fdes,attribfile);
  printf("Writing attribute file.\n");
  nbytes=(natom+nres+1)*sizeof(atomattrib);
  write(fdes,(char *)attribs,(unsigned)nbytes);
  close(fdes);
#endif

/* Write Bond File */
   WriteBondFile(nbond+2,ib);

#ifdef OLD_CODE_REPLACED_BY_LINE_ABOVE
   fdes=open(bondfile,O_WRONLY | O_CREAT | O_TRUNC,00666);
   checkwrite(fdes,bondfile);
   printf("Writing bond file.\n");
   nbytes=(nbond*3+6)*sizeof(int);
   write(fdes,(char *)ib,(unsigned)nbytes);
   close(fdes);
#endif


/* Skip cell and energy stuff */
   /* This one allows i == 4, file too big for buffer... */

   record_count++;
   i = readuf_core(filep,buff,sizeof(buff));
   if (i && i!=4)
      readuf_error(i);

   readuf(filep,buff,sizeof(buff));
   readuf(filep,buff,sizeof(buff));

/* Coordinates (first frame) */
   readuf(filep,buff,sizeof(buff));

   getuffloat(buff,xyz,3*natom);

   for (k=0;k<3;k++)
      ave[k]=0;

   for (i=0;i<natom;i++)
     for (k=0;k<3;k++)
       ave[k] += xyz[3*i+k];

   for (k=0;k<3;k++)
      ave[k] /= natom;
   for (i=0;i<natom;i++)
      for (k=0;k<3;k++) /* BUG: Typecast needs check */

         xx[3*i+k] = (short)(KAPPA*(xyz[3*i+k]-ave[k]));

/* Velocities (skip) */
  readuf(filep,buff,sizeof(buff));

/* Now read remaining frames */
   ifr1 = 1;
   ifr2 = MAXFRAMES;
   ifr=0;

   while (ifr<ifr2)
      {
      record_count++;
      i = readuf_core(filep,buff,sizeof(buff));
      if (i!=0)
         break;

      ifr++;
      printf("%3d\n",ifr);

      /* This one does allow retcd of 4 - so we call readuf_core! */

      record_count++;
      i = readuf_core(filep,buff,sizeof(buff));
      if (i && i!=4)
         readuf_error(i);

      readuf(filep,buff,sizeof(buff));

      readuf(filep,buff,sizeof(buff));
      getuffloat(buff,xyz,3*natmov);
      if (natmov != natom)
         {
         for (i=0;i<natom;i++)
            for (k=0;k<3;k++)
               xx[(ifr*natom+i)*3+k] = xx[3*i+k];
         }

      for (i=0;i<natmov;i++)
         for (k=0;k<3;k++) /* BUG: Typecast needs check */

            xx[(ifr*natom+imove[i]-1)*3+k] = (short)(KAPPA*(xyz[3*i+k]-ave[k]));

      readuf(filep,buff,sizeof(buff));
      }

   fclose_or_exit_if_error(bshFileName,filep);

  /* Write binary coordinate file with all frames */


   frames = ifr-ifr1+1;
   bytesperframe=3*(natom+oneIfInBox)*sizeof(short);

   WriteBinaryCoordinateFileWithAllFrames(
                            xx,
                            natom,
                            bytesperframe,
                            frames);

   free(istmol);
   free(nresmol);
   free(natmol);

/* BUG: CWM returns to caller instead of exiting..... */

/* check out main... is this OK??? */

/* exit(EXIT_SUCCESS); */

}


/*------------------------------------------------------------*/

/* This code was repeated several times in main */

/* so CWM consolidated it.... */

void Copy_optarg_to_defaultFileName(char* defaultFileName,const char* optarg_ptr)
{
   int i=0;

   /* Note that we do _not_ copy over the file extension. */

   /* The global name[] becomes a default for other files. */

   while ((*optarg_ptr)!='.' && (*optarg_ptr) && i<(FILENAME_MAX-5))
      defaultFileName[i++] = *(optarg_ptr++);

   defaultFileName[i]=0;
}

/* BUG - this is not complete documentation. */

/* BUG - We should output this anytime things go wrong on command line */

/* BUG - would be nice to have a more complete help screen. */


void UsageInfo(void)
{
   printf("\nUsage:\n\npreproc -p AmberParmFile -c Coordinates -n movieFileName\\n");
   printf("\n\nOther legacy options include:\n\n");
   printf(" -P pdb_in [-s surf_in][-f #frames][-x][-j]\n\nFlags: -d for Debug\n");
   printf("\n\nSee mddisplay.pdf at\n  http://www.structbio.vanderbilt.edu/~cmoth/mddisplay\nfor more information.\n");
}

void UsageError(void)
{
   UsageInfo();

   // printf("preproc -?\n\nfor more information.\n");
   exit(EXIT_FAILURE);
}

void UsageHelp(void)
{
   UsageInfo();
   exit(EXIT_FAILURE);
}




/*------------------------------------------------------------*/

int main(int argc,char* const* argv)
{
   int i,c,nat,isolv,nsolv;
   extern char *optarg; /* I think this must work with getopt hmmm... */

   float ave[3];
   char *ptr;
   char ch;
   int forceBox = 0;
   int3* ib = (int3*)malloc_or_exit_if_error(MAXDRAW*sizeof(int3),"main - ib");
/* BUG/OLD: int ib[MAXDRAW][3];  Bonds from atom ib[n][0] to ib[n][1] */

   char defaultFileName[FILENAME_MAX];

   memset(defaultFileName,0,sizeof(defaultFileName));

   printf("MD preprocessor version %s\n\n",VERSION);

   /* BUG: need documentation on these things below */

   /* BUG: getopt is not ANSI C */

   while ((c=getopt(argc,argv,"h:s:a:o:p:P:m:c:f:b:n:G:B:Cxjdq")) != -1)
      switch (c)
      {
/*  m - charmm binary psf */
      case 'h':
         sscanf(optarg,"%f",&hbcut);
         break;

      case 'a':
         sprintf(attribFileName,"%s",optarg);
         break;

      case 'd':
         debugFlag = 1; /* Turn on some verbosity for this run. */

         break;

      case 'q':
         forceBox = 1;
         break;

      case 'c':
         coordFileNames[nCoordFiles++]=optarg;
#ifdef DEVELOPLATER
         optind++;
         while (optind < argc && argv[optind][0] != '-')
            {  
            coordFileNames[nCoordFiles++] = argv[optind++];
            }
         optind--;
         optarg = argv[optind];
#endif
         Copy_optarg_to_defaultFileName(defaultFileName,optarg);
         break;

      case 's':
         insurf=optarg;
         break;

      case 'b':
         sprintf(bondFileName,"%s",optarg);
         break;

      case 'o':
         sprintf(binFileName,"%s",optarg);
         break;

      case 'C':
         conectflag=1;
         break;

      case 'p':
         sprintf(parameterFileName,"%s",optarg);
         break;

      case 'P':
         pdbFileName=optarg;
         printf("Pdb file = %s",pdbFileName);
         Copy_optarg_to_defaultFileName(defaultFileName,optarg);
         break;

      case 'm':
         psfFileName=optarg;
         printf("Psf file = %s\n",psfFileName);
         Copy_optarg_to_defaultFileName(defaultFileName,optarg);
         break;

      case 'G':
         gromos_mtFileName=optarg;
         printf("GROMOS molecular topology file = %s\n",gromos_mtFileName);
         Copy_optarg_to_defaultFileName(defaultFileName,optarg);
         break;

      case 'B':
         bshFileName=optarg;
         printf("Biosym .his file = %s\n",bshFileName);
         Copy_optarg_to_defaultFileName(defaultFileName,optarg);
         break;

      case 'f':
         ptr=strstr(optarg,"%");
         if(ptr != NULL)
            {
            sscanf(ptr,"%c%d",&ch,&frinc);
            *ptr=0;
            }
         else
            frinc=1;

         i = sscanf(optarg,"%d%c%d",&fr1,&ch,&fr2);

         if(i<3)
            {
            frames = fr1;
            fr1=0;
            }
         else
            frames = fr2-fr1+1;

         if (frames<=0)
            {
            printf("\n  -f:  invalid frame counts.\n\n");
            exit(EXIT_FAILURE);
            }
         break;

      case 'x':         /* means trajectory file contains box coords */
         conpress=1;
         break;

      case 'n':
         sprintf(binFileName,"%s.cor",optarg);
         sprintf(attribFileName,"%s.attr",optarg);
         sprintf(bondFileName,"%s.bnd",optarg);
         sprintf(surfaceFileName,"%s.srf",optarg);
         break;

      case 'j':
         jumpok=1;
         break;

/* A '?' is returned from getopt() if an unknown command line option */

/* is entered.  Its not the case that the user entered -? for help or something */

      case '?':
         UsageError(); /* Terminates to OS */

         break;
      }

      /* DEFAULT WITH NAME OF PDB OR TRAJECTORY FILE */
      if (defaultFileName[0])
         {
         if (!binFileName[0])
            sprintf(binFileName,"%s.cor",defaultFileName);
         if (!attribFileName[0])
            sprintf(attribFileName,"%s.attr",defaultFileName);
         if (!bondFileName[0])
            sprintf(bondFileName,"%s.bnd",defaultFileName);
         if (!surfaceFileName[0])
            sprintf(surfaceFileName,"%s.srf",defaultFileName);
         }

      /* EXIT IF NOT ENOUGH INFO SPECIFIED */
      if ((coordFileNames[0]==0 || parameterFileName==0) &&
          pdbFileName==0 &&
          psfFileName==0 &&
          bshFileName==0 &&
          gromos_mtFileName==0)
          {
          UsageError(); /* Terminates to OS */

          }

      /* CHECK TO MAKE SURE .cor FILE CAN BE OVERWRITTEN */
/*      i=0; */

      if (nCoordFiles)
         {
         int i;
         int overwrite = 0;
         /* BUG CWM sees no reason to keep this loop */

         for (i=0;i<nCoordFiles;i++)
            if ( ! strcmp(binFileName,coordFileNames[0]) )
               overwrite = 1;

         if (overwrite)
            {
            printf("Output file %s will overwrite input file of same name.\n",
                           binFileName);
            printf("Program aborted.\n\n");
            exit(EXIT_FAILURE);
            }

         {
         FILE* test=fopen(binFileName,"r");
         if (test)
            {
            /* BUG:::: HELLO - don't we need to fclose(test) here....  I think so - so ehre it it! */

            fclose_or_exit_if_error(binFileName,test);
            printf("File %s exists.  OK to overwrite?\n",binFileName);
            if ( ! GetYesOrNo() )
               exit(EXIT_SUCCESS);

            }

         }
   }

   printf("Please wait while the data is loaded.\n");
   if (pdbFileName != 0)
      {
      nat = gencon(ave,ib);
      if (coordFileNames[0])
         dotraj(nat,ave,ib);
      }
   else if (psfFileName)
      {
      nat = readpsf(/*ave,*/  ib); /* CWM notes that readpsf did not modify ave[] */

      dcdtraj(nat);
      }
   else if (gromos_mtFileName)
      {
      nat = read_gromos_mt(&nsolv,&isolv,ib);
      if (coordFileNames[0])
         read_gromos_traj(nat,nsolv,isolv);
      }
   else if (bshFileName)
      {
      /* CWM removed nat - because read_bs_his did not set this variable */

      read_bs_his(ave,ib);
      }
   else
      { // This is what is run for AMBER -p and -c options
      nat = readparm(ib,forceBox);
      dotraj(nat,ave,ib);
      }

   if (insurf)
      surfdude(nat,ave); /*do surfaces if input file is specified*/

   printf("Done.\n\n");
   free(ib);
   return EXIT_SUCCESS;
}

