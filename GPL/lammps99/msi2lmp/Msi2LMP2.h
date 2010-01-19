/******************************** Msi2LMP2.h ********************************
*
* Header file for Msi2LMP2 conversion program.
* Largely uses code from original Msi2LMP conversion program by:
*
*   Steve Lustig, Ph.D.
*   DuPont
*   Central Research                        (E) LUSTIGSR@ESVAX.DNET.DUPONT.COM
*   Experimental Station, Route 141         (T) 302-695-3899
*   Wilmington, DE 19880-0356               (F) 302-695-8207
*
* Modified by Michael Peachey to accomodate the elimination of the use of
* Discover generated .prm files
*
* Last updated 08-05-97
* by:
*
* Michael Peachey
* Cray Intern
* Cray Reasearch, A Silicon Graphics Company
* 655E Lone Oak Drive
* Eagan, MN 55121
*
* n8179@ironwood.cray.com
* peachey@halide.chem.ncsu.edu
* 612-683-3539
*
*/
 
# include <string.h>
# include <stddef.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

#ifdef   MAIN
#define  _EX
#define  _ARG(arg)    = (arg)
#else
#define  _EX           extern
#define  _ARG(arg)
#endif

#define MAX_TYPES    12000
#define MAX_MEMBERS  5
#define MAX_LINE_LENGTH   200
#define MAX_PARAMS        8
#define PI_180  0.01745329252

#define MAX_CONNECTIONS   6
#define MAX_STRING        25


struct Atom
{
  int   molecule;      /* molecule id */
  int   no;            /* atom id */
  char  name[10];      /* atom name */
  float x[3];          /* position vector */
  char  potential[5];  /* atom potential type */
  char  element[2];    /* atom element */
  float q;             /* charge */

  char res_name[8];	/* residue name */
  char res_num[8];      /* residue numer */
  int no_connect;	/* number of connections to atom */
  char connections[MAX_CONNECTIONS][MAX_STRING];
                        /* long form, connection name */
  double bond_order[6]; 
  int conn_no[6];	/* Atom number to which atom is connected */
};

struct BondList
{
  int    no_members;                   /* # members in bond category */
  int    no_params;                    /* # parameters req. to describe bond */
  int    no_types;                     /* # different types of bond category */
  int    count;                        /* # occurences of this category */
  int    found[MAX_TYPES];             /* =1 if params found */
  float  mass[MAX_TYPES];              /* atomic mass (or float storage) */
  char   t[MAX_MEMBERS][5];            /* temporary storage for searches */
  char   a[MAX_TYPES][MAX_MEMBERS][5]; /* unique atom potential types */
  char   e[MAX_TYPES][MAX_MEMBERS][5]; /* atom type equivalences */
  double param[MAX_TYPES][MAX_PARAMS]; /* potential parameters */

  int id[MAX_TYPES][MAX_MEMBERS];	/* atom numbers beloning to item */
  int type_no[MAX_TYPES];	/* unique type number of potential type*/ 
				/* i.e. which 'a' items 'type' belongs to */
  char type[MAX_TYPES][MAX_MEMBERS][5]; /* atom potential type */
  int wc[MAX_TYPES];			/* determies if wild cards have
					   been used in .frc file */
};


_EX  char   rootname[20];
_EX  char   path[20];
_EX  double pbc[9];
_EX  int    periodic   _ARG( 1 ); /* 0= nonperiodic 1= 3-D periodic */
_EX  int    forcefield _ARG( 0 ); /* 0= ClassI      1= ClassII */
_EX  int    *no_atoms;
_EX  int    no_molecules;
_EX  int    replicate[3];
_EX  int    total_no_atoms;
_EX  char   FrcFileName[MAX_LINE_LENGTH];
_EX  FILE   *CarF;
_EX  FILE   *DatF;
_EX  FILE   *FrcF;
_EX  FILE   *PrmF;
_EX  FILE   *MdfF;
_EX  FILE   *RptF;
_EX  struct Atom *atoms;
_EX  struct Molecule *molecules;
_EX  struct BondList vdw, bond, angle, torsion, oop, \
                     bonbon, bonang, angangtor, endbontor, \
                     midbontor, angtor, angang, bonbon13;

#undef   _EX
#undef   _ARG



