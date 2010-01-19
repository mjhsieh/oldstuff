/****************************** MakeLists.c ******************************
*
*  This function builds the BondList structures for atom, bond, angle,
*  torsion, and oop items using data processed from .car and .mdf files.  For 
*  each item the process is done in two main steps:
* 
*  1) The 'id' member is filled in with atom numbers, and the 'type' member 
*     is filled in from potential types.  
*  2) The number of unique potential types is counted while
*      assigning unique types to the 'a' member
*
*  The equivalence types for each item in 'a' are looked up and filled 
*  into the 'e' item in the data structure
*
*  For class II forcefields, the cross term structures are filled using 
*  data from the appropriate diagonal structure
*
*  */

#include "Msi2LMP2.h"
#include "Forcefield.h"

void MakeLists()
{
     char temp[5];
     int i,j,k,m,n; /* counters */
     int match = 0;
     int A,B,C,D; /* counters for torsions*/
     int AtomA, AtomB, AtomC, AtomD;

     void SearchEquivalences(void);
     void AddCrossEquivalences(struct BondList *diag, 
						struct BondList *cross);
     int SortTypes (struct BondList *list, int sub1, int sub2);
					 /* Returns true (1) if any sorting
					    was done, else returns false*/

   /* Zero-out item counters and types counter */

     vdw.count     = vdw.no_types     = 0;
     bond.count    = bond.no_types    = 0;
     angle.count   = angle.no_types   = 0;
     torsion.count = torsion.no_types = 0;
     oop.count     = oop.no_types     = 0;
     angang.count  = angang.no_types  = 0;

   /* Zero-out Assignment Records */

   for( k = 0; k < MAX_TYPES; k++);
   {
      vdw.found[k]       = 0;
      bond.found[k]      = 0;
      angle.found[k]     = 0;
      torsion.found[k]   = 0;
      torsion.wc[k]	 = 0;
      oop.found[k]       = 0;
      bonbon.found[k]    = 0;
      bonbon13.found[k]  = 0;
      angang.found[k]    = 0;
      bonang.found[k]    = 0;
      endbontor.found[k] = 0;
      midbontor.found[k] = 0;
      angangtor.found[k] = 0;
      angtor.found[k]    = 0;
   }

/***************************************************************************/
/*----------------------------------Atoms----------------------------------*/

	/* Atoms are assigned to the vdw structure (non-bond).
		The id of an atom is its atom number.
		*/

   vdw.no_members=1; /* One atom to describe an atom */
   k = 0;
   while(k!=total_no_atoms)
   { /* while */

      vdw.id[k][0] = k; /* for consistency */
      strcpy(vdw.type[k][0], atoms[k].potential);
      vdw.count++; 

	/*----------------Types of Atoms---------------*/
		/* Unique atom potentials are read into the 'a' list.
		   If the type does not exist in the list, it is added */

      match = 0; /* match is true (1) when there is a match (potential type
		 	already counted) */

      for( i = 0; i < vdw.no_types; i++ )
         if( strcmp( vdw.type[k][0], vdw.a[i][0] ) == 0 )
         {
	    match = 1;
            vdw.type_no[k] = i;
         }

      if( match == 0 )
      {
         strcpy(vdw.a[vdw.no_types][0], vdw.type[k][0]);
         vdw.type_no[k] = vdw.no_types;
         vdw.no_types++;
      }

      k++;

   }/*while*/


/***************************************************************************/
/*-------------------------------Bonds-------------------------------------*/
		
	/* Bonds are assigned to the bond structure.  Bonds are found by
	   looping through all of the atoms and assigning a bond to 
	   each connection.  Bonds are only assigned when the atom number
           of the loop is less than the atom number of the connection.  This
	   is done to avoid assigning the same bond twice.  */

   bond.no_members = 2; /* Two atoms to describe a bond */
   for ( i=0; i < total_no_atoms; i++)
   {/* for loop i*/
      for (j=0; j < atoms[i].no_connect; j++)
       if ( atoms[i].no < atoms[i].conn_no[j] )
         {
            bond.id[bond.count][0] = atoms[i].no;
	    bond.id[bond.count][1] = atoms[i].conn_no[j]; 

	    strcpy( bond.type[bond.count][0],atoms[i].potential);
	    strcpy( bond.type[bond.count][1],
                     atoms[bond.id[bond.count][1]].potential);

            SortTypes (&bond, 0, 1); /* Sorts bonds so that lower numbered atom
				     types occupy [0] position.  This makes
				     the unique type searching routines
				     much more efficient */
            bond.count++;
         }
   }/*for loop i*/

/*--------------------------Types of Bonds----------------------------------*/
	/* To get unique types of bonds, next bond type in the bond list
	   is compared against bond types in the 'a' (unique types) list.
	   If bond type not in the 'a' list, it is added and the counter
	   increased. */ 

   k = 0;
   while(k!=bond.count)
   {/*while*/

     match = 0; /* match set to 1 if match is made*/
     for(i=0; i< bond.no_types;i++)
       if(!strcmp(bond.type[k][0], bond.a[i][0])&&
          !strcmp(bond.type[k][1], bond.a[i][1])) 
       {
          match = 1;
          bond.type_no[k] = i;
       }

     if( match == 0 )
     {
      strcpy(bond.a[bond.no_types][0], bond.type[k][0]);
      strcpy(bond.a[bond.no_types][1], bond.type[k][1]);
      bond.type_no[k] = bond.no_types;
      bond.no_types++;
     }
     k++;
   } /* while */

/***************************************************************************/
/*-----------------------------Angles--------------------------------------*/
	/* Angles are assigned to the angle structure.  Angles are found by
		looping through all of the atoms and assinging an angle to
		every permutation of 2 connections with the atom number.
		Because the central atom is unique, there is no 
		accidental repetition of angles */

   angle.no_members = 3;
   for (i=0; i<total_no_atoms; i++)
   { /* for loop i */
      if (atoms[i].no_connect > 1)
         for ( j=0; j < (atoms[i].no_connect - 1); j++)
            for ( k=j+1; k < atoms[i].no_connect; k++)
            { /* for loop k */
               angle.id[angle.count][0] = atoms[atoms[i].conn_no[j]].no;
               angle.id[angle.count][1] = atoms[i].no;
               angle.id[angle.count][2] = atoms[atoms[i].conn_no[k]].no;

               strcpy( angle.type[angle.count][0], 
                  	atoms[angle.id[angle.count][0]].potential);
               strcpy( angle.type[angle.count][1], 
                  	atoms[angle.id[angle.count][1]].potential);
               strcpy( angle.type[angle.count][2], 
                  	atoms[angle.id[angle.count][2]].potential);

               SortTypes (&angle, 0,2); /* Again, simplifies type search */

               angle.count++;
            } /* for loop k */
   } /* for loop i */

/*--------------------------Angles Types------------------------------------*/

   k=0;
   while(k!=angle.count)
   {/*while*/

      match = 0;
      for(i=0;i<angle.no_types;i++)
         if(!strcmp(angle.type[k][0], angle.a[i][0])&&
            !strcmp(angle.type[k][1], angle.a[i][1])&&
            !strcmp(angle.type[k][2], angle.a[i][2])) 
         {
            match = 1;
            angle.type_no[k] = i;
         }

      if( match == 0 )
      {
          strcpy(angle.a[angle.no_types][0], angle.type[k][0]);
          strcpy(angle.a[angle.no_types][1], angle.type[k][1]);
          strcpy(angle.a[angle.no_types][2], angle.type[k][2]);
          angle.type_no[k] = angle.no_types;
          angle.no_types++;
      }
      k++;
   }/* while */

/***************************************************************************/
/*----------------------------Torsions-------------------------------------*/
	/* Torsions are assigned to the torsion structure.  The assignment
	   of torsions involves several loops.  The first loops through all
	   of the atoms.  If the current atom has 2 or more connections, then
	   a second and third loop counts over the number of connections to 
	   the first connection.  If the atom from the second loop has 2 or 
	   more connections (1 back to the first atom), a fourth loop counts 
	   over that atoms connections, creating a torsion 3 1 2 4 or C A B D. 
	   Only atoms where 1<2 are processed to prevent redundant torsions */
   torsion.no_members = 4;
   for (A=0; A < total_no_atoms; A++) /* loop over all atoms */
   { /* A loop */
     AtomA=A;
     if (atoms[AtomA].no_connect > 1) /* process only if A has 2 or more connections*/
       for (B=0; B<atoms[AtomA].no_connect; B++)/* Loop over connections to A (Atom B) */
       { /*B loop */
          AtomB = atoms[A].conn_no[B];
          if (AtomB > AtomA && atoms[AtomB].no_connect > 1)
          /* Continue processing only if B > A and has 2 or more connections */
            for(C=0; C < atoms[A].no_connect; C++) /* C is the other connection to A */
            { /* C for */
               AtomC = atoms[A].conn_no[C];
               if ( AtomC != AtomB )   /* continue if C != B */
                 for (D=0; D < atoms[AtomB].no_connect;D++)
                 { /* D loop */
                    AtomD = atoms[AtomB].conn_no[D];
                    if (AtomD != AtomA && AtomD != AtomC)
                    { /* D if */
                      torsion.id[torsion.count][0] = AtomC;
                      torsion.id[torsion.count][1] = AtomA;
                      torsion.id[torsion.count][2] = AtomB;
                      torsion.id[torsion.count][3] = AtomD;

                      strcpy(torsion.type[torsion.count][0], atoms[AtomC].potential);
                      strcpy(torsion.type[torsion.count][1], atoms[AtomA].potential);
                      strcpy(torsion.type[torsion.count][2], atoms[AtomB].potential);
                      strcpy(torsion.type[torsion.count][3], atoms[AtomD].potential);

                      if( SortTypes(&torsion,1,2)) /* Function returns flag
						      that is set to true when
						      inner members of torsion 
						      were sorted */
		      {
                        m = torsion.id[torsion.count][0];
                        torsion.id[torsion.count][0] = 
						  torsion.id[torsion.count][3];
                        torsion.id[torsion.count][3] = m;

                        strcpy(temp, torsion.type[torsion.count][0]);
                        strcpy(torsion.type[torsion.count][0],
                                               torsion.type[torsion.count][3]);
                        strcpy(torsion.type[torsion.count][3], temp);
                      }
                      else if ( strcmp (torsion.type[torsion.count][1],
                                        torsion.type[torsion.count][2]) == 0 )
                         {
                      if(SortTypes(&torsion,0,3))
                      {
                         /* Must switch innner members to retain connectivity */
                         m=torsion.id[torsion.count][1];
                         torsion.id[torsion.count][1] = torsion.id[torsion.count][2];
                         torsion.id[torsion.count][2] = m;
                      }
}
                      torsion.count++;

                    }	/* D if */
                 } 	/* D loop */
            } 		/* C loop */
       } 		/* B loop */
   } 			/* A loop */

/*------------------Types of torsions----------------*/

   k=0;
   while(k!=torsion.count)
   {/* while */

      match = 0;
      /*Check to see if current torsion is in list*/
      for(i=0;i<torsion.no_types;i++)
         if(!strcmp(torsion.type[k][0], torsion.a[i][0])&&
            !strcmp(torsion.type[k][1], torsion.a[i][1])&&
            !strcmp(torsion.type[k][2], torsion.a[i][2])&&
            !strcmp(torsion.type[k][3], torsion.a[i][3])) 
         {
            match = 1;
            torsion.type_no[k] = i;
         }


      if( match == 0 )
      {
         strcpy(torsion.a[torsion.no_types][0], torsion.type[k][0]);
         strcpy(torsion.a[torsion.no_types][1], torsion.type[k][1]);
         strcpy(torsion.a[torsion.no_types][2], torsion.type[k][2]);
         strcpy(torsion.a[torsion.no_types][3], torsion.type[k][3]);
         torsion.type_no[k] = torsion.no_types;
         torsion.no_types++;
      }
      k++;
   }/* while */

/***************************************************************************/
/*-------------------------------- OOP ------------------------------------*/
	/* Out of Planes are assigned to the oop stucture.  An atom with
	   exactly three connections has an out of plane term.  The central
	   atom is stored in the [1] position */
   oop.no_members = 4;
   for( i = 0; i < total_no_atoms; i++ )
      if (atoms[i].no_connect == 3)
      {
         oop.id[oop.count][1] = atoms[i].no; /* Center atom */
         oop.id[oop.count][0] = atoms[i].conn_no[0];
         oop.id[oop.count][2] = atoms[i].conn_no[1];
         oop.id[oop.count][3] = atoms[i].conn_no[2];

         for(n=0;n<oop.no_members;n++)
            strcpy( oop.type[oop.count][n], 
				atoms[oop.id[oop.count][n]].potential);

         SortTypes(&oop,0,2);  /* Three sorts are performed to get all */
         SortTypes(&oop,2,3);  /* types in the right order */
         SortTypes(&oop,0,2);

         oop.count++;
      }

     /*---------------------False OOP for angang - Class II-----------------*/
	/* Angle-angle cross terms are derived from the oop diagonal.  However,
	   there are angle-angle terms for any atom with three OR MORE
	   connections.  This routine goes back and picks up those atoms with
	   four or more connections, appending them to the oop structure */

   if (forcefield >= 2)    /* Class II additional oop's */
      for ( i = 0; i < total_no_atoms; i++)
         if ( atoms[i].no_connect > 3 )
         {
            for ( j = 0; j < atoms[i].no_connect-2; j++ )
              for ( k = j+1; k < atoms[i].no_connect-1; k++ ) 
                for ( m = k+1; m < atoms[i].no_connect; m++ )
                {
                   oop.id[oop.count][1] = atoms[i].no;
	           oop.id[oop.count][0] = atoms[i].conn_no[j];
	           oop.id[oop.count][2] = atoms[i].conn_no[k];
	           oop.id[oop.count][3] = atoms[i].conn_no[m];

                   for ( n = 0; n < oop.no_members; n++ )
                      strcpy( oop.type[oop.count][n], 
                                        atoms[oop.id[oop.count][n]].potential);

                   SortTypes( &oop, 0, 2 );
                   SortTypes( &oop, 2, 3 );
                   SortTypes( &oop, 0, 2 );

                   oop.count++;
                }   
         }    
  
/*----------------------------OOP Types------------------------------------*/
   k = 0;
   while ( k != oop.count )
   { /* while */
      match = 0;
      for ( i = 0; i < oop.no_types; i++ )
         if ( strcmp(oop.type[k][1], oop.a[i][1]) == 0)
	    if(!strcmp(oop.type[k][0], oop.a[i][0]) &&
               !strcmp(oop.type[k][2], oop.a[i][2]) &&
               !strcmp(oop.type[k][3], oop.a[i][3])) 
            {
               match = 1;
               oop.type_no[k] = i;
            }

         if( match == 0 )
         {
            strcpy(oop.a[oop.no_types][0], oop.type[k][0]);
            strcpy(oop.a[oop.no_types][1], oop.type[k][1]);
            strcpy(oop.a[oop.no_types][2], oop.type[k][2]);
            strcpy(oop.a[oop.no_types][3], oop.type[k][3]);

            oop.type_no[k] = oop.no_types;
            oop.no_types++;
          }
          k++;

   } /* while */
 

/***************************************************************************/
	/* Now that all of the diagonal terms have been built, we will
	   fill in the 'e' items from the equivalence tables in the force
	   field file */

   puts("Searching Equivalence tables...\n");
   SearchEquivalences();
	
	/* Now we can fill in the 2nd order terms because:
	   bonbon, bonang come from angle
	   angtor, angangtor, midbontor, endbontor, bonbon13 come from torsion
	   angang comes from oop */

	/* AddCrossEquivalences copies no_member, no_types, count, a and e
	     to the diagonal structure*/

if ( forcefield >= 2)
{
   AddCrossEquivalences(&angle,&bonbon);
   AddCrossEquivalences(&angle,&bonang);
   AddCrossEquivalences(&torsion,&angangtor);
   AddCrossEquivalences(&torsion,&angtor);
   AddCrossEquivalences(&torsion,&endbontor);
   AddCrossEquivalences(&torsion,&midbontor);
   AddCrossEquivalences(&torsion,&bonbon13);
   AddCrossEquivalences(&oop,&angang);
}
/*Now all that needs to be done is to read the coefficients from the forcefield
   file.  This will be done in the ReadCoefficientsF1() and F2 functions */

}/* End MakeLists() function */



/*****************************SearchEquivalences Function *******************/

void SearchEquivalences(void)
{
   int j;  /* counter */

   void AddEquivalences(int n, int j, struct BondList *list);

   for( j = 0; j < equivalence.entries; j++ )
   {
      AddEquivalences(1,j,&vdw);
      AddEquivalences(2,j,&bond);
      AddEquivalences(3,j,&angle);
      AddEquivalences(4,j,&torsion);
      AddEquivalences(5,j,&oop);
   }
}

void AddEquivalences(int n, int j, struct BondList *list)
{
   int k, m;

   for ( k = 0; k < list->no_types; k++ )
      for ( m = 0; m < list->no_members; m++ )
         if( strcmp(equivalence.data[j].ff_equiv[0], list->a[k][m]) == 0)
            strcpy( list->e[k][m], equivalence.data[j].ff_equiv[n] );
}

/***************************** SortTypes Function ***************************/

int SortTypes (struct BondList *list, int sub1, int sub2)
{
   int return_value=0;
   int r1,r2, temp;
   char tstring[5];
   int FindAtomPotential(char temp[]);

   r1=FindAtomPotential(list->type[list->count][sub1]);
   r2=FindAtomPotential(list->type[list->count][sub2]);

   if (r1>r2)
   {
     temp                        = list->id[list->count][sub1];
     list->id[list->count][sub1] = list->id[list->count][sub2];
     list->id[list->count][sub2] = temp;

     strcpy(tstring, list->type[list->count][sub1]);
     strcpy(list->type[list->count][sub1],list->type[list->count][sub2]);
     strcpy(list->type[list->count][sub2], tstring);

     return_value=1;
   }
   return (return_value);
}

/************************ FindAtomPotential Function **********************/

int FindAtomPotential(char temp[])
{
   int k;
   for(k=0; k < vdw.no_types; k++)
     if (strcmp(temp, vdw.a[k][0])==0) return k;

   printf("Error:  Match not found\n");
return(-1);
}

/************************ AddCrossEquivalences Function **********************/

void AddCrossEquivalences(struct BondList *diag,struct BondList *cross)
{
 int j,k;

 cross->no_members = diag->no_members;
 cross->no_types   = diag->no_types;
 cross->count      = diag->count;

 for(j=0; j < diag->no_types; j++)
   for(k=0; k < diag->no_members; k++)
   {
     strcpy(cross->a[j][k],diag->a[j][k]);
     strcpy(cross->e[j][k],diag->e[j][k]);
   }
}


