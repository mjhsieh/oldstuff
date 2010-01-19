/**************************** ReadCoefficientsF2.c ****************************
*
*  This function uses the BondList structures to find the forcefield
*  coefficients that LAMMPS requires
*
*  This file is for Class 2 forcefields, currently for both
*  CFF91.frc and CFF95.frc
*
*  All equations shown in comments are from the forcefield file
*/

#include "Msi2LMP2.h"
#include "Forcefield.h"

void ReadCoefficientsF2(void)
{
   int i,j,k, temp;	/* counters */
   int c, sw1, sw2; 	/* for expanding angle angle terms, constant, 
			   switch 1, and switch 2 */
   double tparm;	/* used to switch parameters in bonang */
   char tstring[5];	/* string used in switching */

   extern void SearchDynamicFrc(struct BondList *list,struct FrcFieldItem *item,
                                int reverse, int permute);
   extern void SwitchOrderOfParams(struct BondList *list, int sub1, int sub2);
   double CrossSearch(struct BondList *list, int reversible);

/************************* NON-BOND COEFFICIENTS *************************/
	/*  E = eps(ij) [2(r(ij) * / r(ij))**9 - 3(r(ij) * / r(ij))**6]
	    where    r(ij) = [(r(i)**6 + r(j)**6))/2]**(1/6)
        	     eps(ij) = 2 sqrt(eps(i) * eps(j)) *
                               r(i)^3 * r(j)^3/[r(i)^6 + r(j)^6]
	    cffXY.frc supplies r   in col1 and eps in col2
	    LAMMPS    expects  eps in col1 and r   in col2
	    Switch function will be called
     	*/

   printf("   searching nonbonds.....\n");
   vdw.no_params = 2; 
   SearchDynamicFrc(&vdw, &ff_vdw, 0, 0);
   SwitchOrderOfParams(&vdw, 0, 1);

/***************************** BOND COEFFICIENTS *****************************/
	/* E = K2 * (R - R0)^2  +  K3 * (R - R0)^3  +  K4 * (R - R0)^4
           cffXY.frc supplies R0 K2 K3 K4 
           LAMMPS expects same
	*/  

   printf("   searching bonds.....\n");
   bond.no_params = 4;
   SearchDynamicFrc(&bond, &ff_bond, 1, 0);

/**************************** ANGLE COEFFICIENTS *****************************/
	/* Delta = Theta - Theta0
	   E = K2 * Delta^2  +  K3 * Delta^3  +  K4 * Delta^4
           cffXY.frc supplies Theta0 K2 K3 K4
           LAMMPS expects same
        */

   printf("   searching angles.....\n");
   angle.no_params = 4;
   SearchDynamicFrc(&angle, &ff_ang, 1,0);

/************************** TORSION COEFFICIENTS *****************************/
	/* E = SUM(n=1,3) { V(n) * [ 1 + cos(n*Phi - Phi0(n)) ] }
           cffXY.frc supplies V1 Phi0 V2 Phi0 V3 Phi0 
           LAMMPS expects same
        */


   printf("   searching torsions.....\n");
   torsion.no_params = 6;
   SearchDynamicFrc(&torsion, &ff_tor, 1, 0);

/********************** OUT-OF-PLANE COEFFICIENTS ***************************/
	/* E = K * (Chi - Chi0)^2
           cffXY.frc supplies K Chi0
           LAMMPS expects same
	*/

   printf("   searching out-of-planes.....\n");
   oop.no_params = 2;
   SearchDynamicFrc(&oop, &ff_oop, 0, 1);

/************************ BOND-BOND COEFFICIENTS *****************************/
	/* E = K(b,b') * (R - R0) * (R' - R0')
           cffXY.frc supplies K(b,b') 
           LAMMPS expects     K(b,b') R0 R0'
	   Reference bond length must be added for each bond 
        */

   printf("   searching bond-bonds.....\n");
   bonbon.no_params = 1;
   SearchDynamicFrc(&bonbon, &ff_bonbon, 1, 0);

  /* Need to add reference R0 for each bond to parameter list */

   bonbon.no_params = 3;
   for ( k = 0; k < bonbon.no_types; k++ )
   {
      strcpy(bond.t[0], bonbon.e[k][0]);
      strcpy(bond.t[1], bonbon.e[k][1]);
      bonbon.param[k][1] = CrossSearch(&bond, 1);

      strcpy(bond.t[0], bonbon.e[k][1]);
      strcpy(bond.t[1], bonbon.e[k][2]);
      bonbon.param[k][2] = CrossSearch(&bond, 1);
   }

/************************ BOND-ANGLE COEFFICIENTS *****************************/
	/* E = K * (R - R0) * (Theta - Theta0) 
           cffXY.frc supplies K(b,theta) K(b',theta) 
           LAMMPS expects     K(b,b') R0 R0'
	   Reference bond length must be added for each bond 
	   Reference theta0 is not needed

Order of parameters is very important... Can only do forward search.  If not found, switch order of identities and try again.
        */

   printf("   searching bond-angles.....\n");
   bonang.no_params = 2;
   SearchDynamicFrc(&bonang, &ff_bonang, 1,0);

 /* Need to put K's in correct column and add reference bond lengths */
   bonang.no_params = 4;

   for ( k = 0; k < bonang.no_types; k++ )
   {
      if ( bonang.found[k] == 2) /* found switch set to 2 indicates coeff's
				    were found in reverse search, params
				    must be switched */
      {
         tparm = bonang.param[k][0];
         bonang.param[k][0] = bonang.param[k][1];
	 bonang.param[k][1] = tparm;
      }
 
      strcpy(bond.t[0], bonang.e[k][0]);
      strcpy(bond.t[1], bonang.e[k][1]);
      bonang.param[k][2] = CrossSearch(&bond, 1);

      strcpy(bond.t[0], bonang.e[k][1]);
      strcpy(bond.t[1], bonang.e[k][2]);
      bonang.param[k][3] = CrossSearch(&bond, 1);

   }

/************************* ANGLE-ANGLE COEFFICIENTS **************************/
	/* E = K * (Theta - Theta0) * (Theta' - Theta0')
	   cffXY.frc supplies K
	   LAMMPS is expecting K1 K2 K3 Theta0 Theta0 Theta0

	   Each angle angle is stored as an out of plane that will expand 
	   into three angle-angle terms.  To get all of the coefficents LAMMPS 
	   is expecting, a complete list of angle-angle terms is generated by 
	   tripling the angang structure.  All individual coefficients are i
	   looked up, then the list is collapsed.  Finally, the theta0 terms i
	   from the angles are added
	*/

   printf("   searching angle-angles.....\n");
   angang.no_params = 1;

 /* First, add theta0 terms */
       
   for ( k = 0; k < angang.no_types; k++ )
   {
      strcpy( angle.t[0], angang.e[k][0] );     /* A B C is in [3] */
      strcpy( angle.t[1], angang.e[k][1] );
      strcpy( angle.t[2], angang.e[k][2] );
      angang.param[k][3] = CrossSearch(&angle, 1);

      strcpy( angle.t[0], angang.e[k][2] );     /* C B D is in [4] */
      strcpy( angle.t[1], angang.e[k][1] );
      strcpy( angle.t[2], angang.e[k][3] );
      angang.param[k][4] = CrossSearch(&angle, 1);

      strcpy( angle.t[0], angang.e[k][0] );     /* A B D is in [5] */
      strcpy( angle.t[1], angang.e[k][1] );
      strcpy( angle.t[2], angang.e[k][3] );
      angang.param[k][5] = CrossSearch(&angle, 1);
   }

 /* Expand the angleangle list */
   temp = angang.no_types; /* To hold the original number of types */
   /* The fist set of angang is 0 1 2 3, the 0 1 2 - 2 1 3 term */
   c = 0; sw1 = 2; sw2 = 3; /* The second is 0 1 3   3 1 2      */
   for ( i = 0; i < 2; i++)
   {
      for ( j = 0; j < temp; j++)
      {
	 strcpy( angang.e[angang.no_types][1], angang.e[j][1] );
         strcpy( angang.e[angang.no_types][c], angang.e[j][c] );
	 strcpy( angang.e[angang.no_types][sw1], angang.e[j][sw2] );
	 strcpy( angang.e[angang.no_types][sw2], angang.e[j][sw1] );
	 (angang.no_types)++;
      }
      c = 3; sw1 = 0; sw2 = 2; /* the final set is 2 1 0 - 0 1 3 */
   }

 /* Now search the .frc file for terms */

   SearchDynamicFrc(&angang, &ff_angang, 0, 0);

 /* 0 1 2 3 is identical to 3 1 2 0, so if angang.found[i] is zero, switch 
    these and repeat search */

   for ( i = 0; i < angang.no_types; i++ )
      if ( angang.found[i] == 0 )
      {
         strcpy ( tstring, angang.e[i][0] );
         strcpy ( angang.e[i][0], angang.e[i][3] );
         strcpy ( angang.e[i][3], tstring );
      }

   SearchDynamicFrc(&angang, &ff_angang, 0, 0);

 /* Contract angang struct to original size */ 

   for ( i = 0; i < temp; i++ )
   {
      angang.param[i][2] = angang.param[temp+i][0];
      angang.param[i][1] = angang.param[(2*temp)+i][0];
   }
   angang.no_types = temp;
   angang.no_params = 6;

/********************** ANGLE-ANGLE-TORSION COEFFICIENTS *********************/
	/* E = K * (Theta - Theta0) * (Theta' - Theta0') * (Phi - Phi1(0))
	   cffXY.frc supplies K(Ang,Ang,Tor)
	   LAMMPS expects     K(Ang,Ang,Tor) Theta0 Theta0'
	   Reference angles must be added
	*/
 
   printf("   searching angle-angle-torsions.....\n");
   angangtor.no_params = 1;
   SearchDynamicFrc(&angangtor, &ff_angangtor, 1, 0);

 /* Add reference angles */
   angangtor.no_params = 3;
   for ( k = 0; k < angangtor.no_types; k++ )
   {
      strcpy(angle.t[0], angangtor.e[k][0]);
      strcpy(angle.t[1], angangtor.e[k][1]);
      strcpy(angle.t[2], angangtor.e[k][2]);
      angangtor.param[k][1] = CrossSearch(&angle, 1);

      strcpy(angle.t[0], angangtor.e[k][1]);
      strcpy(angle.t[1], angangtor.e[k][2]);
      strcpy(angle.t[2], angangtor.e[k][3]);
      angangtor.param[k][2] = CrossSearch(&angle, 1);      
   }

/************************ END-BOND-TORSION COEFFICIENTS *********************/
        /* E = (R - R0) * SUM { V(n) * cos[n*phi] }
           cffXY.frc supplies F1 F2 F3 F1 F2 F3 
           LAMMPS expects     F1 F2 F3 F1 F2 F3 R0 R0' 
           Reference bonds must be added
        */

   printf("   searching end_bond-torsions.....\n");
   endbontor.no_params = 6;
   SearchDynamicFrc(&endbontor, &ff_endbontor, 1,0);

 /* Switch params if necessary and Add reference bonds */
   endbontor.no_params = 8;
   for ( k = 0; k < endbontor.no_types; k++ )
   {
      if (endbontor.found[k] == 2)
      {
         tparm = endbontor.param[k][0];
	 endbontor.param[k][0] = endbontor.param[k][3];
         endbontor.param[k][3] = tparm;

         tparm = endbontor.param[k][1];
         endbontor.param[k][1] = endbontor.param[k][4];
         endbontor.param[k][4] = tparm;

         tparm = endbontor.param[k][2];
         endbontor.param[k][2] = endbontor.param[k][5];
         endbontor.param[k][5] = tparm;
      }

      strcpy( bond.t[0], endbontor.e[k][0] );
      strcpy( bond.t[1], endbontor.e[k][1] );
      endbontor.param[k][6] = CrossSearch(&bond, 1);

      strcpy( bond.t[0], endbontor.e[k][2] );
      strcpy( bond.t[1], endbontor.e[k][3] );
      endbontor.param[k][7] = CrossSearch(&bond, 1);
   }

/************************ MID-BOND-TORSION COEFFICIENTS *********************/
        /* E = (R - R0) * 
		{ F(1) * cos(phi) + F(2) * cos(2 * phi) + F(3) * cos(3 * phi) }
           cffXY.frc supplies F1 F2 F3 
           LAMMPS expects     F1 F2 F3 R0
           Reference bond must be added
        */

   printf("   searching mid_bond-torsions.....  \n");
   midbontor.no_params = 3;
   SearchDynamicFrc(&midbontor, &ff_midbontor, 1,0);

  /* Add reference bond */
   midbontor.no_params = 4;
   for ( k = 0; k < midbontor.no_types; k++ )
   {
      strcpy( bond.t[0], midbontor.e[k][1] );
      strcpy( bond.t[1], midbontor.e[k][2] );
      midbontor.param[k][3] = CrossSearch(&bond, 1);
   }

/************************ BOND-BOND_1-3 COEFFICIENTS *************************/
	/* E = K(b,b') * (R - R0) * (R' - R0')
           cffXY.frc supplies K 
           LAMMPS expects    (NOT IN CLASS2.TXT) 
	   Reference bond length must be added for each bond 
        */

 printf("   searching bond-bond-1-3.....\n");
 bonbon13.no_params = 1;
 SearchDynamicFrc(&bonbon13, &ff_bonbon13, 1,0);

 /* Add reference bonds */

   bonbon13.no_params = 3;
   for ( k = 0; k < bonbon13.no_types; k++ )
   {
      strcpy(bond.t[0], bonbon13.e[k][0]);
      strcpy(bond.t[1], bonbon13.e[k][1]);
      bonbon13.param[k][1] = CrossSearch(&bond, 1);

      strcpy(bond.t[0], bonbon13.e[k][2]);
      strcpy(bond.t[1], bonbon13.e[k][3]);
      bonbon13.param[k][2] = CrossSearch(&bond, 1);

   }

/************************* ANGLE-TORSION COEFFICIENTS ************************/
	/* E = (Theta - Theta0) *
	    { F(1) * cos(phi)  +  F(2) * cos(2 * phi)  +  F(3) * cos(3 * phi) }
	   cffXY.frc supplies F1 F2 F3 F1 F2 F3 
	   LAMMPS expects     F1 F2 F3 F1 F2 F3 Theta0 Theta0' 
	   Reference angles must be added
	*/

   printf("   searching angle-torsions.....\n");
   angtor.no_params = 6;
   SearchDynamicFrc(&angtor, &ff_angtor, 1, 0);

 /* Switch if found in reverse search and Add reference angles */
   angtor.no_params = 8;
   for ( k = 0; k < angtor.no_types; k++ )
   {
      if (angtor.found[k] == 2)
      {
         tparm = angtor.param[k][0];
         angtor.param[k][0] = angtor.param[k][3];
         angtor.param[k][3] = tparm;

         tparm = angtor.param[k][1];
         angtor.param[k][1] = angtor.param[k][4];
         angtor.param[k][4] = tparm;

         tparm = angtor.param[k][2];
         angtor.param[k][2] = angtor.param[k][5];
         angtor.param[k][5] = tparm;
      }



      strcpy( angle.t[0], angtor.e[k][0] );
      strcpy( angle.t[1], angtor.e[k][1] );
      strcpy( angle.t[2], angtor.e[k][2] );
      angtor.param[k][6] = CrossSearch( &angle, 1);

      strcpy( angle.t[0], angtor.e[k][1] );
      strcpy( angle.t[1], angtor.e[k][2] );
      strcpy( angle.t[2], angtor.e[k][3] );
      angtor.param[k][7] = CrossSearch( &angle, 1);
   }

}/* End function*/


/************************** CrossSearch Function ***************************/

/* CrossSearch is a utility to search structure list->t for
   atoms in list->e and return parameters in vector xparam.
   This is done because elements in list->t come from
   another BondList variable. This allows us to combine data
   from the diagonal potential terms with cross terms.
*/
double CrossSearch(struct BondList *list, int reversible)
{
   int j,k,match,nom1;
   double xparam[MAX_PARAMS];
   nom1 = list->no_members - 1;
   for(k=0; k < list->no_types; k++)
   {
    /* Straight Match */
      for(j=match=0; j < list->no_members; j++)
         if( strcmp(list->e[k][j],list->t[j]) == 0) match++;
       if( match == list->no_members )
       {
          for(j=0; j < list->no_params; j++)
             xparam[j] = list->param[k][j];
          return xparam[0];
       }

    /* Match Reverse, Symmetric Sequence */
       if( reversible )
       {
          for(j=match=0; j < list->no_members; j++)
          {
             if( strcmp(list->e[k][nom1-j],list->t[j]) == 0)
                match++;
          }
          if( match == list->no_members )
          {
             for(j=0; j < list->no_params; j++)
                xparam[j] = list->param[k][j];
             return xparam[0];
          }
       }
   }
   printf("\n\nCould not find term in CrossSearch");

   for(j=0; j < list->no_members; j++)
      printf(" %s ",list->t[j]);

   printf("\n");
   exit(2);
 
   return -1; /* Put here to silence compiler/lint warnings */
}

