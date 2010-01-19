/**************************** ReadCoefficientsF1.c ****************************
*
*  This function uses the BondList structures to find the forcefield
*  coefficients that LAMMPS requires
*
*  This function is specific for Class I forcefields
*
*  All equations shown in comments are from the forcefield file
*/

#include "Msi2LMP2.h"
#include "Forcefield.h"

void ReadCoefficientsF1(void)
{
   int i,j,k; 		/* counters */
   int multiplicity;	/* for getting correct torsion constants */
   double eps, r;	/* Used in transforming vdw constants */

extern	void	SearchDynamicFrc(struct BondList*, struct FrcFieldItem*, int, int);
        void    SwitchOrderOfParams(struct BondList *list, int sub1, int sub2);
	int 	NoConnectionsToAtom(struct BondList *list, int n, int k);
 
/************************* NON-BOND COEFFICIENTS *************************/
        /* E = Aij/r^12 - Bij/r^6
	   where  Aij = sqrt( Ai * Aj )
		  Bij = sqrt( Bi * Bj )
           cvff.frc supplies A   in col1 and B in col2
           LAMMPS expects    eps in col1 and r  in col2
        */

	/* The CVFF forcefield file provides two parameters, A and B, for the 
           nonbond  coefficients.  LAMMPS requires these parameters in terms of
           epsilon and sigma (which it then converts back to A and B....)  
           So this conversion is made using the conversions: epsilon=(B*B)/(4*A)
                                                     and sigma = [(2*A)/B]^(1/6)           from the Discover manual

	In order to be compatible with lammps, we're using the conversions:
	epsilon = (B*B)/(4*A) 
	sigma =   [(A/B)^(1/6)] Note: no factor of 2*/

   printf("   searching nonbonds.....\n");
   vdw.no_params = 2; 
   SearchDynamicFrc(&vdw, &ff_vdw, 0, 0);

   for(j=0;j<vdw.no_types;j++)
      if( ( vdw.param[j][0] != 0.0 ) && ( vdw.param[j][1] != 0.0 ) )
      {
         eps = (vdw.param[j][1]*vdw.param[j][1])/(4.0*vdw.param[j][0]);
         r = pow((vdw.param[j][0]/vdw.param[j][1]),(1.0/6.0));
         vdw.param[j][0] = eps;
         vdw.param[j][1] = r; 
      }
      else
      {
         vdw.param[j][0] = 0.0; vdw.param[j][1] = 0.0;
      }


/***************************** BOND COEFFICIENTS *****************************/
	/* E = K2 * (R - R0)^2
	   cvff.frc supplies R0 in col1 and K2 in col2
 	   LAMMPS expects    K  in col1 and R  in col2
	   Switch function will be invoked
        */

   printf("   searching bonds.....\n");
   bond.no_params = 2; 
   SearchDynamicFrc(&bond, &ff_bond, 1,0);
   SwitchOrderOfParams(&bond, 0, 1);

/**************************** ANGLE COEFFICIENTS *****************************/
	/* E = K2 * (Theta - Theta0)^2
           cvff.frc supplies Theta0 in col1 and K2 in col2
           LAMMPS expects    K  in col1 and Theta0 in col2
           Switch function will be invoked
        */

   printf("   searching angles.....\n");
   angle.no_params = 2; 
   SearchDynamicFrc(&angle, &ff_ang, 1,0);
   SwitchOrderOfParams(&angle, 0, 1);

/************************** TORSION COEFFICIENTS *****************************/
	/* E = Kphi * [ 1 + cos(n*Phi - Phi0) ]
           cvff.frc supplies Kphi in col1, n in col2, and Phi0 in col3
           LAMMPS expects    K    in col1, d in col2, and n    in col3
		d is 1 or -1 corresponding to phi0 of 0 or 180
           Switch function will be invoked to switch col 2 and 3
	   Phi0 will be converted to d 

	   CVFF torsions also have to deal with multiplicity.  If wildcards 
	   were used in finding the frc parameters, the K coefficient in the 
	   torsion must be scaled down by a factor realted to the number of 
	   other connections to the two central atoms.
	*/
 
   printf("   searching torsions.....\n");
   torsion.no_params = 3; 
   SearchDynamicFrc(&torsion, &ff_tor, 1, 0);
   SwitchOrderOfParams(&torsion, 1, 2);

   for ( i = 0; i < torsion.count; i++ )
      if ( torsion.param[i][1] == 0.0 )
         torsion.param[i][1] = 1.0;
      else if ( torsion.param[i][1] == 180.0 )
         torsion.param[i][1] = -1.0;
      else exit (8);

    /* modify torsion for multiplicity */
   for ( k = 0; k < torsion.no_types; k++)
   {
      if (torsion.wc[k])
      {
         multiplicity =( NoConnectionsToAtom(&torsion, 1, k) -1);
         multiplicity *= ( NoConnectionsToAtom(&torsion, 2, k) - 1);
         torsion.param[k][0] /= ((double) multiplicity);
      }
    }


/********************** OUT-OF-PLANE COEFFICIENTS ***************************/
        /* E = Kchi * [ 1 + cos(n*Chi - Chi0) ]
           cvff.frc supplies Kchi in col1, n in col2, and Chi0 in col3
           LAMMPS expects         in col1,   in col2, and      in col3
                d is 1 or -1 corresponding to phi0 of 0 or 180
           Switch function will be invoked to switch col 2 and 3
           Phi0 will be converted to d

	  Out of plane terms are treated as improper torsions
        */

   printf("   searching out-of-planes.....\n");
   oop.no_params = 3;

   SearchDynamicFrc(&oop, &ff_oop, 0, 1);

  /*Assuming OOPS to be like torsion, spaced as k, chi changed to +-1, n*/

   SwitchOrderOfParams(&oop, 1, 2);

   for ( i = 0; i < oop.count; i++ )
      if ( oop.param[i][1] == 0.0 )
         oop.param[i][1] = 1.0;
      else if ( oop.param[i][1] == 180.0 )
         oop.param[i][1] = -1.0;
      else exit (8);


}/* End of ForceFieldStuff() function*/




/*-----------------------Switch Order of Params Function --------------------*/

void SwitchOrderOfParams(struct BondList *list, int sub1, int sub2)
{
   int i;
   double temp;
   for ( i = 0; i < list->count; i++)
   {
      temp = list->param[i][sub1];
      list->param[i][sub1] = list->param[i][sub2];
      list->param[i][sub2] = temp;
   }

}


/*----------------------No Connections to Atom Function ----------------------*/
int NoConnectionsToAtom(struct BondList *list, int n, int k)
{
   int i;

   for ( i = 0; i < list->count; i++)
      if (strcmp(list->type[i][n], list->a[k][n])==0)
         return(atoms[list->id[i][n]].no_connect);
   
   printf("Torsion type not found\n");
   return(-1);
}

