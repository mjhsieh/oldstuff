/****************************** BuildLMP.c ******************************
*
*  This function creates and writes the data file to be used with LAMMPS
*
*/

#include "Msi2LMP2.h"
#include "Forcefield.h"

void BuildLMP()
{
   int j,k;
   char line[MAX_LINE_LENGTH];

   void PrintOutput(struct BondList *list, char tstring[]);

 /* Open file for output */

   sprintf(line, "data.%s", rootname);
   if ((DatF = fopen(line, "w")) == NULL )
   {
      printf("Cannot open %s for output\n", line);
      exit(2);
   }
   printf("Writing LAMMPS molecular data file...\n");

   fprintf(DatF, "LAMMPS data file for %s\n\n", rootname);
   fprintf(DatF, "%d atoms\n", total_no_atoms);
   fprintf(DatF, "%d bonds\n", bond.count);
   fprintf(DatF, "%d angles\n",angle.count);
   fprintf(DatF, "%d dihedrals\n", torsion.count);
   fprintf(DatF, "%d impropers\n", oop.count);

   fprintf(DatF, "\n");

   fprintf(DatF, "%d atom types\n", vdw.no_types);
   if ( bond.count )   fprintf(DatF, "%d bond types\n", bond.no_types);
   if ( angle.count )  fprintf(DatF, "%d angle types\n", angle.no_types); 
   if ( torsion.count) fprintf  (DatF, "%d dihedral types\n", torsion.no_types);
   if ( oop.count )    fprintf  (DatF, "%d improper types\n", oop.no_types);

   fprintf(DatF, "\n");
   fprintf(DatF, "%15f %15f xlo xhi\n", pbc[6], pbc[0]);
   fprintf(DatF, "%15f %15f ylo yhi\n", pbc[7], pbc[1]);
   fprintf(DatF, "%15f %15f zlo zhi\n", pbc[8], pbc[2]);

  
 /* MASSES */

   fprintf(DatF, "\nMasses\n\n");
   for ( k = 0; k < vdw.no_types; k ++)
      for( j = 0; j < atom_types.entries; j++) 
         if( strcmp( vdw.a[k][0], atom_types.data[j].ff_equiv[0]) == 0 ) 
         {
            vdw.mass[k] = atom_types.data[j].ff_param[0];
            break;
         } 

   for(k=0; k < vdw.no_types; k++)
      fprintf(DatF, "%3d %f\n",k+1,vdw.mass[k]);

   fprintf(DatF, "\n");

 /* COEFFICIENTS */

   PrintOutput(&vdw,       "Nonbond");
   PrintOutput(&bond,      "Bond");
   PrintOutput(&angle,     "Angle");
   PrintOutput(&torsion,   "Dihedral");
   PrintOutput(&oop,       "Improper");
   PrintOutput(&bonbon,    "BondBond");
   PrintOutput(&bonang,    "BondAngle");
   PrintOutput(&angang,    "AngleAngle");
   PrintOutput(&angangtor, "AngleAngleTorsion");
   PrintOutput(&endbontor, "EndBondTorsion");
   PrintOutput(&midbontor, "MiddleBondTorsion");
   PrintOutput(&bonbon13,  "BondBond13");
   PrintOutput(&angtor,    "AngleTorsion");

/*----------------------------------------------------------------------*/
 /* ATOMS */

   fprintf(DatF, "Atoms\n\n");
   for(k=0; k < total_no_atoms; k++)
   {

      fprintf(DatF, "%i %i %i %10f %15.9f %15.9f %15.9f\n",
              k+1,
              atoms[k].molecule,
              vdw.type_no[k]+1,
              atoms[k].q,
              atoms[k].x[0],
              atoms[k].x[1],
              atoms[k].x[2]);

   }
   fprintf(DatF, "\n");

 /***** BONDS *****/

   if( bond.count )
   {
      fprintf(DatF, "Bonds\n\n");
      for(k=0; k < bond.count; k++)
         fprintf(DatF, "%i %i %i %i\n",k+1,bond.type_no[k]+1,bond.id[k][0]+1,
							    bond.id[k][1] + 1);
      fprintf(DatF, "\n");
   }

 /***** ANGLES *****/

   if( angle.count )
   {
      fprintf(DatF, "Angles\n\n");
      for(k=0; k < angle.count; k++)
         fprintf(DatF, "%i %i %i %i %i\n",k+1,angle.type_no[k]+1,
                         angle.id[k][0]+1, angle.id[k][1]+1, angle.id[k][2]+1);

   fprintf(DatF, "\n");
   }


 /***** TORSIONS *****/

   if( torsion.count )
   {
      fprintf(DatF, "Dihedrals\n\n");
      for(k=0; k<torsion.count; k++)
         fprintf(DatF, "%i %i %i %i %i %i\n",k+1, torsion.type_no[k]+1,
			torsion.id[k][0]+1, torsion.id[k][1]+1, 
			torsion.id[k][2]+1, torsion.id[k][3]+1);
      fprintf(DatF, "\n");
   }

 /***** OUT-OF-PLANES *****/
   if(oop.count) 
   {
      fprintf(DatF, "Impropers\n\n");
      for (k=0; k<oop.count; k++)
         fprintf(DatF, "%i %i %i %i %i %i \n", k+1, oop.type_no[k]+1, 
                oop.id[k][0]+1, oop.id[k][1]+1,oop.id[k][2]+1, oop.id[k][3]+1);
      fprintf(DatF, "\n");
   }
   fclose(DatF);
}

/*************************** PrintOutput Function ****************************/
void PrintOutput(struct BondList *list, char tstring[])
{
   int j,k;
   if ( list->count)
   {
      fprintf(DatF, "%s Coeffs\n\n", tstring);
      for ( k = 0; k < list->no_types; k++)
      {
         fprintf(DatF, "%i ", k+1);
         for ( j = 0; j < list->no_params; j++)
	    fprintf(DatF, "%20.10f ", list->param[k][j]);
	 fprintf(DatF, "\n");
      }
      fprintf(DatF, "\n");
   }
}
