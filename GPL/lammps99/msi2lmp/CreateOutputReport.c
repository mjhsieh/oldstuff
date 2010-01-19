/******************** CreateOutputReport.c ********************
*
*     Prints an output report listing all parameters and lists
*     any unassigned parameters
*
*/

#include "Msi2LMP2.h"

void CreateOutputReport(void)
{
   int k,x;
   char line[MAX_LINE_LENGTH];

   void PrintParameters(FILE *Fp, char *title, struct BondList *list);
   int MakeUnfoundList(struct BondList*, char *);

 /* Open report file */

   sprintf(line, "report.%s", rootname);
   printf("Report file is named %s\n", line);
   if ( (RptF = fopen(line, "w")) == NULL )
   {
      printf("Cannot open %s\n", line);
      exit(2);
   }

 /***** Parameters *****/

   fprintf(RptF,"There are %i molecules\n",no_molecules);
   for(k=0; k < no_molecules; k++)
      fprintf(RptF,"There are %3i atoms in molecule %i\n",no_atoms[k],k+1);

   fputs("\n------------------------------------------------------\n",RptF);
   fputs("\nFor each type of Interaction:\n",RptF);
   fputs("Column  1   = type number counter\n",RptF);
   fputs("Column  2   = 1 if parameters assigned\n",RptF);
   fputs("              0 otherwise\n",RptF);
   fputs("Column  3   = atom potential types: actual_equivalence\n",RptF);
   fputs("Columns 4-$ = potential parameters assigned\n",RptF);
   fputs("------------------------------------------------------\n\n",RptF);

   PrintParameters(RptF,"NonBonds",&vdw);
   PrintParameters(RptF,"Bonds",&bond);
   PrintParameters(RptF,"Angles",&angle);
   PrintParameters(RptF,"Dihedrals",&torsion);
   PrintParameters(RptF,"Out-Of-Planes",&oop);
   PrintParameters(RptF,"Bond-Bond",&bonbon);
   PrintParameters(RptF,"Bond-Bond_1_3",&bonbon13);
   PrintParameters(RptF,"Bond-Angle",&bonang);
   PrintParameters(RptF,"Angle-Angle",&angang);
   PrintParameters(RptF,"Angle-Torsion",&angtor);
   PrintParameters(RptF,"EndBond-Torsion",&endbontor);
   PrintParameters(RptF,"MidBond-Torsion",&midbontor);
   PrintParameters(RptF,"Angle-Angle-Torsion",&angangtor);

 /***** Unfound Parameters *****/

   MakeUnfoundList(&vdw, "Non-bond");
   MakeUnfoundList(&bond, "Bond");
   MakeUnfoundList(&angle, "Angle");
   MakeUnfoundList(&torsion, "Torsion");
   MakeUnfoundList(&oop, "Out of Plane");
   MakeUnfoundList(&bonbon, "Bond-bond");
   MakeUnfoundList(&bonang, "Bond-Angle");
   MakeUnfoundList(&angangtor, "Angle-Angle-Torsion");
   MakeUnfoundList(&angtor, "Angle-Torsion");
   MakeUnfoundList(&endbontor, "EndBondTorsion");
   MakeUnfoundList(&midbontor, "MidBondTorsion");
   MakeUnfoundList(&bonbon13, "BondBond13");
   x = MakeUnfoundList(&angang, "Angle-angle");
      
   fprintf(RptF, "\nThere are %d unassigned parameters\n", x);

   fclose(RptF);
}

/************************* Print Parameters Function *************************/
void PrintParameters(FILE *Fp, char *title, struct BondList *list)
{
 int j,k;

 fprintf(Fp,"%s   count= %i\n",title,list->count);
 for(k=0; k < list->no_types; k++)
 {
    fprintf(Fp,"%3i %i",k+1,list->found[k]);
    for(j=0; j < list->no_members; j++)
       fprintf(Fp,"%5s_%-5s ",list->a[k][j],list->e[k][j]);
    for(j=0; j < list->no_params; j++)
       fprintf(Fp,"%16.10g ",list->param[k][j]);
    fprintf(Fp,"\n");
 }
 fprintf(Fp,"\n");

}

/************************* MakeUnfoundList Function **************************/

int MakeUnfoundList(struct BondList *list, char *string)
{
   int i;
   static int unfound = 0;

   for ( i = 0; i < list->no_types; i++)
   {
      if ( list->found[i] == 0 )
      {
         fprintf(RptF, "%s %d was not found in the frc file\n", string, i+1);
         unfound++;
      }
   }
return unfound;
}
