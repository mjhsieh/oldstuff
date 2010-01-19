/******************************** Msi2LMP2.c ********************************
*
*  This program takes the .car and .mdf files from MSI/Biosyms's INSIGHT
*  program and produces a file that can be used as the input to 
*  the LAMMPS molecular dynamics program
*
*  The program may be started by supplying information at the command prompt
*  according to the usage described below.  The program may also be started
*  by calling the executable and being prompted for information
*
* USAGE:  msi2lmp2 rootname {f#}
*  -- msi2lmp2 is the name of the executable
*  -- rootname is the base name of the .car and .mdf files
*  -- # is the class of forcefield to use (1 = Class I CVFF)
*					  (2 = Class II CFF91)
*					  (3 = Class II CFF95)
*  -- output of program is a data.root file for LAMMPS input 
*/

#define MAIN

#include "Msi2LMP2.h"

int main (int argc, char *argv[])
{
   int n;		/* Counter */

 /* Functions called from within main */
 /* All code is located in .c file with function name */
   extern void FrcMenu();
   extern void ReadCarFile();
   extern void ReadMdfFile();
   extern void ReadFrcFile();
   extern void MakeLists();
   extern void ReadCoefficientsF1();
   extern void ReadCoefficientsF2();
   extern void BuildLMP();
   extern void CreateOutputReport();

   forcefield = 0;		/* Variable that identifies forcefield to use */ 
   puts("\n\n\nRunning Msi2LMP.....\n\n");

   if (argc < 2) /* If no rootname was supplied, prompt for it */
   {
     printf("Please enter the rootname of the .car and .mdf files: ");
     gets(rootname);
   }
   else /* rootname was supplied as first argument, copy to rootname */
      sprintf(rootname,"%s",argv[1]);

   for ( n = 2; n < argc; n++)	/* Process command line arguments */
   {
      switch ( argv[n][0] )
      {
       case 'f': 
         switch ( argv[n][1] ) 
         {
            case '1':
               forcefield = 1;
               strcpy(FrcFileName, "cvff.frc");
               printf("Using cvff.frc Class I forcefield\n\n");
	       break;

            case '2':
               forcefield = 2;
               strcpy(FrcFileName, "cff91.frc");
               printf("Using cff91.frc Class II forcefield\n\n");	
               break;

	    case '3':
	       forcefield = 3;
               strcpy(FrcFileName, "cff95.frc");
               printf("Using cff95.frc Class II forcefield\n\n");
               break;

            default:
               forcefield = 0; /* menu will be called */
         }
         break;
      }
   }

   if ( forcefield == 0 ) /* If there were no command line arguments except for
		      the rootname, or argument was invalid, give user force 
		      field menu */
      FrcMenu();

   if (forcefield == 0 ) /* Just in case, default forcefield is set */
   {
      forcefield = 3;
      printf("No forcefield specified, using CFF95 Class 2 forcefield\n\n");
      strcpy (FrcFileName, "cff95.frc");
   }

 /* Read in .car file */

   printf("Processing .car file......\n");
   ReadCarFile(); 

 /*Read in .mdf file */

  printf("Processing .mdf file......\n");
  ReadMdfFile();

 /* Read .frc file into memory */

   printf("\nReading %s into memory\n", FrcFileName);
   ReadFrcFile();

 /* Define bonds, angles, etc...*/

   printf("\nBuilding lists......\n\n");
   MakeLists();

 /* Fill in forcefield parameters */

   printf("Finding forcefield parameters\n");
   if ( forcefield >= 2 )
      ReadCoefficientsF2();
   else if ( forcefield == 1)
      ReadCoefficientsF1();

 /* Create ouput files */

   printf("\nCreating LMP input file data.%s\n\n",rootname);
   BuildLMP();

   printf("\nCreating Output Report\n\n");
   CreateOutputReport();

   printf("\nNormal program termination\n");
   return(0);
}
