/****************************** ReadFrcFile.c ******************************
*
*   This routine reads the data from a .frc forcefield file and stores it in
*   dynamically allocated memory.  This allows for fast searches of the
*   file.
*
*/

#define FF_MAIN

#include "Forcefield.h"
#include "Msi2LMP2.h"

void ReadFrcFile(void)
{
   extern void InitializeItems(void);
   extern void SearchAndFill(struct FrcFieldItem *item);

  /* Open Forcefield File */

   if ( (FrcF = fopen(FrcFileName,"r")) == NULL )
   {
      printf("Cannot open %s\n", FrcFileName);
      exit(2);
   }

   InitializeItems(); /* sets keywords, number of members and 
                         number of parameters for each structure */

  /* allocate memory to and search and fill each structure */

   SearchAndFill(&atom_types);
   SearchAndFill(&equivalence);
   SearchAndFill(&ff_vdw);
   SearchAndFill(&ff_bond);
   SearchAndFill(&ff_ang);
   SearchAndFill(&ff_tor);
   SearchAndFill(&ff_oop);

   if (forcefield != 1)  /* Skip cross terms for class I */
   {
      SearchAndFill(&ff_bonbon);
      SearchAndFill(&ff_bonang);
      SearchAndFill(&ff_angtor);
      SearchAndFill(&ff_angangtor);
      SearchAndFill(&ff_endbontor);
      SearchAndFill(&ff_midbontor);
      SearchAndFill(&ff_angang);
      SearchAndFill(&ff_bonbon13);
   }   
   
   fclose(FrcF);
}

