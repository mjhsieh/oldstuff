/****************************** LinkConnections.c *****************************
*
*   This routine takes the longhand connection names processed from .mdf files
*   and converts them to the integer values that match the atom number of the
*   corresponding connection.
*
*/

#include "Msi2LMP2.h"

void LinkConnections()
{
   int mol_ctr;
   int atom_ctr = 0;
   int conn_ctr;
   int begin_atom;		/* marker to beginning atom w/in a molecule*/
   int k;		/* index of atoms structure used in comparison loops*/
   int match = 0; 		/* flag set to true when match assigned */
   char connection[MAX_STRING]; /* connection name token */
   char comp_atom[MAX_STRING];	/* atom to which connection is compared */

/* Outer loop - Each pass completes a molecule */

   for ( mol_ctr = 0; mol_ctr < no_molecules; mol_ctr++)
   {
     begin_atom = atom_ctr; /* Set marker to beginning atom number
				of molecule */

    /* Inner loop - Each pass completes an atom */

     while ( atoms[atom_ctr].molecule == mol_ctr && atom_ctr<total_no_atoms) /* While within atom */
     {
      /* Third nested loop - Each pass completes a connection */

   	for ( conn_ctr = 0; conn_ctr < atoms[atom_ctr].no_connect; conn_ctr++)
	{
	  strcpy ( connection, strtok( atoms[atom_ctr].connections[conn_ctr], "%"));
	  match = 0;
	  k = begin_atom;
	  while (!match)
          {
           strcpy(comp_atom, atoms[k].res_name);
	   strcat(comp_atom, "_");
	   strcat(comp_atom, atoms[k].res_num);
	   strcat(comp_atom, ":");
	   strcat(comp_atom, atoms[k].name);

	/* compare names.  if they match, make assignment and break loop
		going on to next atom.  else increment atom and try again */

	   if (strcmp(connection,comp_atom) == 0)
	   {
	     atoms[atom_ctr].conn_no[conn_ctr] = atoms[k].no;
	     match = 1;
	   }
	   else
	     k++;

	   if ( atoms[atom_ctr].molecule != atoms[k].molecule)
	    {
printf("atom ctr is %d, conn no %d,  while k is %d\n", atom_ctr, conn_ctr, k);
printf("atom %d is %s_%s:%s, connection %d is %s\n", atom_ctr, atoms[atom_ctr].res_name, atoms[atom_ctr].res_num, atoms[atom_ctr].name, conn_ctr, connection);
		printf("Overstepped molecule without finding connection\n");
		exit(2);
	    }

	  } /* End while (comparison) loop */

	} /* End connection loop */

        atom_ctr++;
 
     } /* End atom loop */

   } /* End molecule loop */

} /* End LinkConnections function */
