/****************************** ReadMdfFile.c ******************************
*
*  This function opens the .mdf file and extracts connectivity information
*  into the atoms Atom structure.  It also updates the charge from the .car
*  file because the charge in the .mdf file has more significant figures.
*
*/

#include "Msi2LMP2.h"

/* Prototype for function to process a single atom 
   Returns int that flags end of data file */
int ProcessConnections(char line[], int connect_col_no, 
                       int q_col_no, int *counter);
 
/* Prototype for function that takes connectivty record as stated in 
   .mdf file and fills in any default values */
void MakeConnectFullForm(int *counter);

/* Prototype for function that assigns atom numbers to the connection names */
extern void LinkConnections(void);

void ReadMdfFile(void)
{
   char line[MAX_LINE_LENGTH];	/* Temporary storage for reading lines */
   char *col_no;		/* Pointer to column number stored as char */
   char *col_name;		/* Pointer to column name */
   int connect_col_no = 0;	/* Column number where connection info begins */
   int q_col_no = 0;		/* Column number containg charge information */
   int column_flag=0;		/* Flag for finding connect and q col no */
   int atom_counter=0;		/* Keeps track of current atom number */
   int p_flag = 1; 		/* return value from ProcessConnections() */

/* Open .mdf file for reading */

   sprintf(line,"%s.mdf",rootname);
   printf("   Molecular Data file name is %s\n", line);
   if((MdfF = fopen(line,"r")) == NULL )
   {
      printf("Cannot open %s\n",line);
      exit(2);
   }

/* Search for column containing keyword 'charge' 
     to be able to set loops in processing function*/

   while(!column_flag)
   {
      fgets(line, MAX_LINE_LENGTH, MdfF);
      if (strncmp((strtok(line," ") ),"@column",7) == 0)
      {
         col_no = strtok(NULL, " ");
         col_name = strtok(NULL, " ");
         if ((! strncmp (col_name,"charge",7)))
                        /* Be careful not to get charge_group parameter! */
         {
            q_col_no = atoi(col_no);
            column_flag=1;
         }
      }
   } /* end while loop */

/* Search for column containing keyword 'connect'
     to be able to set loops in processing function*/

   /* Version 3 .mdf files use keyword 'connections' while Ver 4 uses
	'connectivity', although neither is explicitly stated in 
	"Common File Formats".  

	Also Ver 4 is supposed to conatin a 'n_connections' parameter
	that would eliminate the counting loop in the processing function,
	however, there is no guarantee that this column header is used.
   */

   column_flag = 0;
   while(!column_flag)
   {
      fgets(line, MAX_LINE_LENGTH, MdfF);
      if (!strncmp((strtok(line," ")),"@column",7))
      {
         col_no = strtok(NULL, " ");
         col_name = strtok(NULL, " ");
         if ((! strncmp (col_name,"connect",7)))
         {
            connect_col_no = atoi(col_no);
            column_flag=1;
         }
      }
   }   /* end while loop */


/* Process molecules loop.  Each function call processes one molecule
   Loop continues until next section (marked by#) is reached or
    a false return value is returned, indicating end of file*/

   while ( (strncmp(line,"#",1) != 0 ) && (p_flag != 0) )
     p_flag = ProcessConnections(line, connect_col_no, q_col_no, &atom_counter);

/* Assign atom names in connections[] to corresponding atom numbers */

   printf("\nAssigning numbers to connection names.....\n");
   LinkConnections();

/* Close .mdf file */

   if (fclose(MdfF) !=0)
   {
      printf("Error closing %s.car\n", rootname);
      exit(1);
   }

} /* End ReadMdfFile function */




/*---------------------Process Connections Function-----------------------*/

int ProcessConnections(char line[], int connect_col_no, int q_col_no,
			 int *counter)
{
   char *cur_field;	/* For storing current string token */
   int i; 		/* Used in loop counters */
   int connect_no;      /* Connection number within atom */
   int r_val = 1;	/* Return value.  1 = successful
					  0 = EOF encountered */

 /* Search for keyword '@molecule' */

   while( strncmp(line, "@molecule", 9))
     if (fgets(line, MAX_LINE_LENGTH, MdfF) == NULL) 
	exit(3); 

   fgets(line, MAX_LINE_LENGTH, MdfF);
   while (strcmp(line, "\n") == 0) /*  if blank line, skip it and get next */
      fgets(line, MAX_LINE_LENGTH, MdfF);

 /* Processing loop for an atom */

   while(strcmp(line, "\n")) /* signals end of atom record */
   {
      /* Process atom name */

      cur_field = strtok(line,"_");
      sscanf(cur_field, "%s", atoms[*counter].res_name);
      cur_field = strtok(NULL,":");
      sscanf(cur_field, "%s", atoms[*counter].res_num);
      cur_field = strtok(NULL," "); /* Atom name,
					 compare with .car file */
      if (strcmp(atoms[*counter].name, cur_field))
      {
         printf("Names %s from .car file and %s from .mdf file do not match\n",
                                              atoms[*counter].name, cur_field);
         printf("Program Terminating\n");
         exit(4);
      }

  /* Skip unwanted fields until charge column, then update charge */

      for (i=1; i < q_col_no; i++)
        strtok(NULL," ");

      cur_field = strtok(NULL, " ");
      atoms[*counter].q = atof(cur_field);

  /* Continue skipping unwanted fields until connectivity records begin */

      for ( i = (q_col_no + 1); i < connect_col_no; i++)
	strtok(NULL," ");

   /* Process connections */

      connect_no = 0; /* reset connections counter */

      for (i=0; i < MAX_CONNECTIONS; i++)
      {
        cur_field = strtok(NULL, " ");
        sscanf(cur_field, "%s", atoms[*counter].connections[connect_no]);
        if (strcmp(cur_field, "\n")==0)
          break; /* Breaks from atom for loop, does not increase connect_no */
        connect_no++;
      }

   atoms[*counter].no_connect = connect_no;
   MakeConnectFullForm(counter);
   fgets(line, MAX_LINE_LENGTH, MdfF);
   (*counter)++;

   } /* End atom processing loop */
   
/* Read lines until next molecule record is found
   If end of data is encountered, set flag to false and break */

   while( (strncmp(line,"@molecule", 9)) )
     if  (fgets(line, MAX_LINE_LENGTH, MdfF) == NULL) 
     {
       r_val = 0; /* Set flag to false -- EOF has been reached */
       break;
     }

   return r_val;

} /* End ProcessConnections function */




/*-------------------------MakeConnectFullForm Function--------------------*/

void MakeConnectFullForm(int *counter)
{
/* This function processes the connection names after all connections
   for an atom have been read in.
   It replaces any short forms that use implied default values
   with the full form connectivity record */

   int i;  		/* Counter for character array */
   int j;		/* loop counter */
   char tempname[25];   /* name of connection */
   char tempcell[10];   /* Values from connectivity record */
   char tempsym[5];              /* " " */
   char tempbo[6];               /* " " */
   char *charptr;

   for ( j = 0; j < atoms[*counter].no_connect; j++)
   {
  /* If not full name, make name full */
      if (strchr(atoms[*counter].connections[j],':') == NULL)
      {
         strcpy(tempname,atoms[*counter].res_name);
         strcat(tempname,"_");
         strcat(tempname, atoms[*counter].res_num);
         strcat(tempname,":");
         strcat(tempname, atoms[*counter].connections[j]);
         sscanf(tempname, "%s", atoms[*counter].connections[j]);
      }
      else
         sscanf(atoms[*counter].connections[j], "%s", tempname);
  /* Set cell variables */

   i=0;
   charptr = (strchr(tempname,'%'));
   if (charptr != NULL)
   {
      while ( *charptr!='#' && *charptr!='/' && *charptr!='\000')
         tempcell[i++] = *(charptr++);
      tempcell[i] = '\000';
   }
   else
      strcpy(tempcell, "%000");

  /* Set symmetry variables
	-- If not 1, cannot handle at this time */

   i = 0;
   charptr = (strchr(tempname,'#'));
   if (charptr != NULL)
   {
       while ( *charptr != '/' && *charptr !='\000')
       {
          tempsym[i++] = *(charptr++);
          if ((i==2) && (tempsym[1] != '1'))
          {
             printf("Msi2LMP is not equipped to handle symmetry operations\n");
             exit(5);
          }
       }
       tempsym[i] = '\000';
   }
   else
      strcpy(tempsym, "#1");
 
  /* Set bond order and record in data structure */

   i = 0;
   charptr = strchr(tempname,'/');
   if (charptr != NULL)
   {
      charptr++;
      while (*charptr != '\000')
         tempbo[i++] = *(charptr++);
      tempbo[i] = '\000';
   }
   else
      strcpy(tempbo, "1.0");
   atoms[*counter].bond_order[j] = atof(tempbo);
  /* Build connection name and store in atoms data structure */

    strtok( tempname, "%#/");
    strcat( tempname, tempcell);
    strcat( tempname, tempsym);
    strcat( tempname, "/");
    strcat( tempname, tempbo);
    sscanf( tempname, "%s", atoms[*counter].connections[j]);

 }/*End for loop*/

}/* End function MakeNameLong */


