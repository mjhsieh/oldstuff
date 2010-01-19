/***************************** SearchAndFill.c ****************************
*
* This function first allocates memory to the forcefield item
* structures and then reads parameters from the forcefield file into the
* allocated memory
*
*/

#include "Forcefield.h"

void SearchAndFill(struct FrcFieldItem *item)
{
   int i,j; /* counters */
   int ctr = 0;
   char line[MAX_LINE] = "empty";
   char *charptr;
   extern FILE *FrcF;

/*************************ALLOCATE MEMORY FOR STRUCTURE ********************/

 /* Read and discard lines until keyword is found */

   rewind(FrcF);

   while (strncmp(line, item->keyword, strlen(item->keyword)) !=0 )
      fgets( line, MAX_LINE, FrcF );

 /* Count the number of lines until next item is found */

   while( strncmp(fgets(line,MAX_LINE,FrcF), "#", 1) != 0 )
      ctr++;

 /* Allocate the memory using calloc */

   item->data = calloc(ctr, sizeof(struct FrcFieldData));

   if(item->data == NULL)
   {
      printf("Could not allocate memory to %s\n", item->keyword);
      exit(2);
   }

/**********************FILL PARAMETERS AND EQUIVALENCES ********************/

 /* Read lines until keyword is found */

   rewind(FrcF);

   strcpy(line,"empty");
   ctr = 0; /* reset counter to 0 */

   while (strncmp(line, item->keyword, strlen(item->keyword)) != 0 )
     fgets( line, MAX_LINE, FrcF );

 /* Read lines until data starts (when !--- is found) */

   while ( strncmp(line, "!---", 4) != 0 )
      fgets(line, MAX_LINE, FrcF);

   fgets(line, MAX_LINE, FrcF); /* Get first line of data */

 /* Read data into structure */

   while( strncmp( line, "#", 1 ) != 0 )
   {
    /* version number and reference number */

      charptr = strtok(line, " ");
      item->data[ctr].ver = atof(charptr);

      charptr = strtok(NULL, " ");
      item->data[ctr].ref = atoi(charptr);

    /* equivalences */

      for( i = 0; i < item->number_of_members; i++ )
      {
         charptr = strtok(NULL, " ");
         sscanf(charptr, "%s", item->data[ctr].ff_equiv[i]);
      }

    /* parameters -- Because of symmetrical terms, bonang, angtor, and
                     endbontor have to be treated carefully */

         for( i = 0; i < item->number_of_parameters; i++ )
         {
            charptr = strtok(NULL, " ");

            if(charptr == NULL)
            {
               for( j = i; j < item->number_of_parameters; j++ )
                  item->data[ctr].ff_param[j] = item->data[ctr].ff_param[j-i];
               
               break;
            }
            else
               item->data[ctr].ff_param[i] = atof(charptr);
            
         }

      fgets( line, MAX_LINE, FrcF);

      while(strcmp(line,"\n") == 0) /* if blank line encountered, get next */
         fgets( line, MAX_LINE, FrcF );

      ctr++;
   }
   item->entries = ctr;
/*Debugging
 printf("\n%s\n", item->keyword);
 for(i=0;i<ctr;i++)
 {
  for(j=0;j<item->number_of_members;j++)
   printf("%3s ", item->data[i].ff_equiv[j]);
  printf("     ");
  for(j=0;j<item->number_of_parameters;j++)
   printf("%10.5f ", item->data[i].ff_param[j]);
  printf("\n");
 }
*/
}

