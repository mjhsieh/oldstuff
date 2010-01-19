/***************************** SearchDynamicFrc.c *****************************
*
*   This function replaces the search function that searches the frc file.
*   Since the .frc is in memory, this function searches the data structures
*   to fill the parameters in the BondList structures
*
*/

#include "Msi2LMP2.h"
#include "Forcefield.h"

void SearchDynamicFrc(struct BondList *list, struct FrcFieldItem *item, 
						int reverse, int permute )
{
   int i, j, k, m, match;
   int r_match, nom1;

   int Permut023Match(struct BondList *list, struct FrcFieldItem *item, int k, 
						int j,int i0,int i2,int i3);

   for( k = 0; k < list->no_types; k++)
   {
      for ( i = 0; i < item->entries; i++ )
      {

     /* STRAIGHT MATCH */

         list->wc[k]=0;
         for ( j = match = 0; j < list->no_members; j++)
         {
            if( strcmp (list->e[k][j], item->data[i].ff_equiv[j]) == 0) 
	       match++;
            if( strcmp ("*", item->data[i].ff_equiv[j])==0)
	    {
	       match++;
	       list->wc[k]++;
	    }
         } 
         
         if( match == list->no_members )
         {
            for ( m = 0; m < list->no_params; m++ )
               list->param[k][m] = item->data[i].ff_param[m];
            list->found[k]=1;
            break;
         }

     /* REVERSE MATCH */

       if(reverse)
       {
         list->wc[k] = 0;
         nom1 = list->no_members-1;
         for ( j = r_match = 0; j < list->no_members; j++)
         {
            if( strcmp (list->e[k][nom1-j], item->data[i].ff_equiv[j]) == 0) 
               r_match++;
            if( strcmp ("*", item->data[i].ff_equiv[j])==0)
            {
               r_match++; 
	       list->wc[k]++;
	    }
         } 

         if( r_match == list->no_members )
         {
            for ( m = 0; m < list->no_params; m++ )
               list->param[k][m] = item->data[i].ff_param[m];
            list->found[k]=2;
            break;
         }
       }

    /* PERMUTE MATCH */

       if(permute)
       {
         if( (Permut023Match(list,item,k,i,0,2,3) == list->no_members) ||
             (Permut023Match(list,item,k,i,0,3,2) == list->no_members) ||
             (Permut023Match(list,item,k,i,2,0,3) == list->no_members) ||
             (Permut023Match(list,item,k,i,2,3,0) == list->no_members) ||
             (Permut023Match(list,item,k,i,3,0,2) == list->no_members) ||
             (Permut023Match(list,item,k,i,3,2,0) == list->no_members) )
         {
            for ( m = 0; m < list->no_params; m++ )
               list->param[k][m] = item->data[i].ff_param[m];
            list->found[k] = 1;
            break;
         }

       } /* end if(permut) */
      }  
   }

}

/*---------------------- Permute Match Function --------------------------*/

int Permut023Match(struct BondList *list, struct FrcFieldItem *item, 
                                           int k, int j, int i0,int i2,int i3)
{
   int match;
   int temp_int[3];
   char temp_a[3][5];

   match  = ( strcmp(list->e[k][0],item->data[j].ff_equiv[i0]) ==0);
   match += ( strcmp(list->e[k][1],item->data[j].ff_equiv[1])  ==0);
   match += ( strcmp(list->e[k][2],item->data[j].ff_equiv[i2]) ==0);
   match += ( strcmp(list->e[k][3],item->data[j].ff_equiv[i3]) ==0);

   /* In Class I, the oops need to be in 'forcefield' order for LAMMPS to
      choose the correct torsion to calculate.  This still fails, because
      if the types are identical, there is no definite way to obtain the 
      the same torsion that Discover will calculate */

   if ( match == list->no_members && forcefield < 2)
   {
      temp_int[0] = list->id[k][0];
      temp_int[1] = list->id[k][2];
      temp_int[2] = list->id[k][3];
  
      list->id[k][i0] = temp_int[0];
      list->id[k][i2] = temp_int[1];
      list->id[k][i3] = temp_int[2];
  
      strcpy(temp_a[0], list->a[k][0]);
      strcpy(temp_a[1], list->a[k][2]);
      strcpy(temp_a[2], list->a[k][3]);
  
      strcpy(list->a[k][i0], temp_a[0]);
      strcpy(list->a[k][i2], temp_a[1]);
      strcpy(list->a[k][i3], temp_a[2]);
  
   }
   return match;
 }
  
