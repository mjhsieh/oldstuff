/******************************** FrcMenu.c ********************************
*
*   This function currently supplies a menu to the user to choose the
*   forcefield file to be used in the Msi2LMP2 conversion program.
*
*   Created 7-16-97 by Michael Peachey */

#include"Msi2LMP2.h"

void FrcMenu(void)
{
   char menuChoice; /* The value returned by a call to the menu */

   char menu(void);/* Prompts user for a forcefield */
   void menu_option_1(); /* These functions set the forcefield variables */
   void menu_option_2(); /* According to the choices the user made */
   void menu_option_3();
   void menu_option_4();

   menuChoice = menu();  /* Calls to the menu functions */

   switch(menuChoice)
   {
      case '1': menu_option_1(); break;
      case '2': menu_option_2(); break;
      case '3': menu_option_3(); break;
      case '4': menu_option_4(); break;
      case '5': exit(0); break;
      default: printf("This should not happen");
   }

}

/*------------------------- menu() function --------------------------------*/
char menu(void)
{
   char menuChoice[6];

   printf("Select Force Field:\n\n");
   printf("1.  Class I CVFF.frc\n");
   printf("2.  Class II CFF91.frc\n");
   printf("3.  Class II CFF95.frc\n");
   printf("4.  User defined forcefield\n");
   printf("5.  Exit program\n");
   printf("\nEnter a number (1-5): ");

   fgets(menuChoice, 5, stdin);
   printf("\n");

   while ((menuChoice[0] < '1') || (menuChoice[0] > '5'))
   {
      printf("\n%c is not a valid choice, please try again: ", menuChoice[0]);
      fgets(menuChoice, 80, stdin);
   }

   return menuChoice[0];

}
  
/*----------------------- menu_option_1() function -------------------------*/

void menu_option_1()
{
   forcefield = 1;
   strcpy(FrcFileName, "cvff.frc");
   printf("Using cvff.frc Class I forcefield\n\n");
}

/*----------------------- menu_option_2() function -------------------------*/

void menu_option_2()
{
   forcefield = 2;
   strcpy(FrcFileName, "cff91.frc");
   printf("Using cff91.frc Class II forcefield\n\n");
 
}

/*----------------------- menu_option_3() function -------------------------*/

void menu_option_3()
{
   forcefield = 3;
   strcpy(FrcFileName, "cff95.frc");
   printf("Using cff95.frc Class II forcefield\n\n");
}

/*----------------------- menu_option_4() function -------------------------*/

void menu_option_4()
{
   char name[81];

   forcefield = 4;
   puts("\nPlease be careful, this is not recomended!\n");
   puts("\nMust be a Class II forcefield\n");
   puts("Enter the name of your forcefield field with the .frc extension: ");
   gets(name);
   strcpy(FrcFileName, name);
   puts("Using %s Class II forcefield\n\n");
} 

