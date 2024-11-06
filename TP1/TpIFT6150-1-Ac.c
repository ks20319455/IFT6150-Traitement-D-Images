/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Ac.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 06/10/2024                                              */
/* Langage               : C                                                       */
/* Cours                 : IFT6150                                                 */
/*---------------------------------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo1.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/   
/*------------------------------------------------*/

#define NAME_IMG_OUT "image-TpIFT6150-1-Ac"

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/


void LogOperation(float** mtxM, int length, int width) {
  // Appliquer la fonction log
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < width; j++) {
      mtxM[i][j] = log(1+mtxM[i][j]);
    }
  }
}

int main(int argc,char **argv)
 {
  int length, width;
  float** MatriceImgM;

  /* Si l'utilisateur ne donne pas une image comme input au programme,
     on l'initalise par défaut l'image de output de la partie b */
  char *nameImgInput="image-TpIFT6150-1-Ab";
  if(argc>1)
    nameImgInput=argv[1];
    

  /*Allocation memoire pour la matrice image*/
  MatriceImgM=LoadImagePgm(nameImgInput,&length,&width);

  /* Reverse la mise à l'échelle de la matrice M */
  Mult(MatriceImgM,1.0/100.0,length,width);

  /* Operation Log */
  LogOperation(MatriceImgM,length,width);

  Recal(MatriceImgM,length,width);

  /*Sauvegarde de MatriceImgM sous forme d'image pgm*/
  SaveImagePgm(NAME_IMG_OUT,MatriceImgM,length,width);

  /*Liberation memoire pour la matrice */ 
  free_fmatrix_2d(MatriceImgM);    

  /*Commande systeme: visualisation de Ingout.pgm*/
  system("display image-TpIFT6150-1-Ac.pgm&");

  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
