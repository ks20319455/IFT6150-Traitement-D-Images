/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Ab.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 25/09/2024                                              */
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
#define NAME_IMG_IN  "image-TpIFT6150-1-Aa"

#define NAME_IMG_OUT "image-TpIFT6150-1-Ab"

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/

void CenterSpectrum(float** MatriceImgM, int length, int width) {
  int half_length = length / 2;
  int half_width = width/2;
  float temp;

  // Échanger les quadrants
  for (int i = 0; i < half_length; i++) {
    for (int j = 0; j < half_width; j++) {
      // On échange quadrant supérieur gauche (Q1) avec inférieur droit (Q4)
      temp = MatriceImgM[i][j];
      MatriceImgM[i][j] = MatriceImgM[i + half_length][j + half_width];
      MatriceImgM[i + half_length][j + half_width] = temp;

      // On échange quadrant supérieur droit (Q2) avec inférieur gauche (Q3)
      temp = MatriceImgM[i][j + half_width];
      MatriceImgM[i][j + half_width] = MatriceImgM[i + half_length][j];
      MatriceImgM[i + half_length][j] = temp;
    }
  }
}

int main(int argc,char **argv)
 {
  int length, width;
  float** MatriceImgM;

  /*Allocation memoire pour la matrice image*/
  MatriceImgM=LoadImagePgm(NAME_IMG_IN,&length,&width);

  /*Centralization du spectre */
  CenterSpectrum(MatriceImgM, length, width);

  /*Sauvegarde de MatriceImgM sous forme d'image pgm*/
  SaveImagePgm(NAME_IMG_OUT,MatriceImgM,length,width);

  /*Liberation memoire pour la matrice */ 
  free_fmatrix_2d(MatriceImgM);    

  /*Commande systeme: visualisation de Imgout.pgm*/
  system("display image-TpIFT6150-1-Ab.pgm&");

  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
