/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Ca.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 09/10/2024                                              */
/* Langage               : C                                                       */
/* Cours                 : IFT6150                                                 */
/*---------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo1.h"

#define NAME_IMG_OUT_CARRE "image-TpIFT6150-1-Ca_square"
#define NAME_IMG_OUT_CONVOLUTION "image-TpIFT6150-1-Ca_convolution" 

int main(int argc, char **argv)
{
  int i, j;
  int length = 128, width = 128;  
  int square_size;  
  float **MatriceImgR;  
  float **MatriceImgI; 
  float **MatriceImgM;  
  float **MatriceConvR; 
  float **MatriceConvI;

  //Taille carré
  printf("Entrez la taille du carré blanc (entre 1 et 128): ");
  scanf("%d", &square_size);

  if (square_size < 1 || square_size > 128) {
    printf("La taille du carré que vous avez entrée est invalide. Veuillez entrer une taille entre 1 et 128.\n");
    return 1;
  }

   /*Allocation memoire pour les matrices */
  MatriceImgR = fmatrix_allocate_2d(length, width); 
  MatriceImgI = fmatrix_allocate_2d(length, width); 
  MatriceConvR = fmatrix_allocate_2d(length, width); // Matrice pour le résultat réel de l'auto-convolution
  MatriceConvI = fmatrix_allocate_2d(length, width); // Matrice pour le résultat imaginaire de l'auto-convolution

  //Initialisation à zéro
  for (i = 0; i < length; i++)
    for (j = 0; j < width; j++)
      {
    MatriceImgR[i][j] = 0.0;
    MatriceImgI[i][j] = 0.0;     
    MatriceConvR[i][j] = 0.0;
    MatriceConvI[i][j] = 0.0;
      }
  
  int start_x = 0; // Coordonnée X du début du carré
  int start_y = 0; // Coordonnée Y du début du carré

  //Code pour créer notre carré blanc
  for (i = start_y; i < start_y + square_size; i++) 
    for (j = start_x; j < start_x + square_size; j++)
      {
    MatriceImgR[i][j] = 255.0;
      }

  SaveImagePgm(NAME_IMG_OUT_CARRE, MatriceImgR, length, width);
  system("display image-TpIFT6150-1-Ca_square.pgm&");

  FFTDD(MatriceImgR, MatriceImgI, length, width);

  // Auto-convolution dans le domaine fréquentiel
  for (i = 0; i < length; i++)
    for (j = 0; j < width; j++)
      {
    float real = MatriceImgR[i][j];
    float imaginaire = MatriceImgI[i][j];
    MatriceConvR[i][j] = (real * real) - (imaginaire * imaginaire);
    MatriceConvI[i][j] = 2 * real * imaginaire;
      }

  IFFTDD(MatriceConvR, MatriceConvI, length, width); 
  Recal(MatriceConvR, length, width);           
  SaveImagePgm(NAME_IMG_OUT_CONVOLUTION, MatriceConvR, length, width);

  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(MatriceImgR);
  free_fmatrix_2d(MatriceImgI);
  free_fmatrix_2d(MatriceConvR);
  free_fmatrix_2d(MatriceConvI);

  /*Commande systeme: visualisation de Imgout.pgm*/
  system("display image-TpIFT6150-1-Ca_convolution.pgm&");

  printf("\n C'est fini ... \n\n\n");
  return 0;
}
