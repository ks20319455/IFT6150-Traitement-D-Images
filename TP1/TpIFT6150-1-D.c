/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-D.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 09/10/2024                                              */
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
#define ORIGINAL_WIDTH 128
#define ORIGINAL_HEIGHT 128
#define NEW_WIDTH 32
#define NEW_HEIGHT 32
#define NAME_IMG_OUT_DA "image-TpIFT6150-1-Da"
#define NAME_IMG_OUT_DB "image-TpIFT6150-1-Db"

int main()
{
  int x,y;
  float** MatriceImg128;
  float** MatriceImg32;

  MatriceImg128 = fmatrix_allocate_2d(ORIGINAL_HEIGHT,ORIGINAL_WIDTH);
  MatriceImg32 = fmatrix_allocate_2d(NEW_HEIGHT,NEW_WIDTH);

  for (x = 0; x < ORIGINAL_HEIGHT; x++) 
    for (y = 0; y < ORIGINAL_WIDTH; y++)
      {
    int intensity = 128 + 127 * cos(2 * M_PI * (25 * x + 31 * y) / 128);
    if (intensity > 255) intensity = 255;
    if (intensity < 0) intensity = 0;
    MatriceImg128[x][y] = intensity;
      }

  for (int x = 0; x < ORIGINAL_HEIGHT; x += 4)
     for (int y = 0; y < ORIGINAL_WIDTH; y += 4)
       { 
     MatriceImg32[x/4][y/4] = MatriceImg128[x][y];
       }

  /*Sauvegarde de MatriceImg128 sous forme d'image pgm*/
  SaveImagePgm(NAME_IMG_OUT_DA,MatriceImg128,ORIGINAL_HEIGHT,ORIGINAL_WIDTH);
  SaveImagePgm(NAME_IMG_OUT_DB,MatriceImg32,NEW_HEIGHT,NEW_WIDTH);

  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(MatriceImg128);
  free_fmatrix_2d(MatriceImg32);    

  /*Commande systeme: visualisation de Imgout.pgm*/
  system("display image-TpIFT6150-1-Da.pgm&");
  system("display image-TpIFT6150-1-Db.pgm&");
  printf("\n C'est fini ... \n\n\n");
  return 0;
}
