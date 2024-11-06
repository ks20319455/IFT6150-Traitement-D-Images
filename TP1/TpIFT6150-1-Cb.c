/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Cb.c                                        */
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
#define NAME_IMG_IN  "D1r"
#define NAME_IMG_OUT "image-TpIFT6150-1-Cb_filtered" 

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/
/*------------------------------------------------*/
void convolution_in_frequential_domaine(float **MatriceImgR, float **MatriceImgI, float **MatriceFilterR, float **MatriceFilterI, float **MatriceConvR, float **MatriceConvI, int length, int width)
{
  // Multiplication dans le domaine fréquentiel (convolution dans le domaine spatial)
  for (int i = 0; i < length; i++) 
    for (int j = 0; j < width; j++)
      {
   float real_img = MatriceImgR[i][j],image_img = MatriceImgI[i][j];
   float real_filter = MatriceFilterR[i][j],image_filter = MatriceFilterI[i][j];
   MatriceConvR[i][j] = real_img * real_filter - image_img * image_filter;
   MatriceConvI[i][j] = real_img * image_filter + image_img * real_filter;
      }
}

int main(int argc, char **argv)
{
  int i, j;
  int length = 128, width = 128; 
  int filter_size;  // Variable pour stocker la taille du carré
  float **MatriceImgR;  
  float **MatriceImgI;    
  float **MatriceImgFiltreR;
  float **MatriceImgFiltreI;
  float **MatriceConvR; 
  float **MatriceConvI;

  printf("Entrez la taille du filtre (entre 1 et 128): ");
  scanf("%d", &filter_size);

  if (filter_size < 1 || filter_size > 128) {
    printf("La taille du filtre que vous avez entrée est invalide. Veuillez entrer une taille entre 1 et 128.\n");
    return 1;
  }

  MatriceImgR = fmatrix_allocate_2d(length, width);  
  MatriceImgI = fmatrix_allocate_2d(length, width);
  MatriceConvR = fmatrix_allocate_2d(length, width);
  MatriceConvI = fmatrix_allocate_2d(length, width); 
  MatriceImgFiltreR = fmatrix_allocate_2d(length, width); 
  MatriceImgFiltreI = fmatrix_allocate_2d(length, width); 

 
  for (i = 0; i < length; i++)
    for (j = 0; j < width; j++)
      {
    MatriceImgR[i][j] = 0.0;
    MatriceImgI[i][j] = 0.0;
    MatriceConvR[i][j] = 0.0;
    MatriceConvI[i][j] = 0.0;  
    MatriceImgFiltreI[i][j] = 0.0;
    MatriceImgFiltreR[i][j] = 0.0;
      }

  MatriceImgR=LoadImagePgm(NAME_IMG_IN,&length,&width);

  //Code pour créer notre carré blanc
  int start_x = 0;
  int start_y = 0;
  for (i = start_y; i < start_y + filter_size; i++) 
    for (j = start_x; j < start_x + filter_size; j++)
      {
    MatriceImgFiltreR[i][j] = 255.0;
      }
 
  FFTDD(MatriceImgR, MatriceImgI, length, width);  
  FFTDD(MatriceImgFiltreR, MatriceImgFiltreI, length, width);

  convolution_in_frequential_domaine(MatriceImgR, MatriceImgI,MatriceImgFiltreR, MatriceImgFiltreI, MatriceConvR,MatriceConvI,length, width);

  IFFTDD(MatriceConvR, MatriceConvI, length, width); 

  Recal(MatriceConvR, length, width);
  SaveImagePgm(NAME_IMG_OUT, MatriceConvR, length, width); 

  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(MatriceImgR);
  free_fmatrix_2d(MatriceImgI);
  free_fmatrix_2d(MatriceConvR);
  free_fmatrix_2d(MatriceConvI);

  /*Commande systeme: visualisation de Imgout.pgm*/
  system("display image-TpIFT6150-1-Cb_filtered.pgm&");

  printf("\n C'est fini ... \n\n\n");
  return 0;
}
