/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Af.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 08/10/2024                                              */
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
#define NAME_IMG_OUT_ORIGINALE "image-TpIFT6150-1-Af-img"
#define NAME_IMG_OUT_SPECTRE "image-TpIFT6150-1-Af-spc"

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/

void CenterSpectrum(float** mtxM, int length, int width) {
  int half_length = length / 2;
  int half_width = width/2;
  float temp;

  // Échanger les quadrants
  for (int i = 0; i < half_length; i++) {
    for (int j = 0; j < half_width; j++) {
      // Quadrant supérieur gauche (Q1) avec inférieur droit (Q4)
      temp = mtxM[i][j];
      mtxM[i][j] = mtxM[i + half_length][j + half_width];
      mtxM[i + half_length][j + half_width] = temp;

      // Quadrant supérieur droit (Q2) avec inférieur gauche (Q3)
      temp = mtxM[i][j + half_width];
      mtxM[i][j + half_width] = mtxM[i + half_length][j];
      mtxM[i + half_length][j] = temp;
    }
  }
}

void LogOperation(float** mtxM, int length, int width) {
  // Appliquer la fonction log
  CenterSpectrum(mtxM,length,width);
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < width; j++) {
      mtxM[i][j] = log(1+mtxM[i][j]);
    }
  }
}

float Amplitude(float a,float b){
    return sqrt(a*a+b*b);
}


int main(int argc,char **argv)
 {
  int i,j,k;
  int length=128,width=128;
  float** MatriceImgR;
  float** MatriceImgI;
  float** MatriceImgM;
  int size,translation;
  printf("Entrez la taille du carré:");
  scanf("%d",&size);

  printf("Entrez la translation verticale:");
  scanf("%d",&translation);

  /*Allocation memoire pour la matrice image*/
  MatriceImgR=fmatrix_allocate_2d(length,width);
  MatriceImgI=fmatrix_allocate_2d(length,width);
  MatriceImgM=fmatrix_allocate_2d(length,width);

  /*Initialisation a zero de toutes les matrices */
  for(i=0;i<length;i++) 
    for(j=0;j<width;j++) 
      {
	MatriceImgI[i][j]=0.0;
	MatriceImgM[i][j]=0.0;
    MatriceImgR[i][j]=0.0;
      }
  

  /*Initialisation de l'image avec un carré au centre */
  int start_x = (length-size) / 2 - translation;
  int start_y= (width-size) / 2;

  for (int i = start_x; i < start_x + size; i++) {
    for (int j = start_y; j < start_y + size; j++) {
        MatriceImgR[i][j] = 255.0;  // Remplir avec 255 pour le carré blanc
    }
  }

  SaveImagePgm(NAME_IMG_OUT_ORIGINALE,MatriceImgR,length,width);
  system("display image-TpIFT6150-1-Af-img.pgm&");

  /*FFT*/
  FFTDD(MatriceImgR,MatriceImgI,length,width);

  /*Module*/
  Mod(MatriceImgM,MatriceImgR,MatriceImgI,length,width);

  int point_x=1,point_y=1;
  printf("phase[%d][%d]: %f, amplitude:%.2f\n",point_x,point_y,atan2(MatriceImgI[point_x][point_y],MatriceImgR[point_x][point_y]),Amplitude(MatriceImgI[point_x][point_y],MatriceImgR[point_x][point_y]));

  /*Pour visu*/
  Recal(MatriceImgM,length,width);
  LogOperation(MatriceImgM,length,width);
  Mult(MatriceImgM,50.0,length,width);                    

  /*Sauvegarde de MatriceImgM sous forme d'image pgm*/
  SaveImagePgm(NAME_IMG_OUT_SPECTRE,MatriceImgM,length,width);

  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(MatriceImgR);
  free_fmatrix_2d(MatriceImgI); 
  free_fmatrix_2d(MatriceImgM);    

  /*Commande systeme: visualisation de Ingout.pgm*/
  system("display image-TpIFT6150-1-Af-spc.pgm&");

  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
