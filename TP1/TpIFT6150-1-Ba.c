/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Ba.c                                        */
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
#define NAME_IMG_IN  "D1r"

#define NAME_IMG_OUT "image-TpIFT6150-1-Ba"

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
  int i,j,k;
  int length,width;
  float** MatriceImgR;
  float** MatriceImgI;

  /*Allocation memoire pour la matrice image*/
  MatriceImgR=LoadImagePgm(NAME_IMG_IN,&length,&width);
  MatriceImgI=fmatrix_allocate_2d(length,width);

  /*Initialisation a zero de toutes les matrices */
  for(i=0;i<length;i++) 
    for(j=0;j<width;j++) 
      {
	MatriceImgI[i][j]=0.0;
      }

  /*FFT*/
  FFTDD(MatriceImgR,MatriceImgI,length,width);
  CenterSpectrum(MatriceImgR,length,width);
  CenterSpectrum(MatriceImgI,length,width);

  int zones[4][2]={{0,1},{2,5},{6,15},{16,63}};
  float** TempMatriceImgR;
  float** TempMatriceImgI;
  TempMatriceImgI=fmatrix_allocate_2d(length,width);
  TempMatriceImgR=fmatrix_allocate_2d(length,width);
   
  for(int k=0;k<4;k++){
    int lower_limit=zones[k][0];
    int upper_limit=zones[k][1];
    int center_x=length/2;
    int center_y=width/2;
    for(i=0;i<length;i++) 
        for(j=0;j<width;j++) 
          {
       if (((lower_limit <= abs(i - center_x) && abs(i - center_x) <= upper_limit) && (abs(j - center_y) <= upper_limit)) ||
    ((lower_limit <= abs(j - center_y) && abs(j - center_y) <= upper_limit) && (abs(i - center_x) <= upper_limit)))
            {
	    TempMatriceImgI[i][j]=MatriceImgI[i][j];
        TempMatriceImgR[i][j]=MatriceImgR[i][j];
            }
        else
            {
        TempMatriceImgI[i][j]=0.0;
        TempMatriceImgR[i][j]=0.0;
            }
        }
    CenterSpectrum(TempMatriceImgR,length,width);
    CenterSpectrum(TempMatriceImgI,length,width);
    IFFTDD(TempMatriceImgR,TempMatriceImgI,length,width);
    
    Recal(TempMatriceImgR,length,width);

    /*Sauvegarde de MatriceImgM sous forme d'image pgm*/
    char* img_name;
    img_name = malloc(strlen(NAME_IMG_OUT) + 25);
    sprintf(img_name, "%s_%d_%d", NAME_IMG_OUT,lower_limit,upper_limit);
    SaveImagePgm(img_name,TempMatriceImgR,length,width);

    /*Commande systeme: visualisation de l'image*/
    char command[256];
    snprintf(command, sizeof(command), "display %s.pgm&", img_name);
    system(command);
  }

  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(TempMatriceImgR);
  free_fmatrix_2d(TempMatriceImgI);                      
  free_fmatrix_2d(MatriceImgR);
  free_fmatrix_2d(MatriceImgI);   

  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
