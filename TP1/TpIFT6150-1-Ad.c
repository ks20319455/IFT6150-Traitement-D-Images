/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Ad.c                                        */
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
char* NAME_IMG_IN[] = {"D11r", "D46r"};

char* NAME_IMG_OUT[] = {"image-TpIFT6150-1-Ad-D11r-spc","image-TpIFT6150-1-Ad-D46r-spc"};


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
  // Amener l'harmonique (0,0) au centre
  CenterSpectrum(mtxM,length,width);
  // Appliquer la fonction log
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < width; j++) {
      mtxM[i][j] = log(1+mtxM[i][j]);
    }
  }
}


/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/
int main(int argc,char **argv)
 {
  // fais l'opération deux fois pour nos deux images D11r et D46r car D1r et déja fait.
  for(int nIm=0;nIm<2;nIm++){
      int i,j,k;
      int length,width;
      float** MatriceImgR;
      float** MatriceImgI;
      float** MatriceImgM;

      /*Allocation memoire pour la matrice image*/
      MatriceImgR=LoadImagePgm(NAME_IMG_IN[nIm],&length,&width);
      MatriceImgI=fmatrix_allocate_2d(length,width);
      MatriceImgM=fmatrix_allocate_2d(length,width);

      /*Initialisation a zero de toutes les matrices */
      for(i=0;i<length;i++) 
        for(j=0;j<width;j++) 
          {
	    MatriceImgI[i][j]=0.0;
	    MatriceImgM[i][j]=0.0;
          }
  
      /*FFT*/
      FFTDD(MatriceImgR,MatriceImgI,length,width);

      /*Module*/
      Mod(MatriceImgM,MatriceImgR,MatriceImgI,length,width);

      /*Pour visu*/
      Recal(MatriceImgM,length,width);
      LogOperation(MatriceImgM,length,width);
      Mult(MatriceImgM,200.0,length,width);                     
  
      /*Sauvegarde de MatriceImgM sous forme d'image pgm*/
      SaveImagePgm(NAME_IMG_OUT[nIm],MatriceImgM,length,width);

      /*Liberation memoire pour les matrices*/
      free_fmatrix_2d(MatriceImgR);
      free_fmatrix_2d(MatriceImgI); 
      free_fmatrix_2d(MatriceImgM);    

      char command[256];
      snprintf(command, sizeof(command), "display %s.pgm&", NAME_IMG_OUT[nIm]);
      /*Commande systeme: visualisation de l'image*/
      system(command);
  }

  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
