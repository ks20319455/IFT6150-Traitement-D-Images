/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-1-Bb.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 10/10/2024                                              */
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

#define NAME_IMG_OUT "image-TpIFT6150-1-Bb"

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

int main(int argc,char **argv)
 {
  int i,j,k;
  int length,width;
  float** MatriceImgR;
  float** MatriceImgI;
  float** ImagesMatriceImgM[4];
  float** TempMatriceImgR;
  float** TempMatriceImgI;
  float** TempMatriceImgM;

  /*Allocation memoire pour la matrice image*/
  MatriceImgR=LoadImagePgm(NAME_IMG_IN,&length,&width);
  MatriceImgI=fmatrix_allocate_2d(length,width);
  TempMatriceImgM=fmatrix_allocate_2d(length,width);
  TempMatriceImgI=fmatrix_allocate_2d(length,width);
  TempMatriceImgR=fmatrix_allocate_2d(length,width);

  /*Allocation memoire pour les matrices des images tronchées*/
  ImagesMatriceImgM[0]=fmatrix_allocate_2d(length,width);
  ImagesMatriceImgM[1]=fmatrix_allocate_2d(length,width);
  ImagesMatriceImgM[2]=fmatrix_allocate_2d(length,width);
  ImagesMatriceImgM[3]=fmatrix_allocate_2d(length,width);

  /*Initialisation a zero de toutes les matrices */
  for(i=0;i<length;i++) 
    for(j=0;j<width;j++) 
      {
	MatriceImgI[i][j]=0.0;
    TempMatriceImgR[i][j]=0.0;
    TempMatriceImgI[i][j]=0.0;
    TempMatriceImgM[i][j]=0.0;
    ImagesMatriceImgM[0][i][j]=0.0;
    ImagesMatriceImgM[1][i][j]=0.0;
    ImagesMatriceImgM[2][i][j]=0.0;
    ImagesMatriceImgM[3][i][j]=0.0;
      }

  /*FFT*/
  FFTDD(MatriceImgR,MatriceImgI,length,width);
  CenterSpectrum(MatriceImgR,length,width);
  CenterSpectrum(MatriceImgI,length,width);

  int zones[4][2]={{0,1},{2,5},{6,15},{16,63}};

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
    // Copie l'image obtenue par TF inverse dans la matrice correspondante
    for(i=0;i<length;i++) 
        for(j=0;j<width;j++) 
        {
	ImagesMatriceImgM[k][i][j]=TempMatriceImgR[i][j];
        }
  }
  
  for(k=0;k<4;k++){
    //initialise la matrice à zero au début du boucle
    for(i=0;i<length;i++) 
        for(j=0;j<width;j++) 
        {
	TempMatriceImgM[i][j]=0.0;
        }

    //ajoute les images I = I1 + I2 + .... Ik
    for(int nI=0;nI<=k;nI++)
        for(i=0;i<length;i++) 
            for(j=0;j<width;j++) 
            {
	    TempMatriceImgM[i][j] +=ImagesMatriceImgM[nI][i][j];
            }

    Recal(TempMatriceImgM,length,width);

    /*Sauvegarde de TempMatriceImgM sous forme d'image pgm*/
    char* img_name;
    img_name = malloc(strlen(NAME_IMG_OUT) + 25);
    sprintf(img_name, "%s_%d_%d", NAME_IMG_OUT,zones[0][0],zones[k][1]);
    SaveImagePgm(img_name,TempMatriceImgM,length,width);

    /*Commande systeme: visualisation de l'image*/
    char command[256];
    snprintf(command, sizeof(command), "display %s.pgm&", img_name);
    system(command);
  }   

 

  free_fmatrix_2d(TempMatriceImgM);
  free_fmatrix_2d(TempMatriceImgR);
  free_fmatrix_2d(TempMatriceImgI);                      
  free_fmatrix_2d(MatriceImgR);
  free_fmatrix_2d(MatriceImgI); 
  free_fmatrix_2d(ImagesMatriceImgM[0]);                      
  free_fmatrix_2d(ImagesMatriceImgM[1]);
  free_fmatrix_2d(ImagesMatriceImgM[2]); 
  free_fmatrix_2d(ImagesMatriceImgM[3]); 

  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
