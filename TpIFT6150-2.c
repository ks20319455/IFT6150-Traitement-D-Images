/*------------------------------------------------------*/
/* Prog    : TpIFT6150-2.c                              */
/* Auteurs : Kushal Sangwan, Noé Poluain       */
/* Date    :                                            */
/* version :                                            */ 
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo2.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/                              
/*------------------------------------------------*/
#define NAME_IMG_IN  "photograph"
#define NAME_IMG_GRADIENT "TpIFT6150-2-gradient"
#define NAME_IMG_SUPPRESSION "TpIFT6150-2-suppression"
#define NAME_IMG_CANNY "TpIFT6150-2-canny"

#define FOLLOW(i, j)                                                    \
    follow(i, j, tau_l, suppression, gradient_angle, contour, visites, length, width)

void follow(int i, int j, int tau_l, float** suppression, float** gradient_angle, float** contour, float** visites, int length, int width);
/** 
// Define Sobel kernels
float Gx[3][3] = { {-1, 0, 1},
                   {-2, 0, 2},
                   {-1, 0, 1} };

float Gy[3][3] = { { 1, 2, 1},
                   { 0, 0, 0},
                   {-1,-2,-1} };
**/
// Define Sobel kernels
float Gx[3][3] = { {0, 0, 0},
                   {-1, 1, 0},
                   {0, 0, 0} };

float Gy[3][3] = { { 0, 0, 0},
                   { 0, -1, 0},
                   {0,1,0} };



// Function to apply Sobel operator
void apply_sobel(float** image, float** gradient_x, float** gradient_y, int length, int width) {
    int i, j, x, y;
    for (i = 1; i < length - 1; i++) {
        for (j = 1; j < width - 1; j++) {
            float sum_x = 0.0;
            float sum_y = 0.0;
            for (x = -1; x <= 1; x++) {
                for (y = -1; y <= 1; y++) {
                    sum_x += image[i + x][j + y] * Gx[x + 1][y + 1];
                    sum_y += image[i + x][j + y] * Gy[x + 1][y + 1];
                }
            }
            gradient_x[i][j] = sum_x;
            gradient_y[i][j] = sum_y;
        }
    }
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/
int main(int argc,char** argv)
 {
  int i,j,k,l;
  int length,width;
  float tau_l = 50; // entre 0 360
  float tau_h = 100;
  float p_H = 0.9;
  float sigma = 1;

  float** image_r = LoadImagePgm(NAME_IMG_IN, &length, &width);
  float** image_i = fmatrix_allocate_2d(length, width);

  float** gaussienne_r = fmatrix_allocate_2d(length, width);
  float** gaussienne_i = fmatrix_allocate_2d(length, width);

  float** flou_r = fmatrix_allocate_2d(length, width);
  float** flou_i = fmatrix_allocate_2d(length, width);

  float** gradient_norme = fmatrix_allocate_2d(length, width);
  float** gradient_angle = fmatrix_allocate_2d(length, width);

  float** suppression = fmatrix_allocate_2d(length, width);
  float** visites = fmatrix_allocate_2d(length, width);
  float** contour = fmatrix_allocate_2d(length, width);
  
  // Entrer des valeurs
  printf("Entrez la valeur de tau_L: ");
  scanf("%f",&tau_l);
  printf("Entrez la valeur de tau_H: ");
  scanf("%f",&tau_h);
  printf("Entrez l'ecart type du filter Gaussien: ");
  scanf("%f",&sigma);

  float** gradient_x = fmatrix_allocate_2d(length, width);
  float** gradient_y = fmatrix_allocate_2d(length, width);
  

  
  for(i=0; i < length; i++)
      for(j=0; j < width; j++) {
          image_i[i][j] = 0.0;
          float x = (i - length / 2.0);
          float y = (j - width / 2.0);
          gaussienne_r[i][j] = funcgauss2D((i + length/2) % length - length/2,
                                     (j + width/2) % width - width/2, sigma);

          gaussienne_i[i][j] = 0.0;

          contour[i][j] = 0.0;
          
          visites[i][j] = 0;
          gradient_x[i][j]=0;
          gradient_y[i][j]=0;
          flou_r[i][j]=0;
          flou_i[i][j]=0;
      }

  FFTDD(image_r, image_i, length, width);
  FFTDD(gaussienne_r, gaussienne_i, length, width);

  MultMatrix(flou_r, flou_i, image_r, image_i, gaussienne_r, gaussienne_i, length, width);

  IFFTDD(flou_r, flou_i, length, width);

    // Apply the Sobel operator
  apply_sobel(flou_r, gradient_x, gradient_y, length, width);

  for(i=0; i < length; i++)
      for(j=0; j < width; j++) {
          gradient_norme[i][j] = sqrt(gradient_x[i][j] * gradient_x[i][j] + gradient_y[i][j] * gradient_y[i][j]);
          float angleAbsolute = atan2(-gradient_x[i][j], gradient_y[i][j]);
          if(angleAbsolute<0) angleAbsolute += PI;
          angleAbsolute=(angleAbsolute/PI)*180.0;
          // Quantifier l'angle
            if (angleAbsolute < 22.5 || angleAbsolute >= 157.5) {
                gradient_angle[i][j] = 0; // 0°
            } else if (angleAbsolute >= 22.5 && angleAbsolute < 67.5) {
                gradient_angle[i][j] = 45; // 45°
            } else if (angleAbsolute >= 67.5 && angleAbsolute < 112.5) {
                gradient_angle[i][j] = 90; // 90°
            } else {
                gradient_angle[i][j] = 135; // 135°
            }
      }

  Recal(gradient_norme, length, width);
  SaveImagePgm(NAME_IMG_GRADIENT, gradient_norme, length, width);

  int i_ok, j_ok;
  for(i=0; i < length; i++)
      for(j=0; j < width; j++) {
          suppression[i][j] = gradient_norme[i][j];

          i_ok = !(i - 1 < 0 || i + 1 > length - 1);
          j_ok = !(j - 1 < 0 || j + 1 > width - 1);
          
          int angle = (int) gradient_angle[i][j];
          switch(angle) {
          case 0:
              if(i_ok && gradient_norme[i][j] < MAX(gradient_norme[i - 1][j], gradient_norme[i + 1][j]))
                  suppression[i][j] = 0;
              break;
          case 135:
              if(i_ok && j_ok && gradient_norme[i][j] < MAX(gradient_norme[i - 1][j - 1], gradient_norme[i + 1][j + 1]))
                  suppression[i][j] = 0;
              break;
          case 90:
              if(j_ok && gradient_norme[i][j] < MAX(gradient_norme[i][j - 1], gradient_norme[i][j + 1]))
                  suppression[i][j] = 0;
              break;
          case 45:
              if(i_ok && j_ok && gradient_norme[i][j] < MAX(gradient_norme[i + 1][j - 1], gradient_norme[i - 1][j + 1]))
                  suppression[i][j] = 0;
              break;
          }
      }

  Recal(suppression, length, width);
  SaveImagePgm(NAME_IMG_SUPPRESSION, suppression, length, width);

  for(i=0; i < length; i++)
      for(j=0; j < width; j++) {

          if(suppression[i][j] > tau_h)
              FOLLOW(i, j);
      }

  SaveImagePgm(NAME_IMG_CANNY, contour, length, width);    

  printf("Entrez la valeur de p_H: ");
  scanf("%f",&p_H);

  float histogramme[256];
  float cumul = 0;

  compute_histo(gradient_norme, length, width, histogramme);

  tau_h=0;
  
  while(cumul < p_H) {
      cumul += histogramme[(int) tau_h];
      tau_h++;
  }

  tau_l = 0.5 * tau_h;

  for(i=0; i<length; i++) {
      for(j=0; j<width; j++) {
          contour[i][j] = 0;
          visites[i][j] = 0;
      }
  }
  
  
  for(i=0; i < length; i++)
      for(j=0; j < width; j++) {

          if(suppression[i][j] > tau_h)
              FOLLOW(i, j);
      }

  printf("tau_l = %lf, tau_h = %lf\n", tau_l, tau_h);

  SaveImagePgm("TpIFT6150-2-cannySemiAuto", contour, length, width);
  
  return 0; 	 
}

void follow(int i, int j, int tau_l, float** suppression, float** gradient_angle, float** contour, float** visites, int length, int width) {

    if(visites[i][j])
        return;
    
    // Le point est considéré comme contour
    contour[i][j] = 255.0;
    visites[i][j] = 1;

    int angle = gradient_angle[i][j];
    
    switch((angle + 2) % 4) {
    case 0:
        if(i - 1 > 0 && suppression[i - 1][j] > tau_l) {
            FOLLOW(i - 1, j);
        }
            
        if(i + 1 < length && suppression[i + 1][j] > tau_l) {
            FOLLOW(i + 1, j);
        }
        break;
        
    case 1:
        if(i + 1 < length && j - 1 >= 0 && suppression[i + 1][j - 1] > tau_l) {
            FOLLOW(i + 1, j - 1);
        }
        
        if(i - 1 >= 0 && j + 1 < width && suppression[i - 1][j + 1] > tau_l) {
            FOLLOW(i - 1, j + 1);
        }
        break;
        
    case 2:
        if(j - 1 >= 0 && suppression[i][j - 1] > tau_l) {
            FOLLOW(i, j - 1);
        }

        if(j + 1 < width && suppression[i][j + 1] > tau_l) {
            FOLLOW(i, j + 1);
        }
        break;

    case 3:
        if(i - 1 >= 0 && j - 1 >= 0 && suppression[i - 1][j - 1] > tau_l) {
            FOLLOW(i - 1, j - 1);
        }

        if(i + 1 < length && j + 1 < width && suppression[i + 1][j + 1] > tau_l) {
            FOLLOW(i + 1, j + 1);
        }
        
        break;
    }
}
