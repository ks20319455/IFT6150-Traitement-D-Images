/*------------------------------------------------------*/
/* Prog    : ProgDemo11.c                                */
/* Auteur  :                                            */
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

#include "FonctionDemo.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/
/*------------------------------------------------*/
#define NAME_IMG_IN  "lenna"

#define NAME_IMG_OUT0 "ImgOut0"
#define NAME_IMG_OUT1 "ImgOut1" 
#define NAME_IMG_OUT2 "ImgOut2" 
#define NAME_IMG_OUT3 "ImgOut3"  
#define NAME_IMG_OUT4 "ImgOut4"
#define NAME_IMG_OUT5 "ImgOut5" 

#define NB_PROJECTIONS 180
#define LENGTH 128
#define WIDTH  128

#define LENGTH_RADON NB_PROJECTIONS
#define WIDTH_RADON  WIDTH

/*------------------------------------------------*/
/* FONCTIONS -------------------------------------*/                     
/*------------------------------------------------*/

float matrix_val(float** m, int j, int i, int length, int width) {

    if(i < 0 || i >= length)
        i = (i + length) % length;
    
    if(j < 0 || j >= width)
        j = (j + width) % width;

    return m[i][j];
}

void matrix_rotation_ppv(float** src, float** dest, float angle, int l, int w) {
    int x1,y1,x2,y2;
    float xp,yp;

    for(y2=0; y2<l; y2++) {
        for(x2=0; x2<w; x2++) {

            xp =  (x2 - w/2) * cos(-angle) + (y2 - l/2) * sin(-angle) + w/2;
            yp = -(x2 - w/2) * sin(-angle) + (y2 - l/2) * cos(-angle) + l/2;

            x1 = floor(xp);
            y1 = floor(yp);

            dest[y2][x2] = matrix_val(src, x1, y1, l, w);
        }
    }
}

void matrix_rotation_bilin(float** src, float** dest, float angle, int l, int w) {
    int x1,y1,x2,y2;
    float xp,yp,f1,f2,f3;

    for(y2=0; y2<l; y2++) {
        for(x2=0; x2<w; x2++) {

            xp =  (x2 - w/2) * cos(-angle) + (y2 - l/2) * sin(-angle) + w/2;
            yp = -(x2 - w/2) * sin(-angle) + (y2 - l/2) * cos(-angle) + l/2;

            x1 = floor(xp);
            y1 = floor(yp);
                        
            f1 = matrix_val(src, x1, y1, l, w) + (xp - x1)*(matrix_val(src, x1 + 1, y1, l, w) - matrix_val(src, x1, y1, l, w));
            f2 = matrix_val(src, x1, y1 + 1, l, w)
                + (xp - x1) * (matrix_val(src, x1 + 1, y1 + 1, l, w) - matrix_val(src, x1, y1 + 1, l, w));

            dest[y2][x2] = f1 + (yp - y1) * (f2 - f1);
        }
    }
}

void RotateImageBilinear(float** input, float** output, int length, int width, float angle) {
    float centerX = width / 2.0;
    float centerY = length / 2.0;

    float rad = angle * (PI / 180.0);
    float cosTheta = cos(rad);
    float sinTheta = sin(rad);

    // Set the output image to zero
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            output[i][j] = 0.0;
        }
    }

    // Perform rotation with bilinear interpolation
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            float x = j - centerX;
            float y = i - centerY;

            float xRot = x * cosTheta - y * sinTheta + centerX;
            float yRot = x * sinTheta + y * cosTheta + centerY;

            if (xRot >= 0 && xRot < width - 1 && yRot >= 0 && yRot < length - 1) {
                int x0 = (int)xRot;
                int y0 = (int)yRot;
                int x1 = x0 + 1;
                int y1 = y0 + 1;

                float dx = xRot - x0;
                float dy = yRot - y0;

                float topLeft = input[y0][x0];
                float topRight = input[y0][x1];
                float bottomLeft = input[y1][x0];
                float bottomRight = input[y1][x1];

                output[i][j] = (1 - dx) * (1 - dy) * topLeft +
                               dx * (1 - dy) * topRight +
                               (1 - dx) * dy * bottomLeft +
                               dx * dy * bottomRight;
            }
        }
    }
}

// Compute the projection on the x-axis
void ComputeProjection(float** image, float** MatriceRadon,int k) {
    int i,j;
    for(i=0; i<LENGTH; i++)
        for(j=0; j<WIDTH; j++) {
            MatriceRadon[k][j] += image[i][j];
        }   
}

// Compute the Radon transform
void ComputeRadonTransform(float** MatriceImgG, float** MatriceRadon, int length, int width) {
    float** RotatedImage = fmatrix_allocate_2d(length, width);

    for (int k= 0; k< NB_PROJECTIONS; k++) {
        int angle=(k * 180.0 / NB_PROJECTIONS) / 360.0 * 2 * PI;
        // Rotate the image
        matrix_rotation_bilin(MatriceImgG, RotatedImage, angle , LENGTH, WIDTH);

        // Compute the projection along the x-axis
        ComputeProjection(RotatedImage, MatriceRadon, k);

    }

    free_fmatrix_2d(RotatedImage);
}

void CalculProjectionFiltre(float** MatriceRadonRFFT,float** MatriceRadonIFFT,float** MatriceRadon,float* VctR,float* VctI){
    for(int i=0;i<NB_PROJECTIONS;i++){
        for(int j=0; j<WIDTH; j++) {
            VctR[j] = MatriceRadon[i][j];
            VctI[j] = 0.0;
        }

        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);

        FFT1D(VctR, VctI, WIDTH);
        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);

        for(int j=0; j<WIDTH; j++) {
            MatriceRadonRFFT[i][j] = VctR[j];
            MatriceRadonIFFT[i][j] = VctI[j];
        }
    }
}

float AngleDeg(int ptar,int ptac,int ptbr,int ptbc,int ptcr,int ptcc)
{
    float num,den;
    float angle;

    //initialisation
    num=den=angle=0.0;

    //calcul potentiel
    num+=(float)(ptac-ptbc)*(ptcc-ptbc);
    num+=(float)(ptar-ptbr)*(ptcr-ptbr);
    den+=(float)sqrt(CARRE(ptac-ptbc)+CARRE(ptar-ptbr));
    den*=(float)sqrt(CARRE(ptcc-ptbc)+CARRE(ptcr-ptbr));

    if (den!=0.0) angle=acos(num/den);
    else angle=0.0;

    if (ptar>ptcr) angle=2*PI-angle;

    return (angle*(180.0/PI));
}

void radon_to_TF2D(float** tf_r, float** tf_i, float** radon_r, float** radon_i) {
    int i, j, x, y, angle1, rayon1;

    float angle, rayon, f1, f2;

    for(i=0; i<LENGTH/2; i++)
        for(j=0; j<WIDTH; j++) {

            x = j - WIDTH/2;
            y = i - WIDTH/2;

            rayon = fmod(round(sqrt(CARRE(x) + CARRE(y)) + WIDTH/2), WIDTH);

            angle = fmod(180 - round(AngleDeg(x, y, 0, 0, WIDTH, 0)), 180) / (180.0 / NB_PROJECTIONS);
            
            // Interpolation bilinéaire
            angle1 = floor(angle); // x1
            rayon1 = floor(rayon); // y1

            // Calcul de tf_r
            f1 = radon_r[angle1][rayon1] + (angle - angle1)*(radon_r[(angle1 + 1) % NB_PROJECTIONS][rayon1] - radon_r[angle1][rayon1]);
            f2 = radon_r[angle1][(rayon1 + 1) % WIDTH_RADON]
                + (angle - angle1) * (radon_r[(angle1 + 1) % NB_PROJECTIONS][(rayon1 + 1) % WIDTH_RADON] - radon_r[angle1][(rayon1 + 1) % WIDTH_RADON]);
            
            tf_r[i][j] = f1 + (rayon - rayon1) * (f2 - f1);
            tf_r[LENGTH - 1 - i][(WIDTH - j) % WIDTH] = f1 + (rayon - rayon1) * (f2 - f1);
            
            // Calcul de tf_i
            f1 = radon_i[angle1][rayon1] + (angle - angle1)*(radon_i[(angle1 + 1) % NB_PROJECTIONS][rayon1] - radon_i[angle1][rayon1]);
            f2 = radon_i[angle1][(rayon1 + 1) % WIDTH_RADON]
                + (angle - angle1) * (radon_i[(angle1 + 1) % NB_PROJECTIONS][(rayon1 + 1) % WIDTH_RADON] - radon_i[angle1][(rayon1 + 1) % WIDTH_RADON]);
            
            tf_i[i][j] = -(f1 + (rayon - rayon1) * (f2 - f1));
            tf_i[LENGTH - 1 - i][(WIDTH - j) % WIDTH] = f1 + (rayon - rayon1) * (f2 - f1);
        }
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL ---------------------------*/                     
/*------------------------------------------------*/
int main()
 {
  int i,j;
  float rotat;
  
  float** MatriceImgG;
  float** MatriceRadon;
  float** MatriceRadonRFFT;
  float** MatriceRadonIFFT;
  float** MatriceRadonMFFT;
  float** MatRFFT;
  float** MatIFFT;
  float** MatMFFT;
  float** Mat1;
  float** Mat2;
  float*  VctR;
  float*  VctI;
  
  /*Allocation memoire des matrices*/
  MatriceImgG=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatriceRadon=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonRFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonIFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonMFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatRFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatIFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatMFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  Mat1=fmatrix_allocate_2d(LENGTH,WIDTH);
  Mat2=fmatrix_allocate_2d(LENGTH,WIDTH);

  /*Allocation memoire des vecteurs*/
  VctR=fmatrix_allocate_1d(WIDTH);
  VctI=fmatrix_allocate_1d(WIDTH);

  /*Initialisation a zero de toutes les matrices*/
  for(i=0;i<LENGTH;i++) for(j=0;j<WIDTH;j++) 
   { MatriceImgG[i][j]=0.0;
     MatRFFT[i][j]=0.0;
     MatIFFT[i][j]=0.0;
     MatMFFT[i][j]=0.0;
     Mat1[i][j]=0.0;
     Mat2[i][j]=0.0; }

  for(i=0;i<LENGTH_RADON;i++) for(j=0;j<WIDTH_RADON;j++)
     { MatriceRadon[i][j]=0.0;
       MatriceRadonRFFT[i][j]=0.0;
       MatriceRadonIFFT[i][j]=0.0;
       MatriceRadonMFFT[i][j]=0.0; }

  /*Initialisation a zero de tous les vecteurs*/
   for(i=0;i<WIDTH;i++) 
     { VctR[i]=0.0;
       VctI[i]=0.0; }
    
   //On charge l'image dans MatriceImgG
   //----------------------------------
  LoadImagePgm(NAME_IMG_IN,MatriceImgG,LENGTH,WIDTH);

  ComputeRadonTransform(MatriceImgG, MatriceRadon, LENGTH, WIDTH);

  // Normalize the Radon matrix
  
  SaveImagePgm(NAME_IMG_OUT0,MatriceImgG,LENGTH,WIDTH);
  // Save the Radon transform
  Recal(MatriceRadon, LENGTH_RADON, WIDTH_RADON);
  SaveImagePgm(NAME_IMG_OUT1, MatriceRadon, LENGTH_RADON, WIDTH_RADON);

  CalculProjectionFiltre(MatriceRadonRFFT,MatriceRadonIFFT,MatriceRadon,VctR,VctI);
  Mod(MatriceRadonMFFT, MatriceRadonRFFT, MatriceRadonIFFT, LENGTH_RADON, WIDTH_RADON);
  Recal(MatriceRadonMFFT, LENGTH_RADON, WIDTH_RADON);
  Mult(MatriceRadonMFFT, 60, LENGTH_RADON, WIDTH_RADON);
  SaveImagePgm(NAME_IMG_OUT2, MatriceRadonMFFT, LENGTH_RADON, WIDTH_RADON);
  
  radon_to_TF2D(MatRFFT, MatIFFT, MatriceRadonRFFT, MatriceRadonIFFT);

  Mod(MatMFFT, MatRFFT, MatIFFT, LENGTH, WIDTH);
  Mult(MatMFFT, 60, LENGTH, WIDTH);
  SaveImagePgm(NAME_IMG_OUT3, MatMFFT, LENGTH, WIDTH);

  ReMkeImg(MatRFFT, LENGTH, WIDTH);
  ReMkeImg(MatIFFT, LENGTH, WIDTH);

  IFFTDD(MatRFFT, MatIFFT, LENGTH, WIDTH);

  ReMkeImg(MatRFFT, LENGTH, WIDTH);

  Recal(MatRFFT, LENGTH, WIDTH);

  RecalMoy(MatRFFT, MatriceImgG, LENGTH, WIDTH);
  SaveImagePgm(NAME_IMG_OUT4, MatRFFT, LENGTH, WIDTH);


  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(MatriceImgG); 
  free_fmatrix_2d(MatriceRadon);
  free_fmatrix_2d(MatriceRadonRFFT);
  free_fmatrix_2d(MatriceRadonIFFT);
  free_fmatrix_2d(MatriceRadonMFFT);
  free_fmatrix_2d(MatRFFT);
  free_fmatrix_2d(MatIFFT);
  free_fmatrix_2d(MatMFFT);
  free_fmatrix_2d(Mat1); 
  free_fmatrix_2d(Mat2); 

  /*Liberation memoire pour les vecteurs*/
  free(VctR);
  free(VctI);   

  /*Commande systeme: visualisation de Ingout.pgm*/
  //system("gimp ImgOut0.pgm&");
  system("display ImgOut0.pgm&");
  system("display ImgOut1.pgm&");
  
  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
