/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-4.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 24/11/2024                                              */
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

#define LENGTH 128
#define WIDTH  128
// Changer ce nombre à 45 pour avoir 45 projections
#define NB_PROJECTIONS 90

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


/*----------------------------------------------------*/
/* angle en degre entre BA et BC                      */
/* A(ptar,ptac) - B(ptbr,ptbc) - C(ptcr,ptcc)         */
/*----------------------------------------------------*/
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
int main(int argc, char** argv)
{
    int i, j, k;
    float rotat;

    float **MatriceImgG;
    float **MatriceRadon;
    float **MatriceRadonRFFT;
    float **MatriceRadonIFFT;
    float **MatriceRadonMFFT;
    float **MatRFFT;
    float **MatIFFT;
    float **MatMFFT;
    float **Mat1;
    float **Mat2;
    float *VctR;
    float *VctI;

    /*Allocation memoire des matrices */
    MatriceImgG = fmatrix_allocate_2d(LENGTH, WIDTH);
    MatriceRadon = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatriceRadonRFFT = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatriceRadonIFFT = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatriceRadonMFFT = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatRFFT = fmatrix_allocate_2d(LENGTH, WIDTH);
    MatIFFT = fmatrix_allocate_2d(LENGTH, WIDTH);
    MatMFFT = fmatrix_allocate_2d(LENGTH, WIDTH);
    Mat1 = fmatrix_allocate_2d(LENGTH, WIDTH);
    Mat2 = fmatrix_allocate_2d(LENGTH, WIDTH);

    /*Allocation memoire des vecteurs */
    VctR = fmatrix_allocate_1d(WIDTH);
    VctI = fmatrix_allocate_1d(WIDTH);

    /*Initialisation a zero de toutes les matrices */
    for (i = 0; i < LENGTH; i++)
        for (j = 0; j < WIDTH; j++) {
            MatriceImgG[i][j] = 0.0;
            MatRFFT[i][j] = 0.0;
            MatIFFT[i][j] = 0.0;
            MatMFFT[i][j] = 0.0;
            Mat1[i][j] = 0.0;
            Mat2[i][j] = 0.0;
        }

    for (i = 0; i < LENGTH_RADON; i++)
        for (j = 0; j < WIDTH_RADON; j++) {
            MatriceRadon[i][j] = 0.0;
            MatriceRadonRFFT[i][j] = 0.0;
            MatriceRadonIFFT[i][j] = 0.0;
            MatriceRadonMFFT[i][j] = 0.0;
        }

    /*Initialisation a zero de tous les vecteurs */
    for (i = 0; i < WIDTH; i++) {
        VctR[i] = 0.0;
        VctI[i] = 0.0;
    }

    //On charge l'image dans MatriceImgG
    //----------------------------------
    LoadImagePgm(NAME_IMG_IN, MatriceImgG, LENGTH, WIDTH);

    for(k=0; k<NB_PROJECTIONS; k++) {
        
        matrix_rotation_bilin(MatriceImgG, Mat1, (k * 180.0 / NB_PROJECTIONS) / 360.0 * 2 * PI, LENGTH, WIDTH);
        
        for(i=0; i<LENGTH; i++)
            for(j=0; j<WIDTH; j++) {
                MatriceRadon[k][j] += Mat1[i][j];
            }     
        
        for(i=0; i<WIDTH; i++) {
            VctR[i] = MatriceRadon[k][i];
            VctI[i] = 0.0;
        }

        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);
        
        FFT1D(VctR, VctI, WIDTH);
        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);
 
        for(j=0; j<WIDTH; j++) {
            MatriceRadonRFFT[k][j] = VctR[j];
            MatriceRadonIFFT[k][j] = VctI[j];
        }
    }

    /*----------------------------------------------------------*/
    /*Sauvegarde des matrices sous forme d'image pgms */

    Mod(MatriceRadonMFFT, MatriceRadonRFFT, MatriceRadonIFFT, LENGTH_RADON, WIDTH_RADON);
    
    Recal(MatriceRadonMFFT, LENGTH_RADON, WIDTH_RADON);
    Mult(MatriceRadonMFFT, 60, LENGTH_RADON, WIDTH_RADON);
    
    radon_to_TF2D(MatRFFT, MatIFFT, MatriceRadonRFFT, MatriceRadonIFFT);
    
    Mod(MatMFFT, MatRFFT, MatIFFT, LENGTH, WIDTH);

    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    ReMkeImg(MatIFFT, LENGTH, WIDTH);
    
    IFFTDD(MatRFFT, MatIFFT, LENGTH, WIDTH);

    ReMkeImg(MatRFFT, LENGTH, WIDTH);

    Recal(MatRFFT, LENGTH, WIDTH);
    
    RecalMoy(MatRFFT, MatriceImgG, LENGTH, WIDTH);
    
#if NB_PROJECTIONS == 90
      SaveImagePgm(NAME_IMG_OUT4, MatRFFT, LENGTH, WIDTH);
#else
      SaveImagePgm(NAME_IMG_OUT5, MatRFFT, LENGTH, WIDTH);
#endif

    /*Liberation memoire pour les matrices */
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

    /*Liberation memoire pour les vecteurs */
    free(VctR);
    free(VctI);

    /*Commande systeme: visualisation de Imgout.pgm */
    // system("display ImgOut0.pgm&");

    /*retour sans probleme */
    printf("\n C'est fini ... \n\n\n");
    return 0;
}