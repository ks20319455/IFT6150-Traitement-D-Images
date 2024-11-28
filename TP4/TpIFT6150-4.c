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
// Fonction utilitaire pour gérer les indices circulaires

#define IMAGE_WIDTH WIDTH
#define IMAGE_HEIGHT LENGTH

// Fonction utilitaire pour gérer les indices circulaires dans une matrice
float get_matrix_value(float** matrix, int row, int col, int num_rows, int num_cols) {
    if (col < 0 || col >= num_cols)
        col = (col + num_cols) % num_cols;
    if (row < 0 || row >= num_rows)
        row = (row + num_rows) % num_rows;
    return matrix[row][col];
}

// Calcul du rayon (distance à l'origine) à partir de (x, y)
float calculate_radius(int x, int y) {
    return fmod(round(sqrt(CARRE(x) + CARRE(y)) + IMAGE_WIDTH / 2), IMAGE_WIDTH);
}

// Fonction pour calculer l'angle entre trois points dans un plan 2D
float calculate_angle_deg(int A_x, int A_y, int B_x, int B_y, int C_x, int C_y) {
    float numerator, denominator, angle;
    numerator = (A_y - B_y) * (C_y - B_y) + (A_x - B_x) * (C_x - B_x);
    denominator = sqrt(CARRE(A_y - B_y) + CARRE(A_x - B_x)) * sqrt(CARRE(C_y - B_y) + CARRE(C_x - B_x));

    if (denominator != 0.0)
        angle = acos(numerator / denominator);
    else
        angle = 0.0;

    if (A_x > C_x)
        angle = 2 * PI - angle;
    angle=(angle * (180.0 / PI));
    return angle;  // Conversion en degrés
}

// Interpolation bilinéaire générique entre 4 valeurs
float bilinear_interpolation(float** matrix, int lower_x, int lower_y, float rotated_x, float rotated_y, int num_rows, int num_cols) {
    float top_left = get_matrix_value(matrix, lower_y, lower_x, num_rows, num_cols);
    float top_right = get_matrix_value(matrix, lower_y, lower_x + 1, num_rows, num_cols);
    float bottom_left = get_matrix_value(matrix, lower_y + 1, lower_x, num_rows, num_cols);
    float bottom_right = get_matrix_value(matrix, lower_y + 1, lower_x + 1, num_rows, num_cols);

    // Interpolation sur X
    float top = top_left + (rotated_x - lower_x) * (top_right - top_left);
    float bottom = bottom_left + (rotated_x - lower_x) * (bottom_right - bottom_left);

    // Interpolation sur Y
    return top + (rotated_y - lower_y) * (bottom - top);
}

// Fonction de rotation avec interpolation bilinéaire
void rotate_matrix_bilinear(float** input_matrix, float** output_matrix, float angle, int num_rows, int num_cols) {
    int row, col;
    float rotated_x, rotated_y;
    int lower_x, lower_y;

    for (row = 0; row < num_rows; row++) {
        for (col = 0; col < num_cols; col++) {
            // Rotation avec translation du centre
            rotated_x = (col - num_cols / 2) * cos(-angle) + (row - num_rows / 2) * sin(-angle) + num_cols / 2;
            rotated_y = -(col - num_cols / 2) * sin(-angle) + (row - num_rows / 2) * cos(-angle) + num_rows / 2;

            // Indices arrondis pour l'interpolation
            lower_x = (int)floor(rotated_x);
            lower_y = (int)floor(rotated_y);

            // Appel à la fonction d'interpolation bilinéaire
            output_matrix[row][col] = bilinear_interpolation(input_matrix, lower_x, lower_y, rotated_x, rotated_y, num_rows, num_cols);
        }
    }
}

// Calcul du rayon pour un point (x, y)
float compute_radius(int x, int y) {
    return fmod(round(sqrt(CARRE(x) + CARRE(y)) + IMAGE_WIDTH / 2), IMAGE_WIDTH);
}

// Calcul de l'angle entre (x, y) et l'axe de projection
float compute_angle(int x, int y) {
    return fmod(180 - round(calculate_angle_deg(x, y, 0, 0, IMAGE_WIDTH, 0)), 180) / (180.0 / NB_PROJECTIONS);
}


// Fonction de conversion de l'image Radon vers l'espace Fourier 2D
void radon_to_fourier2D(float** fourier_real, float** fourier_imag, float** radon_real, float** radon_imag) {
    int i, j, angle_idx, radius_idx;
    float angle, radius, interpolated_value1, interpolated_value2;

    for (i = 0; i < IMAGE_HEIGHT / 2; i++) {
        for (j = 0; j < IMAGE_WIDTH; j++) {
            int x = j - IMAGE_WIDTH / 2;
            int y = i - IMAGE_HEIGHT / 2;

            radius = compute_radius(x, y);
            angle = compute_angle(x, y);

            // Indices arrondis pour l'interpolation
            angle_idx = (int)floor(angle);
            radius_idx = (int)floor(radius);

            // Vérification des limites des indices
            angle_idx = fmax(0, fmin(angle_idx, NB_PROJECTIONS - 1));
            radius_idx = fmax(0, fmin(radius_idx, WIDTH_RADON - 1));

            // Interpolation bilinéaire pour la composante réelle de Fourier
            interpolated_value1 = bilinear_interpolation(radon_real, radius_idx, angle_idx, radius, angle, NB_PROJECTIONS, WIDTH_RADON);
            interpolated_value2 = bilinear_interpolation(radon_real, radius_idx, (angle_idx + 1) % NB_PROJECTIONS, radius, angle, NB_PROJECTIONS, WIDTH_RADON);
            fourier_real[i][j] = interpolated_value1 + (angle - angle_idx) * (interpolated_value2 - interpolated_value1);

            // Calculer la composante réelle symétrique pour l'image inversée
            int inverse_j = (IMAGE_WIDTH - j) % IMAGE_WIDTH;
            fourier_real[IMAGE_HEIGHT - 1 - i][inverse_j] = interpolated_value1 + (angle - angle_idx) * (interpolated_value2 - interpolated_value1);

            // Interpolation bilinéaire pour la composante imaginaire de Fourier
            interpolated_value1 = bilinear_interpolation(radon_imag, radius_idx, angle_idx, radius, angle, NB_PROJECTIONS, WIDTH_RADON);
            interpolated_value2 = bilinear_interpolation(radon_imag, radius_idx, (angle_idx + 1) % NB_PROJECTIONS, radius, angle, NB_PROJECTIONS, WIDTH_RADON);
            fourier_imag[i][j] = -(interpolated_value1 + (radius - radius_idx) * (interpolated_value2 - interpolated_value1));

            // Calculer la composante imaginaire symétrique pour l'image inversée
            fourier_imag[IMAGE_HEIGHT - 1 - i][inverse_j] = interpolated_value1 + (radius - radius_idx) * (interpolated_value2 - interpolated_value1);
        }
    }
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL ---------------------------*/
/*------------------------------------------------*/
int main(int argc, char** argv)
{
    int i, j, k;

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
    // Chargement et sauvegarde de l'image d'entrée
    LoadImagePgm(NAME_IMG_IN, MatriceImgG, LENGTH, WIDTH);
    SaveImagePgm(NAME_IMG_OUT0, MatriceImgG, LENGTH, WIDTH);

    // Préparation des matrices de Radon et des vecteurs pour les FFT
    for (int k = 0; k < NB_PROJECTIONS; k++) {
        // Calcul de l'angle de rotation pour la projection Radon
        float angle = (k * 180.0 / NB_PROJECTIONS) / 360.0 * 2 * PI;

        // Appliquer la rotation bilinéaire et accumuler les projections
        rotate_matrix_bilinear(MatriceImgG, Mat1, angle, LENGTH, WIDTH);
        
        // Accumuler les projections dans MatriceRadon
        for (int i = 0; i < LENGTH; i++) {
            for (int j = 0; j < WIDTH; j++) {
                MatriceRadon[k][j] += Mat1[i][j];
            }
        }

        // Sauvegarder l'image de la projection à 45° (si nécessaire)
        if (k == 45) {
            SaveImagePgm(NAME_IMG_OUT1, Mat1, LENGTH, WIDTH);
        }

        // Préparer les vecteurs VctR et VctI pour la FFT
        for (int i = 0; i < WIDTH; i++) {
            VctR[i] = MatriceRadon[k][i];
            VctI[i] = 0.0;
        }

        // Appliquer la FFT 1D sur les projections
        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);
        FFT1D(VctR, VctI, WIDTH);
        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);

        // Sauvegarder les résultats de la FFT dans les matrices RadonRFFT et RadonIFFT
        for (int j = 0; j < WIDTH; j++) {
            MatriceRadonRFFT[k][j] = VctR[j];
            MatriceRadonIFFT[k][j] = VctI[j];
        }
    }

    // Recalage et sauvegarde de la matrice Radon
    Recal(MatriceRadon, LENGTH_RADON, WIDTH_RADON);
    SaveImagePgm(NAME_IMG_OUT2, MatriceRadon, LENGTH_RADON, WIDTH_RADON);

    // Calcul et sauvegarde de la matrice Radon après FFT et IFFT
    Mod(MatriceRadonMFFT, MatriceRadonRFFT, MatriceRadonIFFT, LENGTH_RADON, WIDTH_RADON);
    Recal(MatriceRadonMFFT, LENGTH_RADON, WIDTH_RADON);
    Mult(MatriceRadonMFFT, 100, LENGTH_RADON, WIDTH_RADON);
    SaveImagePgm(NAME_IMG_OUT3, MatriceRadonMFFT, LENGTH_RADON, WIDTH_RADON);

    // Transformation Radon vers Fourier 2D
    radon_to_fourier2D(MatRFFT, MatIFFT, MatriceRadonRFFT, MatriceRadonIFFT);
    Mod(MatMFFT, MatRFFT, MatIFFT, LENGTH, WIDTH);
    SaveImagePgm(NAME_IMG_OUT4, MatMFFT, LENGTH, WIDTH);

    // Inverse FFT et sauvegarde de l'image
    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    ReMkeImg(MatIFFT, LENGTH, WIDTH);
    IFFTDD(MatRFFT, MatIFFT, LENGTH, WIDTH);

    // Recalage et normalisation finale
    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    Recal(MatRFFT, LENGTH, WIDTH);
    RecalMoy(MatRFFT, MatriceImgG, LENGTH, WIDTH);

    // Sauvegarde finale de l'image
    SaveImagePgm(NAME_IMG_OUT5, MatRFFT, LENGTH, WIDTH);

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

    /*retour sans probleme */
    printf("\n C'est fini ... \n\n\n");
    return 0;
}