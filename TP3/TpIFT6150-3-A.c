/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-3-A.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 05/11/2024                                              */
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

#include "FonctionDemo3.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/  
/*------------------------------------------------*/
#define INPUT_IMAGE_NAME  "photograph"
#define OUTPUT_IMAGE_NAME_ORIGINAL "photograph_original_A"
#define OUTPUT_IMAGE_NAME_DEGRADED_NOISELESS "photograph_degraded_withoutNoise_A"
#define OUTPUT_IMAGE_NAME_RESTORED_NOISELESS "photograph_restored_withoutNoise_A"
#define OUTPUT_IMAGE_NAME_DEGRADED_NOISY "photograph_degraded_withNoise_A"
#define OUTPUT_IMAGE_NAME_RESTORED_NOISY "photograph_restored_withNoise_A"

void landweber_restore(float**, float**, float**, float**, float**, float**, float, int, int, int);
void apply_clamp(float**,int,int);

int main(int argc, char** argv) {
    int num_iterations ,filter_size;  // Nombre d'itérations pour l'algorithme de Landweber et taille du filtre 
    int img_length, img_width, i, j; // Dimensions de l'image
    float noise_variance, alpha = 1.0; 
    
    // Déclaration des matrices pour les différentes étapes du traitement d'image
    float** input_real;   // Partie réelle de l'image d'entrée
    float** input_imag;   // Partie imaginaire de l'image d'entrée
    float** low_pass_filter_real; // Filtre passe-bas pour le flou
    float** low_pass_filter_imag;
    float** degraded_real; // Image dégradée (avec flou et/ou bruit)
    float** degraded_imag;
    float** restored_real; // Image restaurée (après traitement)
    float** restored_imag;
   
    printf("Entrez la largeur du filtre passe-bas : ");
    scanf("%d", &filter_size);  // Lecture de la taille du filtre

    printf("Entrez le nombre d'itérations pour la méthode de Landweber : ");
    scanf("%d", &num_iterations);  // Lecture du nombre d'itérations

    // Chargement de l'image d'entrée
    input_real = LoadImagePgm(INPUT_IMAGE_NAME, &img_length, &img_width);
    SaveImagePgm(OUTPUT_IMAGE_NAME_ORIGINAL, input_real ,img_length, img_width);
    
    // Allocation de mémoire pour les matrices
    input_imag = fmatrix_allocate_2d(img_length, img_width);
    low_pass_filter_real = fmatrix_allocate_2d(img_length, img_width);
    low_pass_filter_imag = fmatrix_allocate_2d(img_length, img_width);
    degraded_real = fmatrix_allocate_2d(img_length, img_width);
    degraded_imag = fmatrix_allocate_2d(img_length, img_width);
    restored_real = fmatrix_allocate_2d(img_length, img_width);
    restored_imag = fmatrix_allocate_2d(img_length, img_width);

    // Initialisation des matrices à zéro
    for (i = 0; i < img_length; i++) {
        for (j = 0; j < img_width; j++) {
            float limit = filter_size / 2.0;
            // Création d'un filtre passe-bas uniforme pour appliquer le flou
            if ((i <  limit || i >= img_length - limit) &&
                (j < limit || j >= img_width - limit)) {
                low_pass_filter_real[i][j] = 1.0 / (filter_size * filter_size);
            } else {
                low_pass_filter_real[i][j] = 0.0;
            }
            low_pass_filter_imag[i][j] = 0.0;
            restored_imag[i][j] = 0.0;
        }
    }

    // Application du flou à l'image d'entrée : g = image + flou
    FFTDD(input_real, input_imag, img_length, img_width);
    FFTDD(low_pass_filter_real, low_pass_filter_imag, img_length, img_width);
    
    MultMatrix(degraded_real, degraded_imag, input_real, input_imag, low_pass_filter_real, low_pass_filter_imag, img_length, img_width);

    // Conversion de l'image dégradée dans le domaine spatial
    IFFTDD(degraded_real, degraded_imag, img_length, img_width);
    IFFTDD(input_real, input_imag, img_length, img_width);
    
    SaveImagePgm(OUTPUT_IMAGE_NAME_DEGRADED_NOISELESS, degraded_real, img_length, img_width);
    
    /*******************************************************/
    /* Restaurer l'image g (floue mais sans bruit)        */
    /*******************************************************/
    landweber_restore(restored_real, restored_imag, degraded_real, low_pass_filter_real, low_pass_filter_imag,
                      input_real, alpha, num_iterations, img_length, img_width);

    apply_clamp(restored_real,img_length,img_width);
    
    SaveImagePgm(OUTPUT_IMAGE_NAME_RESTORED_NOISELESS, restored_real, img_length, img_width);

    printf("Entrez la variance du bruit : ");
    scanf("%f", &noise_variance); // Lecture de la variance du bruit
    
    // Ajout de bruit gaussien à l'image floue : g = g + bruit =  image + flou + bruit
    add_gaussian_noise(degraded_real, img_length, img_width, noise_variance);

    SaveImagePgm(OUTPUT_IMAGE_NAME_DEGRADED_NOISY, degraded_real, img_length, img_width);

    /*******************************************************/
    /* Restaurer l'image g (floue et bruitée)             */
    /*******************************************************/
    landweber_restore(restored_real, restored_imag, degraded_real, low_pass_filter_real, low_pass_filter_imag,
                      input_real, alpha, num_iterations, img_length, img_width);
    
    apply_clamp(restored_real,img_length,img_width);
    
    SaveImagePgm(OUTPUT_IMAGE_NAME_RESTORED_NOISY, restored_real, img_length, img_width);

    // Libération de la mémoire allouée
    free_fmatrix_2d(input_real);
    free_fmatrix_2d(input_imag);
    free_fmatrix_2d(degraded_real);
    free_fmatrix_2d(degraded_imag);
    free_fmatrix_2d(restored_real);
    free_fmatrix_2d(restored_imag);
    free_fmatrix_2d(low_pass_filter_real);
    free_fmatrix_2d(low_pass_filter_imag);
    
    printf("\nLe traitement est terminé avec succès.\n\n");
    return 0;
}

// Fonction pour soustraire deux matrices
void SubtractMatrices(float** result, float** mat1, float** mat2, int length, int width) {
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            result[i][j] = mat2[i][j] - mat1[i][j];
        }
    }
}

// Fonction pour additionner deux matrices
void AddMatrices(float** result, float** mat1, float** mat2, int length, int width) {
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            result[i][j] = mat2[i][j] + mat1[i][j];
        }
    }
}

// Fonction pour calculer la différence au carré entre deux matrices
float ComputeSquaredDifference(float** mat1, float** mat2, int length, int width) {
    float total_diff = 0.0;
    
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            total_diff += SQUARE(mat1[i][j] - mat2[i][j]);
        }
    }

    return total_diff;
}

// Fonction pour réinitialiser manuellement une matrice à zéro
void reset_matrix(float** matrix, int length, int width) {
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            matrix[i][j] = 0.0f;  // Réinitialisation de chaque élément à zéro
        }
    }
}

void apply_clamp(float** f, int height, int width) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            f[i][j] = fmax(0, fmin(f[i][j], 255));
        }
    }
}

/**
 * Fonction de restauration d'image utilisant l'algorithme de Landweber
**/
void landweber_restore(float** restored_real, float** restored_imag, 
                       float** degraded_real, float** filter_real, float** filter_imag,
                       float** input_real, float alpha, 
                       int num_iterations, int length, int width) {

    // Allocation de matrices temporaires pour les calculs intermédiaires
    float** intermediate_real1 = fmatrix_allocate_2d(length, width);  // Matrice pour stocker les résultats intermédiaires
    float** intermediate_imag1 = fmatrix_allocate_2d(length, width);  // Matrice imaginaire intermédiaire
    float** intermediate_real2 = fmatrix_allocate_2d(length, width);  // Matrice pour les résultats intermédiaires après la soustraction
    float** intermediate_imag2 = fmatrix_allocate_2d(length, width);  // Matrice imaginaire intermédiaire après soustraction

    // Calcul de la différence initiale au carré entre l'image d'entrée et l'image dégradée
    float initial_diff = ComputeSquaredDifference(input_real, degraded_real, length, width);
    float isnr;  // Signal-to-noise ratio (ISNR)

    // Initialisation de l'image restaurée avec l'image dégradée
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            restored_real[i][j] = degraded_real[i][j];  // L'image restaurée commence comme l'image dégradée
        }
    }

    // Début des itérations de l'algorithme de Landweber
    for (int iteration = 0; iteration < num_iterations; iteration++) {

        // Calcul du ISNR (Image Signal-to-Noise Ratio) pour mesurer la qualité de la restauration
        isnr = 10 * log10(initial_diff / ComputeSquaredDifference(input_real, restored_real, length, width));
        printf("Itération: %d , ISNR: %f\n", iteration, isnr);  // Affichage de l'ISNR à chaque itération

        // Réinitialisation des matrices intermédiaires à zéro pour la prochaine itération
        reset_matrix(intermediate_real1,length,width);
        reset_matrix(intermediate_imag1, length, width);
        reset_matrix(intermediate_real2,length, width );
        reset_matrix(intermediate_imag2,length, width );

        // Étape 1: Calcul de la transformée de Fourier de f
        FFTDD(restored_real, restored_imag, length, width);

        // Étape 2: Appliquer le filtre passe-bas en multipliant dans le domaine fréquentiel
        MultMatrix(intermediate_real1, intermediate_imag1, filter_real, filter_imag, restored_real, restored_imag, length, width);

        // Étape 3: Transformation inverse de Fourier pour revenir au domaine spatial
        IFFTDD(intermediate_real1, intermediate_imag1, length, width);

        // Étape 4: Transformation inverse de Fourier de l'image restaurée
        IFFTDD(restored_real, restored_imag, length, width);
        
        // Étape 5: Soustraction 
        SubtractMatrices(intermediate_real2, intermediate_real1, degraded_real, length, width);  // calcul de g - f
        
        // Étape 6: Calcul de la transformée de Fourier de la différence
        FFTDD(intermediate_real2, intermediate_imag2, length, width);

        // Étape 7: Appliquer le filtre passe-bas à la différence dans le domaine fréquentiel
        MultMatrix(intermediate_real1, intermediate_imag1, intermediate_real2, intermediate_imag2, filter_real, filter_imag, length, width);

        // Étape 8: Transformation inverse de Fourier de la matrice résultante
        IFFTDD(intermediate_real1, intermediate_imag1, length, width);

        // Étape 9: Multiplication par le facteur d'apprentissage alpha
        Mult(intermediate_real1, alpha, length, width);

        // Étape 10: Ajout du résultat à l'image restaurée
        AddMatrices(restored_real, intermediate_real1, restored_real, length, width);

    }

    // Finalisation du calcul de l'ISNR après toutes les itérations
    isnr = 10 * log10(initial_diff / ComputeSquaredDifference(input_real, restored_real, length, width));
    printf("Final - ISNR: %f\n", isnr);  // Affichage de l'ISNR final

    // Libération de la mémoire pour les matrices intermédiaires
    free_fmatrix_2d(intermediate_real1);
    free_fmatrix_2d(intermediate_imag1);
    free_fmatrix_2d(intermediate_real2);
    free_fmatrix_2d(intermediate_imag2);
}
