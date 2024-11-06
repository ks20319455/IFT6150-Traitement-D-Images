/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-3-B.c                                        */
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
#define INPUT_IMAGE   "photograph"
#define NOISED_IMAGE   "photograph_bruite_B"
#define DENOISED_IMAGE "photograph_debruite_B"


/*------------------------------------------------*/
/* Fonction pour calculer la différence carrée entre deux images */
/*------------------------------------------------*/
float squared_difference(float** image1, float** image2, int height, int width) {
    float total_diff = 0.0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            total_diff += SQUARE(image1[i][j] - image2[i][j]);
        }
    }
    return total_diff;
}

void apply_clamp(float** image, int height, int width) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            image[i][j] = fmax(0, fmin(image[i][j], GREY_LEVEL));
        }
    }
}

int main(int argc, char** argv) {
    int height, width, num_levels;
    float noise_variance, threshold;

    // Demander à l'utilisateur la variance du bruit, le nombre de niveaux et le seuil
    printf("Entrez la variance du bruit: ");
    scanf("%f", &noise_variance);
    
    printf("Entrez le nombre de niveaux à traiter : ");
    scanf("%d", &num_levels);

    printf("Entrez le seuil : ");
    scanf("%f", &threshold);

    // Charger l'image d'origine
    float** original_image = LoadImagePgm(INPUT_IMAGE, &height, &width);
    float** noisy_image = LoadImagePgm(INPUT_IMAGE, &height, &width);
    float** haar_coeffs = fmatrix_allocate_2d(height, width);
    float** tmp_image = fmatrix_allocate_2d(height, width);

    // Ajouter du bruit gaussien à l'image
    add_gaussian_noise(noisy_image, height, width, noise_variance);
    SaveImagePgm(NOISED_IMAGE, noisy_image, height, width);

    // Appliquer la transformée de Haar à l'image bruitée
    haar2D_complete(noisy_image, haar_coeffs, num_levels, height, width);

    // Calcul de la taille de l'image moyenne après réduction (dépend du nombre de niveaux)
    int reduced_image_size = height / powf(2, num_levels);

    // Seuillage des coefficients de Haar
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (i < reduced_image_size && j < reduced_image_size)
                continue;

            if (fabs(haar_coeffs[i][j]) >= threshold) 
                haar_coeffs[i][j] = GREY_LEVEL * haar_coeffs[i][j] / fabs(haar_coeffs[i][j])
                    * (fabs(haar_coeffs[i][j]) - threshold) / (GREY_LEVEL - threshold);
            else 
                // Seuillage proportionnel
                haar_coeffs[i][j] = 0;
        }
    }

    // Appliquer la transformée inverse de Haar
    ihaar2D_complete(haar_coeffs, tmp_image, num_levels, height, width);

    // Calcul du ISNR (Image Signal-to-Noise Ratio)
    float isnr = 10 * log10(squared_difference(original_image, noisy_image, height, width) /
                            squared_difference(original_image, tmp_image, height, width));

    printf("ISNR : %f\n", isnr);

    // Sauvegarder l'image débruitée
    apply_clamp(tmp_image, height, width);
    SaveImagePgm(DENOISED_IMAGE, tmp_image, height, width);

    // Libérer la mémoire allouée
    free_fmatrix_2d(original_image);
    free_fmatrix_2d(noisy_image);
    free_fmatrix_2d(haar_coeffs);
    free_fmatrix_2d(tmp_image);

    // Fin du programme
    printf("\nC'est fini ...\n\n");
    return 0;
}