/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-3-C.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 05/11/2024                                              */
/* Langage               : C                                                       */
/* Cours                 : IFT6150                                                 */
/*---------------------------------------------------------------------------------*/

/*------------------------------------------------*/
/* Fichiers inclus -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo3.h"

/*------------------------------------------------*/
/* Définitions -----------------------------------*/  
/*------------------------------------------------*/
#define IMAGE_IN_NAME  "photograph"
#define IMAGE_OUT_NAME_ORIGINAL "photograph_original_C"
#define IMAGE_OUT_NAME_DEGRADED "photograph_degraded_C"
#define IMAGE_OUT_NAME_RESTORED "photograph_restaured_C"  



void landweber_restore(float** restoredImage, float** restoredImageImag, 
                             float** noisyImageReal,
                             float** filterReal, float** filterImag, 
                             float** originalImageReal, float alpha, 
                             int iterations, int height, int width, int printStatus);
float ComputeSquaredDifference(float** mat1, float** mat2, int height, int width);
void apply_clamp(float** image, int height, int width);

int main(int argc, char** argv) {
    int iterations,height, width;
    float noiseVariance;
    int filterSize,nbLevels = 3;

    float** inputImageReal;    /* Image d'entrée réelle */
    float** inputImageImag;    /* Image d'entrée imaginaire */
    float** filterReal;        /* Filtre réel */
    float** filterImag;        /* Filtre imaginaire */
    float** noisyImageReal;    /* Image bruitée réelle */
    float** noisyImageImag;    /* Image bruitée imaginaire */
    float** restoredImageReal; /* Image restaurée réelle */
    float** restoredImageImag; /* Image restaurée imaginaire */
    float** previousRestoredImageReal; /* Image restaurée précédente */
    float** zeroMatrix;        /* Matrice de zéros */
    float** haarWavelet;       /* Coefficients de transformation Haar */


    /* Demander les paramètres à l'utilisateur */
    printf("Entrez la taille du filtre passe-bas : ");
    scanf("%d", &filterSize);  

    printf("Entrez la variance du bruit : ");
    scanf("%f", &noiseVariance);

    printf("Entrez le nombre d'itérations (Landweber) : ");
    scanf("%d", &iterations);

    /* Charger l'image d'entrée */
    inputImageReal = LoadImagePgm(IMAGE_IN_NAME, &height, &width);
    SaveImagePgm(IMAGE_OUT_NAME_ORIGINAL, inputImageReal, height, width);
    
    /* Allocation des matrices */
    inputImageImag = fmatrix_allocate_2d(height, width);
    filterReal = fmatrix_allocate_2d(height, width);
    filterImag = fmatrix_allocate_2d(height, width);
    noisyImageReal = fmatrix_allocate_2d(height, width);
    noisyImageImag = fmatrix_allocate_2d(height, width);
    restoredImageReal = fmatrix_allocate_2d(height, width);
    restoredImageImag = fmatrix_allocate_2d(height, width);
    previousRestoredImageReal = fmatrix_allocate_2d(height, width);
    haarWavelet = fmatrix_allocate_2d(height, width);
    zeroMatrix = fmatrix_allocate_2d(height, width);

    /* Initialisation des matrices */
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            zeroMatrix[i][j] = 0.0;
            inputImageImag[i][j] = 0.0;
            noisyImageImag[i][j] = 0.0;
            filterImag[i][j] = 0.0;
            restoredImageImag[i][j]= 0.0;
            haarWavelet[i][j]=0.0;
        }
    }

    /* Initialisation du filtre passe-bas */
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float limit = filterSize / 2.0;
            if ((i < limit || i >= height - limit) && 
                (j < limit || j >= width - limit)) {
                filterReal[i][j] = 1.0 / (filterSize * filterSize);
            } else {
                filterReal[i][j] = 0.0;
            }
        }
    }
    
    /* Ajouter du flou à l'image d'entrée */
    FFTDD(inputImageReal, inputImageImag, height, width);
    FFTDD(filterReal, filterImag, height, width);
    
    MultMatrix(noisyImageReal, noisyImageImag, inputImageReal, inputImageImag, filterReal, filterImag, height, width);

    /* Retourner dans le domaine spatial */
    IFFTDD(noisyImageReal, noisyImageImag, height, width);
    IFFTDD(inputImageReal, inputImageImag, height, width);

    /* Ajouter du bruit gaussien */
    add_gaussian_noise(noisyImageReal, height, width, noiseVariance);
    SaveImagePgm(IMAGE_OUT_NAME_DEGRADED, noisyImageReal, height, width);

    /* Initialiser l'image restaurée avec l'image bruitée */
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            restoredImageReal[i][j] = noisyImageReal[i][j];
        }
    }

    /* Application de l'algorithme de Landweber */
    landweber_restore(restoredImageReal, restoredImageImag, noisyImageReal, 
                           filterReal, filterImag, inputImageReal, 1.0, iterations, height, width, 1);

    /* Application de la transformation et du filtrage en ondelettes */
    int iterationCount = 0, shouldContinue = 1;
    int imageSizeAfterWavelet = height / powf(2, nbLevels);
    
    float isnr;
    while (shouldContinue) {
        isnr= 10 * log10(ComputeSquaredDifference(inputImageReal, noisyImageReal, height, width) /
                                ComputeSquaredDifference(inputImageReal, restoredImageReal, height, width));
        
        printf("%d ISNR : %f\n", iterationCount, isnr);

        /* Sauvegarder l'image restaurée précédente */
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                previousRestoredImageReal[i][j] = restoredImageReal[i][j];
            }
        }

        /* Landweber itération supplémentaire */
        landweber_restore(restoredImageReal, restoredImageImag, noisyImageReal,
                               filterReal, filterImag, inputImageReal, 1, 1, height, width, 0);
        
        /* Transformation Haar et filtrage */
        haar2D_complete(restoredImageReal, haarWavelet, nbLevels, height, width);

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (i < imageSizeAfterWavelet && j < imageSizeAfterWavelet)
                    continue;
                if (haarWavelet[i][j] != 0) {
                    haarWavelet[i][j] = MAX(0, SQUARE(haarWavelet[i][j]) - 3 * noiseVariance) / haarWavelet[i][j];
                }
            }
        }

        /* Retourner l'image à l'état restauré */
        ihaar2D_complete(haarWavelet, restoredImageReal, nbLevels, height, width);
        
        /* Vérification de la condition d'arrêt */
        float diffPrevRestored = ComputeSquaredDifference(restoredImageReal, previousRestoredImageReal, height, width);
        float diffPrevRestoredVsZero = ComputeSquaredDifference(previousRestoredImageReal, zeroMatrix, height, width);
        shouldContinue = (diffPrevRestored / diffPrevRestoredVsZero) > 0.00003;

        iterationCount++;
    }

    /* Calculer l'ISNR final */
    isnr = 10 * log10(ComputeSquaredDifference(inputImageReal, noisyImageReal, height, width) /
                     ComputeSquaredDifference(inputImageReal, restoredImageReal, height, width));
    printf("Final - ISNR : %f\n", isnr);

    /* Appliquer la fonction porte et sauvegarder l'image restaurée */
    apply_clamp(restoredImageReal, height, width);
    SaveImagePgm(IMAGE_OUT_NAME_RESTORED, restoredImageReal, height, width);
    
    /* Libération de la mémoire */
    free_fmatrix_2d(inputImageReal);
    free_fmatrix_2d(inputImageImag);
    free_fmatrix_2d(noisyImageReal);
    free_fmatrix_2d(noisyImageImag);
    free_fmatrix_2d(restoredImageReal);
    free_fmatrix_2d(restoredImageImag);
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
                       int num_iterations, int length, int width,int should_print) {

    // Allocation de matrices temporaires pour les calculs intermédiaires
    float** intermediate_real1 = fmatrix_allocate_2d(length, width);  // Matrice pour stocker les résultats intermédiaires
    float** intermediate_imag1 = fmatrix_allocate_2d(length, width);  // Matrice imaginaire intermédiaire
    float** intermediate_real2 = fmatrix_allocate_2d(length, width);  // Matrice pour les résultats intermédiaires après la soustraction
    float** intermediate_imag2 = fmatrix_allocate_2d(length, width);  // Matrice imaginaire intermédiaire après soustraction

    // Calcul de la différence initiale au carré entre l'image d'entrée et l'image dégradée
    float initial_diff = ComputeSquaredDifference(input_real, degraded_real, length, width);
    float isnr;  // Signal-to-noise ratio (ISNR)


    // Début des itérations de l'algorithme de Landweber
    for (int iteration = 0; iteration < num_iterations; iteration++) {

        // Calcul du ISNR (Image Signal-to-Noise Ratio) pour mesurer la qualité de la restauration
        isnr = 10 * log10(initial_diff / ComputeSquaredDifference(input_real, restored_real, length, width));
        if(should_print) printf("Itération: %d , ISNR: %f\n", iteration, isnr);  // Affichage de l'ISNR à chaque itération

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
    if(should_print) printf("Final - ISNR : %f\n", isnr);  // Affichage de l'ISNR final

    // Libération de la mémoire pour les matrices intermédiaires
    free_fmatrix_2d(intermediate_real1);
    free_fmatrix_2d(intermediate_imag1);
    free_fmatrix_2d(intermediate_real2);
    free_fmatrix_2d(intermediate_imag2);
}


