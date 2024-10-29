/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-2.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 26/10/2024                                              */
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
#include <float.h>

#include "FonctionDemo2.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/                              
/*------------------------------------------------*/
#define INPUT_IMAGE_NAME  "photograph"
#define GRADIENT_IMAGE_NAME "TpIFT6150-2-gradient"
#define SUPPRESSION_IMAGE_NAME "TpIFT6150-2-suppression"
#define CANNY_IMAGE_NAME "TpIFT6150-2-canny"
#define CANNY_SEMIAUTO_IMAGE_NAME "TpIFT6150-2-cannySemiAuto"
#define IS_CONNECTED(i, j) isWeaklyConnected(i, j, suppressionMap, gradientAngles, contourMap, visitedMap, lowThreshold, imgHeight, imgWidth)


// Déclaration en avance de isWeaklyConnected
void isWeaklyConnected(int i, int j, float** suppressionMap, 
                       float** gradientAngles, float** contourMap, 
                       float** visitedMap, int lowThreshold, 
                       int imgHeight, int imgWidth);

/* Vérifier si le sommet a des indexes valides ou pas et s'il a une valeur moins que tau_l 
et puis appeler la fonction isWeaklyConnected qui explore le reste du composant connexe*/
void connectIfValid(int ni, int nj, float** suppressionMap, 
                    float** gradientAngles, float** contourMap, 
                    float** visitedMap, int lowThreshold, 
                    int imgHeight, int imgWidth) {
    if (ni >= 0 && ni < imgHeight && nj >= 0 && nj < imgWidth && 
        suppressionMap[ni][nj] > lowThreshold) {
        isWeaklyConnected(ni, nj, suppressionMap, gradientAngles, 
                          contourMap, visitedMap, lowThreshold, imgHeight, imgWidth);
    }
}

/*Dans cette fonction, on essaie de trouver un composant connexe i.e les sommets 
qui sont connectés à notre sommet principale( ayant valeur> tau_h) et qui ont une valeur> tau_l*/
void isWeaklyConnected(int i, int j, float** suppressionMap, 
                       float** gradientAngles, float** contourMap, 
                       float** visitedMap, int lowThreshold, 
                       int imgHeight, int imgWidth) {

    // Vérifier si le point a déjà été visité
    if (visitedMap[i][j]) return;

    // Obtenir l'angleTangent du gradient au pixel courant
    int angleTangent = (int) gradientAngles[i][j];

    // Marquer le point comme contour et comme visité
    contourMap[i][j] = 255.0;
    visitedMap[i][j] = 1;

    // Switch en fonction de l'angleTangent du gradient
    switch (angleTangent) {
        case 0:
            connectIfValid(i, j - 1, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Gauche
            connectIfValid(i, j + 1, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Droite
            break;

        case 45:
            connectIfValid(i - 1, j - 1, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Haut-Gauche
            connectIfValid(i + 1, j + 1, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Bas-Droite
            break;

        case 90:
            connectIfValid(i - 1, j, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Haut
            connectIfValid(i + 1, j, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Bas
            break;

        case 135:
            connectIfValid(i + 1, j - 1, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Bas-Gauche
            connectIfValid(i - 1, j + 1, suppressionMap, gradientAngles, 
                           contourMap, visitedMap, lowThreshold, 
                           imgHeight, imgWidth); // Haut-Droite
            break;
    }
}

/** 
// On a aussi essayé d'obtenir des résultats  avec sobel kernels
float Gx[3][3] = { {-1, 0, 1},
                   {-2, 0, 2},
                   {-1, 0, 1} };

float Gy[3][3] = { { 1, 2, 1},
                   { 0, 0, 0},
                   {-1,-2,-1} };
**/

float Gx[3][3] = { {0, 0, 0},
                       {0, -1, 1},
                       {0, 0, 0} };

float Gy[3][3] = { { 0, 0, 0},
                       { 0, -1, 0},
                       {0, 1, 0} };

// Fonction pour appliquer l'opérateur du gradient
void applyGradientFunction(float** image, float** gradientX, float** gradientY, int height, int width) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float sumX = 0.0;
            float sumY = 0.0;
            for (int x = -1; x <= 1; x++) {
                for (int y = -1; y <= 1; y++) {
                    int imageX = i + x;
                    int imageY = j + y;

                    // Fais répétition si on est au bord de l'image
                    imageX = imageX < 0 ? 0 : (imageX >= height ? height - 1 : imageX);
                    imageY = imageY < 0 ? 0 : (imageY >= width ? width - 1 : imageY);

                    int pixelValue = image[imageX][imageY];
                    sumX += pixelValue * Gx[x + 1][y + 1];
                    sumY += pixelValue * Gy[x + 1][y + 1];
                }
            }
            gradientX[i][j] = sumX;
            gradientY[i][j] = sumY;
        }
    }
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/
int main(int argc, char** argv) {
    int imgHeight, imgWidth;
    float gaussianSigma = 1.0; // variance du filtre gaussien
    float lowThreshold = 0.0; // tau_l
    float highThreshold = 100.0; // tau_h
    float pH = 0.70; // entre 0.70 et 0.95

    // Chargement de l'image d'entrée
    float** inputImageReal = LoadImagePgm(INPUT_IMAGE_NAME, &imgHeight, &imgWidth);

    // Allocation des matrices pour les différentes étapes de traitement
    float** inputImageImaginary = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** gaussianFilterReal = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** gaussianFilterImaginary = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** blurredImageReal = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** blurredImageImaginary = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** gradientX = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** gradientY = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** gradientMagnitude = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** gradientAngles = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** suppressionMap = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** visitedMap = fmatrix_allocate_2d(imgHeight, imgWidth);
    float** contourMap = fmatrix_allocate_2d(imgHeight, imgWidth);
  
    // Entrée utilisateur pour les seuils
    printf("Veuillez entrer la valuer de tau_L (seuil bas): ");
    scanf("%f", &lowThreshold);
    printf("Veuillez entrer la valuer de tau_H (seuil haut): ");
    scanf("%f", &highThreshold);
    printf("Veillez entrer l'écart type du filtre gaussien: ");
    scanf("%f", &gaussianSigma);

    // Initialiser les matrices
    for (int i = 0; i < imgHeight; i++) {
        for (int j = 0; j < imgWidth; j++) {
            int centerX = imgHeight / 2;
            int centerY = imgWidth / 2;
            int adjustedX = (i + centerX) % imgHeight - centerX;
            int adjustedY = (j + centerY) % imgWidth - centerY;
            gaussianFilterReal[i][j] = funcgauss2D(adjustedX, adjustedY, gaussianSigma);
            inputImageImaginary[i][j] = 0.0;
            gaussianFilterImaginary[i][j] = 0.0;
            contourMap[i][j] = 0.0;
            visitedMap[i][j] = 0;
            gradientX[i][j] = 0;
            gradientY[i][j] = 0;
            blurredImageReal[i][j] = 0;
            blurredImageImaginary[i][j] = 0;
        }
    }

    // Effectuer FFT et convolution avec filtre gaussien pour obtenir l'image flou
    FFTDD(inputImageReal, inputImageImaginary, imgHeight, imgWidth);
    FFTDD(gaussianFilterReal, gaussianFilterImaginary, imgHeight, imgWidth);
    MultMatrix(blurredImageReal, blurredImageImaginary, inputImageReal, inputImageImaginary, gaussianFilterReal, gaussianFilterImaginary, imgHeight, imgWidth);
    IFFTDD(blurredImageReal, blurredImageImaginary, imgHeight, imgWidth);

    // Appliquer l'opérateur du gradient
    applyGradientFunction(blurredImageReal, gradientX, gradientY, imgHeight, imgWidth);

    // Calculer l'amplitude et l'angleTangent du gradient
    for (int i = 0; i < imgHeight; i++) {
        for (int j = 0; j < imgWidth; j++) {
            gradientMagnitude[i][j] = sqrt(gradientX[i][j] * gradientX[i][j] + gradientY[i][j] * gradientY[i][j]);
            float angleTangent = atan2(-gradientX[i][j], gradientY[i][j]);
            if (angleTangent < 0) angleTangent += PI;
            angleTangent = (angleTangent / PI) * 180.0;

            // Quantifier l'angleTangent
            if (angleTangent < 22.5 || angleTangent >= 157.5) {
                gradientAngles[i][j] = 0; // 0°
            } else if (angleTangent >= 22.5 && angleTangent < 67.5) {
                gradientAngles[i][j] = 45; // 45°
            } else if (angleTangent >= 67.5 && angleTangent < 112.5) {
                gradientAngles[i][j] = 90; // 90°
            } else {
                gradientAngles[i][j] = 135; // 135°
            }
        }
    }

    // Normalisation et sauvegarde de l'image du gradient
    Recal(gradientMagnitude, imgHeight, imgWidth);
    SaveImagePgm(GRADIENT_IMAGE_NAME, gradientMagnitude, imgHeight, imgWidth);

    // Suppression des non-maximums
    for (int i = 0; i < imgHeight; i++) {
        for (int j = 0; j < imgWidth; j++) {
            suppressionMap[i][j] = gradientMagnitude[i][j];

            // Vérifier les limites
            int iGreaterZero = (i - 1 >= 0);
            int iLessHeight = (i + 1 < imgHeight);
            int jGreaterZero = (j - 1 >= 0);
            int jLessWidth = (j + 1 < imgWidth);

            // Obtenir l'angleTangent
            int angleTangent = (int) gradientAngles[i][j];

            // Obtenir les voisins en fonction de l'angleTangent
            float neighbor1 = FLT_MAX, neighbor2 = FLT_MAX;
            switch (angleTangent) {
                case 0:
                    neighbor1 = iGreaterZero ? gradientMagnitude[i - 1][j] : FLT_MAX; // Voisin du haut
                    neighbor2 = iLessHeight ? gradientMagnitude[i + 1][j] : FLT_MAX; // Voisin du bas
                    break;
                case 45:
                    neighbor1 = iLessHeight && jGreaterZero ? gradientMagnitude[i + 1][j - 1] : FLT_MAX; // Voisin du bas-gauche
                    neighbor2 = iGreaterZero && jLessWidth ? gradientMagnitude[i - 1][j + 1] : FLT_MAX; // Voisin du haut-droite
                    break;
                case 90:
                    neighbor1 = jGreaterZero ? gradientMagnitude[i][j - 1] : FLT_MAX; // Voisin de gauche
                    neighbor2 = jLessWidth ? gradientMagnitude[i][j + 1] : FLT_MAX; // Voisin de droite
                    break;
                case 135:
                    neighbor1 = iGreaterZero && jGreaterZero ? gradientMagnitude[i - 1][j - 1] : FLT_MAX; // Voisin du haut-gauche
                    neighbor2 = iLessHeight && jLessWidth ? gradientMagnitude[i + 1][j + 1] : FLT_MAX; // Voisin du bas-droite
                    break;
            }

            // Suppression si c'est un non-maximum
            if (gradientMagnitude[i][j] < MAX(neighbor1, neighbor2)) {
                suppressionMap[i][j] = 0; // Mettre à zéro si le seuil est dépassé
            }
        }
    }

    // Normaliser et sauvegarder l'image de suppression
    Recal(suppressionMap, imgHeight, imgWidth);
    SaveImagePgm(SUPPRESSION_IMAGE_NAME, suppressionMap, imgHeight, imgWidth);

    // pour faire le seuillage par hystérésis, on traverse l'image obtenu, essayant de trouver les composants connexes pour les sommets où valeur > tau_h
    for(int i=0; i < imgHeight; i++)
      for(int j=0; j < imgWidth; j++) {
            if(suppressionMap[i][j] > highThreshold)
              //essayer de trouver le composant connexes i.e les sommets qui sont voisins avec valeur> tau_l
              IS_CONNECTED(i, j);
      }

    SaveImagePgm(CANNY_IMAGE_NAME, contourMap, imgHeight, imgWidth);    

    printf("Veillez entrer la valeur de pH: ");
    scanf("%f",&pH);
    
    //réinitialisatio des matrices et paramètres 
    for(int i=0; i<imgHeight; i++) {
        for(int j=0; j<imgWidth; j++) {
            contourMap[i][j] = 0;
            visitedMap[i][j] = 0;
        }
    }
    highThreshold=0;
    lowThreshold=0;

    float histogramme[256];
    float cumulative = 0.00;

    compute_histo(gradientMagnitude, imgHeight, imgWidth, histogramme);
    
    // computation de fonction de répartition 
    while(cumulative < pH) {
        cumulative = cumulative+ histogramme[(int) highThreshold];
        highThreshold++;
    }

    lowThreshold = 0.5 * highThreshold;

    printf("Seuil Haut Tau_h = %f, Seuil Bas Tau_l = %f\n", highThreshold, lowThreshold);

    //seuillage par hystérésis
    for(int i=0; i < imgHeight; i++)
        for(int j=0; j < imgWidth; j++) 
            if(suppressionMap[i][j] > highThreshold)
                IS_CONNECTED(i, j);

    SaveImagePgm(CANNY_SEMIAUTO_IMAGE_NAME, contourMap, imgHeight, imgWidth);
  

    // Libérer la mémoire
    free_fmatrix_2d(inputImageImaginary);
    free_fmatrix_2d(gaussianFilterReal);
    free_fmatrix_2d(gaussianFilterImaginary);
    free_fmatrix_2d(blurredImageReal);
    free_fmatrix_2d(blurredImageImaginary);
    free_fmatrix_2d(gradientMagnitude);
    free_fmatrix_2d(gradientX);
    free_fmatrix_2d(gradientY);
    free_fmatrix_2d(gradientAngles);
    free_fmatrix_2d(suppressionMap);
    free_fmatrix_2d(visitedMap);
    free_fmatrix_2d(contourMap);

    return 0;
}