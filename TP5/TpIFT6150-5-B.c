/*---------------------------------------------------------------------------------*/
/* Prog                  : TpIFT6150-5-B.c                                        */
/* Auteurs               : Kushal Sangwan & Noé Poulain                            */
/* Courrier Électroniques: kushal.sangwan@umontreal.ca & noe.poulain@umontreal.ca  */
/* Date                  : 02/12/2024                                              */
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
#include <time.h>

#include "FonctionDemo5.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/      
/*------------------------------------------------*/
#define NAME_IMG_OUT "B-segmentation"

// fonction pour êchanger deux variables
void swap(float* a, float* b) {
    float temp = *a;  
    *a = *b;        
    *b = temp;      
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/

int main(int argc, char *argv[])
{
    int i,j;
    int seuilMV;
    float mean0,mean1,var0,var1; /* parametres statistiques des deux classes : a estimer  */
    int length,width;

    float** y;  /* image a segmenter */
    float** x;  /* champs d'etiquettes a estimer */
    float*  VectHist;
	
 	if(argc<2){
		printf("Usage :\n\t TpIFT6150-4-B image_a_segmenter\n\n");
		return 0;
	}
	
    /* Ouvrir l'image d'entree */
    char* pathImage = argv[argc - 1];
    y = LoadImagePgm(pathImage, &length, &width);
    x = fmatrix_allocate_2d(length, width);

	/* Lancer l'algorithme des K-Moyennes */
    float* histogramme = malloc(sizeof(float) * GREY_LEVEL);
    compute_histo(y, length, width, histogramme);

    srand(14);
    // la probabilité que mean0 et mean1 soit pareil est presque zero
    mean0 = randomize() * GREY_LEVEL;
    mean1=  randomize() * GREY_LEVEL;

     // On échange le mean0  et mean1 pour assure que mean0 reste le centre de la classe des pixels noirs
    if(mean0 > mean1)
        swap(&mean0,&mean1);

    float sum0, sum1, next_mean0, next_mean1, total0, total1;
    int  seuil_classe1, it = 0, isConverged = 0;
    
    while(!isConverged) {
        total0 = 0;
        total1 = 0;
        sum0 = 0;
        sum1 = 0;
        
        seuil_classe1 = (mean1 - mean0) / 2.0 + mean0;
        
        for(i=0; i<seuil_classe1; i++) {
            sum0 += i*histogramme[i];
            total0 += histogramme[i];
        }

       
        for(i=seuil_classe1; i<=GREY_LEVEL; i++) {
            sum1 += i*histogramme[i];
            total1 += histogramme[i];
        }

        next_mean0 = sum0/total0;
        next_mean1 = sum1/total1;
        it++;
        
        isConverged = next_mean0 == mean0 && next_mean1 == mean1;

        printf("%d : %.1f %.1f   Vs   %.1f  %.1f\n", it, mean0, mean1, next_mean0, next_mean1);

        mean0 = next_mean0;
        mean1 = next_mean1;
    }

    total0 = 0;
    total1 = 0;
    sum0=0;
    sum1=0;
   
    for(i=0; i<seuil_classe1; i++) {
        total0 += histogramme[i];
        sum0 += histogramme[i] * SQUARE(i - mean0);
    }
  
    for(i=seuil_classe1; i<=GREY_LEVEL; i++) {
        total1 += histogramme[i];
        sum1 += histogramme[i] * SQUARE(i - mean1);
    }

    var0 = sum0 / total0;
    var1 = sum1 / total1;

    printf("(%.1f  %.1f)    (%.1f %.1f)\n", mean0, mean1, var0, var1);
    
    for(i=0; i<length; i++)
        for(j=0; j<width; j++) {
            int depasseClasse1 = y[i][j] >= seuil_classe1;
            x[i][j] = (depasseClasse1) * GREY_LEVEL;
        }

    /* Sauvegarde du champ d'etiquettes */
    SaveImagePgm(NAME_IMG_OUT, x, length, width);
    
    /* Liberer la memoire des images */
    free_fmatrix_2d(x);
    free_fmatrix_2d(y);

    /*Commande systeme: visualisation de Ingout.pgm*/
    system("display B-segmentation.pgm&"); 
  
    /*retour sans probleme*/ 
    printf("\n C'est fini ... \n\n\n");
    return 0; 	 
}