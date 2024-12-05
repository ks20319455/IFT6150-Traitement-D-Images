/* Prog                  : TpIFT6150-5-D.c                                        */
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

#include "FonctionDemo5.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/      
/*------------------------------------------------*/
#define NAME_IMG_OUT "D-segmentation"

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/
int main(int argc, char *argv[])
{
    int i,j,k,l;
    float alpha;
    float mean0,mean1,var0,var1; /* parametres statistiques des deux classes : a estimer  */
    int length,width;

    float** y;  /* image a segmenter */
    float** x;  /* champs d'etiquettes a estimer */
    float*  VectHist;
	
 	if(argc<3){
		printf("Usage :\n\t TpIFT6150-4-D image_a_segmenter alpha\n\n");
		return 0;
	}

    alpha = atof(argv[argc - 1]);
	
    /* Ouvrir l'image d'entree */
    y = LoadImagePgm(argv[argc - 2], &length, &width);
    x = fmatrix_allocate_2d(length, width);

	/* Lancer l'algorithme des K-Moyennes */

    // Centres initialisés au hasard
    float next_mean0 = -1, next_mean1 = -1;

    float* histogramme = malloc(sizeof(float) * GREY_LEVEL);
    compute_histo(y, length, width, histogramme);

    srand(time(NULL));
    mean0 = randomize() * 255;

    // On ne doit pas avoir deux centres égaux
    while((mean1 = randomize() * 255) == mean0);

    float tmp;
    
    if(mean0 > mean1) {
        // Échange les valeurs pour assurer que mean0 == centre de la classe des pixels noirs
        tmp = mean0;
        mean0 = mean1;
        mean1 = tmp;
    }

    int convergence = 0, debut_classe_1;
    float somme0, somme1, total0, total1, new_mean0, new_mean1;

    while(!convergence) {

        somme0 = 0;
        somme1 = 0;
        total0 = 0;
        total1 = 0;
        
        debut_classe_1 = (mean1 - mean0) / 2.0 + mean0;

        // Classe 0
        for(i=0; i<debut_classe_1; i++) {
            somme0 += i*histogramme[i];
            total0 += histogramme[i];
        }

        // Classe 1
        for(i=debut_classe_1; i<=GREY_LEVEL; i++) {
            somme1 += i*histogramme[i];
            total1 += histogramme[i];
        }

        new_mean0 = somme0/total0;
        new_mean1 = somme1/total1;
        
        convergence = new_mean0 == mean0 && new_mean1 == mean1;

        mean0 = new_mean0;
        mean1 = new_mean1;
    }

    total0 = total1 = 0;
    
    // Classe 0
    for(i=0; i<debut_classe_1; i++) {
        somme0 += histogramme[i] * SQUARE(i - mean0);
        total0 += histogramme[i];
    }
    var0 = somme0 / total0;

    // Classe 1
    for(i=debut_classe_1; i<=GREY_LEVEL; i++) {
        somme1 += histogramme[i] * SQUARE(i - mean1);
        total1 += histogramme[i];
    }
    var1 = somme1 / total1;
    
    for(i=0; i<length; i++)
        for(j=0; j<width; j++) {
            x[i][j] = (y[i][j] > debut_classe_1);
        }

	/* Lancer l'algorithme de Recuit Simulé */
    
    float** vraissemblance0 = fmatrix_allocate_2d(length, width);
    float** vraissemblance1 = fmatrix_allocate_2d(length, width);

    // Pré-calcul des vraissemblances
    for(i=0; i<length; i++)
        for(j=0; j<width; j++) {
            vraissemblance0[i][j] = SQUARE(y[i][j]-mean0)/(2*var0) + log(sqrt(2*PI*var0));
            vraissemblance1[i][j] = SQUARE(y[i][j]-mean1)/(2*var1) + log(sqrt(2*PI*var1));
        }

    
    int changes = 1, iteration = 1;
    int voisins0, voisins1;
    float u0, u1, energie_globale;
    int classe;

    double t0 = 50, tmin = 0.1, tp = t0,
        p0, p1, p_total;
    
	while(tp > tmin) {
        changes = 0;
        energie_globale = 0;
        for(i=0; i<length; i++)
            for(j=0; j<width; j++) {

                voisins0 = 0;
                voisins1 = 0;

                for(k=fmax(0, i-1); k<=fmin(i+1, length - 1); k++)
                    for(l=fmax(0, j-1); l<=fmin(j+1, width - 1); l++) {
                        if(k == i && j == l)
                            continue;
                        
                        voisins0 += x[k][l] == 0;
                        voisins1 += x[k][l] == 1;
                    }

                u0 = vraissemblance0[i][j] + alpha * voisins1;
                u1 = vraissemblance1[i][j] + alpha * voisins0;

                p0 = exp(-1/tp * u0);
                p1 = exp(-1/tp * u1);

                p_total = p0 + p1;

                p0 /= p_total;
                p1 /= p_total;

                if(randomize() < p0) {
                    classe = 0;
                    energie_globale += u0;
                } else {
                    classe = 1;
                    energie_globale += u1;
                }

                x[i][j] = classe;
            }
        
        if(iteration % 10 == 0)
            printf("#%03d : énergie=%f; température=%f\n", iteration, energie_globale, tp);
        
        tp = 0.99 * tp;
        iteration++;
    }
    
    /* Sauvegarde du champ d'etiquettes */
    Recal(x, length, width);
    SaveImagePgm(NAME_IMG_OUT, x, length, width);
    
    /* Liberer la memoire des images */
    free_fmatrix_2d(x);
    free_fmatrix_2d(y);
    free_fmatrix_2d(vraissemblance0);
    free_fmatrix_2d(vraissemblance1);
    
    /*Commande systeme: visualisation de Ingout.pgm*/
    system("display D-segmentation.pgm&"); 
  
    /*retour sans probleme*/ 
    printf("\n C'est fini ... \n\n\n");

    return 0; 	 
}
