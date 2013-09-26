/*
IPOL SIFT
Copyright (C) 2013, Ives Rey Otero, CMLA ENS Cachan 
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20130318 (March 18, 2013)

== Patent Warning and Licence =================================================

The SIFT method is patented 

    [3] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89
  
 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.

*/
/**
 @file sift_keypoint.c
 * @brief data structures to store information relative to keypoint
 *
 * @li struct kypt      keypoint data structure.
 * @li struct kypt_lst  list of keypoints with variable length.
 * @li print,save, read  for lists of keypoints.
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */










#include "sift_keypoint.h"



/** @brief Memory allocation for keypoint structure
 *
 * Allocate memory and set to default values.
 *
 * 
 * @param n_ori   : number of bins per histogram
 * @param n_hist  : number of histogram per dimension 
 *
 * the feature vector is n_ori*n_hist*n_hist long.
 */
struct kypt* malloc_kypt(int n_ori, int n_hist){
    // malloc with ONE element
    struct kypt* key =  (struct kypt*)malloc(sizeof(struct kypt));  
    
    key->n_ori  = n_ori;
    key->n_hist = n_hist; 
    key->x=0.;
    key->y=0.;
    key->sigma=0.;
    key->theta=0.;
    key->val=0;
    key->o=0;
    key->s=0;
    key->i=0;
    key->j=0;
    key->descr = (double*)malloc(n_ori*n_hist*n_hist*sizeof(double));
    for(int n=0;n<n_ori*n_hist*n_hist;n++){
        key->descr[n]=0.;
    }  
   
    
    
    return key;
}

 
/** @brief Memory allocation for keypoint structure with EXTRA information
 *
 * The extra information is :
 *    - extrema interpolation.
 *    - orientation assignment (histogram).
 *    - on edge unstable keypoints filtering.
 *
 * extra parameters :
 * @param n_bins  : number of histogram bins for orientation attribution
 * @param nMaxIntrp  : max number of successive interpolation (5 in practice)
 */
struct kypt* malloc_kypt_extra(int n_ori, int n_hist, int n_bins, int nMaxIntrp){
    
    struct kypt* key = malloc_kypt(n_ori, n_hist);
    
    key->n_bins = n_bins;
    key->nMaxIntrp = nMaxIntrp; 
    key->pos_intrp    =    (int*)malloc(3*nMaxIntrp*sizeof(int));
    key->offset_intrp = (double*)malloc(3*nMaxIntrp*sizeof(double));  
    for(int i=0;i<3*nMaxIntrp;i++){
        key->pos_intrp[i]=0;
        key->offset_intrp[i]=0.;
    }
    key->hsEdgeResp=0;
    key->hist_prOri=(double*)malloc(n_bins*sizeof(double));
    key->n_prOri=0;
    key->ori_prOri=(double*)malloc(n_bins*sizeof(double));
    for(int i=0;i<n_bins;i++){
        key->hist_prOri[i]=0.;
        key->ori_prOri[i]=0.;
    }
    key->nsamples_descr=0;
    return key;
}

struct kypt* malloc_kypt_from_model_and_copy(struct kypt* key){
    
    int n_hist  = key->n_hist;   /* number of histograms per direction in the feature vector */
    int n_ori   = key->n_ori;    /* number of bins in each orientation histogram */
    int n_bins  = key->n_bins;   /* number of bins in the orientation histogram for principal orientation attribution */
    int nMaxIntrp  = key->nMaxIntrp;
    
    struct kypt* keycp = malloc_kypt_extra(n_ori,n_hist,n_bins,nMaxIntrp);
    
    keycp->i = key->i;
    keycp->j = key->j;
    keycp->s = key->s;
    keycp->o = key->o;
    keycp->x = key->x;
    keycp->y = key->y;
    keycp->sigma = key->sigma;
    keycp->val   = key->val;
    keycp->theta = key->theta;
    keycp->hsEdgeResp = key->hsEdgeResp;
    keycp->n_hist = n_hist;
    keycp->n_ori  = n_ori;
    keycp->n_bins  = n_bins;
    keycp->nMaxIntrp  = nMaxIntrp;
    for(int n=0;n<n_hist*n_hist*n_ori;n++){
        keycp->descr[n] = key->descr[n];
    }
    keycp->nsamples_descr = key->nsamples_descr;
    for(int n=0;n<n_bins;n++){
        keycp->hist_prOri[n] = key->hist_prOri[n];  
    }
    keycp->n_prOri = key->n_prOri;
    keycp->nIntrp = key->nIntrp;
    return keycp;
}










struct kypt_lst* malloc_kypt_lst(){
    struct kypt_lst * keys = (struct kypt_lst*)malloc(sizeof(struct kypt_lst));
    keys->lst = (struct kypt **)malloc(20*sizeof(struct kypt*));
    keys->msize = 20;
    keys->card=0;
    return keys;
}
    
void realloc_kypt_lst(struct kypt_lst* keys){
    keys->lst = (struct kypt**)realloc(keys->lst, 2*keys->msize*sizeof(struct kypt*));
    keys->msize = 2*keys->msize;
}


void add_kypt_to_lst(struct kypt* key, struct kypt_lst* keys){
    if(keys->card > keys->msize-5){ /* checks if there is enough memory allocated */
        realloc_kypt_lst(keys);
    }
    keys->lst[keys->card]=key;
    keys->card+=1;
}



void free_kypt(struct kypt* key){
    free(key->descr);   
    free(key);    
}

void free_kypt_extra(struct kypt* key){
    free(key->pos_intrp);
    free(key->offset_intrp);
    free(key->hist_prOri);
    free(key->ori_prOri);
    free_kypt(key);
}

// void free_kypt_lst(struct kypt_lst* keys){
//     for(int k=0;k<keys->card;k++){
//         free_kypt(keys->lst[k]);
//     }
//     free(keys->lst);
//     free(keys);
// }

void free_kypt_lst(struct kypt_lst* keys){
    for(int k=0;k<keys->card;k++){
        free_kypt_extra(keys->lst[k]);
    }
    free(keys->lst);
    free(keys);
}


/** @brief Save in ASCII file keypoint detection position.
 *
 * 
 * in  @param keys : list of keypoint structure.
 *     @param dim  : dimension of the feature vector
 * out @param in ASCII file (name) 
 *  The format is
 *    x  y  sigma  theta  fvec[1] fvec[2] ... fvec[dim]
 * 
 * 
 */
void save_kypt_lst_detection(struct kypt_lst* keys, char* name){
    FILE* file = fopen(name,"w");
    struct kypt * key;
    for(int k=0;k<keys->card;k++){
        key = keys->lst[k];
        fprintf(file,"%12.5g%12.5g%12.5g%12.5g%7i%7i\n",key->x,
                                                        key->y,
                                                        key->sigma,
                                                        key->theta,
                                                        key->o,
                                                        key->s);
    }
    fclose(file);
}



/** @brief Save in ASCII file keypoint standard information.
 *
 * 
 * in  @param keys : list of keypoint structure.
 *     @param dim  : dimension of the feature vector
 * out @param in ASCII file (name)
 * 
 *  The format is 
 *    x  y  sigma  theta  fvec[1] fvec[2] ... fvec[dim]
 * 
 * 
 */
void save_kypt_lst_description(struct kypt_lst* keys, int dim, char* name){
    FILE* file = fopen(name,"w");
    struct kypt * key;
    for(int k=0;k<keys->card;k++){
        key = keys->lst[k];
        fprintf(file,"%12.5f%12.5f%12.5f%12.5f",key->x,
                                                key->y,
                                                key->sigma,
                                                key->theta);
        for(int n=0;n<dim;n++){
            fprintf(file,"%10.3i", (int)key->descr[n]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}



/** @brief Save in ASCII file keypoint standard and extra information.
 * 
 * in @param keys : list of keypoint structure.
 *    @param name : name of ASCII file.
 * 
 *  The format is
 *    x  y  sigma  theta  fvec[1] ... fvec[dim] oct sca  orihist[1] ... orihist[n_bins]
 * 
 * 
 * @param dim    : dimension of the feature vector,
 * @param n_bins : number of bins in the orientation histogram.
 *
 * 
 */
void save_kypt_lst_extra(struct kypt_lst* keys, int dim, int n_bins, char* name){
    FILE* file = fopen(name,"w");
    struct kypt * key;
    for(int k=0;k<keys->card;k++){
        key = keys->lst[k];
        // the standard output...
        fprintf(file,"%12.5f%12.5f%12.5f%12.5f",key->x,
                                                key->y,
                                                key->sigma,
                                                key->theta);
        for(int n=0;n<dim;n++){
            fprintf(file,"%10.3i", (int)key->descr[n]);
        }
        fprintf(file,"                           ");
        
        // ...plus extra information.
        fprintf(file,"%10i%10i",key->o,key->s);
        for(int n=0;n<n_bins;n++){
            fprintf(file,"%10.3f", key->hist_prOri[n]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}





/** @brief Print in standard output keypoints list
 * 
 *     @param dim : dimension of the feature vector
 * 
 * outputs 
 *    x  y  sigma theta  fvec[1]...fvec[dim]
 * 
 */
void print_kypt_lst(struct kypt_lst* keys, int dim){
    for(int k=0;k<keys->card;k++){
        printf("%12.5f%12.5f%12.5f%12.5f",keys->lst[k]->x,
                                          keys->lst[k]->y,
                                          keys->lst[k]->sigma,
                                          keys->lst[k]->theta);
        for(int n=0;n<dim;n++){
            printf("%10.3i",(int)keys->lst[k]->descr[n]);
        }
        printf("\n");
    }
}



/** @brief Print in standard output keypoints list
 * 
 * @param dim    : dimension of the feature vector
 * @param n_bins : number of bins in the orientation histogram
 * 
 * outputs
 *    x  y  sigma theta  fvec[1]...fvec[dim]  o  s  orihist[1]...orihist[dim]
 * 
 */
void print_kypt_lst_extra(struct kypt_lst* keys, int dim, int n_bins){
    for(int k=0;k<keys->card;k++){
        printf("%12.5f%12.5f%12.5f%12.5f",keys->lst[k]->x,
                                          keys->lst[k]->y,
                                          keys->lst[k]->sigma,
                                          keys->lst[k]->theta);
        for(int n=0;n<dim;n++){
            printf("%10.3i",(int)keys->lst[k]->descr[n]);
        }
        // and extra information
        printf("     %10i %10i    ",keys->lst[k]->o,keys->lst[k]->s);
        for(int n=0;n<n_bins;n++){
            printf(" %10.3f ", keys->lst[k]->hist_prOri[n]);
        }
        printf("\n");
    }
}



/** @brief read keypoint list from txt file
 *
 * @param keys    =  output list of keypoints
 * @param name    =  input filename
 * 
 * @param n_hist  the descriptor has (n_hist*n_hist) weighted coefficients.
 * @param n_ori   the descriptor has (n_hist*n_hist*n_ori) coefficients
 * 
 */
void read_kypt_lst(struct kypt_lst* keys, char* label, int n_hist, int n_ori){
    
    char text[4096];
    FILE* stream;
    stream = fopen(label,"r");
    int pos, read;
  
    while(fgets(text, 4096, stream) != NULL){
        pos = 0; read = 0;
//         struct kypt* key = malloc_kypt(n_ori,n_hist);
        struct kypt* key = malloc_kypt_extra(n_ori, n_hist, 0, 0);
        
        sscanf(text+pos,"%lf  %lf  %lf  %lf %n",&key->x    /* read coordinates */
                                               ,&key->y
                                               ,&key->sigma
                                               ,&key->theta
                                               ,&read);  pos+=read;
        
        for(int i=0;i<n_hist*n_hist*n_ori;i++){            /* read descriptor */
            sscanf(text+pos,"%lf  %n",&(key->descr[i]),&read);
            pos +=read;
        }
        add_kypt_to_lst(key,keys);
    }
    save_kypt_lst_description(keys,n_hist*n_hist*n_ori, "test_read.txt");
}




/** @brief read keypoint list from txt file
 *
 * keys    =  output list of keypoints
 * name    =  input filename
 * 
 * n_hist  
 * n_ori   ,  descriptor has (n_hist*n_hist*n_ori) coefficients
 * 
 * n_bins  ,  number of bins in the orientation histogram.
 * 
 * 
 * EXPECTED data format 
 * 
 *  x y sigma theta {fv_i}_i\in[1,n_hist^2*n_ori] oct sca {ori_i}_i\in[1,n_bins]
 * 
 * 
 */
void read_kypt_lst_extra(struct kypt_lst* keys,
                         char* label,
                         int n_hist,
                         int n_ori,
                         int n_bins){

                       /* to avoid buffer overflow */
    char text[64096];  /*TEMP  BIG - high values of n_hist, n_ori and n_bins lead to long txt lines */
    FILE* stream;
    stream = fopen(label,"r");
    int pos, read;
    
    /* Loop over lines (each line is a keypoint) */ 
    while(fgets(text, 64096, stream) != NULL){
        pos = 0; read = 0;
        struct kypt* key = malloc_kypt_extra(n_ori,n_hist,n_bins,5);
        
        /* Read keypoint coordinates */
        sscanf(text+pos,"%lf  %lf  %lf  %lf  %n",&key->x,
                                                 &key->y,
                                                 &key->sigma,
                                                 &key->theta,
                                                 &read); 
        pos+=read;
        
        /* Read keypoint feature vector */
        for(int i=0;i<n_hist*n_hist*n_ori;i++){  /* descriptor */
            sscanf(text+pos,"%lf  %n",&(key->descr[i]),&read);
            pos +=read;
        }
        
        /* Read keypoint octave and scale index */
        sscanf(text+pos," %i %i %n",&key->o,&key->s,&read); pos+=read; 
        
        /* Read keypoint orientation histogram */
        for(int i=0;i<n_bins;i++){     /* orientation histogram */ 
            sscanf(text+pos,"  %lf  %n",&(key->hist_prOri[i]),&read);
            pos +=read;
        }
        
        /* Add keypoint to the list of keypoints */
        add_kypt_to_lst(key,keys);
    }
}












struct intTrio_lst* malloc_intTrio_lst(){
    struct intTrio_lst* coords = (struct intTrio_lst*)malloc(sizeof(struct intTrio_lst));
    coords->lst = (struct intTrio**)malloc(200*sizeof(struct intTrio*));
    coords->card=0;
    coords->msize=200;
    return coords;
}

void realloc_intTrio_lst(struct intTrio_lst* coords){
    coords->lst = (struct intTrio**)realloc(coords->lst, (2*coords->msize+10)*sizeof(struct intTrio*));
    coords->msize = 2*coords->msize;
}

void add_intTrio_to_lst(struct intTrio* coord, struct intTrio_lst * coords){
    if (coords->card > coords->msize-5){
        realloc_intTrio_lst(coords);
    }
    coords->lst[coords->card]=coord;
    coords->card +=1;
}

void free_intTrio_lst(struct intTrio_lst * coords){
    int k;
    for(k=0;k<coords->card;k++){
        free(coords->lst[k]);
    }
    free(coords->lst);
    free(coords);
}
