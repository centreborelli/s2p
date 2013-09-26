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
 * @file sift_scalespace.c
 * @brief data structures to store the scale-space
 *
 * @li struct kypt      keypoint data structure.
 * @li struct kypt_lst  list of keypoints with variable length.
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */






#include "sift_scalespace.h"



/** ************************************ ALLOCATION *******************************************/


/** @brief Octave structure memory allocation
 *
 *  An octave is a stack of images sharing same width, height and intersample distance.
 *  Each image of the stack has a level of blur.
 * 
 * @param nSca   : number of images in the stack.
 * @param delta  : intersample distance
 * @param width 
 * @param height
 * @param sigma  : levels of blur for each image in the stack
 *  
 *
 */
struct octa *malloc_octa(double delta, int width, int height, int nSca, double* sigmas){
    
    /* Memory allocation */
    struct octa* octave = (struct octa*)malloc(sizeof(struct octa));   /* pointer to a (structure octave) */

    /* Modifying structure members */
    octave->delta  = delta;
    octave->width  = width;
    octave->height = height;
    octave->nSca   = nSca;
    octave->sigmas  = (double*)malloc(nSca*sizeof(double));                /* pointer to the array of level of simulated blurs */
    octave->imStack = (double*)malloc(nSca*width*height*sizeof(double));   /* pointer to the array of pixels  */
    for(int i=0;i<nSca;i++){
        octave->sigmas[i] = sigmas[i];
    }
    return octave;
}

/** @brief allocation and attributes modification
 *
 *       Allocate memory for a scalespace structure 
 *       and modify attributes that will be used to actually build the scalespace
 * 
 * @param nOct     : number of octaves
 * @param deltas   : intersample distance in each octave of the scalespace
 * @param widths   : image width in each octave of the scalespace
 * @param heights  :
 * @param nScas    : number of images in each octave of the scalespace
 * @param sigmas   : array of blurs for each image in each octave of the scalespace
 * 
 */
struct scsp* malloc_scsp(int nOct, double* deltas, int* widths, int* heights, int* nScas, double** sigmas){
    // a malloc for ONE element
    struct scsp* scalespace = (struct scsp *)malloc(sizeof(struct scsp));   // pointer to a (struct scalespace)
    scalespace->nOct = nOct;
    for(int o=0;o<nOct;o++){
        scalespace->octaves[o] = malloc_octa(deltas[o], widths[o], heights[o], nScas[o], sigmas[o]);
    }   
    return scalespace;
}



void free_octa(struct octa *octave){
    free(octave->sigmas);    /* the memory allocated to store simulated values */
    free(octave->imStack);   /* the memory allocated to store the octave samples */
    free(octave);          
}


//return value de free
void free_scsp(struct scsp *scalespace){
    int nOct = scalespace->nOct;
    for(int o=0;o<nOct;o++)
        free_octa(scalespace->octaves[o]);
    free(scalespace);
}



/** **********************************  ALLOCATION WRAPPERS ****************************************/


/** @brief Allocation via copy a scalespace structure
 * 
 * - Copy structure (number of octaves, dimensions, sampling rates, levels of blur)
 * - Doesn't copy the image stacks
 * 
 */
struct scsp* malloc_scsp_from_model(struct scsp* model_scsp){
    
    struct scsp * copy_scsp;
    
    int nOct = model_scsp->nOct;
    
    int* nScas   = (int*)malloc(nOct*sizeof(int));
    int* widths  = (int*)malloc(nOct*sizeof(int));
    int* heights = (int*)malloc(nOct*sizeof(int));
    double* deltas  = (double*)malloc(nOct*sizeof(double));
    double** sigmas = (double**)malloc(nOct*sizeof(double*));
    
    for(int o=0;o<nOct;o++){
        nScas[o]   = model_scsp->octaves[o]->nSca;
        widths[o]  = model_scsp->octaves[o]->width;
        heights[o] = model_scsp->octaves[o]->height;
        deltas[o]  = model_scsp->octaves[o]->delta;
        
        sigmas[o]  = (double*)malloc(nScas[o]*sizeof(double));
        for(int i=0;i<nScas[o];i++)
            sigmas[o][i] = model_scsp->octaves[o]->sigmas[i];
    }
    
    copy_scsp = malloc_scsp(nOct, deltas, widths, heights, nScas, sigmas);
    
    free(deltas);
    free(widths);
    free(heights);
    free(nScas);
    for(int o=0;o<nOct;o++)
        free(sigmas[o]);
    free(sigmas);    
    
    return copy_scsp;
}



/** @brief Allocation of a scalespace with David Lowe's original settings
 *
 */
struct scsp* malloc_scsp_lowe_scalespace(int nOct,            /* # of octaves            */
                                         int nSca,                      /* # of scales of detection (WARNING different from the # number of allocated and computed scale) */
                                         int im_width, int im_height,   /* # input image dimension */
                                         bool dozoomin,                 /* bool flag for prior zoom     */
                                         double sigma_min){             /* minimal scale in each octave (relatively to the sampling rate) */
    
                                         
    int* widths     = (int*)malloc(nOct*sizeof(int));
    int* heights    = (int*)malloc(nOct*sizeof(int));
    int* nScas      = (int*)malloc(nOct*sizeof(int));
    double* deltas  = (double*)malloc(nOct*sizeof(double));
    double** sigmas = (double**)malloc(nOct*sizeof(double*));    /* nSca+3 = nSca + 2 extra scale to search 3d extrema + 1 extra scale to compute DoG */
    
    if(dozoomin){
        widths[0] = 2*im_width;
        heights[0] = 2*im_height;
        deltas[0] = 0.5;
    }
    else{
        widths[0] = im_width;
        heights[0] = im_height;
        deltas[0] = 1;
    }
    
    for(int o=1;o<nOct;o++){
            widths[o] = widths[o-1]/2;     /*integer division*/
            heights[o] = heights[o-1]/2;
            deltas[o] = deltas[o-1]*2.0;
    }
    
    for(int o=0;o<nOct;o++){
        nScas[o] = nSca+3;  /* 3 extra images in the stack, 1 for dog computation and 2 for 3d discrete extrema definition*/
        sigmas[o] = (double*)malloc(nScas[o]*sizeof(double));
        for(int s=0;s<nSca+3;s++){
            
            sigmas[o][s] = deltas[o]/deltas[0]*sigma_min*pow(2.0,(double)s/(double)nSca);
            
        }
    }
    
    struct scsp* scalespace = malloc_scsp(nOct, deltas, widths, heights, nScas, sigmas);   
    
    free(deltas);
    free(widths);
    free(heights);
    free(nScas);
    for(int o=0;o<nOct;o++)
        free(sigmas[o]);
    free(sigmas);    
    return scalespace;
}








/** @brief Allocation of a scalespace with David Lowe's original structure.
 * 
 * 
 * The set of simulated level of blur
 *    sigmas[o][s] = deltas[o]*sigma_min*pow(2.0,(double)s/(double)nSca);
 * 
 *
 */
struct scsp* malloc_scsp_lowe_scalespace_bis(int nOct,                 /* # of octaves            */
                                         int nSca,                     /* # of scales of detection (excluding the 3 auxiliary scales */
                                         int im_width, int im_height,  /* # input image dimension */
                                         double delta_min,             /* minimal inter distance sample     */
                                         double sigma_min){            /* minimal scale in each octave (relatively to the sampling rate) */
    
                                         
    int* widths     = (int*)malloc(nOct*sizeof(int));
    int* heights    = (int*)malloc(nOct*sizeof(int));
    int* nScas      = (int*)malloc(nOct*sizeof(int));
    double* deltas  = (double*)malloc(nOct*sizeof(double));
    double** sigmas = (double**)malloc(nOct*sizeof(double*));    /* nSca+3 = nSca + 2 extra scales to search 3d extrema + 1 extra scale to compute DoG */
    
    
    assert(delta_min <=1);
    deltas[0] = delta_min;
    heights[0] = (int)(im_height/delta_min);
    widths[0] = (int)(im_width/delta_min);
    
    
    for(int o=1;o<nOct;o++){
        widths[o] = widths[o-1]/2;     /*integer division*/
        heights[o] = heights[o-1]/2;
        deltas[o] = deltas[o-1]*2.0;
    }
    
    
    for(int o=0;o<nOct;o++){
        nScas[o] = nSca+3;  /* 3 extra images in the stack, 1 for dog computation and 2 for 3d discrete extrema definition*/
        sigmas[o] = (double*)malloc(nScas[o]*sizeof(double));
        for(int s=0;s<nSca+3;s++){ /* nSca images + 3 auxiliary images*/
            
            
            sigmas[o][s] = deltas[o]/deltas[0]*sigma_min*pow(2.0,(double)s/(double)nSca);
                     
            
        }
    }
    
    struct scsp* scalespace = malloc_scsp(nOct, deltas, widths, heights, nScas, sigmas);   
    
    free(deltas);
    free(widths);
    free(heights);
    free(nScas);
    for(int o=0;o<nOct;o++)
        free(sigmas[o]);
    free(sigmas);    
    return scalespace;
}












/** @brief allocation for Difference of Gaussian scalespace.
 * 
 *   Characteristics of a DoG scalespace
 *      - one image less in each octave than the scalespace it is based on
 *      - same level of blurs than the scalespace it is based on (this is an approximation)
 * 
 */
struct scsp* malloc_scsp_dog_from_scalespace(struct scsp* scalespace){           /* minimal scale in each octave (relatively to the sampling rate) */
    
    /* copy and allocate using a model */
    struct scsp* dog = malloc_scsp_from_model(scalespace);
    
    /* modify some allocations */
    int nOct = scalespace->nOct;
    int nSca;
    int width;
    int height;
    
    for(int o=0;o<nOct;o++){
        
        free(dog->octaves[o]->sigmas);
        free(dog->octaves[o]->imStack);
        
        nSca   = scalespace->octaves[o]->nSca-1;
        width  = scalespace->octaves[o]->width;
        height = scalespace->octaves[o]->height;
        
        dog->octaves[o]->nSca = nSca;
        dog->octaves[o]->sigmas  = (double*)malloc(nSca*sizeof(double));
        dog->octaves[o]->imStack = (double*)malloc(nSca*width*height*sizeof(double));
        
        for(int i=0;i<nSca;i++){
            dog->octaves[o]->sigmas[i] = scalespace->octaves[o]->sigmas[i];
        }
    }
    return dog;

}




void print_scsp_gray(struct scsp* scalespace, char* basename){
    
    struct octa* octave;
    double* image;
    char name[FILENAME_MAX];
    
    int nOct = scalespace->nOct;
    int nSca, width, height;
    double sigma;
    for(int o=0;o<nOct;o++){
        octave = scalespace->octaves[o];
        nSca   = octave->nSca;
        width  = octave->width;
        height = octave->height;
        for(int s=0;s<nSca;s++){
            image = &octave->imStack[s*width*height];
            sigma = octave->sigmas[s];
            sprintf(name,"%s_oct%iof%i_scale%iof%i_sigma%f.png",basename,o,nOct,s,nSca,sigma);
            printImage_LinearConversion(image, width,height, name);            
        }
    }
}


void nn_interp(double* in,int win, int hin, double* out, int wout,int hout){
    int k,l;
    for(int i=0;i<hout;i++){
        k = (int)((double)i/(double)hout*(double)hin);
        assert(k>-1); assert(k<hin);
        for(int j=0;j<wout;j++){
            l = (int)((double)j/(double)wout*(double)win);
            assert(l>-1); assert(l<win);
            out[i*wout+j] = in[k*win+l];
        }
    }       
}


void print_scsp_gray_nearestneighor(struct scsp* scalespace, char* basename){
    
    /* pointers for dereferencing */
    struct octa* octave;
    double* image;
    
    char name[FILENAME_MAX];
    
    int nOct = scalespace->nOct;
    int nSca, width, height;
    double sigma;
    int widthALL = scalespace->octaves[0]->width;
    int heightALL = scalespace->octaves[0]->height;
    double* imtemp = (double*)malloc(widthALL*heightALL*sizeof(double));
    double* imrgb = (double*)malloc(3*widthALL*heightALL*sizeof(double));
    for(int o=0;o<nOct;o++){
        octave = scalespace->octaves[o];
        nSca   = octave->nSca;
        width  = octave->width;
        height = octave->height;
        for(int s=0;s<nSca;s++){
            image = &octave->imStack[s*width*height];
            nn_interp(image,width,height,
                      imtemp,widthALL,heightALL);
            sigma = octave->sigmas[s];
            sprintf(name,"%s_o%i_s%i_delta%f_sigma%f.png",basename,o+1,s,octave->delta,sigma);
            sprintf(name,"%s_o%03i_s%03i.png",basename,o,s);
            printImage_LinearConversion(imtemp, widthALL,heightALL, name);        
        }
    }
    free(imtemp);
    free(imrgb);
}


void print_scsp_rgb(struct scsp* scalespace, char* basename){
    
    /* pointers for dereferencing */
    struct octa* octave;
    double* image;
    
    char name[FILENAME_MAX];

    int nOct = scalespace->nOct;    
    int nSca, width, height;
    double sigma;
    
    int widthALL = scalespace->octaves[0]->width;
    int heightALL = scalespace->octaves[0]->height;
    double* imtemp = (double*)malloc(widthALL*heightALL*sizeof(double));
    double* imrgb = (double*)malloc(3*widthALL*heightALL*sizeof(double));
    
    double octamin, octamax;
    
    /* measuring max and min value in the entire scalespace */
    octamax = -100000;
    octamin = +100000;
    for(int o=0;o<nOct;o++){
        octave = scalespace->octaves[o];
        nSca   = octave->nSca;
        width  = octave->width;
        height = octave->height;
        for(int k=0;k<nSca*width*height;k++){
            if (octamax < octave->imStack[k]){octamax = octave->imStack[k];}
            if (octamin > octave->imStack[k]){octamin = octave->imStack[k];}
        }   
    }
        
    for(int o=0;o<nOct;o++){
        octave = scalespace->octaves[o];
        nSca   = octave->nSca;
        width  = octave->width;
        height = octave->height;
        for(int s=1;s<nSca-1;s++){ // temp : excluding the auxiliary images          
            image = &octave->imStack[s*width*height];
            nn_interp(image,width,height,imtemp,widthALL,heightALL);
            
            gray2hsv(imtemp,imrgb,widthALL,heightALL,octamin,octamax);
            
            sigma = octave->sigmas[s];
            sprintf(name,"%s_o%03i_s%03i_sigma_%07f.png",basename,o,s,sigma);
            sprintf(name,"%s_o%03i_s%03i.png",basename,o,s);
            printrgb(imrgb, widthALL,heightALL, name);        
        }
    }
    
    free(imrgb);
    free(imtemp); 
}
