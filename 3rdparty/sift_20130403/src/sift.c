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
 * @file sift.c
 * @brief [[MAIN]] The SIFT method 
 *
 * @li basic SIFT transform applied to one image
 * @li verbose SIFT transform 
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */














#include "sift_main.h"
#include "sift_matching.h"
#include "io_png.h"


/** @brief Main SIFT routine
 * 
 * takes one image as input.        
 * outputs the extracted keypoints in the standard output.
 * 
 * @param flag for the SIFT transform
 * 
 * 0 = one image -> one txt file 
 * 1 = one image -> all txt files
 * 2 = one image -> all txt files + scalespace and DoG
 * 
 */
int main(int argc, char **argv){

    
    int verbflag;
    char* label;
    
    int n_oct;
    int n_spo;
    double sigma_min;
    double delta_min;
    double sigma_in;
    double C_DoG;
    double C_edge;
    int n_bins;         
    double lambda_ori;  
    double t;           
    int n_hist;
    int n_ori;
    double lambda_descr;
    
    
    switch (argc){
        case 2:
            /* STANDARD SIFT : standard parameters, single output */
            verbflag = 0;
            /* standard parameters */
            n_oct        = 8;          
            n_spo        = 3;
            sigma_min    = 0.8;
            delta_min    = 0.5;
            sigma_in     = 0.5;
            C_DoG        = 0.03;
            C_edge       = 10;
            n_bins       = 36;
            lambda_ori   = 1.5;
            t            = 0.8;           
            n_hist       = 4;
            n_ori        = 8;
            lambda_descr = 6;
            break;
        case 16:
            /* CUSTOMIZABLE SIFT : lists and scalespace */
            verbflag = 2;
            label        = argv[2];
            /* user's choice of parameters */
            n_oct        = atoi(argv[3]);      
            n_spo        = atoi(argv[4]);      
            sigma_min    = atof(argv[5]);      
            delta_min    = atof(argv[6]);      
            sigma_in     = atof(argv[7]);      
            C_DoG        = atof(argv[8]);      
            C_edge       = atof(argv[9]);      
            n_bins       = atoi(argv[10]);      
            lambda_ori   = atof(argv[11]);      
            t            = atof(argv[12]);      
            n_hist       = atoi(argv[13]);      
            n_ori        = atoi(argv[14]);      
            lambda_descr = atof(argv[15]); 
            break;
        default:
            printf("\n SIFT transform of an image \n\n");
            printf("  sift image \n");
            printf("            algorithm with standard parameters\n");
            printf("            output list : x y sigma theta featurevector[]\n\n");
            printf("  sift image label  n_oct  n_spo  sigma_min  delta_min  sigma_in"); 
            printf(" C_DoG  C_edge  n_bins  lambda_ori  t   n_hist n_ori  lambda_descr\n");
            printf("                    8      3      0.8        0.5        0.5     "); 
            printf(" 0.03   10      36      1.5         0.8 4      8      6\n");
            printf("            algorithm with customizable parameters\n");
            printf("            output list : x y sigma theta featurevector[] histogram[] octa sca\n\n ");
            return -1;
    }
    
    if(sigma_min<sigma_in){
        sigma_min = sigma_in;
    }
    

    /** Loading image */
    size_t width, height;
    float* imf32 = io_png_read_f32_rgb(argv[1], &width, &height);
    double* image = (double*)malloc(width*height*sizeof(double));
    for(unsigned int i = 0;i<width*height;i++)
        image[i] = (double)((imf32[i]+imf32[width*height+i]+imf32[2*width*height+i])/3./256.); //RGB2GRAY
    
    /** List of keypoints in image */
    struct kypt_lst* keys = malloc_kypt_lst();   /* keypoints */
    
    /** Adapting the value of the C_DoG to be equivalent to what it would be for n_spo = 3 */
    double C_DoG_adapt = (exp(log(2)/(double)n_spo)-1)/(exp(log(2)/(double)3)-1)*C_DoG;
    
    /** Running the SCALE INVARIANT FEATURE TRANSFORM on the image 1 */
    sift_transform(image, width, height,
                   keys,label, verbflag,
                   n_oct,n_spo,sigma_min,delta_min,sigma_in,
                   C_DoG_adapt, C_edge,
                   n_bins,lambda_ori,t,n_hist,n_ori,lambda_descr); 
    
    free(image);
    free_kypt_lst(keys);
    return 0;
}
