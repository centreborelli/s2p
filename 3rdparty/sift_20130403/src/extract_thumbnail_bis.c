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
 * @file extract_thumbnail_bis.c
 * @brief Illustration of the descriptor computation
 *
 * @li extract the patch used for reference orientation attribution
 * @li save the corresponding  gradient field into an ASCII file
 *     
 * @li extract the patch used for descriptor computation
 * @li save the corresponding gradient field into an ASCII file
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "discrete_representation.h"
#include "io_png.h"
#include "std_image_lib.h"

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )




void thumbnail_ori(double x_key, double y_key, double sigma_key,
                 double* im, int width, int height,
                 double lambda_ori,
                 char* name_thumbnail){     // saved as an image (interpolated bilinearly)
    
   
    
    int radius = (int)(3*lambda_ori*sigma_key);
    double* thumb = (double*)malloc((2*radius+1)*(2*radius+1)*sizeof(double));
    for(int k = 0 ;k<(2*radius+1)*(2*radius+1);k++){
        thumb[k]=0.0;
    }
    
    int i,j;
    for(int si = (int)x_key-radius; si<= (int)x_key+radius; si++){
        for(int sj = (int)y_key-radius; sj <= (int)y_key+radius; sj++){
            if((si>=0)&&(si<height)&&(sj>=0)&&(sj<width)){
                i = si-((int)x_key-radius);
                j = sj-((int)y_key-radius);
                thumb[i*(2*radius+1)+j] = im[si*width+sj];
            }
        }    
    }
    
    printImage_LinearConversion(thumb,(2*radius+1), (2*radius+1),name_thumbnail);
    free(thumb);
}



/** @brief Extract the normalized descriptor in the form of a thumbnail.
 * 
 * 
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key
 * @param theta_key      keypoint coordinates (as perceived in the octave).
 * @param lambda_descr   (standard value lambda_descr=6)
 * @param Nhist          (standard value Nhist = 4)
 *                       - The descriptor width is 2*lambda_descr*sigma_key*(1+1/Nhist);
 * 
 * 
 * also requires Nhist because of the width of the descriptor. 
 * 
 * ----OUTPUT----
 * @output name_out      png image file with the thumbnail.
 * 
 * 
 * L HYPOTHESE  :  LA MEME DISTANCE INTER-PIXEL, C'EST ÇA QUI DETERMINE LA GRILLE.
 * 
 */
void thumbnail_descr(double x_key, double y_key, double sigma_key, double theta_key,
                                  double* image, int width, int height,
                                  double lambda_descr, int Nhist,
                                  char* name_out){
        
//     int i,j;             // thumbnail samples coordinates.
    double Xi,Yj;        // thumbnail coordinates.
    double x,y;          // image coordinates.
    int im,jm,ip,jp;     // auxiliary, for the bilinear interpolation.
//     int ime,jme,ipe,jpe;  // auxiliary for the image extension 
    
    
    int  w_thumb = (int)(2*lambda_descr*(1+1/(double)Nhist)*sigma_key);
    int  h_thumb = w_thumb;
    double* thumb = (double*)malloc(w_thumb*h_thumb*sizeof(double));
    
    for(int k=0;k<w_thumb*h_thumb;k++){thumb[k]=0.0;}
    
    for(int i=0;i<h_thumb;i++){       
        for(int j=0;j<w_thumb;j++){
            
            Xi = (double)i-(double)w_thumb/2.0;
            Yj = (double)j-(double)h_thumb/2.0;
            
            x = x_key+Xi*cos(theta_key)-Yj*sin(theta_key);
            y = y_key+Xi*sin(theta_key)+Yj*cos(theta_key);
                        
            /// Compute the image value according to a bilinear interpolation model.
            im=(int)x;   ip=(int)x+1;
            jm=(int)y;   jp=(int)y+1;
            
           //image extension by symmetrization
            
//             ime = im;    ipe = ip;
//             jme = jm;    jpe = jp;
//             
//             if(ime<0){ime=-1-ime;}    if(ipe<0){ipe=-1-ipe;}
//             if(jme<0){jme=-1-jme;}    if(jpe<0){jpe=-1-jpe;}
//                        
//             if(ime>height-1){ime=2*height-1-ime;}    if(jme>width-1) {jme=2*width-1-jme;}
//             if(ipe>height-1){ipe=2*height-1-ipe;}    if(jpe>width-1) {jpe=2*width-1-jpe;}
            
            //image extension with zeros
            if((im>=0)&&(im<height)&&(jm>0)&&(jm<width)){
                  thumb[i*w_thumb+j] = (x-im)*(y-jm)*image[ip*width+jp]
                                     + (x-im)*(jp-y)*image[ip*width+jm]
                                     + (ip-x)*(y-jm)*image[im*width+jp]
                                     + (ip-x)*(jp-y)*image[im*width+jm];
            
            }
            
               
                               
        }
    }
    printImage_LinearConversion(thumb,w_thumb,h_thumb,name_out);
    free(thumb);
}



/** @brief Extract normalized patch P_ori, patch for the principal orientation attribution.
 * 
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key      keypoint coordinates.
 * 
 * @param gradX
 * @param gradX          precomputed gradient relative to the nearest scale.
 * 
 * @param name_output_file
 * 
 * ----PARAM----
 * @param sigma_window   (= 1.5)   std dev of Gaussian window, the analysis of the gradient orientation is local
 * @param K              (= 3) The histogram cover a square patch of width 2*K*sigma_window (=3*1.5*sigma_key)
 * 
 * ----OUTPUT----
 * @output text          all samples in the normalized patch P_ori are stored in a txt file. (input for gnuplot)
 *
 * ----RETURN----
 * @return count         number of pixels contributing to the histogram.
 *
 */
void gradient_field_ori(double x_key, double y_key, double sigma_key,
                 double* gradX, double* gradY, int width, int height,
                 double lambda_ori,
                 char* name_grad_field){     // saved as a data txt file
    

    double M;      // gradient magnitude
    int    si,sj;  // sample coordinates in the image referential
    double sX,sY;  // sample coordinates in the keypoint referential
    double sOri;
    
    
    /** output file */
    FILE* file = fopen(name_grad_field,"w");
    
   
    /// Contributing pixels are inside a patch [siMin;siMax] X [sjMin;sjMax] of K*sigma_window*sigma_key (=4.5*sigma_key) 
    int siMin = MAX(0, (int)(x_key-3*lambda_ori*sigma_key+0.5));
    int sjMin = MAX(0, (int)(y_key-3*lambda_ori*sigma_key+0.5));
    int siMax = MIN((int)(x_key+3*lambda_ori*sigma_key+0.5), height-1);
    int sjMax = MIN((int)(y_key+3*lambda_ori*sigma_key+0.5), width-1);
    
    /// For each pixel inside the patch.
    for(si=siMin;si<=siMax;si++){
        for(sj=sjMin;sj<=sjMax;sj++){
            
            /// Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            sX = ((double)si-x_key)/sigma_key;
            sY = ((double)sj-y_key)/sigma_key;
            
            /// Compute the gradient orientation (theta) on keypoint's invariant referential
            sOri = fmod(atan2(gradY[si*width+sj], gradX[si*width+sj]) + 10*2*PI, 2*PI);
                
            /// Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to gradients far from the keypoint.
            M = sqrt(gradY[si*width+sj]*gradY[si*width+sj]+gradX[si*width+sj]*gradX[si*width+sj])
                *exp(-(sX*sX+sY*sY)/(2*lambda_ori*lambda_ori));

            /// Save in file
            fprintf(file,"%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n", sX, sY, cos(sOri), sin(sOri), M);
            
        }
    }
    fclose(file);
}




/** @brief 
 *  
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key
 * @param theta_key      keypoint coordinates.
 * 
 * also requires Nhist because of the actual width of the descriptor. 
 * 
 * ----OUTPUT----
 * @output name_out      name of a text file where the gradient is sampled 
 * 
 */
void gradient_field_descr(double x_key, double y_key, double sigma_key, double theta_key,
                                  double* gradX, double* gradY, int width, int height,
                                  double lambda_descr, int Nhist,
                                  char* name_out){
   
    
    double M;      // gradient magnitude
    double sX,sY,sOri;        // pixel coordinates in the descriptor referential
    
    /** output file */
    FILE* file = fopen(name_out,"w");
    
    /// Contributing pixels are inside a patch [siMin;siMax]X[sjMin;sjMax] of 2*lambda_descr*sigma_key*(nhist+1)/nhist
    int siMin = MAX(0, (int)(x_key-sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5));
    int sjMin = MAX(0, (int)(y_key-sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5));
    int siMax = MIN((int)(x_key+sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5), height-1);
    int sjMax = MIN((int)(y_key+sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5), width-1);
        
    
    /// For each pixel inside the patch.
    for(int si=siMin;si<siMax;si++){
        for(int sj=sjMin;sj<sjMax;sj++){
            
            /// Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            sX = (cos(-theta_key)*((double)si-x_key)-sin(-theta_key)*((double)sj-y_key))/sigma_key;
            sY = (sin(-theta_key)*((double)si-x_key)+cos(-theta_key)*((double)sj-y_key))/sigma_key;
            
            /// Does this sample fall inside the descriptor area ?
            if((fabs(sX)<lambda_descr*(1+1/(double)Nhist))&(fabs(sY)<lambda_descr*(1+1/(double)Nhist))){                

                /// Compute the gradient orientation (theta) on keypoint's invariant referential.
                sOri = fmod(atan2(gradY[si*width+sj], gradX[si*width+sj]) - theta_key + 10*2*PI,2*PI);                
                
                /// Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to gradients far from the keypoint.
                M = sqrt(gradY[si*width+sj]*gradY[si*width+sj]+gradX[si*width+sj]*gradX[si*width+sj])
                    *exp(-(sX*sX+sY*sY)/(2*lambda_descr*lambda_descr));
                
                /// Save in output file
                fprintf(file,"%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n", sX, sY, cos(sOri), sin(sOri), M);                                                             
            }
        }
    }
    fclose(file);
}






int main(int argc, char **argv){

    if(argc != 16) {
        printf(" Illustrates: - the reference orientation computation \n");
        printf("              - the descriptor computation \n");
        printf("Usage: extract image, \n");
        printf("                   x, y, sigma, theta, o, s, \n");
        printf("                               delta_min, sigma_min, sigma_in, n_spo, \n");
        printf("                                               lambda_ori, lambda_descr, nhist, Name \n");
        printf("\n");
        printf("     outputs : - [Name]_thumbnail_ori_hist.png\n");
        printf("               - [Name]_thumbnail_weighted_hists.png\n");
        printf("               - [Name]_gradient_field_ori.t \n");
        printf("               - [Name]_gradient_field_descr.t \n");
        printf("\n");
        printf("extract image x y sigma theta oct sca 0.5 0.8 0.5 3 1.5 6 4 \n");
        return false;
    }
    
    char filename[FILENAME_MAX];
    
    /** Read input image */
    size_t width, height;
    float* imf32 = io_png_read_f32_rgb(argv[1], &width, &height);
    double* image = (double*)malloc(width*height*sizeof(double));
    for (int i = 0;i<width*height;i++)
        image[i] = (double)((30 * imf32[i] + 59* imf32[width*height+i] + 11* imf32[2*width*height+i])/100);           
    
    /** Read parameters 
     * (keypoint coordinates and method parameters */
    double x         = atof(argv[2]);   
    double y         = atof(argv[3]);
    double sigma     = atof(argv[4]);   
    double theta     = atof(argv[5]);
    int oct          = atoi(argv[6]);
    int sca          = atoi(argv[7]);
    double delta_min = atof(argv[8]);
    double sigma_min = atof(argv[9]);
    double sigma_in  = atof(argv[10]);
    int nspo         = atoi(argv[11]);
    double lambda_ori   = atof(argv[12]);
    double lambda_descr = atof(argv[13]);
    int nhist           = atoi(argv[14]);
    
    /** Compute the seed image of the scale-space */    // We assume here that delta_min = 2^{-iter}
    int width_seed  = (int)(width/delta_min);
    int height_seed = (int)(height/delta_min);
    double* seed    = (double*)malloc(width_seed*height_seed*sizeof(double));
    assert(delta_min<=1);
    oversample_bilin(image , width, height, seed, width_seed, height_seed, delta_min);
        
    /** Compute the image (o,s) */
    double delta_o   = delta_min*pow(2,oct);
    int width_o      = (int)((double)width_seed/pow(2,oct));
    int height_o     = (int)((double)height_seed/pow(2,oct));
    double* im_o_s   = (double*)malloc(width_o*height_o*sizeof(double));
    double sigma_o_s = delta_o*sigma_min/delta_min*pow(2,(double)sca/(double)nspo);
    //
    double* im_temp  = (double*)malloc(width_seed*height_seed*sizeof(double)); 
    /* Gaussian blur */
    double sig_blur  = sqrt(sigma_o_s*sigma_o_s - sigma_in*sigma_in);
    add_gaussian_blur(seed, im_temp, width_seed, height_seed, sig_blur);     //TODO maybe, for speed - replace by Fourier
    /* Sub-sampling */
    subsample_by_intfactor(im_temp, im_o_s, width_seed, height_seed, pow(2,oct));
    
    /** Compute gradient of image (o,s) */
    double* gradX = (double*)malloc(width_o*height_o*sizeof(double));
    double* gradY = (double*)malloc(width_o*height_o*sizeof(double));
    compute_gradient(im_o_s,gradX,gradY,width_o,height_o);
    
    /** Gradient field in P_ori - Normalized patch for reference orientation */
    sprintf(filename,"%s_gradient_field_ori.t",argv[15]);
    gradient_field_ori(x/delta_o, y/delta_o, sigma/delta_o,
                       gradX, gradY, width_o, height_o, lambda_ori,
                       filename);
    
    /** Gradient field in P_descr - Normalized patch for descriptor */
    sprintf(filename,"%s_gradient_field_descr.t",argv[15]);
    gradient_field_descr(x/delta_o, y/delta_o, sigma/delta_o, theta,
                         gradX, gradY, width_o, height_o, lambda_descr, nhist,
                         filename);
    
    /** Thumbnail of P_ori */
    sprintf(filename,"%s_thumbnail_ori_hist.png",argv[15]);
    thumbnail_ori(x/delta_o, y/delta_o, sigma/delta_o,
                  im_o_s, width_o, height_o, lambda_ori,
                  filename);
    
    /** Thumbnail of P_descr */
    sprintf(filename,"%s_thumbnail_weighted_hists.png",argv[15]);
    thumbnail_descr(x/delta_o, y/delta_o, sigma/delta_o, theta,
                    im_o_s, width_o, height_o, lambda_descr, nhist,
                    filename);
    
}
