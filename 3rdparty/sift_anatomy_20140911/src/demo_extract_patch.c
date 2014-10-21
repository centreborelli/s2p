/*
IPOL SIFT
Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20140911 (September 11th, 2014)

This C ANSI source code is related to the IPOL publication

    [1] "Anatomy of the SIFT Method."
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2013.
        http://www.ipol.im/pub/algo/rd_anatomy_sift/

An IPOL demo is available at
        http://www.ipol.im/pub/demo/rd_anatomy_sift/



== Patent Warning and License =================================================

The SIFT method is patented 

    [2] "Method and apparatus for identifying scale invariant features
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
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "lib_discrete.h"
#include "lib_util.h"
#include "io_png.h"

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define ABS(x) ((x)<0?-(x):(x))



/** @brief Image subsampling by an integer factor k
 * 
 *  Computes the gradient via the centered finite difference scheme
 *  [-1/2,0,+1/2]
 * 
 * \param in [i,j] , (0 <= i <= h-1) (0 <= j <= w-1)  . 
 * \param out [i,j]= in[k*i,k*j] ,  (0 <= i <= int(h/k)-1) (0 <= j <= int(w/k)-1)
 * 
 */
void subsample_by_intfactor(const float* in, float* out, int wi, int hi, int factor)
{
    int wo = wi/factor;
    int ho= hi/factor;
    for(int i=0;i<ho;i++){
        int i_p=factor*i;
        for(int j=0;j<wo;j++){
            int j_p=factor*j;
            out[i*wo+j] = in[i_p*wi+j_p];
        }
    }
}



void printImage_LinearConversion(const float *im, int w, int h, const char *name)
{
    float *imtmp = xmalloc(w*h*sizeof(float));
    linear_conversion(im, imtmp, w * h);
    for(int i =0; i<w*h; i++){ imtmp[i] = 255*imtmp[i]; }
    io_png_write_f32(name, imtmp, w, h, 1);
    xfree(imtmp);
}

void thumbnail(float x_key, float y_key, float sigma_key, float theta_key,
               const float* image, int width, int height,
               float lambda_descr,
               const char* name_out)
{

    float Xi,Yj;        // thumbnail coordinates.
    float x,y;          // image coordinates.
    int im,jm,ip,jp;     // auxiliary, for the bilinear interpolation.


    int  w_thumb = (int)(2*lambda_descr*sigma_key);
    int  h_thumb = w_thumb;
    float* thumb = xmalloc(w_thumb*h_thumb*sizeof(float));

    for(int k=0;k<w_thumb*h_thumb;k++){thumb[k]=0.0;}

    float ct = cos(theta_key);
    float st = sin(theta_key);

    for(int i=0;i<h_thumb;i++){
        for(int j=0;j<w_thumb;j++){

            Xi = (float)i-(float)w_thumb/2.0;
            Yj = (float)j-(float)h_thumb/2.0;

            x = x_key + Xi * ct - Yj * st;
            y = y_key + Xi * st + Yj * ct;

            /// Compute the image value according to a bilinear interpolation model.
            im=(int)x;   ip=(int)x+1;
            jm=(int)y;   jp=(int)y+1;
            //image extension with zeros
            if((im>=0) && (im<height) && (jm>0) && (jm<width)){
                  thumb[i*w_thumb+j] = (x-im)*(y-jm)*image[ip*width+jp]
                                     + (x-im)*(jp-y)*image[ip*width+jm]
                                     + (ip-x)*(y-jm)*image[im*width+jp]
                                     + (ip-x)*(jp-y)*image[im*width+jm];
            }
        }
    }
    printImage_LinearConversion(thumb,w_thumb,h_thumb,name_out);
    xfree(thumb);
}




int main(int argc, char **argv)
{

    if(argc != 14) {
        fprintf(stderr," Illustrates: - the reference orientation computation \n");
        fprintf(stderr,"              - the descriptor computation \n");
        fprintf(stderr,"Usage: extract image,    \n");
        fprintf(stderr,"                   x, y, sigma, theta, \n");
        fprintf(stderr,"                               delta_min, sigma_min, sigma_in, n_spo, \n");
        fprintf(stderr,"                                               lambda_ori, lambda_descr, nhist, Name \n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     outputs : - [Name]_thumbnail_ori_hist.png\n");
        fprintf(stderr,"               - [Name]_thumbnail_weighted_hists.png\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"extract image x y sigma theta 0.5 0.8 0.5 3 1.5 6 4 Name \n");
        return -1;
    }

    char filename[FILENAME_MAX];

    /** Read input image */
    size_t width, height;
    float* image = io_png_read_f32_rgb(argv[1], &width, &height);

    /** Read parameters 
     * (keypoint coordinates and method parameters */
    float x             = atof(argv[2]);
    float y             = atof(argv[3]);
    float sigma         = atof(argv[4]);
    float theta         = atof(argv[5]);
    float delta_min     = atof(argv[6]);
    float sigma_min     = atof(argv[7]);
    float sigma_in      = atof(argv[8]);
    int nspo            = atoi(argv[9]);
    float lambda_ori    = atof(argv[10]);
    float lambda_descr  = atof(argv[11]);
    int nhist           = atoi(argv[12]);

    // compute octave and scale
    int oct, sca;
    int a = (int)(round( nspo * log( sigma / sigma_min) /M_LN2  ));
    oct = (a-1)/nspo;
    if (oct > -1){
        sca = (a-1)%nspo + 1;
    }
    else{
        oct = 0;
        sca = 0;
    }



    /** Compute the seed image of the scale-space */    // We assume here that delta_min = 2^{-iter}
    int width_seed  = (int)(width/delta_min);
    int height_seed = (int)(height/delta_min);
    float* seed    = xmalloc(width_seed*height_seed*sizeof(float));
    assert(delta_min<=1);
    sift_oversample_bilin(image , width, height, seed, width_seed, height_seed, delta_min);

    /** Compute the image (o,s) */
    float delta_o = delta_min*pow(2,oct);
    int width_o = (int)((float)width_seed/pow(2,oct));
    int height_o     = (int)((float)height_seed/pow(2,oct));
    float* im_o_s   = xmalloc(width_o*height_o*sizeof(float));
    float sigma_o_s = delta_o*sigma_min/delta_min*pow(2,(float)sca/(float)nspo);
    //
    float* im_temp  = xmalloc(width_seed*height_seed*sizeof(float));
    /* Gaussian blur */
    float sig_blur  = sqrt(sigma_o_s*sigma_o_s - sigma_in*sigma_in);
    sift_add_gaussian_blur(seed, im_temp, width_seed, height_seed, sig_blur);
    /* Sub-sampling */
    subsample_by_intfactor(im_temp, im_o_s, width_seed, height_seed, pow(2,oct));

    snprintf(filename, FILENAME_MAX, "%s_thumbnail_ori_hist.png", argv[13]);

    /** Thumbnail of P_ori */
    thumbnail(x/delta_o, y/delta_o, sigma/delta_o, 0.0,
              im_o_s, width_o, height_o, 3*lambda_ori,
              filename);

    /** Thumbnail of P_descr */
    sprintf(filename,"%s_thumbnail_weighted_hists.png",argv[13]);
    thumbnail(x/delta_o, y/delta_o, sigma/delta_o, theta,
              im_o_s, width_o, height_o, (nhist+1)*lambda_descr/nhist,
              filename);

}








