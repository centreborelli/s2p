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
 * @file discrete_representation.c
 * @brief simple image transformations 
 *
 * This is a front-end to libpng, with routines to:
 *      @li Separable discrete convolutions
 *      @li Gaussian blur via discrete convolution with truncated kernel
 *      @li 2D Gradient
 *      @li Subsampling by integer factor
 *      @li bilinear interpolation
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */






#include "discrete_representation.h"






/** @brief Compute image gradient via symmetric finite difference schemes
 *
 *   image extension :  symmetrization at x=-1/2 
 * 
 */
void compute_gradient(double* im, double* im_x, double* im_y, int width, int height){

    double* p_in;
    double* p_in_p;
    double* p_in_m;
    double* p_out;
    
    /** Computing im_y      */
    /* pixels in the center */
    p_in   = &im[1];
    p_out  = &im_y[1];
    p_in_p = p_in+1;
    p_in_m = p_in-1;
    for(int i=1;i<height*width-1;i++){   /* produces false value on borders */
        *p_out = (*p_in_p-*p_in_m)*0.5;
        p_out++;
        p_in_m++;
        p_in_p++;
    }
    /*    pixels on borders - symmetrization at y=-1/2  */
    for(int i=0;i<height;i++){
        im_y[i*width]           = im[i*width+1]- im[i*width];
        im_y[i*width+(width-1)] = im[i*width+(width-1)] - im[i*width+(width-2)];
    }
    
    
    /** Computing im_x      */
    /* pixels in the center */
    p_in   = &im[width];
    p_out  = &im_x[width];
    p_in_p = p_in+width;
    p_in_m = p_in-width;
    for(int i=width;i<height*width-width;i++){   /* produces false value on borders */
        *p_out = (*p_in_p-*p_in_m)*0.5;
        p_out++;
        p_in_m++;
        p_in_p++;
    }
    /*    pixels on borders - symmetrization at x=-1/2  */
    for(int i=0;i<width;i++){
        im_x[i]                   = im[i+width]-im[i];
        im_x[(height-1)*width+i]  = im[(height-1)*width+i] - im[(height-1-1)*width+i];
    }
}



/** @brief Build mono-dimensional sampled Gaussian kernel
 *   of size 2*rad+1 and standard deviation sigma
 */
double* malloc_gaussian_kernel(double sigma, int rad){
    double* gker = (double*)malloc((2*rad+1)*sizeof(double));
    double tmp = 0;
    
    assert(sigma>=0);
    
    if(sigma>0){
        for(int i=-rad;i<=rad;i++){
            gker[rad+i] = exp(-0.5*(double)i*(double)i/sigma/sigma);
            tmp += gker[rad+i];
        }
        for(int i=-rad;i<=rad;i++)
            gker[rad+i]/=tmp;
    }else{
       for(int i=-rad;i<=rad;i++){ gker[rad+i]=0.0;}
       gker[rad]=1.;
    }
    return gker;
}       


void add_gaussian_blur(double* in, double* out, int width, int height, double sigma){
    int r_gker = (int)(4*sigma);
    double* gker = malloc_gaussian_kernel(sigma, r_gker);
    convolve(in,out,width,height,gker,r_gker,gker,r_gker);
    free(gker);
}

/** @brief sub sampling by factor 2, keeping sample (0,0) */
void subsample_by2(double* in, double* out, int in_width, int in_height){
    int out_width = in_width/2;
    int out_height= in_height/2;
    int i,j,i_p,j_p;
    for(i=0;i<out_height;i++){
        i_p=2*i;
        for(j=0;j<out_width;j++){
            j_p=2*j;
            out[i*out_width+j] = in[i_p*in_width+j_p];
        }
    }
}

/** @brief sub sampling by an integer factor, keeping sample (0,0) */
void subsample_by_intfactor(double* in, double* out, int in_width, int in_height, int factor){
    int out_width = in_width/factor;
    int out_height= in_height/factor;
    for(int i=0;i<out_height;i++){
        int i_p=factor*i;
        for(int j=0;j<out_width;j++){
            int j_p=factor*j;
            out[i*out_width+j] = in[i_p*in_width+j_p];
        }
    }
}


/** @brief Interpolate the image with a bilinear model
 * 
 *  the inter-pixel distance in the output image is delta_out
 * 
 * @param in  : input digital image with (in_width X in_height) samples.
 * @param out : output digita image with (out_width X out_height) samples,
 *            with out_width  = 2 X in_width;
 *             and out_height = 2 X out_height
 * 
 * 
 * 
 */
void oversample_by2_bilin(double* in, double* out, int in_width, int in_height){
        
    int out_width = 2*in_width;
    int out_height= 2*in_height;
    
    for(int i=0;i<out_width*out_height;i++)
        out[i] = 0.0;
    
    for(int i=0;i<in_height-1;i++){
        for(int j=0;j<in_width-1;j++){
            out[(2*i)*out_width+(2*j)]     = in[i*in_width+j];
            out[(2*i)*out_width+(2*j+1)]   = 0.5*(in[i*in_width+j]+in[i*in_width+j+1]);
            out[(2*i+1)*out_width+(2*j)]   = 0.5*(in[i*in_width+j]+in[(i+1)*in_width+j]);
            out[(2*i+1)*out_width+(2*j+1)] = 0.25*(in[i*in_width+j] + in[i*in_width+(j+1)] +in[(i+1)*in_width+j] + in[(i+1)*in_width+(j+1)]);
        }
    }
    
    /* last column */
    for(int i=0;i<in_height-1;i++){
        out[(2*i)*out_width+(2*in_width-2)]   =  in[i*in_width+(in_width-1)];
        out[(2*i)*out_width+(2*in_width-1)]   =  in[i*in_width+(in_width-1)];
        out[(2*i+1)*out_width+(2*in_width-2)] = (in[i*in_width+(in_width-1)]+in[(i+1)*in_width+(in_width-1)])/2.0;
        out[(2*i+1)*out_width+(2*in_width-1)] = (in[i*in_width+(in_width-1)]+in[(i+1)*in_width+(in_width-1)])/2.0;
    }
    
    /* last line */
    for(int j=0;j<in_width-1;j++){
        out[(2*in_height-2)*out_width+2*j]   = in[(in_height-1)*in_width+j];
        out[(2*in_height-1)*out_width+2*j]   = in[(in_height-1)*in_width+j];
        out[(2*in_height-2)*out_width+2*j+1] = (in[(in_height-1)*in_width+j]+in[(in_height-1)*in_width+j+1])/2.0;
        out[(2*in_height-1)*out_width+2*j+1] = (in[(in_height-1)*in_width+j]+in[(in_height-1)*in_width+j+1])/2.0;
    }
    
    
    /* down left corner */
    out[(2*in_height-2)*out_width+(2*in_width-2)] = in[(in_height-1)*in_width+(in_width-1)];
    out[(2*in_height-2)*out_width+(2*in_width-1)] = in[(in_height-1)*in_width+(in_width-1)];
    out[(2*in_height-1)*out_width+(2*in_width-1)] = in[(in_height-1)*in_width+(in_width-1)];
    out[(2*in_height-2)*out_width+(2*in_width-2)] = in[(in_height-1)*in_width+(in_width-1)];
    
                                           
}



/** @brief Interpolate the image with a bilinear model
 * 
 *  the inter-pixel distance in the output image is delta_out
 * 
 *  in  : input digital image with (in_width X in_height) samples.
 *  out : output digita image with (out_width X out_height) samples,
 *        with out_width  = \lfloor in_width  / delta_out \rfloor
 *         and out_height = \lfloor in_height / delta_out \rfloor
 * 
 * 
 */
void oversample_bilin(double* in , int in_width,  int in_height,
                      double* out, int out_width, int out_height,
                      double delta_out){
    
    assert(delta_out<=1);

    
    int im, jm, ip, jp;    
    double x, y;
    
    for(int i=0;i<out_height;i++){
        for(int j=0;j<out_width;j++){
            
            x= i*delta_out;
            y= j*delta_out;
            im = (int)x;
            jm = (int)y;
            ip = (int)x+1;
            jp = (int)y+1;

            //image extension by symmetrization
            if(ip>in_height-1){ip=2*in_height-1-ip;}
            if(im>in_height-1){im=2*in_height-1-im;}
            if(jp>in_width-1){jp=2*in_width-1-jp;}
            if(jm>in_width-1){jm=2*in_width-1-jm;}

                               
            out[i*out_width+j] = (x-(int)x)  *(y-(int)y)  * in[ip*in_width+jp]
                               + (x-(int)x)  *((int)y+1-y)* in[ip*in_width+jm]
                               + ((int)x+1-x)*(y-(int)y)  * in[im*in_width+jp]
                               + ((int)x+1-x)*((int)y+1-y)* in[im*in_width+jm];                               
                              
        }
    }
}































/** @brief Apply a convolution with a separable kernel 
 *
 * @param in             Input image of width X height samples
 * @param out            Output image (same dimension)
 *
 * @param xker           Kernel applied along direction x
 *                          radius = r_xker
 *                          width = 2*r_xker+1
 *
 * @param yker           Kernel applied along direction y
 *                          radius = r_yker
 *                          width = 2*r_yker+1
 *
 *  Compute:
 *
 *   out[i,j] = \sum_{-r_xker \leq k \leq +r_xker} 
 *                    xker[k]
 *                    . \sum_{-r_xker \leq l \leq +r_xker}
 *                              yker[l].in[i-k,j-l]
 *
 *  Border condition:  symmetrization at border
 *                     (at x = -1./2 and x = width-1./2)
 *
 */
void convolve(double* in, double* out, int width, int height,
              double* xker, int r_xker,
              double* yker, int r_yker){
    int i,j;
    int i_p,j_p;
    int k;
    double tmp;
    
    double* im_tmp = (double*)malloc(width*height*sizeof(double));
                
    /* convolution along x coordinates */
    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            tmp=0;
            for(k=-r_xker;k<=r_xker;k++){
                i_p = (i-k+2*height)%(2*height);    // TODO to manage module for negative value
                if(i_p>height-1){i_p=2*height-1-i_p;}
              
                tmp += xker[r_xker+k]*in[i_p*width+j];
            }
            im_tmp[i*width+j]=tmp;
        }
    }
    
    /* convolution along y coordinates */
    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            tmp=0;
            for(k=-r_yker;k<=r_yker;k++){
                j_p = (j-k+2*width)%(2*width);    // TODO to manage module for negative value
                if(j_p>width-1){j_p=2*width-1-j_p;}
                tmp += yker[r_yker+k]*im_tmp[i*width+j_p];
            }
            out[i*width+j]=tmp;
        }
    }
    
    free(im_tmp);
}
