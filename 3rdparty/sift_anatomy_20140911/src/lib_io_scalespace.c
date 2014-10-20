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




This file also implements the colomap published in the paper
   [2] "Diverging Color Maps for Scientific Visualization."
        Kenneth Moreland
        International Symposium on Visual Computing 2009




*/
/**
 * @file sift_scalespace.c
 * @brief data structures to store the scale-space
 *
 * @li struct keypoint      keypoint data structure.
 * @li struct sift_keypoints  list of keypoints with variable length.
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "lib_discrete.h"
#include "lib_io_scalespace.h"
#include "lib_util.h"
#include "io_png.h"






/**
 *
 *  Assumes the values of imageIn are ranging from 0 to 1
 *
 *
 * */
void printImage(const float *imageIn, int w, int h, const char *name)
{

    float *imf32 = (float *) malloc(w * h * sizeof(float));
    for (int i = 0; i < w * h; i++) {
        imf32[i] = 255*(float) imageIn[i];
    }
    io_png_write_f32(name, imf32, w, h, 1);
    xfree(imf32);
}


void printImage_LinearConversion(const float *imageIn, int w, int h, const char *name)
{
    float *imTemp = xmalloc(w*h*sizeof(float));
    linear_conversion(imageIn, imTemp, w * h); // imTemp values are in [0,1]
    printImage(imTemp, w, h, name);
    xfree(imTemp);
}


void printrgb(const float *rgb, int w, int h, const char *name)
{
    io_png_write_f32(name, rgb, w, h, 3);
}




/** @brief Gray to false color (using the HSV colormap)
 *
 *   for DoG illustration in new code
 *
 */
void gray2hsv(const float *gray, float *rgb, int w, int h, float vmin, float vmax)
{
    // The color wheel used as a colormap here is a restriction of the hsv
    // color space to the circle defined by
    float saturation = 1.0;
    float value = 1.0;
    // and hue \in [0,360].

    for (int i = 0; i < 3 * w * h; i++)  rgb[i] = 0;

    float max = array_max(gray, w*h);
    float min = array_min(gray, w*h);

    for (int i = 0; i < w * h; i++) {
        float hue = (gray[i] - min) / (max - min) * 359.0;

        int t = (int) (hue / 60.0);
        float red, green, blue;
        /*
         * The color map is defined by a set of 3 piecewise linear
         * functions
         */
        switch (t) {
        case 0:
            red = value;
            green =  value * (1.0 - (1.0 - (hue / 60.0 - (float) t)) * saturation);
            blue = value * (1.0 - saturation);
            break;
        case 1:
            red = value * (1.0 - (hue / 60.0 - (float) t) * saturation);
            green = value;
            blue = value * (1.0 - saturation);
            break;
        case 2:
            red = value * (1.0 - saturation);
            green = value;
            blue =
                value * (1.0 - (1.0 - (hue / 60.0 - (float) t)) * saturation);
            break;
        case 3:
            red = value * (1.0 - saturation);
            green = value * (1.0 - (hue / 60.0 - (float) t) * saturation);
            blue = value;
            break;
        case 4:
            red =  value * (1.0 - (1.0 - (hue / 60.0 - (float) t)) * saturation);
            green = value * (1.0 - saturation);
            blue = value;
            break;
        case 5:
            red = value;
            green = value * (1.0 - saturation);
            blue = value * (1.0 - (hue / 60.0 - (float) t) * saturation);
            break;
        default:
            red =0;
            green = 0;
            blue = 0;
            break;
        }
        rgb[i] = 250 * red;
        rgb[w * h + i] = 250 * green;
        rgb[2 * w * h + i] = 250 * blue;
    }

}







void print_sift_scalespace_gray(const struct sift_scalespace* scalespace, const char* basename)
{
    char name[FILENAME_MAX];
    int nOct = scalespace->nOct;
    for(int o = 0; o < nOct; o++){
        const struct octa* octave = scalespace->octaves[o];
        int nSca = octave->nSca;
        int w = octave->w;
        int h = octave->h;
        for(int s = 0; s < nSca; s++){
            const float* image = &octave->imStack[s*w*h];
            sprintf(name,"%s_o%03i_s%03i.png",basename,o,s);
            printImage(image, w,h, name);
        }
    }
}


void nearestneighbor_interp(const float* in,int win, int hin, float* out, int wout,int hout)
{
    if ((hin != hout) || (win != wout)){
        float fh = (float)hin/(float)hout;
        float fw = (float)win/(float)wout;
        for(int i = 0; i < hout; i++){
            int k = floor( fh * (float)i );
            assert(k>-1); assert(k < hin);
            for(int j = 0; j < wout; j++){
                int l = floor( fw * (float)j );
                assert(l > -1); assert(l < win);
                out[i*wout+j] = in[k*win+l];
            }
        }
    }else{
        for(int i = 0; i < wout*hout; i++){
            out[i] = in[i];
        }
    }
}






void Msh2Lab(float M, float s, float h, float* L, float* a, float* b)
{
    *L = M*cos(s);
    *a = M*sin(s)*cos(h);
    *b = M*sin(s)*sin(h);
}


void Lab2xyz(float L, float a, float b, float* x, float* y, float* z)
{
    float vY = (L + 16.0)/116.0;
    float vX = a/500.0 + vY;
    float vZ = vY - b/200.0;

    if (vY*vY*vY > 0.008856){
        vY = vY*vY*vY;
    }else{
        vY = (vY - 16.0/116.0) / 7.787;
    }
    if (vX*vX*vX > 0.008856){
        vX = vX*vX*vX;
    }
    else{
        vX = (vX - 16.0/116.0) / 7.787;
    }
    if (vZ*vZ*vZ > 0.008856){
        vZ = vZ*vZ*vZ;
    }else{
        vZ = (vZ - 16.0/116.0) / 7.787;
    }

    *x =  95.047*vX;  //ref_X =  95.047  Observer= 2 degrees, Illuminant= D65
    *y = 100.000*vY;  //ref_Y = 100.000
    *z = 108.883*vZ;  //ref_Z = 108.883
}



void xyz2rgb(float x, float y, float z, float *r, float* g, float* b)
{
    float var_x = x/100;  //x from 0 to  95.047  (observer = 2°, illuminant = d65)
    float var_y = y/100;  //y from 0 to 100.000
    float var_z = z/100;  //z from 0 to 108.883

    float var_r = var_x *  3.2406 + var_y * -1.5372 + var_z * -0.4986;
    float var_g = var_x * -0.9689 + var_y *  1.8758 + var_z *  0.0415;
    float var_b = var_x *  0.0557 + var_y * -0.2040 + var_z *  1.0570;

    if ( var_r > 0.0031308 ){
        var_r = 1.055 * ( pow(var_r,  ( 1.0 / 2.4 )) ) - 0.055;
    }else{
        var_r = 12.92 * var_r;
    }
    if ( var_g > 0.0031308 ){
        var_g = 1.055 * ( pow(var_g,  ( 1.0 / 2.4 )) ) - 0.055;
    }else{
        var_g = 12.92 * var_g;
    }
    if ( var_b > 0.0031308 ){
        var_b = 1.055 * ( pow(var_b,  ( 1.0 / 2.4 )) ) - 0.055;
    }else{
        var_b = 12.92 * var_b;
    }
    *r = var_r*255;
    *g = var_g*255;
    *b = var_b*255;
}


/** @brief gray 2 Msh
 *
 *      @param in : array of w*h float  
 *      @param out : array of w*h color points in Msh color space (polar Lab)
 *                     ... of 3*w*h float 
 * 
 *
 */
void gray2Msh2rgb(const float* in, float* out, int w, int h)
{
    float M, s, hue,    L, a, b,   x, y, z ;
    float max = array_max(in, w*h);
    float min = array_min(in, w*h);
    float mid = (max + min)/2.0;

    for(int i = 0; i < w*h; i++){
        if ( in[i] < mid ){
            float a = (in[i] - min) / (mid - min);
            M = 80.0 + (88.0 - 80.0)*a;
            s = 1.08 - 1.08*a;
            hue = 0.50 + (1.061 - 0.5)*a; 
        }else{
            float a = (in[i] - mid) / (max - mid);
            M = 88.0 + (80.0 - 88.0)*a;
            s = 1.08*a;
            hue = 1.061 + (-1.1 - 1.061)*a;
        }
        Msh2Lab(M, s, hue, &L, &a, &b);
        Lab2xyz(L, a, b, &x, &y, &z);
        xyz2rgb(x, y, z, &out[i], &out[w*h+i], &out[2*w*h+i]);
    }
}


void print_sift_scalespace_rgb(const struct sift_scalespace* scalespace, const char* basename)
{
    char name[FILENAME_MAX];
    int wALL = scalespace->octaves[0]->w;
    int hALL = scalespace->octaves[0]->h;
    float* imtemp = xmalloc(wALL*hALL*sizeof(float));
    float* imrgb = xmalloc(3*wALL*hALL*sizeof(float));

    int nOct = scalespace->nOct;
    for(int o = 0; o < nOct; o++){
        const struct octa* octave = scalespace->octaves[o];
        int nSca = octave->nSca;
        int w = octave->w;
        int h = octave->h;
        for(int s = 1; s < nSca-1; s++){
            const float* image = &octave->imStack[s*w*h];
            nearestneighbor_interp(image,w,h,imtemp,wALL,hALL);
            gray2Msh2rgb(imtemp, imrgb, wALL, hALL);
            sprintf(name,"%s_o%03i_s%03i.png",basename,o,s);
            printrgb(imrgb, wALL,hALL, name);
        }
    }
    xfree(imrgb);
    xfree(imtemp);
}


void print_sift_scalespace_gray_nearestneighbor(const struct sift_scalespace* scalespace, const char* basename)
{
    char name[FILENAME_MAX];
    int wALL = scalespace->octaves[0]->w;
    int hALL = scalespace->octaves[0]->h;
    float* imtemp = xmalloc(wALL*hALL*sizeof(float));

    int nOct = scalespace->nOct;
    for(int o = 0; o < nOct; o++){
        const struct octa* octave = scalespace->octaves[o];
        int nSca = octave->nSca;
        int w = octave->w;
        int h = octave->h;
        for(int s=0;s<nSca;s++){
            const float* image = &octave->imStack[s*w*h];
            nearestneighbor_interp(image,w,h,
                      imtemp,wALL,hALL);
            sprintf(name,"%s_o%03i_s%03i.png",basename,o,s);
            printImage(imtemp, wALL,hALL, name);
        }
    }
    xfree(imtemp);
}


void print_sift_scalespace_rgb_nointerp(const struct sift_scalespace* scalespace, const char* basename)
{
    char name[FILENAME_MAX];
    int nOct = scalespace->nOct;
    for(int o = 0; o < nOct; o++){
        struct octa* octave = scalespace->octaves[o];
        int nSca   = octave->nSca;
        int w  = octave->w;
        int h = octave->h;

        for(int s = 1; s < nSca - 1; s++){
            const float* image = &octave->imStack[s*w*h];
            float* imrgb = xmalloc(3*w*h*sizeof(float));
            for(int i = 0; i < 3*w*h; i++){
                imrgb[i] = 0.0;
            }
            gray2Msh2rgb(image, imrgb, w, h);
            sprintf(name,"%s_o%03i_s%03i.png",basename,o,s);
            printrgb(imrgb, w,h, name);
            xfree(imrgb);
        }
    }
}
