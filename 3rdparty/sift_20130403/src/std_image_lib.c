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
 * @file std_image_lib.c
 * @brief  Save images
 *
 * @li gray to false color routines
 * @li contrast changes
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */







#include "std_image_lib.h"










void linear_conversion(double *imIn, double *imOut, int length){
    double a, b;
    double minVal = +10000000000;
    double maxVal = -10000000000;
    for (int i = 0; i < length; i++) {
        if (imIn[i] > maxVal)
        maxVal = imIn[i];
        if (imIn[i] < minVal)
            minVal = imIn[i];
    }
    a = (250.0 - 0.0) / (maxVal - minVal);
    b = -a * minVal;
    for (int i = 0; i < length; i++)
        imOut[i] = a * imIn[i] + b;
}

void printImage(double *imageIn, int width, int height, char *name){
    float *imf32 = (float *) malloc(width * height * sizeof(float));
    for (int i = 0; i < width * height; i++) {
        imf32[i] = (float) imageIn[i];
    }
    io_png_write_f32(name, imf32, width, height, 1);
    free(imf32);
}



void printImage_LinearConversion(double *imageIn, int width, int height, char *name){
    double *imTemp = (double*)malloc(width*height*sizeof(double));
    linear_conversion(imageIn, imTemp, width * height);
    printImage(imTemp, width, height, name);
    free(imTemp);
}


void printColorImage(double *red, double *green, double *blue, int width, int height, char *name){
    float *imf32 = (float*)malloc(3*width*height*sizeof(float));
    for (int i = 0; i < width * height; i++) {
        imf32[i] = red[i];
        imf32[i + width * height] = green[i];
        imf32[i + 2 * width * height] = blue[i];
    }
    io_png_write_f32(name, imf32, width, height, 3);
}

void printrgb(double *rgb, int width, int height, char *name){
    float *imf32 = (float*)malloc(3*width*height*sizeof(float));    
    for(int i=0;i<3*width*height;i++){
        imf32[i] = rgb[i];
    }
    io_png_write_f32(name, imf32, width, height, 3);
}






/** @brief Gray to false color (using the HSV colormap)
 *
 *   for DoG illustration in new code
 *
 */
void gray2hsv(double *gray, double *rgb, int width, int height, double vmin, double vmax){

    // The color wheel used as a colormap here is a restriction of the hsv 
    // color space to the circle defined by
    double saturation = 1.0;
    double value = 1.0;
    // and hue \in [0,360].

    for (int pix = 0; pix < 3 * width * height; pix++) {
        rgb[pix] = 0;
    }


    double monmax = -1000000.0;
    double monmin = +1000000.0;
    for (int pix = 0; pix < width * height; pix++) {
        if (gray[pix] > monmax) {
            monmax = gray[pix];
        }
        if (gray[pix] < monmin) {
            monmin = gray[pix];
        }
    }
    
    double max, min;
    max = monmax;
    min = monmin;

    for (int pix = 0; pix < width * height; pix++) {
        //double hue = (gray[pix] - min) / (max - min) * 360.0;
        double hue = (gray[pix] - min) / (max - min) * 359.0;
        
        int t = (int) (hue / 60.0);
        double red, green, blue;
        /*
         * The color map is defined by a set of 3 piecewise linear
         * functions 
         */
        switch (t) {
        case 0:
            red = value;
            green =  value * (1.0 - (1.0 - (hue / 60.0 - (double) t)) * saturation);
            blue = value * (1.0 - saturation);
            break;
        case 1:
            red = value * (1.0 - (hue / 60.0 - (double) t) * saturation);
            green = value;
            blue = value * (1.0 - saturation);
            break;
        case 2:
            red = value * (1.0 - saturation);
            green = value;
            blue =
                value * (1.0 - (1.0 - (hue / 60.0 - (double) t)) * saturation);
            break;
        case 3:
            red = value * (1.0 - saturation);
            green = value * (1.0 - (hue / 60.0 - (double) t) * saturation);
            blue = value;
            break;
        case 4:
            red =  value * (1.0 - (1.0 - (hue / 60.0 - (double) t)) * saturation);
            green = value * (1.0 - saturation);
            blue = value;
            break;
        case 5:
            red = value;
            green = value * (1.0 - saturation);
            blue = value * (1.0 - (hue / 60.0 - (double) t) * saturation);
            break;
        default:
            red =0;
            green = 0;
            blue = 0;
            break;
        }
        rgb[pix] = 250 * red;
        rgb[width * height + pix] = 250 * green;
        rgb[2 * width * height + pix] = 250 * blue;
    }

}
