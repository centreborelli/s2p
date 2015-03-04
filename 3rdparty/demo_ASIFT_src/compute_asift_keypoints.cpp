// Copyright (c) 2008-2011, Guoshen Yu <yu@cmap.polytechnique.fr>
// Copyright (c) 2008-2011, Jean-Michel Morel <morel@cmla.ens-cachan.fr>
//
// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// Jean-Michel Morel and Guoshen Yu, Method and device for the invariant 
// affine recognition recognition of shapes (WO/2009/150361), patent pending. 
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.
//
// 
//*------------------------ compute_asift_keypoints -------------------------*/
// Compute the ASIFT keypoints on the input image. 
// 
// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
// 
// Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image 
//            Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009. 
// Reference: ASIFT online demo (You can try ASIFT with your own images online.) 
//			  http://www.ipol.im/pub/algo/my_affine_sift/
/*---------------------------------------------------------------------------*/


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "compute_asift_keypoints.h"

#ifdef _OPENMP
#include <omp.h>
#endif


#define ABS(x)    (((x) > 0) ? (x) : (-(x)))


/* InitSigma gives the amount of smoothing applied to the image at the
first level of each octave.  In effect, this determines the sampling
needed in the image domain relative to amount of smoothing.  Good
values determined experimentally are in the range 1.2 to 1.8.
*/
/* float InitSigma_aa = 1.0;*/
static float InitSigma_aa = 1.6;

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/* Gaussian convolution kernels are truncated at this many sigmas from
the center.  While it is more efficient to keep this value small,
experiments show that for consistent scale-space analysis it needs
a value of about 3.0, at which point the Gaussian has fallen to
only 1% of its central value.  A value of 2.0 greatly reduces
keypoint consistency, and a value of 4.0 is better than 3.0.
*/
const float GaussTruncate1 = 4.0;


/* --------------------------- Blur image --------------------------- */


/* Same as ConvBuffer, but implemented with loop unrolling for increased
speed.  This is the most time intensive routine in keypoint detection,
so deserves careful attention to efficiency.  Loop unrolling simply
sums 5 multiplications at a time to allow the compiler to schedule
operations better and avoid loop overhead.  This almost triples
speed of previous version on a Pentium with gcc.
*/
void ConvBufferFast(float *buffer, float *kernel, int rsize, int ksize)
{
  int i;
  float *bp, *kp, *endkp;
  float sum;

  for (i = 0; i < rsize; i++) {
    sum = 0.0;
    bp = &buffer[i];
    kp = &kernel[0];
    endkp = &kernel[ksize];

    /* Loop unrolling: do 5 multiplications at a time. */
    //      while (kp + 4 < endkp) {
    //      sum += (double) bp[0] * (double) kp[0] + (double)  bp[1] * (double) kp[1] + (double) bp[2] * (double) kp[2] +
    //             (double) bp[3] * (double) kp[3] + (double)  bp[4] * (double) kp[4];
    //      bp += 5;
    //      kp += 5;		  
    //     }
    //      /* Do 2 multiplications at a time on remaining items. */
    //     while (kp + 1 < endkp) {
    //       sum += (double) bp[0] * (double) kp[0] + (double)  bp[1] * (double) kp[1];
    //	   bp += 2;
    //	   kp += 2;
    //	 }
    //      /* Finish last one if needed. */
    //		if (kp < endkp) {
    //		sum += (double) *bp * (double) *kp;
    //		}

    while (kp < endkp) {
      sum += *bp++ * *kp++;			  
    }

    buffer[i] = sum;
  }
}

/* Convolve image with the 1-D kernel vector along image rows.  This
is designed to be as efficient as possible.  Pixels outside the
image are set to the value of the closest image pixel.
*/
void ConvHorizontal(vector<float>& image, int width, int height, float *kernel, int ksize)
{
  int rows, cols, r, c, i, halfsize;
  float buffer[4000];
  vector<float> pixels(width*height);


  rows = height;
  cols = width;

  halfsize = ksize / 2;
  pixels = image;
  assert(cols + ksize < 4000);

  for (r = 0; r < rows; r++) {
    /* Copy the row into buffer with pixels at ends replicated for
    half the mask size.  This avoids need to check for ends
    within inner loop. */
    for (i = 0; i < halfsize; i++)
      buffer[i] = pixels[r*cols];
    for (i = 0; i < cols; i++)
      buffer[halfsize + i] = pixels[r*cols+i];
    for (i = 0; i < halfsize; i++)
      buffer[halfsize + cols + i] = pixels[r*cols+cols-1];

    ConvBufferFast(buffer, kernel, cols, ksize);
    for (c = 0; c < cols; c++)
      pixels[r*cols+c] = buffer[c];
  }
  image = pixels;	 
}


/* Same as ConvHorizontal, but apply to vertical columns of image.
*/
void ConvVertical(vector<float>& image, int width, int height, float *kernel, int ksize)
{
  int rows, cols, r, c, i, halfsize;
  float buffer[4000];
  vector<float> pixels(width*height);

  rows = height;
  cols = width;

  halfsize = ksize / 2;
  pixels = image;
  assert(rows + ksize < 4000);

  for (c = 0; c < cols; c++) {
    for (i = 0; i < halfsize; i++)
      buffer[i] = pixels[c];
    for (i = 0; i < rows; i++)
      buffer[halfsize + i] = pixels[i*cols+c];
    for (i = 0; i < halfsize; i++)
      buffer[halfsize + rows + i] = pixels[(rows - 1)*cols+c];

    ConvBufferFast(buffer, kernel, rows, ksize);
    for (r = 0; r < rows; r++)
      pixels[r*cols+c] = buffer[r];
  }

  image = pixels;	 
}



/* 1D Convolve image with a Gaussian of width sigma and store result back
in image.   This routine creates the Gaussian kernel, and then applies
it in horizontal (flag_dir=0) OR vertical directions (flag_dir!=0).
*/
void GaussianBlur1D(vector<float>& image, int width, int height, float sigma, int flag_dir)
{
  float x, kernel[100], sum = 0.0;
  int ksize, i;

  /* The Gaussian kernel is truncated at GaussTruncate sigmas from
  center.  The kernel size should be odd.
  */
  ksize = (int)(2.0 * GaussTruncate1 * sigma + 1.0);
  ksize = MAX(3, ksize);    /* Kernel must be at least 3. */
  if (ksize % 2 == 0)       /* Make kernel size odd. */
    ksize++;
  assert(ksize < 100);

  /* Fill in kernel values. */
  for (i = 0; i <= ksize; i++) {
    x = i - ksize / 2;
    kernel[i] = exp(- x * x / (2.0 * sigma * sigma));
    sum += kernel[i];
  }
  /* Normalize kernel values to sum to 1.0. */
  for (i = 0; i < ksize; i++)
    kernel[i] /= sum;

  if (flag_dir == 0)
  {
    ConvHorizontal(image, width, height, kernel, ksize);
  }
  else
  {
    ConvVertical(image, width, height, kernel, ksize);
  }
}


void compensate_affine_coor1(float *x0, float *y0, int w1, int h1, float t1, float t2, float Rtheta)
{
  float x_ori, y_ori;	
  float x_tmp, y_tmp;

  float x1 = *x0;
  float y1 = *y0;


  Rtheta = Rtheta*PI/180;

  if ( Rtheta <= PI/2 )
  {
    x_ori = 0;
    y_ori = w1 * sin(Rtheta) / t1;
  }
  else
  {
    x_ori = -w1 * cos(Rtheta) / t2;
    y_ori = ( w1 * sin(Rtheta) + h1 * sin(Rtheta-PI/2) ) / t1;
  }

  float sin_Rtheta = sin(Rtheta);
  float cos_Rtheta = cos(Rtheta);


  /* project the coordinates of im1 to original image before tilt-rotation transform */
  /* Get the coordinates with respect to the 'origin' of the original image before transform */
  x1 = x1 - x_ori;
  y1 = y1 - y_ori;
  /* Invert tilt */
  x1 = x1 * t2;
  y1 = y1 * t1;
  /* Invert rotation (Note that the y direction (vertical) is inverse to the usual concention. Hence Rtheta instead of -Rtheta to inverse the rotation.) */
  x_tmp = cos_Rtheta*x1 - sin_Rtheta*y1;
  y_tmp = sin_Rtheta*x1 + cos_Rtheta*y1;
  x1 = x_tmp;
  y1 = y_tmp;		

  *x0 = x1;
  *y0 = y1;
}


/* -------------- MAIN FUNCTION ---------------------- */

int compute_asift_keypoints(vector<float>& image, int width, int height, int num_of_tilts, int verb, vector< vector< keypointslist > >& keys_all, siftPar &siftparameters)
// Compute ASIFT keypoints in the input image.
// Input:
// image: input image
// width, height: width and height of the input image.
// num_of_tilts: number of tilts to simulate.
// verb: 1/0 --> show/don not show verbose messages. (1 for debugging) 
// keys_all (output): ASIFT keypoints. It is a 2D matrix with varying rows and columns. Each entry keys_all[tt][rr] 
//	stores the SIFT keypoints calculated on the image with the simulated tilt index tt and simulated rotation index rr (see the code below). In the coordinates of the keypoints,
//	the affine distortions have been compensated. 
// siftparameters: SIFT parameters.
//
// Output: the number of keypoints 
{	
  vector<float> image_t, image_tmp1, image_tmp;	

  float t_min, t_k;
  int num_tilt, tt, num_rot_t2, rr;
  int fproj_o;
  float fproj_p, fproj_bg;
  char fproj_i;
  float *fproj_x4, *fproj_y4;
  //  float frot_b=0;
  float frot_b=128;
  char *frot_k;
  int  counter_sim=0, num_sim;
  int flag_dir = 1;
  float BorderFact=6*sqrt(2.);

  int num_keys_total=0;


  fproj_o = 3;
  fproj_p = 0;
  fproj_i = 0;
  fproj_bg = 0;
  fproj_x4 = 0;
  fproj_y4 = 0;

  frot_k = 0;

  num_rot_t2 = 10;

  t_min = 1;
  t_k = sqrt(2.);


  num_tilt = num_of_tilts;


  if ( num_tilt < 1)
  {
    printf("Number of tilts num_tilt should be equal or larger than 1. \n");
    exit(-1);	
  }

  image_tmp1 = image;


  /* Calculate the number of simulations, and initialize keys_all */
  keys_all = std::vector< vector< keypointslist > >(num_tilt);	
  for (tt = 1; tt <= num_tilt; tt++)
  {
    float t = t_min * pow(t_k, tt-1);

    if ( t == 1 )
    {
      counter_sim ++;

      keys_all[tt-1] = std::vector< keypointslist >(1);
    }
    else
    {
      int num_rot1 = round(num_rot_t2*t/2);        
      if ( num_rot1%2 == 1 )
      {
        num_rot1 = num_rot1 + 1;
      }
      num_rot1 = num_rot1 / 2;
      counter_sim +=  num_rot1;

      keys_all[tt-1] = std::vector< keypointslist >(num_rot1);	
    }         		
  }

  num_sim = counter_sim;

  if ( verb )
  {
    printf("%d affine simulations will be performed. \n", num_sim);
  }

  counter_sim = 0;



  /* Affine simulation (rotation+tilt simulation) */
  // Loop on tilts. 
#ifdef _OPENMP
  omp_set_nested(1);
#endif
#pragma omp parallel for private(tt)
  for (tt = 1; tt <= num_tilt; tt++)
  {
    float t = t_min * pow(t_k, tt-1);

    float t1 = 1;
    float t2 = 1/t;

    // If tilt t = 1, do not simulate rotation. 	
    if ( t == 1 )
    {					
      // copy the image from vector to array as compute_sift_keypoints uses only array.				
      float *image_tmp1_float = new float[width*height];
      for (int cc = 0; cc < width*height; cc++)
        image_tmp1_float[cc] = image_tmp1[cc];

      compute_sift_keypoints(image_tmp1_float,keys_all[tt-1][0],width,height,siftparameters);

      delete[] image_tmp1_float;

    }
    else
    {
      // The number of rotations to simulate under the current tilt.   
      int num_rot1 = round(num_rot_t2*t/2);        

      if ( num_rot1%2 == 1 )
      {
        num_rot1 = num_rot1 + 1;
      }
      num_rot1 = num_rot1 / 2;
      float delta_theta = PI/num_rot1;		

      // Loop on rotations.    
#pragma omp parallel for private(rr)
      for ( int rr = 1; rr <= num_rot1; rr++ ) 
      {
        float theta = delta_theta * (rr-1);
        theta = theta * 180 / PI;

        vector<float> image_t;
        int width_r, height_r;

        // simulate a rotation: rotate the image with an angle theta. (the outside of the rotated image are padded with the value frot_b)
        frot(image, image_t, width, height, &width_r, &height_r, &theta, &frot_b , frot_k);

        /* Tilt */			 
        int width_t = (int) (width_r * t1);
        int height_t = (int) (height_r * t2);  

        int fproj_sx = width_t;
        int fproj_sy = height_t;     

        float fproj_x1 = 0;
        float fproj_y1 = 0;
        float fproj_x2 = width_t;
        float fproj_y2 = 0;
        float fproj_x3 = 0;	     
        float fproj_y3 = height_t;

        /* Anti-aliasing filtering along vertical direction */
        /* sigma_aa = InitSigma_aa * log2(t);*/
        float sigma_aa = InitSigma_aa * t / 2;
        GaussianBlur1D(image_t,width_r,height_r,sigma_aa,flag_dir);


        // simulate a tilt: subsample the image along the vertical axis by a factor of t.
        vector<float> image_tmp(width_t*height_t);			 
        fproj (image_t, image_tmp, width_r, height_r, &fproj_sx, &fproj_sy, &fproj_bg, &fproj_o, &fproj_p, &fproj_i , fproj_x1 , fproj_y1 , fproj_x2 , fproj_y2 , fproj_x3 , fproj_y3, fproj_x4, fproj_y4); 

        vector<float> image_tmp1 = image_tmp;	 

        if ( verb )
        {
          printf("Rotation theta = %.2f, Tilt t = %.2f. w=%d, h=%d, sigma_aa=%.2f, \n", theta, t, width_t, height_t, sigma_aa);
        }


        float *image_tmp1_float = new float[width_t*height_t];
        for (int cc = 0; cc < width_t*height_t; cc++)
          image_tmp1_float[cc] = image_tmp1[cc];	 

        // compute SIFT keypoints on simulated image. 	 
        keypointslist keypoints;
        keypointslist keypoints_filtered;
        compute_sift_keypoints(image_tmp1_float,keypoints,width_t,height_t,siftparameters);

        delete[] image_tmp1_float;		

        /* check if the keypoint is located on the boundary of the parallelogram (i.e., the boundary of the distorted input image). If so, remove it to avoid boundary artifacts. */
        if ( keypoints.size() != 0 )
        {
          for ( int cc = 0; cc < (int) keypoints.size(); cc++ )
          {		      

            float x0, y0, x1, y1, x2, y2, x3, y3 ,x4, y4, d1, d2, d3, d4, scale1, theta1, sin_theta1, cos_theta1, BorderTh;

            x0 = keypoints[cc].x;
            y0 = keypoints[cc].y;
            scale1= keypoints[cc].scale;

            theta1 = theta * PI / 180;
            sin_theta1 = sin(theta1);
            cos_theta1 = cos(theta1);

            /* the coordinates of the 4 submits of the parallelogram */
            if ( theta <= 90 )
            {
              x1 = height * sin_theta1;
              y1 = 0;			 
              y2 = width * sin_theta1;
              x3 = width * cos_theta1;
              x4 = 0;
              y4 = height * cos_theta1;
              x2 = x1 + x3;
              y3 = y2 + y4;

              /* note that the vertical direction goes from top to bottom!!! 
              The calculation above assumes that the vertical direction goes from the bottom to top. Thus the vertical coordinates need to be reversed!!! */
              y1 = y3 - y1;
              y2 = y3 - y2;
              y4 = y3 - y4;
              y3 = 0;

              y1 = y1 * t2;
              y2 = y2 * t2;
              y3 = y3 * t2;
              y4 = y4 * t2;
            }
            else
            {
              y1 = -height * cos_theta1;
              x2 = height * sin_theta1;
              x3 = 0;
              y3 = width * sin_theta1;				 
              x4 = -width * cos_theta1;
              y4 = 0;
              x1 = x2 + x4;
              y2 = y1 + y3;

              /* note that the vertical direction goes from top to bottom!!! 
              The calculation above assumes that the vertical direction goes from the bottom to top. Thus the vertical coordinates need to be reversed!!! */
              y1 = y2 - y1;
              y3 = y2 - y3;
              y4 = y2 - y4;
              y2 = 0;

              y1 = y1 * t2;
              y2 = y2 * t2;
              y3 = y3 * t2;
              y4 = y4 * t2;
            }		       		    

            /* the distances from the keypoint to the 4 sides of the parallelogram */
            d1 = ABS((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1)) / sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
            d2 = ABS((x3-x2)*(y2-y0)-(x2-x0)*(y3-y2)) / sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
            d3 = ABS((x4-x3)*(y3-y0)-(x3-x0)*(y4-y3)) / sqrt((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3));
            d4 = ABS((x1-x4)*(y4-y0)-(x4-x0)*(y1-y4)) / sqrt((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4));

            BorderTh = BorderFact*scale1;

            if (!((d1<BorderTh) || (d2<BorderTh) || (d3<BorderTh) || (d4<BorderTh) ))
            {				 					   
              // Normalize the coordinates of the matched points by compensate the simulate affine transformations
              compensate_affine_coor1(&x0, &y0, width, height, 1/t2, t1, theta);
              keypoints[cc].x = x0;
              keypoints[cc].y = y0;

              keypoints_filtered.push_back(keypoints[cc]);	 
            }				   
          }
        }			 
        keys_all[tt-1][rr-1] = keypoints_filtered;
      }		 
    }         
  }

  {		
    for (tt = 0; tt < (int) keys_all.size(); tt++)
      for (rr = 0; rr < (int) keys_all[tt].size(); rr++)
      {
        num_keys_total += (int) keys_all[tt][rr].size();
      }				
      printf("%d ASIFT keypoints are detected. \n", num_keys_total);
  }

  return num_keys_total;
}
