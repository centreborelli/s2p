/*
IPOL SIFT
Copyright (C) 2013, Ives Rey Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20130318 (March 18, 2013)

== Patent Warning and Licence =================================================

The SIFT method is pattented 

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
 * @file sift_description.c
 * @brief Computation the SIFT feature vector
 *
 * @li Attribution of a principal orientation
 * @li Computation of the SIFT feature vector
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */





#ifndef _SIFT_DESCRIPTION_H_
#define _SIFT_DESCRIPTION_H_

#include <math.h>
#include <assert.h>
#include <stdbool.h>

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

//#ifndef PI
#define PI 3.141592
//#endif



// TEMP
// #include <stdio.h>
// #include <stdlib.h>

int accumulate_orientation_histogram(double x_key, double y_key, double sigma_key,
                                     double* gradX, double* gradY, int width, int height,
                                     int nbins, double sigma_window, double K,
                                     double* hist);

/** with the same parameterization than in the pdf */
int accumulate_orientation_histogram_bis(double x_key, double y_key, double sigma_key,
                                         double* gradX, double* gradY, int width, int height,
                                         int nbins, double lambda_ori,
                                         double* hist);


int extract_principal_orientations(double* hist, int nbins,
                                   double threshold, bool flagsmooth,
                                   double* tmp_oris);


int extract_sift_feature_vector(double x_key, double y_key, double sigma_key, double theta_key,
                                double* gradX, double* gradY, int width, int height,
                                int Nhist, int Nbins, double sigma_window, double K,
                                double* descr);

/** with the same parameterization than in the pdf */
int extract_sift_feature_vector_bis(double x_key, double y_key, double sigma_key, double theta_key,
                                double* gradX, double* gradY, int width, int height,
                                int Nhist, int Nbins, double lambda_descr,
                                double* descr);


void threshold_and_quantize_sift_descriptor(double* descr, int Nhist, int Nbins,
                                            bool flagthreshold, double threshold);
                                           // int* qtdescr);


/** Alternative computation of the feature vector :
 *  - bilinear spatial distribution of each sample contribution
 */
int extract_sift_feature_vector_2(double x_key, double y_key, double sigma_key, double theta_key,
                                  double* gradX, double* gradY, int width, int height,
                                  int Nhist, int Nbins, double sigma_window, double K,
                                  double* descr);

 
#endif
