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
 * @file sift_description.c
 * @brief Computation the SIFT feature vector
 *
 * @li Attribution of a principal orientation
 * @li Computation of the SIFT feature vector
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */





#include "sift_description.h"



/** @brief Accumulate gradient orientation histogram around a keypoint
 * 
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key      keypoint coordinates.
 * 
 * @param gradX
 * @param gradX          precomputed gradient relative to the nearest scale.
 * 
 * ----PARAM----
 * @param sigma_window   (= 1.5)   std dev of Gaussian window, the analysis of the gradient orientation is local
 * @param nbins          (= 36)  number of bins covering the range [0,2pi]
 * @param K              (= 3) The histogram cover a square patch of width 2*K*sigma_window (=3*1.5*sigma_key)
 * 
 * ----OUTPUT----
 * @output hist          the output is vec of nbins valus
 *
 * ----RETURN----
 * @return count         number of pixels contributing to the histogram.
 *
 */
int accumulate_orientation_histogram(double x_key, double y_key, double sigma_key,
                                     double* gradX, double* gradY, int width, int height,
                                     int nbins, double sigma_window, double K,
                                     double* hist){

    
    int count=0;   // number of pixel contributing to the histogram
    double M;      // gradient magnitude
    int    si,sj;  // sample coordinates in the image referential
    double sX,sY;  // sample coordinates in the keypoint referential
    double sOri;
    int    gamma;  // index of nearest bin
    
    /// Initialise output vector
    for(int i= 0;i<nbins;i++){hist[i] = 0.0;}
   
    /// Contributing pixels are inside a patch [siMin;siMax] X [sjMin;sjMax] of K*sigma_window*sigma_key (=4.5*sigma_key) 
    int siMin = MAX(0, (int)(x_key-K*sigma_window*sigma_key+0.5)); 
    int sjMin = MAX(0, (int)(y_key-K*sigma_window*sigma_key+0.5));  
    int siMax = MIN((int)(x_key+K*sigma_window*sigma_key+0.5), height-1);
    int sjMax = MIN((int)(y_key+K*sigma_window*sigma_key+0.5), width-1);
    
    
    
    

    /// For each pixel inside the patch.
    for(si=siMin;si<=siMax;si++){
        for(sj=sjMin;sj<=sjMax;sj++){

            /// Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            sX = ((double)si-x_key)/sigma_key;
            sY = ((double)sj-y_key)/sigma_key;
            
            /// Compute the gradient orientation (theta) on keypoint's invariant referential
            sOri = fmod(atan2(gradY[si*width+sj], gradX[si*width+sj]) + 10*2*PI, 2*PI);
                
            /// Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to gradients far from the keypoint.
            M = sqrt(gradY[si*width+sj]*gradY[si*width+sj]+gradX[si*width+sj]*gradX[si*width+sj])*exp(-(sX*sX+sY*sY)/(2*sigma_window*sigma_window));

            /// Determine the bin index in the circular histogram
            gamma = (int)(sOri/(2*PI)*(double)nbins+0.5)%nbins;
                
            /// Add the contribution to the orientation histogram
            hist[gamma] += M;
            
        }
    }
    count = (siMax-siMin)*(sjMax-sjMin);
    return count;
}







/** @brief Accumulate gradient orientation histogram around a keypoint
 * 
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key      keypoint coordinates.
 * 
 * @param gradX
 * @param gradX          precomputed gradient relative to the nearest scale.
 * 
 * ----PARAM----
 * @param lambda_ori     (= 1.5)
 *                       - The patch P^ori is ( 6 X lambda_ori X sigma_key )
 *                       - The Gaussian window has a standard deviation of lambda_ori X sigma_key
 * 
 * @param nbins          (= 36)  number of bins covering the range [0,2pi]
 *
 * 
 * ----OUTPUT----
 * @output hist          the output is vec of nbins valus
 *
 * ----RETURN----
 * @return count         number of pixels contributing to the histogram.
 *
 */
int accumulate_orientation_histogram_bis(double x_key, double y_key, double sigma_key,
                                         double* gradX, double* gradY, int width, int height,
                                         int nbins, double lambda_ori,
                                         double* hist){

    int count=0;   /* number of contributing pixels. */
    double M;      // gradient magnitude
    int    si,sj;  // sample coordinates in the image referential
    double sX,sY;  // sample coordinates in the keypoint referential
    double sOri;
    int    gamma;  // index of nearest bin
    
    /// Initialise output vector
    for(int i= 0;i<nbins;i++){hist[i] = 0.0;}
   
    /// Contributing pixels are inside a patch [siMin;siMax] X [sjMin;sjMax] of width  6*lambda_ori*sigma_key (=9*sigma_key) 
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
            M = sqrt(gradY[si*width+sj]*gradY[si*width+sj]+gradX[si*width+sj]*gradX[si*width+sj])*exp(-(sX*sX+sY*sY)/(2*lambda_ori*lambda_ori));

            /// Determine the bin index in the circular histogram
            gamma = (int)(sOri/(2*PI)*(double)nbins+0.5)%nbins;
            //printf("Gamma: %i , M: %f , sOri : %f \n ",gamma,M,sOri);
                
            /// Add the contribution to the orientation histogram
            hist[gamma] += M;
            
        }
    }
    count = (siMax-siMin)*(sjMax-sjMin);
    return count;
}



















void smooth_circular_histogram(int niter, double* hist, int nbins);
    


/** @brief Extract principal orientations from gradient orientation histogram
 *
 * ----INPUT----
 * @param hist           histogram of gradient orientation.
 * @param nbins          (=36) number of bins covering range [0,2pi].
 *
 * ----PARAM----
 * @param threshold      (0.8) local maxima are considered secondary
 *                       principal orientation if they exceed
 *                       threshold times the absolute maximum.
 *
 * @param flagsmooth     bool, if=1 then the histogram is smoothed
 *
 * ----OUTPUT----
 * @output tmp_oris      pointer to tab storing principal orientations.
 *
 * ----RETURN----
 * @return nori          number of principal orientation (>=1)
 */
int extract_principal_orientations(double* hist, int nbins,
                                   double threshold, bool flagsmooth,
                                   double* tmp_oris){

    int i,i_prev,i_next;   // bin indices
    int o;                 // index of principal orientation, and index tmp_oris.
    double tmp, offset;    // for normalization and interpolation

    if(flagsmooth){
        /// Smooth histogram : 3 iterated box filters
        smooth_circular_histogram(3,hist,nbins);
    }

    /// Linear normalization, so the absolute maximum equals 1
    tmp =-1.0; for(i=0;i<nbins;i++){if(hist[i]>tmp){tmp=hist[i];}}
    for(i=0;i<nbins;i++){hist[i]=hist[i]/tmp;}


    /// DEBUG DEBUG check linear normalisation
    tmp =-1.0; for(i=0;i<nbins;i++){if(hist[i]>tmp){tmp=hist[i];}}
  
    /// Search for local extrema in the histogram
    o=0; // at this stage, no principal orientation
    for(i=0;i<nbins;i++){
        i_prev=(i-1+nbins)%nbins;
        i_next=(i+1)%nbins;
        if((hist[i]>threshold)&(hist[i]>hist[i_prev])&(hist[i]>hist[i_next])){

            /// Quadratic interpolation of the position of each local maximum
            offset = (hist[i_prev]-hist[i_next])/(2*(hist[i_prev]+hist[i_next]+2*hist[i]));
            
            /// Add to vector of principal orientations (expressed in [0,2pi]
            tmp_oris[o] = ((double)i+0.5+offset)*2*PI/(double)nbins;
            o++; // new principal orientation
        }
    }
    
    // return the number of principal orientations
    return o;
}









/** @brief Extract keypoint feature vector
 * 
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key
 * @param theta_key      keypoint coordinates.
 * 
 * @param gradX
 * @param gradX          precomputed gradient relative to the nearest scale.
 * 
 * ----PARAM----
 * @param sigma_window   (=6) standard deviation of Gaussian window to make the representation local
 * @param Nhist          (=4) number of histograms in each of the two directions,
 * @param Nbins          (=8) number of bins covering the range [0,2pi]
 * @param K              (=2), Nhist histograms cover uniformly a patch of width K*sigma_window. (in the normalized referential)
 *                            i.e. the descriptor covers a (rotated) window of width = K*sigma_window*sigma (in the original image referential)
 * 
 * ----OUTPUT----
 * @output descr         the output is a vector of double here but it should be a vector of integers, descr being quantized. 
 * 
 * ----RETURN----
 * @return count         number of sample contributing to the descriptor
 */
int extract_sift_feature_vector(double x_key, double y_key, double sigma_key, double theta_key,
                                double* gradX, double* gradY, int width, int height,
                                int Nhist, int Nbins, double sigma_window, double K,
                                double* descr){

    
    int count=0;   // number of contributing samples
    double M;      // gradient magnitude
    int    si,sj;             // pixel coordinates in the image referential
    double sX,sY,sOri;        // pixel coordinates in the descriptor referential
    double alpha, beta,gamma; // used to compute histogram indices and weighting of trilinear distribution
    
    /// Contributing pixels are inside a patch [siMin;siMax]X[sjMin;sjMax] of K*(1+1/Nhist)*sigma_window*sigma_key  TODO 2>2.5 TODO TODO
    int siMin = MAX(0, (int)(x_key-sqrt(2.0)*K*sigma_window*sigma_key*(1+1./(2*(double)Nhist))+0.5));  int siMax = MIN((int)(x_key+sqrt(2.0)*K*sigma_window*sigma_key*(1+1./(2*(double)Nhist))+0.5), height-1);
    int sjMin = MAX(0, (int)(y_key-sqrt(2.0)*K*sigma_window*sigma_key*(1+1./(2*(double)Nhist))+0.5));  int sjMax = MIN((int)(y_key+sqrt(2.0)*K*sigma_window*sigma_key*(1+1./(2*(double)Nhist))+0.5), width-1);
        
    /// Initialise descr tab
    for(int i=0;i<Nhist*Nhist*Nbins;i++){descr[i] = 0.0;}
    
    /// For each pixel inside the patch.
    for(si=siMin;si<siMax;si++){
        for(sj=sjMin;sj<sjMax;sj++){
            
            /// Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            sX = (cos(-theta_key)*((double)si-x_key)-sin(-theta_key)*((double)sj-y_key))/sigma_key;
            sY = (sin(-theta_key)*((double)si-x_key)+cos(-theta_key)*((double)sj-y_key))/sigma_key;
            
            /// Does this sample fall inside the descriptor area ?
            if((fabs(sX)<0.5*(K*sigma_window*(((double)Nhist+1)/(double)Nhist)))&(fabs(sY)<0.5*(K*sigma_window*(((double)Nhist+1)/(double)Nhist)))){
                count++;

                /// Compute the gradient orientation (theta) on keypoint's invariant referential.
                sOri = fmod(atan2(gradY[si*width+sj], gradX[si*width+sj]) - theta_key + 10*2*PI,2*PI);                
                assert(sOri>0.0); assert(sOri<2*PI);
                
                /// Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to gradients far from the keypoint.
                M = sqrt(gradY[si*width+sj]*gradY[si*width+sj]+gradX[si*width+sj]*gradX[si*width+sj])
                    *exp(-(sX*sX+sY*sY)/(2*sigma_window*sigma_window));
                
                /// Determine the histogram and bin indices, Compute the (tri)linear weightings ...
                alpha = sX/(K*sigma_window/(double)Nhist) + ((double)Nhist-1.0)/2.0;
                beta  = sY/(K*sigma_window/(double)Nhist) + ((double)Nhist-1.0)/2.0;
                gamma = sOri/(2*PI)*(double)Nbins;
                
                ///    ...and add contributions to respective bins in different histograms.
                for(int i=MAX(0,(int)floor(alpha));i<=MIN((int)ceil(alpha),Nhist-1);i++){
                    for(int j=MAX(0,(int)floor(beta));j<=MIN((int)ceil(beta),Nhist-1);j++){ // looping through all surrounding histograms.
                        
                        
                        assert(i*Nhist*Nbins+j*Nbins+ ((int)ceil(gamma)+2*Nbins)%Nbins>=0);
                        assert(i*Nhist*Nbins+j*Nbins+ ((int)ceil(gamma)+2*Nbins)%Nbins<Nhist*Nhist*Nbins);
                        assert(i*Nhist*Nbins+j*Nbins+ (int)floor(gamma)>=0);
                        assert(i*Nhist*Nbins+j*Nbins+ (int)floor(gamma)<Nhist*Nhist*Nbins);

                        
                        // Contribution to left bin.                        
//                         printf("descripteur p   (i,j,k)=(%i,%i,%i)   \n",i,j,((int)ceil(gamma))%Nbins);
//                         printf("descripteur m   (i,j,k)=(%i,%i,%i)   \n",i,j,((int)floor(gamma)+Nbins)%Nbins);
                         
                        descr[i*Nhist*Nbins+j*Nbins+ (int)floor(gamma)] += (1.-(gamma-floor(gamma)))
                                                                          *(1.0-fabs((double)i-alpha))
                                                                          *(1.0-fabs((double)j-beta))
                                                                          *M;
                        // Contribution to right bin.
                        descr[i*Nhist*Nbins+j*Nbins+ ((int)ceil(gamma)+2*Nbins)%Nbins] += (1.0-(ceil(gamma)-gamma))
                                                                               *(1.0-fabs((double)i-alpha))
                                                                               *(1.0-fabs((double)j-beta))
                                                                               *M;
                                                                               
                                                                               
                                                                               
                    }
                }
            }
        }
    }
    // return number of samples used
    return count;
}








/** @brief Extract keypoint feature vector
 * 
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key
 * @param theta_key      keypoint coordinates.
 * 
 * @param gradX
 * @param gradX          precomputed gradient relative to the nearest scale.
 * 
 * ----PARAM----
 * @param lambda_descr   (=6)
 *                       - The gaussian window has a standard deviation of lambda_descr X=* sigma_key
 *                       - The patch P^descr is ( 2 * lambda_descr * sigma_key * (1+1/Nhist) wide
 *                       
 * @param Nhist          (=4) number of histograms in each of the two directions,
 * @param Nbins          (=8) number of bins covering the range [0,2pi]
 *
 * ----OUTPUT----
 * @output descr         the output is a vector of double here but it should be a vector of integers, descr being quantized. 
 * 
 * ----RETURN----
 * @return count         number of sample contributing to the descriptor
 */
int extract_sift_feature_vector_bis(double x_key, double y_key, double sigma_key, double theta_key,
                                double* gradX, double* gradY, int width, int height,
                                int Nhist, int Nbins, double lambda_descr,
                                double* descr){

    
    int count=0;   // number of contributing samples
    double M;      // gradient magnitude
    int    si,sj;             // pixel coordinates in the image referential
    double sX,sY,sOri;        // pixel coordinates in the descriptor referential
    double alpha, beta,gamma; // used to compute histogram indices and weighting of trilinear distribution
    
    /// Contributing pixels are inside a patch [siMin;siMax]X[sjMin;sjMax] of 2*lambda_descr*sigma_key*(nhist+1)/nhist
    int siMin = MAX(0, (int)(x_key-sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5));
    int sjMin = MAX(0, (int)(y_key-sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5));
    int siMax = MIN((int)(x_key+sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5), height-1);
    int sjMax = MIN((int)(y_key+sqrt(2.0)*lambda_descr*sigma_key*(1+1/(double)Nhist)+0.5), width-1);
        
    /// Initialise descr tab
    for(int i=0;i<Nhist*Nhist*Nbins;i++){descr[i] = 0.0;}
    
    /// For each pixel inside the patch.
    for(si=siMin;si<siMax;si++){
        for(sj=sjMin;sj<sjMax;sj++){
            
            /// Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            sX = (cos(-theta_key)*((double)si-x_key)-sin(-theta_key)*((double)sj-y_key))/sigma_key;
            sY = (sin(-theta_key)*((double)si-x_key)+cos(-theta_key)*((double)sj-y_key))/sigma_key;
            
            /// Does this sample fall inside the descriptor area ?
            if((fabs(sX)<lambda_descr*(((double)Nhist+1)/(double)Nhist))&(fabs(sY)<lambda_descr*(((double)Nhist+1)/(double)Nhist))){
                count++;

                /// Compute the gradient orientation (theta) on keypoint's invariant referential.
                sOri = fmod(atan2(gradY[si*width+sj], gradX[si*width+sj]) - theta_key + 10*2*PI,2*PI);                
                assert(sOri>0.0); assert(sOri<2*PI);
                
                /// Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to gradients far from the keypoint.
                M = sqrt(gradY[si*width+sj]*gradY[si*width+sj]+gradX[si*width+sj]*gradX[si*width+sj])
                    *exp(-(sX*sX+sY*sY)/(2*lambda_descr*lambda_descr));
                
                /// Determine the histogram and bin indices, Compute the (tri)linear weightings ...
                alpha = sX/(2*lambda_descr/(double)Nhist) + ((double)Nhist-1.0)/2.0;
                beta  = sY/(2*lambda_descr/(double)Nhist) + ((double)Nhist-1.0)/2.0;
                gamma = sOri/(2*PI)*(double)Nbins;
                
                ///    ...and add contributions to respective bins in different histograms.
                for(int i=MAX(0,(int)floor(alpha));i<=MIN((int)ceil(alpha),Nhist-1);i++){
                    for(int j=MAX(0,(int)floor(beta));j<=MIN((int)ceil(beta),Nhist-1);j++){ // looping through all surrounding histograms.
                        
                        
                        assert(i*Nhist*Nbins+j*Nbins+ ((int)ceil(gamma)+2*Nbins)%Nbins>=0);
                        assert(i*Nhist*Nbins+j*Nbins+ ((int)ceil(gamma)+2*Nbins)%Nbins<Nhist*Nhist*Nbins);
                        assert(i*Nhist*Nbins+j*Nbins+ (int)floor(gamma)>=0);
                        assert(i*Nhist*Nbins+j*Nbins+ (int)floor(gamma)<Nhist*Nhist*Nbins);

                        
                        // Contribution to left bin.                        
                        descr[i*Nhist*Nbins+j*Nbins+ (int)floor(gamma)] += (1.-(gamma-floor(gamma)))
                                                                          *(1.0-fabs((double)i-alpha))
                                                                          *(1.0-fabs((double)j-beta))
                                                                          *M;
                        // Contribution to right bin.
                        descr[i*Nhist*Nbins+j*Nbins+ ((int)ceil(gamma)+2*Nbins)%Nbins] += (1.0-(ceil(gamma)-gamma))
                                                                               *(1.0-fabs((double)i-alpha))
                                                                               *(1.0-fabs((double)j-beta))
                                                                               *M;
                                                                               
                                                                               
                                                                               
                    }
                }
            }
        }
    }
    // return number of samples used
    return count;
}











/** @brief Modify descriptor for robustness and speed
 *
 *  - Threshold bins exceeding the value of threshold times the l2 norm
 *  (to increase robustness to saturation effects)
 *
 *  - Quantize descrtiptor values
 *  (to increase the speed of distance computations)
 *
 * ----INPUT----
 * @param descr                 double histogram of length Nhist*Nhist*Nbins (=128)
 * 
 * ----PARAM----
 * @param flagthreshold         bool (if true, threshold the descriptor)
 * @param threshold             threshold applied to the l2-norm of the descriptor
 *
 * TODO change (now it's only double tab)
 * ----OUTPUT----
 * @param qtdescr               int histogram of length Nhist*Nhist*Nbins and values in [0..255]
 *
 */
void threshold_and_quantize_sift_descriptor(double* descr, int Nhist, int Nbins,
                                            bool flagthreshold, double threshold){
                                            //int* qtdescr){
    int i;
    double l2norm, maxVal;
 
    
    maxVal=-1.0; 
    for(i=0;i<Nhist*Nhist*Nbins;i++){if(descr[i]>maxVal){maxVal=descr[i];}}
    //for(i = 0;i<Nhist*Nhist*Nbins;i++){qtdescr[i] = (int)(descr[i]/maxVal*255);}   //TODO TODO TODO
    for(i=0;i<Nhist*Nhist*Nbins;i++){descr[i] = descr[i]/maxVal*255;}
    
    
    /// Threshold bins
    if(flagthreshold){
        l2norm=0;
        for(i=0;i<Nhist*Nhist*Nbins;i++){
            l2norm += descr[i]*descr[i];
            
        } l2norm = sqrt(l2norm);
        for(i=0;i<Nhist*Nhist*Nbins;i++){
            descr[i] = MIN(descr[i],threshold*l2norm);
        }
    }
    
     /// Quantize the descriptor to get an integer vector
     maxVal=-1.0; 
     for(i=0;i<Nhist*Nhist*Nbins;i++){if(descr[i]>maxVal){maxVal=descr[i];}}
     //for(i = 0;i<Nhist*Nhist*Nbins;i++){qtdescr[i] = (int)(descr[i]/maxVal*255);}   //TODO TODO TODO
     for(i=0;i<Nhist*Nhist*Nbins;i++){descr[i] = floor(descr[i]/maxVal*255);}
//    
}













/** @brief Iterative box filter of width 3 bins
 *
 * ----INPUT----
 * @param hist           : histogram
 * @param nbins          : number of bins
 *
 * ----PARAM----
 * @param niter          : number of iteration
 *
 * equivalent to gaussian smoothing of width TODO
 */
void smooth_circular_histogram(int niter, double* hist, int nbins){
    
    int i,i_prev,i_next;    // bin indices
    double hist_tmp[nbins]; // CHECK TODO validity
    
    /// Initialization
    for(i=0;i<nbins;i++){hist_tmp[i] = hist[i];}
     
    /// Convolution with box filters
    for(;niter>0;niter--){
        for(i=0;i<nbins;i++){hist_tmp[i] = hist[i];}
        for(i=0;i<nbins;i++){
            i_prev = (i-1+nbins)%nbins;
            i_next = (i+1)%nbins;
            hist[i] = (hist_tmp[i_prev]+hist_tmp[i]+hist_tmp[i_next])/3.;
        }
    }
}










/** @brief Extract keypoint feature vector WITHOUT trilinear distribution
 *
 * without spatial bilinear distribution
 * but with linear attribution in the orientation dimension
 * 
 */
int extract_sift_feature_vector_2_COARSE(double x_key, double y_key, double sigma_key, double theta_key,
                                  double* gradX, double* gradY, int width, int height,
                                  int Nhist, int Nbins, double sigma_window, double K,
                                  double* descr){

    int i;   
    int count=0;         /* number of samples in normalized patch P_descr */
    double M;            /* gradient magnitude */
    int    si,sj;        /* pixel coordinates in the image referential */
    double sX,sY,sOri;   /* pixel coordinates in the descriptor referential */
    int alpha, beta;     /* histogram indices */
    double gamma;        /* used to compute histogram indices and weighting of trilinear distribution */
    
    /// contributing pixels are inside a patch [siMin;siMax]X[sjMin;sjMax] of K*sigma_window*sigma_key
    int siMin = MAX(0, (int)(x_key-sqrt(2.0)*K*sigma_window*sigma_key+0.5));  int siMax = MIN((int)(x_key+sqrt(2.0)*K*sigma_window*sigma_key+0.5), height-1);
    int sjMin = MAX(0, (int)(y_key-sqrt(2.0)*K*sigma_window*sigma_key+0.5));  int sjMax = MIN((int)(y_key+sqrt(2.0)*K*sigma_window*sigma_key+0.5), width-1);
        
    // initialise descr tab
    for(i = 0;i<Nhist*Nhist*Nbins;i++){descr[i] = 0.0;}
    
    /// For each pixel inside the patch.
    for(si=siMin;si<siMax;si++){
        for(sj=sjMin;sj<sjMax;sj++){
            
            /// Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            sX = (cos(-theta_key)*((double)si-x_key)-sin(-theta_key)*((double)sj-y_key))/sigma_key;
            sY = (sin(-theta_key)*((double)si-x_key)+cos(-theta_key)*((double)sj-y_key))/sigma_key;
            
            /// Does this sample fall inside the descriptor area ?
            if((fabs(sX)<0.5*K*sigma_window)&(fabs(sY)<0.5*K*sigma_window)){
                count++;

                /// Compute the gradient orientation (theta) on keypoint's invariant referential.
                sOri = fmod(atan2(gradY[si*width+sj], gradX[si*width+sj]) - theta_key + 10*2*PI,2*PI);                
                
                /// Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to gradients far from the keypoint.
                M = sqrt(gradY[si*width+sj]*gradY[si*width+sj]+gradX[si*width+sj]*gradX[si*width+sj])*exp(-(sX*sX+sY*sY)/(2*sigma_window*sigma_window));
                /// Determine the histogram and bin indices...
                alpha = (int)((sX+0.5*K*sigma_window)/(K*sigma_window/(double)Nhist));
                beta  = (int)((sY+0.5*K*sigma_window)/(K*sigma_window/(double)Nhist));
                gamma = sOri/(2*PI)*(double)Nbins;
                
                ///    ...and add contributions to respective bins.
                assert((int)floor(gamma)>=0);
                assert((int)floor(gamma)<Nbins);
                assert((int)ceil(gamma)%Nbins>=0);
                assert((int)ceil(gamma)%Nbins<Nbins);
                
                /// ONLY two bins are affected - No spatial bilinear dispatching
                descr[alpha*Nhist*Nbins+beta*Nbins+ (int)floor(gamma)]      += (1.0-(gamma-floor(gamma)))*M;
                descr[alpha*Nhist*Nbins+beta*Nbins+ (int)ceil(gamma)%Nbins] += (1.0-(ceil(gamma)-gamma))*M;
            }
        }
    }
    return count; /* returns the number of samples in the normalized patch */
}
