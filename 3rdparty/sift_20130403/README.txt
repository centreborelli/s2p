IPOL SIFT
Copyright (C) 2013, Ives Rey Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20130327 (March 27, 2013)


===============================================================================
== Overview ===================================================================

This C ANSI source code implements the SIFT method originaly described in
    [1] "Distinctive Image Features from Scale Invariant"
        D.G Lowe
        International Journal of Computer Vision, vol 60, pp 91-110, 2004.
    
The SIFT method is also described in IPOL publication
    [2] "An Anatomy of the SIFT Method." 
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2013.
        http://www.ipol.im/pub/algo/rd_sift_anatomy/
        
It is featured in the IPOL demo page
        http://www.ipol.im/pub/demo/rd_sift_anatomy/
  
  
  

===============================================================================
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
 
 



===============================================================================
== Compiling (Linux) ==========================================================

To compile IPOL_SIFT

Type

    make

in the directory where the Makefile is located. The compilation of the source
 code provides two executables:

1) IPOL_SIFT  applies the SIFT method to a PNG image. Its uses either standard
              parameters (as documented in [1]) user selected parameters.
2) MATCHING   matches the SIFT keypoints extracted from two image.

     
===============================================================================     
=== The standard IPOL_SIFT ====================================================


The standard IPOL_SIFT command is
     
     ./ipol_sift image

The output is a list of keypoints (one keypoint per line) with the formatting 

          x y sigma theta fv[1] fv[2]  ...   fv[128]

where (x,y) is the position, (sigma) the detection scale, (theta) the reference 
orientation and (fv) the 128-D feature vector. The output can be redirected
into an ASCII file with

     ./ipol_sift image > keypoints.txt
  
     

===============================================================================
=== The customizable IPOL_SIFT ================================================


The customizable IPOL_SIFT command is 
   
   
     ./ipol_sift image label n_oct   n_spo   sigma_min
                                delta_min   sigma_in   C_DoG
                                     C_edge    n_bins    lambda_ori    t
                                             n_hist   n_ori   lambda_descr  
     ./ipol_sift image label 8 3 0.8 0.5 0.5 0.03 10 36 1.5 0.8 4 8 6.0 

     
The output is a list of keypoints with the following formatting
     
      x  y  sigma  theta  octa  sca fv[1] fv[2] ... fv[d]
                                               orihist[1] ... orihist[n_bins]

where (octa) and (sca) are octave and scale indices, fv the d-dimensional
feature vector, and orihist the orientation histogram.     
        
     
The following table summarizes the parameters of the customizable IPOL_SIFT.
(For details, please refer to the IPOL publicatio  "An anatomy of 
the SIFT method.")


-------------------------------------------------------------------------------
Parameter    std value  Definition
-------------------------------------------------------------------------------
n_oct           8       number of octaves in the scale-space,
n_spo           3       number of scales per octave,
sigma_min       0.8     minimal level of blur featured in the scale-space,
delta_min       0.5     minimal inter-pixel distance  featured in the
                        scale-space,
sigma_in        0.5     assumed level of blur in the input image, 

C_DoG           0.03    threshold on DoG operator (expressed for n_spo = 3),
C_edge         10       threshold on the ratio of principal curvatures,

n_bins         36       number of bins in the orientation histogram,
lambda_ori      1.5     sets the width of the orientation patch,
t               0.8     reference orientation minimal value in the histogram, 

n_hist          4       the descriptor is composed of n_histXnhist weighted
                        histograms,
n_ori           8       each weighted histogram is composed of n_ori bins,
lambda_descr    6.0     sets the width of the descriptor patch.

label                   string used to label all the extra output files.
-------------------------------------------------------------------------------

     

     
     
     
The following tables lists the extra files produced
-------------------------------------------------------------------------------
  Extra output files                      Content
-------------------------------------------------------------------------------
 1) extra_NES_[label].txt                 discrete 3D extrema of DoG,
 2) extra_DoGSoftThresh_[label].txt       discrete 3D extrema passing a
                                          conservative threshold on DoG, 
 3) extra_ExtrInterp_[label].txt          interpolated 3D extrema,
 4) extra_ExtrInterpREJ_[label].txt       detectop, discarded during
                                          interpolation,
 5) extra_DoGThresh_[label].txt           interpolated extrema passing the
                                          threshold on DoG,
 6) extra_OnEdgeResp_[label].txt          interpolated extrema passing the
                                          Harris-Stephen edgeness test,
 7) extra_OnEdgeRespREJ_[label].txt       interpolated extrema failling the
                                          Harris-Stephen edgeness test,
 8) extra_OriAssigned_[label].txt         keypoints with reference orientation,
 9) extra_OriAssignedMULT_[label].txt     subset of keypoints with multiple
                                          reference orientation.

       Each line of theses files follows the data formatting
           x  y  sigma  theta  octa  sca
       where (octa) is the octave index and (sca) is the scale index.

10) extra_keypoints_[label].txt           keypoints with descriptors
     
        Each line of this file follows the data formatting
            x  y  sigma  theta  octa  sca fv[1] fv[2] ... fv[d] ...
                                         ... orihist[1] ... orihist[n_bins] 
        where (fv) is the feature vector of dimension d=n_hist*n_hist*n_ori 
        and (orihist) is the orientation histogram of n_bins bins.
     



===============================================================================     
=== The IPOL_SIFT_MATCHING ====================================================

The MATCHING command is

      ./ipol_sift_matching keys1.txt keys2.txt methodFlag C n_hist n_ori n_bins
      ./ipol_sift_matching keys1.txt keys2.txt 1 0.6 4 8 36
    
The output is a list of matches with the following formatting
      x1  y1  sigma1  theta1   x2  y2  sigma2  theta 2       
     
      
The algorithm parameters are summarized in the following table
-------------------------------------------------------------------------------
Parameter    std value   Definition
-------------------------------------------------------------------------------
methodFlag    (1 / 0)    Selects matching method.
                         - relative(1) : threshold on the distance ratio
                                         (distance to nearest/
                                           distance to second nearest distance)
                         - absolute(0) : threshold on the distance to the
                                         nearest keypoint.

C_match      (0.6/ 100)  threshold value,

n_hist          4        the descriptor is composed of n_histXnhist weighted
                         histograms,
n_ori           8        each weighted histogram is composed of n_ori bins,
n_bins         36        Number of bins in the orientation histogram.
-------------------------------------------------------------------------------
      
     

Extra information relative to the matching algorithm is stored in:

-------------------------------------------------------------------------------
Filename                     Content 
-------------------------------------------------------------------------------
1) OUTpairsAll.txt           All pair of keypoints,
2) OUTmatches.txt            Only pairs satisfying the matching criteria,
3) OUTkeys1Matching.txt      Subset of matching keypoints in image 1,
4) OUTkeys2Matching.txt      Subset of matching keypoints in image 2,
5) OUTkeys1NotMatching.txt   Subset of non matching keypoints in image 1,
-------------------------------------------------------------------------------


Files 1) and 2 ) have the following formatting:
      
      key1  key2a  key2b
      
where (key1) designates a keypoint in image1, (key2a) and (key2b) designate
respectively the nearest and the second nearest neighbors in image 2.
The data relative to each keypoint is formatted as follows

      x  y  sigma  theta  fv[1] fv[2] ... fv[d]
                        octa sca   orihist[1] ... orihist[n_bins] 

where (fv) is the feature vector of dimension d=n_hist*n_hist*n_ori  and
(orihist) is the orientation histogram of n_bins bins.    
Files 3) 4) and 5) have the same formatting:

      x  y  sigma  theta  fv[1] fv[2] ... fv[d]
                        octa sca   orihist[1] ... orihist[n_bins]  




===============================================================================
== Generating the doxygen documentation =======================================

In the src/ directory type :
  
   doxygen -g 
   doxygen Doxyfile

doxygen documentation in directory ./html/ 


      
===============================================================================
== Acknowledgements ===========================================================

a
