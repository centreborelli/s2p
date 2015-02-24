#ifndef _LIBSTEREO_H
#define _LIBSTEREO_H


#include "../libraryStable/libraryBasic.h"

#include <iostream>
//#include <omp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>


#define L2 0
#define L2M 1
#define STEREOVERBOSE 0



namespace libIIPStable {
	
	
	
	//! Pixel Precision Parameters
	struct strParameters
	{
	
        
       //! GF. input image segmentation used to regularize the result 
       // THIS OPTION IS NOT WORKING YET
       flimage segments;
        
        //! Prolate Window
        flimage prolate;
        int flagWinWeighted;     //! If it is flat distances can be accelerated or even computed by integral image
                                 //! A non square window and constant is non weighted.
        
        int flagListKernels;     //! A weighted kernel can have any shape, the ListKernel is a list of offsets with weights
                                 //! If this flag is activated the kernel is interpteted by rows as : {offset_x, offset_y, kern_values}
        
        //! For the moment SSD, SSD-mean of the window
        int itypeDist;
		
        
		//! Number of precisions used for pixelian
		int inPrecisions;
      //! factor of disparity refinement by cubic interpolation
      int cubicRefinementFactor;
		
        
		//! Number of scales and current
		int nScales;
        int currentScale;
        
        //! (dmin, dmax) minimal and maximal for current scale.  Computed as the initial / 2^currentscale.
		// Necessary when a point is detected as non correctly correlated. These values are passed to the next scale
		float dmin0, dmax0, idmin0, idmax0;
        
        
        //! Subpixel reciprocity flag. In order to avoid computation if not
        //! used afterwards. Only activated if inPrecision > 1.
        int flagSubRecip;
        float valueSubRecip;

        
        //! If activated left-right and right-left consistency is checked. Usually activated
		int flagRecip;
        float valueRecip;
        
        
        //! Alternative left-right and right-left checkings.
        //! Applied only with parameter 1
        int flagRecipDual;
        int flagiRecip;
        
        
        //! Self similarity. Flag already used for pixelian in order
        //! to avoid computation of best match if not used afterwards
        int flagSelfSimilarity;
        float valueSelfSimilarity;
        
        //! Grain Filter
        int flagGrainFilter;
        float valueGrainFilter;

        //! Cost Filter
        int flagCostFilter;
        float valueCostFilter;
        
        //! Horizontal variance. To be used with squared windows.
        //! Being replaced by flagHorInfo
        int flagVariance;
        float valueVariance;
        
        
        //! Bernard precision term.
        //! Rewrites as the sum of horizontal derivatives if kernel is constant
        int flagHorInfo;
        float valueHorInfo;
        
        
        
        //! Remove isolated points
        //! Compares number of known neighbours to percent of the window
        int flagRemoveIsolated;
        float valueRemoveIsolated;
        
        
        
        //! Noise std
		float fSigma;
        
        
        //! Min dist check and value comparison
        int flagMinDist;
        int flagMinDistDilatation;
        float valueMinDist;
        int flagRegressionMinDist;
        
		
        //! Multiple windows flag (only used by:
        int flagMultiWin;


        //! Flag to allow combining restul of last scale with the previous one 
        //  This densifies the result.
        int flagCombineLastScale;


	};
    
    
    
    
    
    
    
    
    
    
    
	//! SubPixel Parameters
	
    struct strSubPixelPar
	{
		
    
		flimage prolate;
		int sflagWinWeighted;
		
		int inScales;
		int fnScales;
        
        int sflagRafa;
        float sValueRafa;

		int sflagSameDepth;
        
        float sfSigma;
		
    
	 /*
		
		
		int sflagPrecision;
		flimage sPrecision;
		
		
		int  sflagRecip;
		flimage idispr;
		flimage imaskr;
		flimage sReciprocity;
	*/
        
	};
	
    
    
	
    
    
    
    
    
    
    
    
    
    
    
    
    
    //************** Not here yet ****************//
    
    
    //! Window
    void  stereo_build_gaussian_kernel(flimage &prolate, float alpha, float beta);
    
    
    //! Pixelian correlation
    void stereo_pixel_chain_one_direction(cflimage &input, cflimage &input2, flimage &imask, flimage &Dmin, flimage &Dmax, flimage &odisp, flimage &omask, flimage &odist, flimage &odist2, flimage &oselfdist, flimage &olrfdist, flimage &ocdif, strParameters &strPar);
    
    
    //! Multi Scale Pixelian correlation
    
    
    void set_strParameters_for_current_scale(strParameters & strIn, strParameters & strOut, int iiScale);
    void update_dmin_dmax(flimage &Dmin, flimage &Dmax, flimage &out, flimage &omask, float minim, float maxim, int pwidth, int pheight);
    void update_dmin_dmax_window(flimage &Dmin, flimage &Dmax, flimage &out, flimage &omask, float minim, float maxim, flimage &prolate);
    void stereo_pixel_chain(cflimage &input, cflimage &input2, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &odisp, flimage &odisp2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar);

    void stereo_pixel_chain_multi_window(cflimage &input, cflimage &input2, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &odisp, flimage &odisp2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar);
    
    void stereo_pixel_multiscale_chain(cflimage &input, cflimage &input2, int &in, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &out, flimage &out2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar);
    
    
    
    
    //! Strobe effect
    void stereo_find_local_minima(float *v, int n, float &dif1, float &dif0);
    
    
    
    
    void stereo_detect_occlusions(flimage &idisparity, flimage &imask, flimage &omask);

    
    
    
    //////////////////////////////    //////////////////////////////    //////////////////////////////
    //////////////////////////////    //////////////////////////////    //////////////////////////////
    //////////////////////////////    //////////////////////////////    //////////////////////////////
    
    void stereo_diagonal_flat_window(int win, int d, int inverse, flimage &kernel);
    
    void stereo_adapt_window_rafa(int ipx,int ipy, flimage &corrwindow, flimage &grad, float fRafaThresh);
    void stereo_adapt_window_to_depth(int ipx, int ipy, flimage &corrwindow, flimage &idisparity, flimage &imask, float fprec);
    
    void stereo_close_segments(flimage &imask, flimage &omask, int iLength);
    
    //! Checks Reciprocity
    void stereo_check_pixelian_reciprocity(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold, flimage &oDif);
    void stereo_check_pixelian_reciprocity_remove_only_if_inconsistent(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold, flimage &oDif);
    
    void stereo_check_pixelian_reciprocity_dual(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold);
	
    void stereo_check_inverse_pixelian_reciprocity(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold);
    
    
    
    //! Check Enough information
    void dilate_mask_by_window_support(flimage &mask, flimage &corrwindow) ;
    void stereo_check_integral_derivatives(cflimage &input, flimage &omask, float fThreshold, float fSigma, flimage &corrwindow);
    
    
    void stereo_compute_horizontal_variance(cflimage & input, cflimage & variance, int pwidth, int pheight);
    void stereo_check_horizontal_variance(cflimage &input, flimage &mask, float fThreshold, float fSigma, int flagWholePatch, int pwidth, int pheight, flimage &oVAR);
    
    
    //! Check Strobe effect and self similarity
    void stereo_check_strobe_and_self_simililarity_effect(cflimage &input, flimage &odist1, flimage &odist2, flimage &imask, flimage &omask,  int iflagTrans, int itypeDist, float fTrans, float fSigma, float fMult, flimage &prolate, flimage &oDif, int flagListKernels);
    
    
    //! Check MinDist
    void stereo_check_min_dist(flimage &idisparity, flimage &idistance, flimage &imask, flimage &omask, float fThreshold, flimage &prolate, flimage &oDif);
    void stereo_check_regression_min_dist(flimage &idisparity, flimage &idistance, flimage &imask, flimage &omask, float fThreshold, flimage &prolate, flimage &oDif);

    
    
    //! Correlation
    void stereo_subpixel_computation(cflimage &input, cflimage &input2, flimage &idisparity, flimage &imask, flimage &out,flimage &omask, flimage &odist, strSubPixelPar &strPar);
    
    
    //! Stereo remove isolated
    void stereo_remove_isolated_points(flimage &imask, int iradius, flimage &omask, float fPercent);
    void stereo_remove_isolated_points_window(flimage &imask, flimage &prolate, flimage &omask, float fPercent);
    void stereo_grain_filter(flimage &imask, int grainarea, flimage &omask);
    float stereo_cost_filter(flimage &cost, flimage &imask, float discard_statistic , flimage &omask);

    
    //! Apply Mask
    void stereo_apply_mask(flimage &flDisparity,flimage &flBin, cflimage &out, float minF, float maxF);
    
    
    
    
    
    
}


#endif
