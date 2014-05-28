/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

#include "core.hpp"



namespace cv
{

/*!
 Semi-Global Block Matching Stereo Correspondence Algorithm

 The class implements the original SGBM stereo correspondence algorithm by H. Hirschmuller and some its modification.
 */
class CV_EXPORTS_W StereoSGBM
{
public:
    enum { DISP_SHIFT=4, DISP_SCALE = (1<<DISP_SHIFT) };

    //! the default constructor
    CV_WRAP StereoSGBM();

    //! the full constructor taking all the necessary algorithm parameters
    CV_WRAP StereoSGBM(int minDisparity, int numDisparities, int SADWindowSize,
               int P1=0, int P2=0, int disp12MaxDiff=0,
               int preFilterCap=0, int uniquenessRatio=0,
               int speckleWindowSize=0, int speckleRange=0,
               bool fullDP=false);
    //! the destructor
    virtual ~StereoSGBM();

    //! the stereo correspondence operator that computes disparity map for the specified rectified stereo pair
    CV_WRAP_AS(compute) virtual void operator()(InputArray left, InputArray right,
                             OutputArray _disp , OutputArray _cost);

    CV_PROP_RW int minDisparity;
    CV_PROP_RW int numberOfDisparities;
    CV_PROP_RW int SADWindowSize;
    CV_PROP_RW int preFilterCap;
    CV_PROP_RW int uniquenessRatio;
    CV_PROP_RW int P1;
    CV_PROP_RW int P2;
    CV_PROP_RW int speckleWindowSize;
    CV_PROP_RW int speckleRange;
    CV_PROP_RW int disp12MaxDiff;
    CV_PROP_RW bool fullDP;

protected:
    Mat buffer;
};

//! filters off speckles (small regions of incorrectly computed disparity)
CV_EXPORTS_W void filterSpeckles( InputOutputArray img, double newVal, int maxSpeckleSize, double maxDiff,
                                  InputOutputArray buf=noArray() );

//! computes valid disparity ROI from the valid ROIs of the rectified images (that are returned by cv::stereoRectify())
CV_EXPORTS_W Rect getValidDisparityROI( Rect roi1, Rect roi2,
                                        int minDisparity, int numberOfDisparities,
                                        int SADWindowSize );

//! validates disparity using the left-right check. The matrix "cost" should be computed by the stereo correspondence algorithm
CV_EXPORTS_W void validateDisparity( InputOutputArray disparity, InputArray cost,
                                     int minDisparity, int numberOfDisparities,
                                     int disp12MaxDisp=1 );

//! reprojects disparity image to 3D: (x,y,d)->(X,Y,Z) using the matrix Q returned by cv::stereoRectify
CV_EXPORTS_W void reprojectImageTo3D( InputArray disparity,
                                      OutputArray _3dImage, InputArray Q,
                                      bool handleMissingValues=false,
                                      int ddepth=-1 );

CV_EXPORTS_W  int estimateAffine3D(InputArray src, InputArray dst,
                                   OutputArray out, OutputArray inliers,
                                   double ransacThreshold=3, double confidence=0.99);

}

