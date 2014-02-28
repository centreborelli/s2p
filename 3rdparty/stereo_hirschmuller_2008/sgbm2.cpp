/**
 * @file sgbm.cpp
 * @brief wrapper OpenCV's implementation of the SGBM stereo algorithm (Hirschmuller'08)
 * @author Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
 *
 * Copyright (c) 2012-2013, Gabriele Facciolo
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or redistribute it
 * under the terms of the simplified BSD License. You should have received a
 * copy of this license along this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <cstdio>
extern "C" {
#include "iio.h"
}

using namespace cv;

void paste(Mat &dest, Mat &overlay, int px, int py)
{
    overlay.copyTo(dest(Range(py, py + overlay.rows),
        Range(px, px + overlay.cols)));
}

// paste a rectangular region of source into target.
// x0,y0 provides the alignement of target over the source, it is the position
// over target of the 0,0 pixel of source
// relative_x/ystart is the starting point of the copy,
// while (w,h) is the size of the area to be copied
void paste(Mat &dest, Mat &overlay, 
		 int x0, int y0, 
       int relative_xstart, int relative_ystart,
       int w, int h)
{
   Mat overlaysub = overlay(Range(relative_ystart, relative_ystart+h), Range(relative_xstart, relative_xstart+w));
   overlaysub.copyTo( dest(Range(y0+relative_ystart, y0+relative_ystart + overlaysub.rows),
            Range(x0+relative_xstart, x0+relative_xstart + overlaysub.cols)));
}

// FROM OPENCV's DOCUMENTATION
//
//
// The class implements the modified H. Hirschmuller algorithm [HH08] that
// differs from the original one as follows:
//
// By default, the algorithm is single-pass, which means that you consider only
// 5 directions instead of 8. Set fullDP=true to run the full variant of the
// algorithm but beware that it may consume a lot of memory.  The algorithm
// matches blocks, not individual pixels. Though, setting SADWindowSize=1
// reduces the blocks to single pixels.  Mutual information cost function is
// not implemented. Instead, a simpler Birchfield-Tomasi sub-pixel metric from
// [BT98] is used. Though, the color images are supported as well. Some pre-
// and post- processing steps from K. Konolige algorithm StereoBM::operator()
// are included, for example: pre-filtering (CV_STEREO_BM_XSOBEL type) and
// post-filtering (uniqueness check, quadratic interpolation and speckle
// filtering).
//
// PARAMETERS :
//
// minDisparity – Minimum possible disparity value. Normally, it is zero but
// sometimes rectification algorithms can shift images, so this parameter needs
// to be adjusted accordingly.
// numDisparities – Maximum disparity minus minimum disparity. The value is
// always greater than zero. In the current implementation, this parameter must
// be divisible by 16.
// SADWindowSize – Matched block size. It must be an odd number >=1 .
// Normally, it should be somewhere in the 3..11 range.
// P1 – The first parameter controlling the disparity smoothness. See below.
// P2 – The second parameter controlling the disparity smoothness.  The larger
// the values are, the smoother the disparity is. P1 is the penalty on the
// disparity change by plus or minus 1 between neighbor pixels. P2 is the
// penalty on the disparity change by more than 1 between neighbor pixels. The
// algorithm requires P2 > P1 . See stereo_match.cpp sample where some
// reasonably good P1 and P2 values are shown (like
// 8*number_of_image_channels*SADWindowSize*SADWindowSize and
// 32*number_of_image_channels*SADWindowSize*SADWindowSize , respectively).
// disp12MaxDiff – Maximum allowed difference (in integer pixel units) in the
// left-right disparity check. Set it to a non-positive value to disable the
// check.
// preFilterCap – Truncation value for the prefiltered image pixels.  The
// algorithm first computes x-derivative at each pixel and clips its value by
// [-preFilterCap, preFilterCap] interval. The result values are passed to the
// Birchfield-Tomasi pixel cost function.
// uniquenessRatio – Margin in percentage by which the best (minimum) computed
// cost function value should “win” the second best value to consider the found
// match correct. Normally, a value within the 5-15 range is good enough.
// speckleWindowSize – Maximum size of smooth disparity regions to consider
// their noise speckles and invalidate. Set it to 0 to disable speckle
// filtering. Otherwise, set it somewhere in the 50-200 range.
// speckleRange – Maximum disparity variation within each connected component.
// If you do speckle filtering, set the parameter to a positive value, it will
// be implicitly multiplied by 16. Normally, 1 or 2 is good enough.  fullDP –
// Set it to true to run the full-scale two-pass dynamic programming algorithm.
// It will consume O(W*H*numDisparities) bytes, which is large for 640x480
// stereo and huge for HD-size pictures. By default, it is set to false .




Mat sgbm_stereo(Mat &u1, Mat &u2, int mindisp, int maxdisp, int SADwin, int P1, int P2, int LRdiff) {
    // swap max and min
    int tmp=-maxdisp;
    maxdisp=-mindisp;
    mindisp=tmp;

    // prepare the StereoSGBM object
    // the number of disparities has to be a multiple of 16
    StereoSGBM sgbm;
    sgbm.numberOfDisparities = (int) 16 * ceil((maxdisp - mindisp) / 16.0);
    sgbm.minDisparity = mindisp;
    sgbm.SADWindowSize = SADwin;
    sgbm.P1 = P1;
    sgbm.P2 = P2;
    sgbm.disp12MaxDiff = LRdiff;

    sgbm.fullDP = 1; // to run the full-scale two-pass dynamic programming algorithm
    sgbm.preFilterCap = 63;
    sgbm.uniquenessRatio = 10;
    sgbm.speckleWindowSize = 50;
    sgbm.speckleRange = 1;

    fprintf (stderr, "mind:%d maxd:%d win:%d P1:%d P2:%d LR:%d\n", mindisp,
        mindisp + sgbm.numberOfDisparities, SADwin, P1, P2, LRdiff);
    fprintf (stderr, "fullDP:%d preFilterCap:%d uniquenessRatio:%d"
        "speckleWin:%d SpeckleRange:%d\n", sgbm.fullDP, sgbm.preFilterCap,
        sgbm.uniquenessRatio, sgbm.speckleWindowSize, sgbm.speckleRange);

    // CROP TRICK
    // To avoid losing two vertical bands on the left and right sides of the
    // reference image
    // the stereo matching is done by the sgbm() operator
    Mat uu1(u1.rows, (int) (u1.cols + max(-mindisp, 0) ), u1.type());
    Mat uu2(u2.rows, (int) (u2.cols + max(-mindisp, 0) ), u2.type());
    paste(uu1, u1,  0, 0);
    paste(uu2, u2,  0, 0);
    Mat disp, ddisp;
    sgbm(uu1, uu2, ddisp);
    disp = ddisp(Range::all(), Range( 0, u1.cols));

//    // WITHOUT CROP TRICK
//    Mat disp, disp_8u;
//    sgbm(u1, u2, disp);
//

    Mat fdisp(u1.rows, u1.cols,DataType<float>::type);
    for (int y=0; y<u1.rows; y++)
        for (int x=0; x<u1.cols; x++)
            if (disp.at<int16_t>(y, x) == -((-mindisp+1)*16))
                fdisp.at<float>(y,x) = NAN;
            else
                fdisp.at<float>(y,x) = -((float) disp.at<int16_t>(y, x)) / 16.0;
                // sgbm output disparity map is a 16-bit signed single-channel
                // image of the same size as the input image. It contains
                // disparity values scaled by 16. So, to get the floating-point
                // disparity map, you need to divide each disp element by 16.
   return fdisp;
}



Mat sgbm_stereo_slice(Mat &u1, Mat &u2, int mindisp, int maxdisp, int SADwin, int P1, int P2, int LRdiff) 
{
   int nc1 = u1.cols, nr1 = u1.rows;
   int nc2 = u2.cols;

   //// TODO TUNE THESE PARAMETERS FOR SGBM
   //
   int MIN_SLICEH = 100;   // SLICE HEIGHT
   int h = 100;            // OVERLAP HEIGHT (2x)
   int sliceh = fmax(2 * h, MIN_SLICEH);
   int halfh_up = h / 2;
   int halfh_down = h - (int) (h / 2) - 1;

   Mat disp (nr1, nc1, DataType<float>::type);

   // NO NEED FOR SLICING: pass directly 
   if (sliceh >= nr1) {
      disp = sgbm_stereo(u1, u2, mindisp, maxdisp, SADwin, P1, P2, LRdiff );
      return disp;
   }

#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (int s = 0; s < nr1; s += sliceh) {

      // compute slice and padding parameterr
      int ymin = fmax(0, s - halfh_up);
      int ymax = fmin(s + sliceh + halfh_down, nr1);
      int top_padding = s - ymin;
      int bottom_padding = ymax - fmin(s + sliceh, nr1);

      Mat sliceu1 = u1(Range( ymin, ymax ), Range::all());
      Mat sliceu2 = u2(Range( ymin, ymax ), Range::all());
      Mat slicedisp  = sgbm_stereo(sliceu1, sliceu2, mindisp, maxdisp, SADwin, P1, P2, LRdiff );

#ifdef _OPENMP
#pragma omp critical
#endif
      paste(disp, slicedisp,
		  0, ymin, 0, top_padding, nc1,
		  ymax - ymin - top_padding - bottom_padding);

      // ALTERNATIVE WAY TO DO THE PASTE
//      for(int y=0;y<ymax - ymin - top_padding - bottom_padding;y++) {
//      for(int x=0;x<nc1;x++) {
//         disp.at<float>(ymin+top_padding+y,x) = slicedisp.at<uint8_t>(top_padding+y,x);
//      }
//      }

   }
   return disp;
}




void stereoLRRL(Mat &dL, Mat &dR, Mat *out, float tolerance)
{
  int nc = dL.cols,      nr =dL.rows;
  int rightnc = dR.cols, rightnr =dR.rows;

//  printf("%f\n",tolerance);

  /* explore the left image */
  for(int y=0;y<nr;y++) 
     for(int x=0;x<nc;x++) {
        out->at<float>(y,x)=0;
        int L = (int) dL.at<float>(y,x);
        if((x+L)<rightnc && (x+L)>=0 && fabs(dR.at<float>(y,x+L) +(float)L) <= tolerance)
           out->at<float>(y,x) = 255;
     }
}



int main(int c, char** v)
{
    if (c < 4) {
        fprintf(stderr, "\tusage: %s im1 im2 out [mindisp(0) maxdisp(64)"
            "SADwindow(1) P1(0) P2(0) LRdiff(1)]\n", v[0]);
        exit(1);
    }
    // TODO : READ WITH IIO AND CONVERT TO INT16 or INT8? INTERNALLY
    Mat u1 = imread(v[1], 0);
    Mat u2 = imread(v[2], 0);
    int nch = u1.channels();

    int i = 4;
    // Warning: here we convert between openCV (and Hirschmuller, and
    // Middlebury) disparity convention and ours. What we call disparity is the
    // opposite of what they call disparity.
    // the min and max disparities given as inputs of this binary follow our
    // convention, while in the code the disparities follow openCV convention.
    int mindisp = (c>i) ? atoi(v[i]) : 0;  i++;
    int maxdisp = (c>i) ? atoi(v[i]) : 64; i++;
    int SADwin  = (c>i) ? atoi(v[i]) : 1;  i++;
    int P1      = (c>i) ? atoi(v[i]) :  8*nch*SADwin*SADwin; i++;
    int P2      = (c>i) ? atoi(v[i]) : 32*nch*SADwin*SADwin; i++;
    int LRdiff  = (c>i) ? atoi(v[i]) : 1;  i++;

    if (mindisp>=maxdisp){
        fprintf(stderr, "\terror: mindisp(%d) > maxdisp(%d)\n", mindisp, maxdisp); 
        exit(1);
    }

    Mat disp  = sgbm_stereo_slice(u1, u2, mindisp, maxdisp, SADwin, P1, P2, LRdiff );
    Mat disp2 = sgbm_stereo_slice(u2, u1, -maxdisp , -mindisp, SADwin, P1, P2, LRdiff );

    // left right test
    Mat LRmask (disp.rows, disp.cols, DataType<float>::type);
    stereoLRRL(disp, disp2, &LRmask, 2);

    // convert back the disparities to our convention before saving
    float *odisp = (float*) malloc(disp.rows*disp.cols*sizeof(float));
    for (int y=0; y<disp.rows; y++)
       for (int x=0; x<disp.cols; x++){
          if(LRmask.at<float>(y, x) )
             odisp[x+y*disp.cols] = disp.at<float>(y, x);
          else
             odisp[x+y*disp.cols] = NAN;
       }
    

    // save
    iio_save_image_float(v[3], odisp, disp.cols, disp.rows);
    free(odisp);
    return 0;
}
