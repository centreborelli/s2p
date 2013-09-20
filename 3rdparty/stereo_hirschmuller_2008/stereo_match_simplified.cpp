#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <cstdio>

using namespace cv;

static void print_help(char *bin_name)
{
    printf("\nOpenCV implementation of Heiko Hirschmuller Semi-Global Matching (SGM)\n");
    printf("\tusage: %s im1 im2 out [mindisp(0) maxdisp(64) SADwindow(1) P1(0) P2(0) LRdiff(1)]\n", bin_name);
    //                0  1   2   3      4          5            6          7    8        9
}


int main(int c, char **v)
{
    if(c != 4 && c != 10)
    {
        print_help(*v);
        return 0;
    }

    // mandatory arguments
    const char* im1_filename = v[1];
    const char* im2_filename = v[2];
    const char* disparity_filename = v[3];
    Mat im1 = imread(im1_filename, -1); // -1: no color conversion
    Mat im2 = imread(im2_filename, -1);
    int cn = im1.channels();

    // default values for optional parameters
    int mindisp = 0;
    int maxdisp = 64;
    int SADwindow = 1;
    int P1 =  8*cn*SADwindow*SADwindow;
    int P2 = 32*cn*SADwindow*SADwindow;
    int LRdiff = 1;

    // read optional parameters
    if (c == 10)
    {
        mindisp = atoi(v[4]);
        maxdisp = atoi(v[5]);
        SADwindow = atoi(v[6]);
        P1 = atoi(v[7]);
        P2 = atoi(v[8]);
        LRdiff = atoi(v[9]);
    }

    // prepare the StereoSGBM object
    StereoSGBM sgbm;
    sgbm.minDisparity = mindisp;
    // the number of disparities has to be a multiple of 16
    sgbm.numberOfDisparities = (int) 16 * ceil((maxdisp - mindisp)/16.0);
    sgbm.SADWindowSize = SADwindow;
    sgbm.P1 = P1;
    sgbm.P2 = P2;
    sgbm.disp12MaxDiff = LRdiff;

    sgbm.fullDP = 1; // to run the full-scale two-pass dynamic programming algorithm
    sgbm.preFilterCap = 63;
    sgbm.uniquenessRatio = 10;
    sgbm.speckleWindowSize = 0;
    sgbm.speckleRange = 1;
// from openCV documentation:
//  preFilterCap – Truncation value for the prefiltered image pixels. The
//  algorithm first computes x-derivative at each pixel and clips its value
//  by [-preFilterCap, preFilterCap] interval. The result values are passed
//  to the Birchfield-Tomasi pixel cost function.
//  uniquenessRatio – Margin in percentage by which the best (minimum)
//  computed cost function value should “win” the second best value to
//  consider the found match correct. Normally, a value within the 5-15 range
//  is good enough.
//  speckleWindowSize – Maximum size of smooth disparity regions to consider
//  their noise speckles and invalidate. Set it to 0 to disable speckle
//  filtering. Otherwise, set it somewhere in the 50-200 range.
//  speckleRange – Maximum disparity variation within each connected
//  component. If you do speckle filtering, set the parameter to a positive
//  value, it will be implicitly multiplied by 16. Normally, 1 or 2 is good
//  enough.

    // run the sgbm
    Mat disp, disp8;
    sgbm(im1, im2, disp);

    // save the disparity map
    disp.convertTo(disp8, CV_8U);
    imwrite(disparity_filename, disp8);

    return 0;
}
