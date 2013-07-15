// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// David Lowe  "Method and apparatus for identifying scale invariant 
// features in an image and use of same for locating an object in an 
// image",  U.S. Patent 6,711,293.
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

#define _USE_MATH_DEFINES
#include <cmath>

#include "demo_lib_sift.h"
#include "library.h"
#include "filter.h"
#include "domain.h"
#include "splines.h"

#include "libLWImage/LWImage.h"
#include "libNumerics/numerics.h"
typedef LWImage<float> flimage;

#include <vector>
#include <cassert>
#include <cstdio>

#define DEBUG 0

void default_sift_parameters(siftPar &par)
{
	par.OctaveMax=3;
	par.DoubleImSize = 1;
	par.order = 7;
	par.InitSigma = 1.6f;
	par.BorderDist = 5; /*15 by Zhongwei*/
	par.Scales = 3;
	par.PeakThresh = 255.0f * 0.06f / 3.0f;
	par.EdgeThresh = 0.14f /*0.06*/ /*0.14*/; /* 0.14 by Zhongwei*/
	par.EdgeThresh1 = 0.16f /*0.08*/ /*0.16*/; /* 0.16 by Zhongwei*/
	par.OriBins  = 36;
	par.OriSigma = 1.5;
	par.OriHistThresh = 0.8f;
	par.MaxIndexVal = 0.2f;
	par.MagFactor  = 3;
	par.IndexSigma  = 1.0;
	par.IgnoreGradSign = 0;
	par.MatchRatio = 0.6f;
	par.MatchXradius = 1000000.0f;
	par.MatchYradius = 1000000.0f;

	par.noncorrectlylocalized = 0;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// SIFT Keypoint detection 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OctaveKeypoints(flimage& image, float octSize, keypointslist& keys,siftPar &par);

void FindMaxMin(const flimage* dogs,  const flimage& blur, int s, float octSize ,keypointslist& keys,siftPar &par);

bool LocalMaxMin(float val, const flimage& dog, int y0, int x0);

int NotOnEdge(const flimage& dog, int r, int c, float octSize,siftPar &par);

float FitQuadratic(std::vector<float>& offset,  const flimage* dogs, int r, int c);

void InterpKeyPoint(
	const flimage* dogs, int s, int r, int c,
	const flimage& grad, LWImage<bool>& map,
	float octSize, keypointslist& keys, int movesRemain,siftPar &par);

void AssignOriHist(
	const flimage& grad, float octSize,
	float octScale, float octRow, float octCol, keypointslist& keys,siftPar &par);

void SmoothHistogram(
	float* hist, int bins);
	
float InterpPeak(
	float a, float b, float c);

void MakeKeypoint(
	const flimage& grad, float octSize, float octScale,
	float octRow, float octCol, float angle, keypointslist& keys,siftPar &par);

void MakeKeypointSample(
	keypoint& key, const flimage& grad,
	float scale, float row, float col,siftPar &par);
	
void NormalizeVec(
	float* vec);
	
void KeySampleVec(
	keypoint& key, const flimage& grad,
	float scale, float row, float col,siftPar &par);
	
void KeySample(
	float index[IndexSize][IndexSize][OriSize], keypoint& key,
	const flimage& grad,
	float scale, float row, float col,siftPar &par);

void AddSample(
	float index[IndexSize][IndexSize][OriSize], keypoint& key,
	const flimage& grad,
	int r, int c, float rpos, float cpos, float rx, float cx,siftPar &par);

void PlaceInIndex(
	float index[IndexSize][IndexSize][OriSize],
	float mag, float ori, float rx, float cx,siftPar &par);

void compute_sift_keypoints(float *input, keypointslist& keys, int width, int height, siftPar &par)
{

	flimage image;

	/// Make zoom of image if necessary
	float octSize = 1.0f;
	if (par.DoubleImSize){
		//printf("... compute_sift_keypoints :: applying zoom\n");
		image = alloc_image<float>(2*width, 2*height);
		apply_zoom(input,image.data,2.0,par.order,width,height);
		octSize *= 0.5f;
	} else
        image = alloc_image( make_image(input,width,height) );

    /// Apply initial smoothing to input image to raise its smoothing to par.InitSigma.  
    /// We assume image from camera has smoothing of sigma = 0.5, which becomes sigma = 1.0 if image has been doubled. 
    /// increase = sqrt(Init^2 - Current^2)
    float curSigma=0.5f;
    if (par.DoubleImSize) curSigma *= 2.0f;
    if (par.InitSigma > curSigma ) {
		if (DEBUG) printf("Convolving initial image to achieve std: %f \n", par.InitSigma);
		float sigma = (float)sqrt(par.InitSigma*par.InitSigma -
                                  curSigma*curSigma);
		gaussian_convolution(image.data, image.data, image.w, image.h, sigma);
	}

	/// Convolve by par.InitSigma at each step inside OctaveKeypoints by steps of 
	/// Subsample of factor 2 while reasonable image size 

	/// Keep reducing image by factors of 2 until one dimension is
	/// smaller than minimum size at which a feature could be detected.
	int 	minsize = 2 * par.BorderDist + 2;
	//printf("... compute_sift_keypoints :: maximum number of scales : %d\n", par.OctaveMax);

    for(int i=(par.DoubleImSize?-1:0);
        i<=par.OctaveMax && image.w>minsize && image.h>minsize;
        i++) {
		if(DEBUG) printf("Calling OctaveKeypoints \n");
		OctaveKeypoints(image, octSize, keys,par);

		// image is blurred inside OctaveKeypoints and therefore can be sampled
		flimage aux = alloc_image<float>(image.w/2, image.h/2);
		if(DEBUG) printf("Sampling initial image \n");
		sample(image.data, aux.data, 2.0f, image.w, image.h);
        free(image.data);
		image = aux;

		octSize *= 2.0f;
	}

    free(image.data);
/*	printf("sift::  %d keypoints\n", keys.size());
	printf("sift::  plus non correctly localized: %d \n", 	par.noncorrectlylocalized);*/

}

static void permute(flimage im[3])
{
    float* tmp = im[0].data;
    im[0].data = im[1].data;
    im[1].data = im[2].data;
    im[2].data = tmp;
}

/////////////////////////////////////////////////
/// EXTREMA DETECTION IN ONE SCALE-SPACE OCTAVE:
/////////////////////////////////////////////////

/// par.Scales determine how many steps we perform to pass from one scale to the next one:  sigI --> 2*sigI
/// At each step we pass from    sigma_0 --> sigma0 * (1 + R)
/// At the last step  sigI * (1 + R)^par.Scales = 2 * sigI
/// (1+R) = 2 ^(1 / par.Scales)  it is called sigmaRatio

/// It seems that blur[par.Scales+1] is compared in two succesive iterations
void OctaveKeypoints(flimage & image, float octSize, keypointslist& keys,siftPar &par)
{
	flimage* blur = new flimage[3];
	flimage* dogs = new flimage[3];
 	float sigmaRatio = (float) pow(2.0, 1.0 / (double) par.Scales);

	blur[0] = image; // First level is input to this routine
	float prevSigma = par.InitSigma; // Input image has par.InitSigma smoothing

    for(int i=1; i < 3; i++)
        blur[i] = alloc_image<float>(image.w, image.h);
    for(int i=0; i < 3; i++)
        dogs[i] = alloc_image<float>(image.w, image.h);
    
    int iDog = 0;
    for(int i=1; i < par.Scales+3; i++) {
		if (DEBUG) printf("Convolving scale: %d \n", i);
        if(i > 2)
            permute(blur);
        int j = (i==1)? 1: 2;
		float sigma = prevSigma*sqrt(sigmaRatio*sigmaRatio-1.0f);
		gaussian_convolution(blur[j-1].data, blur[j].data, blur[j].w, blur[j].h,
                             sigma);
		prevSigma *= sigmaRatio;
        
		/// dogs[i] = dogs[i] - blur[i+1]
		combine(blur[j-1].data,1.0f, blur[j].data,-1.0f,
                dogs[iDog].data, dogs[iDog].w*dogs[iDog].h); 
        if(iDog < 2)
            ++iDog;
        else {
            if (DEBUG) printf("************************scale: %d\n", i-2);
            FindMaxMin(dogs, blur[0], i-2, octSize, keys,par);
            permute(dogs);
        }
    }

    std::swap(image.data,blur[0].data);

    for(int i=1; i < 3; i++)
        free(blur[i].data);
    for(int i=0; i < 3; i++)
        free(dogs[i].data);
	delete [] blur;
	delete [] dogs;
}
	

/////////////////////////////////////////////////
///Find the local maxima and minima of the DOG images in scale space.  Return the keypoints for these locations.
/////////////////////////////////////////////////

/// For each point at each scale we decide if it is a local maxima:
/// - |dogs(x,y,s)| > 0.8 * par.PeakThresh
/// - Local maximum or minimum in s-1,s,s+1
/// - NotonEdge:  ratio of the two principle curvatures of the DOG function at this point be below a threshold.


/// blur[par.Scales+1] is not used in order to look for extrema
/// while these could be computed using avalaible blur and dogs
void FindMaxMin(const flimage* dogs,  const flimage& blur, int s,
                float octSize, keypointslist& keys,siftPar &par)
{

	int width = dogs[0].w, height = dogs[0].h;

	/* Create an image map in which locations that have a keypoint are
	marked with value 1.0, to prevent two keypoints being located at
	same position.  This may seem an inefficient data structure, but
	does not add significant overhead.
	*/
	LWImage<bool> map  = alloc_image<bool>(width,height);
	flimage grad = alloc_image<float>(width,height,2);
    grad.planar = false; // Contiguous norm and dir
    for(int i=map.sizeBuffer()-1; i>=0; i--)
        map.data[i]=false;
    for(int i=grad.sizeBuffer()-1; i>=0; i--)
        grad.data[i]=0.0f;
	
    /* For each intermediate image, compute gradient and orientation
    images to be used for keypoint description.  */
    compute_gradient_orientation(blur.data, grad.data, blur.w, blur.h);
	
    /* Only find peaks at least par.BorderDist samples from image border, as
    peaks centered close to the border will lack stability. */
    assert(par.BorderDist >= 2);
    float val;
    int partialcounter = 0;
    for (int r = par.BorderDist; r < height - par.BorderDist; r++) 
        for (int c = par.BorderDist; c < width - par.BorderDist; c++) {
            /* Pixel value at (c,r) position. */
            val = *dogs[1].pixel(c,r);	

            /* DOG magnitude must be above 0.8 * par.PeakThresh threshold
            (precise threshold check will be done once peak
            interpolation is performed).  Then check whether this
            point is a peak in 3x3 region at each level, and is not
            on an elongated edge.
            */
            if (fabs(val) > 0.8 * par.PeakThresh) {
                if(LocalMaxMin(val, dogs[0], r, c) &&
                   LocalMaxMin(val, dogs[1], r, c) &&
                   LocalMaxMin(val, dogs[2], r, c) &&
                   NotOnEdge(dogs[1], r, c, octSize,par)) {
                    partialcounter++;
                    if (DEBUG) printf("%d:  (%d,%d,%d)  val: %f\n",partialcounter, s,r,c,val);
                    InterpKeyPoint(dogs, s, r, c, grad,
                                   map, octSize, keys, 5,par);	
                }
            }
		}
    free(map.data);
    free(grad.data);
}

/* Return TRUE iff val is a local maximum (positive value) or
   minimum (negative value) compared to the 3x3 neighbourhood that
   is centered at (row,col).
*/
bool LocalMaxMin(float val, const flimage& dog, int y0, int x0)
{
	// For efficiency, use separate cases for maxima or minima, and
	// return as soon as possible
	if (val > 0.0) {
		for (int x = x0 - 1; x <= x0 + 1; x++)
			for (int y = y0  - 1; y <= y0 + 1; y++){
				if (*dog.pixel(x,y) > val) return false;
			}
	} else {
		for (int x = x0 - 1; x <= x0 + 1; x++)
			for (int y = y0  - 1; y <= y0 + 1; y++){
				if (*dog.pixel(x,y) < val) return false;
			}
	}

	return true;
}



/* Returns FALSE if this point on the DOG function lies on an edge.
   This test is done early because it is very efficient and eliminates
   many points.  It requires that the ratio of the two principle
   curvatures of the DOG function at this point be below a threshold.

   Edge threshold is higher on the first scale where SNR is small in 
   order to reduce the number of unstable keypoints. 
*/
int NotOnEdge(const flimage& dog, int r, int c, float octSize,siftPar &par)
{
	/* Compute 2x2 Hessian values from pixel differences. */
	float	H00 = *dog.pixel(c,r-1) - 2.0f * *dog.pixel(c,r) + *dog.pixel(c,r+1), /* AMIR: div by ? */
		H11 = *dog.pixel(c-1,r) - 2.0f * *dog.pixel(c,r) + *dog.pixel(c+1,r),
		H01 = ( (*dog.pixel(c+1,r+1) - *dog.pixel(c-1,r+1)) - (*dog.pixel(c+1,r-1) - *dog.pixel(c-1,r-1)) ) / 4.0f;

	/* Compute determinant and trace of the Hessian. */
	float	det = H00 * H11 - H01 * H01,	/// Det H = \prod l_i
		trace = H00 + H11;		/// tr H = \sum l_i	

	/// As we do not desire edges but only corners we demand l_max / l_min less than a threshold
	/// In practice if A = k B,     A*B = k B^2
	///				(A + B)^2 = (k+1)^2 * B^2
	///				k B^2 >  t * (k+1)^2 * B^2 sii   k  / (k+1)^2 > t
	/// This is a decreasing function for k > 1 and value 0.3 at k=1.
	/// Setting t = 0.08, means k<=10 
	 
	/* To detect an edge response, we require the ratio of smallest
	   to largest principle curvatures of the DOG function
	   (eigenvalues of the Hessian) to be below a threshold.  For
	   efficiency, we use Harris' idea of requiring the determinant to
	   be above par.EdgeThresh times the squared trace, as for eigenvalues
	   A and B, det = AB, trace = A+B.  So if A = 10B, then det = 10B**2,
	   and trace**2 = (11B)**2 = 121B**2, so par.EdgeThresh = 10/121 =
	   0.08 to require ratio of eigenvalues less than 10.
	*/
	if (octSize <= 1)
		return (det > par.EdgeThresh1 * trace * trace);
	else
		return (det > par.EdgeThresh * trace * trace);

}


/* Create a keypoint at a peak near scale space location (s,r,c), where
   s is scale (index of DOGs image), and (r,c) is (row, col) location.
   Add to the list of keys with any new keys added.
*/
void InterpKeyPoint(
    const flimage* dogs, int s, int r, int c,
	const flimage& grad, LWImage<bool>& map,
	float octSize, keypointslist& keys, int movesRemain,siftPar &par)
{
	
	/* Fit quadratic to determine offset and peak value. */
  std::vector<float> offset(3);
	float peakval = FitQuadratic(offset, dogs, r, c);
	if (DEBUG) printf("peakval: %f, of[0]: %f  of[1]: %f  of[2]: %f\n", peakval, offset[0], offset[1], offset[2]);

	/* Move to an adjacent (row,col) location if quadratic interpolation
	   is larger than 0.6 units in some direction (we use 0.6 instead of
	   0.5 to avoid jumping back and forth near boundary).  We do not
	   perform move to adjacent scales, as it is seldom useful and we
	   do not have easy access to adjacent scale structures.  The
	   movesRemain counter allows only a fixed number of moves to
	   prevent possibility of infinite loops.
	*/
	int newr = r, newc = c;
	if (offset[1] > 0.6 && r < dogs[0].h - 3)
		newr++;
	else if (offset[1] < -0.6 && r > 3)
		newr--;

	if (offset[2] > 0.6 && c < dogs[0].w - 3)
		newc++;
	else if (offset[2] < -0.6 && c > 3)
		newc--;

	if (movesRemain > 0  &&  (newr != r || newc != c)) {
		InterpKeyPoint(dogs, s, newr, newc, grad, map,
                       octSize, keys,movesRemain - 1,par);
		return;
	}

	/* Do not create a keypoint if interpolation still remains far
	   outside expected limits, or if magnitude of peak value is below
	   threshold (i.e., contrast is too low). */
	if (fabs(offset[0]) > 1.5 || fabs(offset[1]) > 1.5 ||
		fabs(offset[2]) > 1.5 || fabs(peakval) < par.PeakThresh)		
		{
			if (DEBUG) printf("Point not well localized by FitQuadratic\n"); 	
			par.noncorrectlylocalized++;
			return;
		}
	
	/* Check that no keypoint has been created at this location (to avoid
	   duplicates).  Otherwise, mark this map location.
	*/
	if (*map.pixel(c,r)) return;
	*map.pixel(c,r) = true;

	/* The scale relative to this octave is given by octScale.  The scale
	   units are in terms of sigma for the smallest of the Gaussians in the
	   DOG used to identify that scale.
	*/
	float octScale = par.InitSigma * pow(2.0f, (s + offset[0]) / (float) par.Scales);

	/// always use histogram of orientations
	//if (UseHistogramOri)
    AssignOriHist(grad, octSize, octScale,
                  r + offset[1], c + offset[2], keys, par);
	//else
	//	AssignOriAvg(
	//		grad, ori, octSize, octScale,
	//		r + offset[1], c + offset[2], keys);
}



/* Apply the method developed by Matthew Brown (see BMVC 02 paper) to
   fit a 3D quadratic function through the DOG function values around
   the location (s,r,c), i.e., (scale,row,col), at which a peak has
   been detected.  Return the interpolated peak position as a vector
   in "offset", which gives offset from position (s,r,c).  The
   returned value is the interpolated DOG magnitude at this peak.
*/
float FitQuadratic(std::vector<float>& offset, const flimage* dogs, int r, int c)
{
  typedef libNumerics::flnum flnum;
  libNumerics::matrix<flnum> H(3,3);

	/* Select the dog images at peak scale, dog1, as well as the scale
	   below, dog0, and scale above, dog2. */
	const flimage& dog0 = dogs[0];
	const flimage& dog1 = dogs[1];
	const flimage& dog2 = dogs[2];

	/* Fill in the values of the gradient from pixel differences. */
    libNumerics::vector<flnum> g(3);
	g(0) = (*dog2.pixel(c,r)   - *dog0.pixel(c,r))   / 2.0f;
	g(1) = (*dog1.pixel(c,r+1) - *dog1.pixel(c,r-1)) / 2.0f;
	g(2) = (*dog1.pixel(c+1,r) - *dog1.pixel(c-1,r)) / 2.0f;

	/* Fill in the values of the Hessian from pixel differences. */
	H(0,0) = *dog0.pixel(c,r)   - 2.0f * *dog1.pixel(c,r) + *dog2.pixel(c,r);
	H(1,1) = *dog1.pixel(c,r-1) - 2.0f * *dog1.pixel(c,r) + *dog1.pixel(c,r+1);
	H(2,2) = *dog1.pixel(c-1,r) - 2.0f * *dog1.pixel(c,r) + *dog1.pixel(c+1,r);
	H(0,1) = H(1,0) = ( (*dog2.pixel(c,r+1) - *dog2.pixel(c,r-1)) -
                        (*dog0.pixel(c,r+1) - *dog0.pixel(c,r-1)) ) /4.0f;
	H(0,2) = H(2,0) = ( (*dog2.pixel(c+1,r) - *dog2.pixel(c-1,r)) -
                        (*dog0.pixel(c+1,r) - *dog0.pixel(c-1,r)) ) /4.0f;
	H(1,2) = H(2,1) = ( (*dog1.pixel(c+1,r+1) - *dog1.pixel(c-1,r+1)) -
                        (*dog1.pixel(c+1,r-1) - *dog1.pixel(c-1,r-1)) ) /4.0f;

	/* Solve the 3x3 linear sytem, Hx = -g. Result, x, gives peak offset.
	   Note that SolveLinearSystem destroys contents of H. */
    libNumerics::vector<flnum> sol(-g);
    libNumerics::solveLU(H, sol);
    for(int i=0; i < 3; i++)
      offset[i] = static_cast<float>(sol(i));

	/* Also return value of DOG at peak location using initial value plus
	   0.5 times linear interpolation with gradient to peak position
	   (this is correct for a quadratic approximation).
	*/
	return static_cast<float>(*dog1.pixel(c,r) + 0.5f*dot(sol,g));
}

/// - Compute histogram of orientation in a neighborhood weighted by gradient and distance to center
/// - Look for local (3-neighborhood) maximum with valuer larger or equal than par.OriHistThresh * maxval


/* Assign an orientation to this keypoint.  This is done by creating a
   Gaussian weighted histogram of the gradient directions in the
   region.  The histogram is smoothed and the largest peak selected.
   The results are in the range of -PI to PI.
*/
void AssignOriHist(
	const flimage& grad, float octSize,
	float octScale, float octRow, float octCol,keypointslist& keys,siftPar &par)
{
	int	bin, prev, next;
	float* hist = new float[par.OriBins];
	float	distsq, dif, weight, angle, interp;
	float radius2, sigma2;		

	int	row = (int) (octRow+0.5),
		col = (int) (octCol+0.5),
		rows = grad.h,
		cols = grad.w;

	for (int i = 0; i < par.OriBins; i++) hist[i] = 0.0;

	/* Look at pixels within 3 sigma around the point and sum their
	  Gaussian weighted gradient magnitudes into the histogram. */
	float	sigma = par.OriSigma * octScale;
	int	radius = (int) (sigma * 3.0);
	int rmin = std::max(0,row-radius);
	int cmin = std::max(0,col-radius);
	int rmax = std::min(row+radius,rows-2);
	int cmax = std::min(col+radius,cols-2);
	radius2 = (float)(radius * radius);
	sigma2 = 2.0f*sigma*sigma;

	for (int r = rmin; r <= rmax; r++) {
		for (int c = cmin; c <= cmax; c++) {
            dif = (r - octRow);	distsq = dif*dif;
            dif = (c - octCol);	distsq += dif*dif;
			const float* g=grad.pixel(c,r);

            if (g[0] > 0.0  &&  distsq < radius2 + 0.5) {
                weight = exp(- distsq / sigma2);
					
                /* Ori is in range of -PI to PI. */
                bin = (int) (par.OriBins * (g[1] + M_PI + 0.001) / (2.0 * M_PI));
                assert(bin >= 0 && bin <= par.OriBins);
                bin = std::min(bin, par.OriBins - 1);
                hist[bin] += weight * g[0];
            }
		}
	}


	/* Apply smoothing 6 times for accurate Gaussian approximation. */
	for (int i = 0; i < 6; i++)
		SmoothHistogram(hist, par.OriBins);

	/* Find maximum value in histogram. */
	float maxval = 0.0;
	for (int i = 0; i < par.OriBins; i++) 
		if (hist[i] > maxval) maxval = hist[i];

	/* Look for each local peak in histogram.  If value is within
	  par.OriHistThresh of maximum value, then generate a keypoint. */
	for (int i = 0; i < par.OriBins; i++) {
		prev = (i == 0 ? par.OriBins - 1 : i - 1);
		next = (i == par.OriBins - 1 ? 0 : i + 1);

		if (	hist[i] > hist[prev]  &&  hist[i] > hist[next]  &&
			hist[i] >= par.OriHistThresh * maxval ) {
	
			/* Use parabolic fit to interpolate peak location from 3 samples.
			  Set angle in range -PI to PI. */
			interp = InterpPeak(hist[prev], hist[i], hist[next]);
			angle = 2.0f * M_PI * (i + 0.5f + interp) / (float)par.OriBins - M_PI;
			assert(angle >= -M_PI  &&  angle <= M_PI);
		
			if (DEBUG) printf("angle selected: %f \t location: (%f,%f)\n", angle, octRow, octCol);

			/* Create a keypoint with this orientation. */
			MakeKeypoint(
				grad, octSize, octScale,
				octRow, octCol, angle, keys,par);
		}

	}
	delete [] hist;
}



/* Smooth a histogram by using a [1/3 1/3 1/3] kernel.  Assume the histogram
   is connected in a circular buffer.
*/
void SmoothHistogram(float* hist, int bins)
{
	float prev, temp;

	prev = hist[bins - 1];
	for (int i = 0; i < bins; i++) {
		temp = hist[i];
		hist[i] = ( prev + hist[i] + hist[(i + 1 == bins) ? 0 : i + 1] ) / 3.0f;
		prev = temp;
	}
}


/* Return a number in the range [-0.5, 0.5] that represents the
   location of the peak of a parabola passing through the 3 evenly
   spaced samples.  The center value is assumed to be greater than or
   equal to the other values if positive, or less than if negative.
*/
float InterpPeak(float a, float b, float c)
{
	if (b < 0.0) {
		a = -a; b = -b; c = -c;
	}
	assert(b >= a  &&  b >= c);
	return 0.5f * (a - c) / (a - 2.0f * b + c);
}





/* Joan Pau: Add a new keypoint to a vector of keypoints
   Create a new keypoint and return list of keypoints with new one added.
*/
void MakeKeypoint(
	const flimage& grad, float octSize, float octScale,
	float octRow, float octCol, float angle, keypointslist& keys,siftPar &par)
{
	keypoint newkeypoint;
	newkeypoint.x = octSize * octCol;	/*x coordinate */
	newkeypoint.y = octSize * octRow;	/*y coordinate */
	newkeypoint.scale = octSize * octScale;	/* scale */
	newkeypoint.angle = angle;		/* orientation */
	MakeKeypointSample(newkeypoint,grad,octScale,octRow,octCol,par);
	keys.push_back(newkeypoint);
}



/* Use the parameters of this keypoint to sample the gradient images
     at a set of locations within a circular region around the keypoint.
     The (scale,row,col) values are relative to current octave sampling.
     The resulting vector is stored in the key.
*/
void MakeKeypointSample(
	keypoint& key, const flimage& grad,
	float scale, float row, float col,siftPar &par)
{
	/* Produce sample vector. */
	KeySampleVec(key, grad, scale, row, col,par);


	/* Normalize vector.  This should provide illumination invariance
	for planar lambertian surfaces (except for saturation effects).
	Normalization also improves nearest-neighbor metric by
	increasing relative distance for vectors with few features.
	It is also useful to implement a distance threshold and to
	allow conversion to integer format.
	*/
	NormalizeVec(key.vec);

	/* Now that normalization has been done, threshold elements of
	index vector to decrease emphasis on large gradient magnitudes.
	Admittedly, the large magnitude values will have affected the
	normalization, and therefore the threshold, so this is of
	limited value.
	*/
	bool changed = false;
	for (int i = 0; i < VecLength; i++)
		if (key.vec[i] > par.MaxIndexVal) {
			key.vec[i] = par.MaxIndexVal;
			changed = true;
		}

	if (changed) NormalizeVec(key.vec);

	/* Convert float vector to integer. Assume largest value in normalized
	vector is likely to be less than 0.5. */
	/// QUESTION: why is the vector quantized to integer
	int intval;
	for (int i = 0; i < VecLength; i++) {
		intval =  static_cast<int>(512.0 * key.vec[i]);
		key.vec[i] = static_cast<float>(std::min(255, intval));
	}
}

/* Normalize length of vec to 1.0.
*/
void NormalizeVec(float* vec)
{
	float val, fac;
	
	float sqlen = 0.0;
	for (int i = 0; i < VecLength; i++) {
		val = vec[i];
		sqlen += val * val;
	}
	fac = 1.0f / sqrt(sqlen);

	for (int i = 0; i < VecLength; i++)
		vec[i] *= fac;
}


/* Create a 3D index array into which gradient values are accumulated.
   After filling array, copy values back into vec.
*/
void KeySampleVec(
	keypoint& key, const flimage& grad,
	float scale, float row, float col,siftPar &par)
{
	
	float index[IndexSize][IndexSize][OriSize];

	/* Initialize index array. */
	for (int i = 0; i < IndexSize; i++)
		for (int j = 0; j < IndexSize; j++)
			for (int k = 0; k < OriSize; k++)
				index[i][j][k] = 0.0;


	KeySample(index, key, grad, scale, row, col, par);


	/* Unwrap the 3D index values into 1D vec. */
	int v = 0;
	for (int i = 0; i < IndexSize; i++)
		for (int j = 0; j < IndexSize; j++)
			for (int k = 0; k < OriSize; k++)
				key.vec[v++] = index[i][j][k];
}



/* Add features to vec obtained from sampling the grad and ori images
   for a particular scale.  Location of key is (scale,row,col) with respect
   to images at this scale.  We examine each pixel within a circular
   region containing the keypoint, and distribute the gradient for that
   pixel into the appropriate bins of the index array.
*/
void KeySample(
	float index[IndexSize][IndexSize][OriSize], keypoint& key,
	const flimage& grad, float scale, float row, float col,siftPar &par)
{
	float rpos, cpos, rx, cx;

	int	irow = (int) (row + 0.5),
		icol = (int) (col + 0.5);
	float	sine   = (float) sin(key.angle),
		cosine = (float) cos(key.angle);

	/* The spacing of index samples in terms of pixels at this scale. */
	float	spacing = scale * par.MagFactor;

	/* Radius of index sample region must extend to diagonal corner of
	index patch plus half sample for interpolation. */
	float	radius = 1.414f * spacing * (IndexSize + 1) / 2.0f;
	int	iradius = (int) (radius + 0.5);

	/* Examine all points from the gradient image that could lie within the
	index square. */
	for (int i = -iradius; i <= iradius; i++) {
		for (int j = -iradius; j <= iradius; j++) {

			/* Rotate sample offset to make it relative to key orientation.
			 Uses (row,col) instead of (x,y) coords.  Also, make subpixel
			 correction as later image offset must be an integer.  Divide
			 by spacing to put in index units.
			*/

			 /* Guoshen Yu, inverse the rotation */
			 rpos = ((cosine * i - sine * j) - (row - irow)) / spacing;
			 cpos = ((sine * i + cosine * j) - (col - icol)) / spacing;

			 /*
			 rpos = ((cosine * i + sine * j) - (row - irow)) / spacing;
			 cpos = ((- sine * i + cosine * j) - (col - icol)) / spacing;*/

			 /* Compute location of sample in terms of real-valued index array
			 coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
			 weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. */
			 rx = rpos + IndexSize / 2.0f - 0.5f;
			 cx = cpos + IndexSize / 2.0f - 0.5f;
	
			/* Test whether this sample falls within boundary of index patch. */
			if (	rx > -1.0 && rx < (float) IndexSize  &&
				cx > -1.0 && cx < (float) IndexSize )
				AddSample(
					index, key, grad,
					irow + i, icol + j, rpos, cpos, rx, cx,par);
		}
	}
}


/* Given a sample from the image gradient, place it in the index array.
*/
void AddSample(
	float index[IndexSize][IndexSize][OriSize], keypoint& key,
	const flimage& grad,
	int r, int c, float rpos, float cpos, float rx, float cx,siftPar &par)
{
	/* Clip at image boundaries. */
	if (r < 0  ||  r >= grad.h  ||  c < 0  ||  c >= grad.w)
		return;

    const float* g = grad.pixel(c,r);
	/* Compute Gaussian weight for sample, as function of radial distance
	   from center.  Sigma is relative to half-width of index. */
	float	sigma  = par.IndexSigma * 0.5f * IndexSize,
		weight = exp(- (rpos * rpos + cpos * cpos) / (2.0f * sigma * sigma)),
		mag    = weight * g[0];

	/* Subtract keypoint orientation to give ori relative to keypoint. */
	float	ori = g[1] - key.angle;

	/* Put orientation in range [0, 2*PI].  If sign of gradient is to
	   be ignored, then put in range [0, PI]. */
	if (par.IgnoreGradSign) {
		while (ori > M_PI ) ori -= M_PI;
		while (ori < 0.0f) ori += M_PI;
	} else {
		while (ori > 2.0f*M_PI) ori -= 2.0f*M_PI;
		while (ori < 0.0f   ) ori += 2.0f*M_PI;
	}
	PlaceInIndex(index, mag, ori, rx, cx,par);
}


/* Increment the appropriate locations in the index to incorporate
   this image sample.  The location of the sample in the index is (rx,cx).
*/
void PlaceInIndex(
	float index[IndexSize][IndexSize][OriSize],
	float mag, float ori, float rx, float cx,siftPar &par)
{
	int	orr, rindex, cindex, oindex;
	float	rweight, cweight, oweight;
	float  *ivec;

	float	oval = OriSize * ori / (par.IgnoreGradSign ? M_PI : 2.0f*M_PI);

	int	ri = static_cast<int>( (rx >= 0.0) ? rx : rx - 1.0),	/* Round down to next integer. */
		  ci = static_cast<int>( (cx >= 0.0) ? cx : cx - 1.0),
		  oi = static_cast<int>( (oval >= 0.0) ? oval : oval - 1.0);
	float	rfrac = rx - ri,			/* Fractional part of location. */
		    cfrac = cx - ci,
		    ofrac = oval - oi;
	assert(
		ri >= -1  &&  ri < IndexSize  &&
		oi >=  0  &&  oi <= OriSize  &&
		rfrac >= 0.0  &&  rfrac <= 1.0);

	/* Put appropriate fraction in each of 8 buckets around this point
		in the (row,col,ori) dimensions.  This loop is written for
		efficiency, as it is the inner loop of key sampling. */
	for (int r = 0; r < 2; r++) {
		rindex = ri + r;
		if (rindex >=0 && rindex < IndexSize) {
			rweight = mag * ((r == 0) ? 1.0f - rfrac : rfrac);

			for (int c = 0; c < 2; c++) {
				cindex = ci + c;
				if (cindex >=0 && cindex < IndexSize) {
					cweight = rweight * ((c == 0) ? 1.0f - cfrac : cfrac);
					ivec = index[rindex][cindex];
					for (orr = 0; orr < 2; orr++) {
						oindex = oi + orr;
						if (oindex >= OriSize)  /* Orientation wraps around at PI. */
							oindex = 0;
						oweight = cweight * ((orr == 0) ? 1.0f - ofrac : ofrac);
						ivec[oindex] += oweight;
					}
				}
			}
		}
	}
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// SIFT keypoint matching 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float DistSquared(keypoint &k1,keypoint &k2, float tdist, siftPar &par)
{
	
	float dif;
	float distsq = 0.0;

	if (fabsf(k1.x - k2.x) > par.MatchXradius || fabsf(k1.y - k2.y) > par.MatchYradius) return tdist;

	float *ik1 = k1.vec;
	float *ik2 = k2.vec;

	for (int i = 0; i < VecLength && distsq <= tdist; i++) {
		dif = ik1[i] - ik2[i];
		distsq += dif * dif;
	}

	return distsq;
}


/* This searches through the keypoints in klist for the two closest
   matches to key.  It returns the ratio of the distance to key of the
   closest and next to closest keypoints in klist, while bestindex is the index
   of the closest keypoint.
*/
float CheckForMatch(
keypoint& key, keypointslist& klist, keypointslist::size_type& min,siftPar &par)
{
	float	dsq, distsq1, distsq2;
	distsq1 = distsq2 = 1000000000000.0f;

	for (keypointslist::size_type j=0; j< klist.size(); j++){
	
		dsq = DistSquared(key, klist[j], distsq2,par);
		
		if (dsq < distsq1) {
			distsq2 = distsq1;
			distsq1 = dsq;
			min = j;
		} else if (dsq < distsq2)
			distsq2 = dsq;
	}

	return distsq1/distsq2 ;
}


void compute_sift_matches(
	 keypointslist& keys1,  keypointslist& keys2,
	matchingslist& matchings,siftPar &par)
{
	float sqminratio = par.MatchRatio * par.MatchRatio;
		
	for (keypointslist::size_type i=0; i< keys1.size(); i++) {
        keypointslist::size_type imatch=0;
		float sqratio = CheckForMatch(keys1[i], keys2, imatch,par);
		if(sqratio < sqminratio)
			matchings.push_back( Match(keys1[i].x,keys1[i].y,
                                       keys2[imatch].x, keys2[imatch].y) );
	}
}
