// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

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


#ifndef _CLIBSIFT_H_
#define _CLIBSIFT_H_



///////////// Description
/// For each octave:
///    	- Divide in par.Scales scales
/// 	- Convolve and compute differences of convolved scales
///	- Look for a 3x3 multiscale extrema and contraste enough and with no predominant direction (no 1d edge)

/// For each extrema
///	- Compute orientation histogram in neighborhood.
///	- Generate a keypoint for each mode with this orientation


/// For each keypoint
///	 - Create vector 



///////////// Possible differences with MW
/// Gaussian convolution


#include <stdlib.h>
#include <assert.h>

#include "numerics1.h"
#include "library.h"
#include "filter.h"
#include "domain.h"
#include "splines.h"
#include "flimage.h"


#include <vector>

// BASIC STRUCTURES:

// Keypoints:
#define OriSize  8
#define IndexSize  4
#define VecLength  IndexSize * IndexSize * OriSize


/* Keypoint structure:
	position:	x,y
	scale:		s
	orientation:	angle
	descriptor:	array of gradient orientation histograms in a neighbors */
struct keypoint {
	float	x,y,
		scale,
		angle;
	float	vec[VecLength];
};


/* Keypoint structure:
 position:	x,y
 scale:		s
 orientation:	angle
 descriptor:	array of gradient orientation histograms in a neighbors */
struct keypoint_char {
	float	x,y,
	scale,
	angle;
	unsigned char	vec[VecLength];
};

/* Keypoint structure:
 position:	x,y
 scale:		s
 orientation:	angle
 descriptor:	array of gradient orientation histograms in a neighbors */
struct keypoint_short {
	float	x,y,
	scale,
	angle;
	unsigned short	vec[VecLength];
};

/* Keypoint structure:
 position:	x,y
 scale:		s
 orientation:	angle
 descriptor:	array of gradient orientation histograms in a neighbors */
struct keypoint_int {
	float	x,y,
	scale,
	angle;
	unsigned int	vec[VecLength];
};

/* List of keypoints: just use the standard class vector: */
typedef std::vector<keypoint> keypointslist;

/* List of keypoints: just use the standard class vector: */
typedef std::vector<keypoint_char> keypointslist_char;
typedef std::vector<keypoint_short> keypointslist_short;
typedef std::vector<keypoint_int> keypointslist_int;



/* Matching: just use the standard class pair: */
typedef std::pair<keypoint,keypoint> matching;


/* List of matchings: just use the standard class vector: */
typedef std::vector<matching> matchingslist;


struct siftPar
{

int OctaveMax;

int DoubleImSize;

int order;


/* InitSigma gives the amount of smoothing applied to the image at the
   first level of each octave.  In effect, this determines the sampling
   needed in the image domain relative to amount of smoothing.  Good
   values determined experimentally are in the range 1.2 to 1.8.
*/
float  InitSigma /*= 1.6*/;
 

/* Peaks in the DOG function must be at least BorderDist samples away
   from the image border, at whatever sampling is used for that scale.
   Keypoints close to the border (BorderDist < about 15) will have part
   of the descriptor landing outside the image, which is approximated by
   having the closest image pixel replicated.  However, to perform as much
   matching as possible close to the edge, use BorderDist of 4.
*/
int BorderDist /*= 5*/;


/* Scales gives the number of discrete smoothing levels within each octave.
   For example, Scales = 2 implies dividing octave into 2 intervals, so
   smoothing for each scale sample is sqrt(2) more than previous level.
   Value of 2 works well, but higher values find somewhat more keypoints.
*/

int Scales /*= 3*/;


/// Decreasing PeakThresh allows more non contrasted keypoints
/* Magnitude of difference-of-Gaussian value at a keypoint must be above
   this threshold.  This avoids considering points with very low contrast
   that are dominated by noise.  It is divided by Scales because more
   closely spaced scale samples produce smaller DOG values.  A value of
   0.08 considers only the most stable keypoints, but applications may
   wish to use lower values such as 0.02 to find keypoints from low-contast
   regions.
*/

//#define  PeakThreshInit  255*0.04 
//#define  PeakThresh      PeakThreshInit / Scales
float PeakThresh  /*255.0 * 0.04 / 3.0*/;

/// Decreasing EdgeThresh allows more edge points
/* This threshold eliminates responses at edges.  A value of 0.08 means
   that the ratio of the largest to smallest eigenvalues (principle
   curvatures) is below 10.  A value of 0.14 means ratio is less than 5.
   A value of 0.0 does not eliminate any responses.
   Threshold at first octave is different.
*/
float  EdgeThresh  /*0.06*/;
float  EdgeThresh1 /*0.08*/;


/* OriBins gives the number of bins in the histogram (36 gives 10
   degree spacing of bins).
*/
int OriBins  /*36*/;


/* Size of Gaussian used to select orientations as multiple of scale
     of smaller Gaussian in DOG function used to find keypoint.
     Best values: 1.0 for UseHistogramOri = FALSE; 1.5 for TRUE.
*/
float OriSigma  /*1.5*/;


/// Look for local (3-neighborhood) maximum with valuer larger or equal than OriHistThresh * maxval
///  Setting one returns a single peak
/* All local peaks in the orientation histogram are used to generate
   keypoints, as long as the local peak is within OriHistThresh of
   the maximum peak.  A value of 1.0 only selects a single orientation
   at each location.
*/
float OriHistThresh  /*0.8*/;


/// Feature vector is normalized to has euclidean norm 1.
/// This threshold avoid the excessive concentration of information on single peaks
/* Index values are thresholded at this value so that regions with
   high gradients do not need to match precisely in magnitude.
   Best value should be determined experimentally.  Value of 1.0
   has no effect.  Value of 0.2 is significantly better.
*/
float  MaxIndexVal  /*0.2*/;


/* This constant specifies how large a region is covered by each index
   vector bin.  It gives the spacing of index samples in terms of
   pixels at this scale (which is then multiplied by the scale of a
   keypoint).  It should be set experimentally to as small a value as
   possible to keep features local (good values are in range 3 to 5).
*/
int  MagFactor   /*3*/;


/* Width of Gaussian weighting window for index vector values.  It is
   given relative to half-width of index, so value of 1.0 means that
   weight has fallen to about half near corners of index patch.  A
   value of 1.0 works slightly better than large values (which are
   equivalent to not using weighting).  Value of 0.5 is considerably
   worse.
*/
float   IndexSigma  /*1.0*/;

/* If this is TRUE, then treat gradients with opposite signs as being
   the same.  In theory, this could create more illumination invariance,
   but generally harms performance in practice.
*/
int  IgnoreGradSign  /*0*/;



float MatchRatio  /*0.6*/;

/*
   In order to constrain the research zone for matches.
   Useful for example when looking only at epipolar lines
*/

float MatchXradius /*= 1000000.0f*/;
float MatchYradius /*= 1000000.0f*/;

int noncorrectlylocalized;

};

//////////////////////////////////////////////////////////
/// SIFT 
//////////////////////////////////////////////////////////

void default_sift_parameters(siftPar &par);

void compute_sift_keypoints(float *input,  keypointslist& keypoints,int width, int height, siftPar &par);




// MATCHING DETECTION FUNCTION:
void compute_sift_matches( keypointslist& keys1, keypointslist& keys2, matchingslist& matchings, siftPar &par);

#endif // _LIBSIFT_H_



