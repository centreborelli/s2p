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
//*----------------------------- demo_ASIFT  --------------------------------*/
// Detect corresponding points in two images with the ASIFT method. 

// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
// 
// Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image 
//            Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009. 
// Reference: ASIFT online demo (You can try ASIFT with your own images online.) 
//			  http://www.ipol.im/pub/algo/my_affine_sift/
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "demo_lib_sift.h"
#include "io_png/io_png.h"

#include "library.h"
#include "frot.h"
#include "fproj.h"
#include "compute_asift_keypoints.h"
#include "compute_asift_matches.h"

# define IM_X 800
# define IM_Y 600

int main(int argc, char **argv)
{			
	
    if ((argc != 8) && (argc != 9)) {
        std::cerr << " ******************************************************************************* " << std::endl
				  << " ***************************  ASIFT image matching  **************************** " << std::endl
				  << " ******************************************************************************* " << std::endl
				  << "Usage: " << argv[0] << " imgIn1.png imgIn2.png imgOutVert.png imgOutHori.png " << std::endl
										  << "           matchings.txt keys1.txt keys2.txt [Resize option: 0/1] " << std::endl
									      << "- imgIn1.png, imgIn2.png: input images (in PNG format). " << std::endl
										  << "- imgOutVert.png, imgOutHori.png: output images (vertical/horizontal concatenated, " << std::endl
				                          << "  in PNG format.) The detected matchings are connected by write lines." << std::endl
										  << "- matchings.txt: coordinates of matched points (col1, row1, col2, row2). " << std::endl
										  << "- keys1.txt keys2.txt: ASIFT keypoints of the two images." << std::endl
										  << "- [optional 0/1]. 1: input images resize to 800x600 (default). 0: no resize. " << std::endl 
   				  << " ******************************************************************************* " << std::endl
				  << " *********************  Jean-Michel Morel, Guoshen Yu, 2010 ******************** " << std::endl
				  << " ******************************************************************************* " << std::endl;
        return 1;
    }
	
	//////////////////////////////////////////////// Input
	// Read image1
    float * iarr1;
    size_t w1, h1;
    if (NULL == (iarr1 = read_png_f32_gray(argv[1], &w1, &h1))) {
        std::cerr << "Unable to load image file " << argv[1] << std::endl;
        return 1;
    }
    std::vector<float> ipixels1(iarr1, iarr1 + w1 * h1);
	free(iarr1); /*memcheck*/
	
	// Read image2
    float * iarr2;
    size_t w2, h2;
    if (NULL == (iarr2 = read_png_f32_gray(argv[2], &w2, &h2))) {
        std::cerr << "Unable to load image file " << argv[2] << std::endl;
        return 1;
    }
    std::vector<float> ipixels2(iarr2, iarr2 + w2 * h2);
	free(iarr2); /*memcheck*/	
	
	///// Resize the images to area wS*hW in remaining the apsect-ratio	
	///// Resize if the resize flag is not set or if the flag is set unequal to 0
	float wS = IM_X;
	float hS = IM_Y;
	
	float zoom1=0, zoom2=0;	
	int wS1=0, hS1=0, wS2=0, hS2=0;
	vector<float> ipixels1_zoom, ipixels2_zoom;	
		
	int flag_resize = 1;
	if (argc == 9)
	{	
		flag_resize = atoi(argv[8]);
	}
	
	if ((argc == 8) || (flag_resize != 0))
	{
		cout << "WARNING: The input images are resized to " << wS << "x" << hS << " for ASIFT. " << endl 
		<< "         But the results will be normalized to the original image size." << endl << endl;
		
		float InitSigma_aa = 1.6;
		
		float fproj_p, fproj_bg;
		char fproj_i;
		float *fproj_x4, *fproj_y4;
		int fproj_o;

		fproj_o = 3;
		fproj_p = 0;
		fproj_i = 0;
		fproj_bg = 0;
		fproj_x4 = 0;
		fproj_y4 = 0;
				
		float areaS = wS * hS;

		// Resize image 1 
		float area1 = w1 * h1;
		zoom1 = sqrt(area1/areaS);
		
		wS1 = (int) (w1 / zoom1);
		hS1 = (int) (h1 / zoom1);
		
		int fproj_sx = wS1;
		int fproj_sy = hS1;     
		
		float fproj_x1 = 0;
		float fproj_y1 = 0;
		float fproj_x2 = wS1;
		float fproj_y2 = 0;
		float fproj_x3 = 0;	     
		float fproj_y3 = hS1;
		
		/* Anti-aliasing filtering along vertical direction */
		if ( zoom1 > 1 )
		{
			float sigma_aa = InitSigma_aa * zoom1 / 2;
			GaussianBlur1D(ipixels1,w1,h1,sigma_aa,1);
			GaussianBlur1D(ipixels1,w1,h1,sigma_aa,0);
		}
			
		// simulate a tilt: subsample the image along the vertical axis by a factor of t.
		ipixels1_zoom.resize(wS1*hS1);
		fproj (ipixels1, ipixels1_zoom, w1, h1, &fproj_sx, &fproj_sy, &fproj_bg, &fproj_o, &fproj_p, 
			   &fproj_i , fproj_x1 , fproj_y1 , fproj_x2 , fproj_y2 , fproj_x3 , fproj_y3, fproj_x4, fproj_y4); 
		
		
		// Resize image 2 
		float area2 = w2 * h2;
		zoom2 = sqrt(area2/areaS);
				
		wS2 = (int) (w2 / zoom2);
		hS2 = (int) (h2 / zoom2);
		
		fproj_sx = wS2;
		fproj_sy = hS2;     
		
		fproj_x2 = wS2;
		fproj_y3 = hS2;
		
		/* Anti-aliasing filtering along vertical direction */
		if ( zoom1 > 1 )
		{
			float sigma_aa = InitSigma_aa * zoom2 / 2;
			GaussianBlur1D(ipixels2,w2,h2,sigma_aa,1);
			GaussianBlur1D(ipixels2,w2,h2,sigma_aa,0);
		}
			
		// simulate a tilt: subsample the image along the vertical axis by a factor of t.
		ipixels2_zoom.resize(wS2*hS2);
		fproj (ipixels2, ipixels2_zoom, w2, h2, &fproj_sx, &fproj_sy, &fproj_bg, &fproj_o, &fproj_p, 
			   &fproj_i , fproj_x1 , fproj_y1 , fproj_x2 , fproj_y2 , fproj_x3 , fproj_y3, fproj_x4, fproj_y4); 
	}
	else 
	{
		ipixels1_zoom.resize(w1*h1);	
		ipixels1_zoom = ipixels1;
		wS1 = w1;
		hS1 = h1;
		zoom1 = 1;
		
		ipixels2_zoom.resize(w2*h2);	
		ipixels2_zoom = ipixels2;
		wS2 = w2;
		hS2 = h2;
		zoom2 = 1;
	}

	
	///// Compute ASIFT keypoints
	// number N of tilts to simulate t = 1, \sqrt{2}, (\sqrt{2})^2, ..., {\sqrt{2}}^(N-1)
	int num_of_tilts1 = 7;
	int num_of_tilts2 = 7;
//	int num_of_tilts1 = 1;
//	int num_of_tilts2 = 1;
	int verb = 0;
	// Define the SIFT parameters
	siftPar siftparameters;	
	default_sift_parameters(siftparameters);

	vector< vector< keypointslist > > keys1;		
	vector< vector< keypointslist > > keys2;	
	
	int num_keys1=0, num_keys2=0;
	
	
	cout << "Computing keypoints on the two images..." << endl;
	time_t tstart, tend;	
	tstart = time(0);

	num_keys1 = compute_asift_keypoints(ipixels1_zoom, wS1, hS1, num_of_tilts1, verb, keys1, siftparameters);
	num_keys2 = compute_asift_keypoints(ipixels2_zoom, wS2, hS2, num_of_tilts2, verb, keys2, siftparameters);
	
	tend = time(0);
	cout << "Keypoints computation accomplished in " << difftime(tend, tstart) << " seconds." << endl;
	
	//// Match ASIFT keypoints
	int num_matchings;
	matchingslist matchings;	
	cout << "Matching the keypoints..." << endl;
	tstart = time(0);
	num_matchings = compute_asift_matches(num_of_tilts1, num_of_tilts2, wS1, hS1, wS2, 
										  hS2, verb, keys1, keys2, matchings, siftparameters);
	tend = time(0);
	cout << "Keypoints matching accomplished in " << difftime(tend, tstart) << " seconds." << endl;
	
	///////////////// Output image containing line matches (the two images are concatenated one above the other)
	int band_w = 20; // insert a black band of width band_w between the two images for better visibility
	
	int wo =  MAX(w1,w2);
	int ho = h1+h2+band_w;
	
	float *opixelsASIFT = new float[wo*ho];
	
	for(int j = 0; j < (int) ho; j++)
		for(int i = 0; i < (int) wo; i++)  opixelsASIFT[j*wo+i] = 255;		
	
	/////////////////////////////////////////////////////////////////// Copy both images to output
	for(int j = 0; j < (int) h1; j++)
		for(int i = 0; i < (int) w1; i++)  opixelsASIFT[j*wo+i] = ipixels1[j*w1+i];				
	
	for(int j = 0; j < (int) h2; j++)
		for(int i = 0; i < (int) (int)w2; i++)  opixelsASIFT[(h1 + band_w + j)*wo + i] = ipixels2[j*w2 + i];	
	
	//////////////////////////////////////////////////////////////////// Draw matches
	matchingslist::iterator ptr = matchings.begin();
	for(int i=0; i < (int) matchings.size(); i++, ptr++)
	{		
		draw_line(opixelsASIFT, (int) (zoom1*ptr->first.x), (int) (zoom1*ptr->first.y), 
				  (int) (zoom2*ptr->second.x), (int) (zoom2*ptr->second.y) + h1 + band_w, 255.0f, wo, ho);		
	}
			
	///////////////////////////////////////////////////////////////// Save imgOut	
	write_png_f32(argv[3], opixelsASIFT, wo, ho, 1);
	
	delete[] opixelsASIFT; /*memcheck*/
	
	/////////// Output image containing line matches (the two images are concatenated one aside the other)
	int woH =  w1+w2+band_w;
	int hoH = MAX(h1,h2);
	
	float *opixelsASIFT_H = new float[woH*hoH];
	
	for(int j = 0; j < (int) hoH; j++)
		for(int i = 0; i < (int) woH; i++)  opixelsASIFT_H[j*woH+i] = 255;
	
	/////////////////////////////////////////////////////////////////// Copy both images to output
	for(int j = 0; j < (int) h1; j++)
		for(int i = 0; i < (int) w1; i++)  opixelsASIFT_H[j*woH+i] = ipixels1[j*w1+i];				
	
	for(int j = 0; j < (int) h2; j++)
		for(int i = 0; i < (int) w2; i++)  opixelsASIFT_H[j*woH + w1 + band_w + i] = ipixels2[j*w2 + i];	

	
	//////////////////////////////////////////////////////////////////// Draw matches
	matchingslist::iterator ptrH = matchings.begin();
	for(int i=0; i < (int) matchings.size(); i++, ptrH++)
	{		
		draw_line(opixelsASIFT_H, (int) (zoom1*ptrH->first.x), (int) (zoom1*ptrH->first.y), 
				  (int) (zoom2*ptrH->second.x) + w1 + band_w, (int) (zoom2*ptrH->second.y), 255.0f, woH, hoH);		
	}
	
	///////////////////////////////////////////////////////////////// Save imgOut	
	write_png_f32(argv[4], opixelsASIFT_H, woH, hoH, 1);
	
	delete[] opixelsASIFT_H; /*memcheck*/
	
	////// Write the coordinates of the matched points (row1, col1, row2, col2) to the file argv[5]
	std::ofstream file(argv[5]);
	if (file.is_open())
	{		
		// Write the number of matchings in the first line
		file << num_matchings << std::endl;
		
		matchingslist::iterator ptr = matchings.begin();
		for(int i=0; i < (int) matchings.size(); i++, ptr++)		
		{
			file << zoom1*ptr->first.x << "  " << zoom1*ptr->first.y << "  " <<  zoom2*ptr->second.x << 
			"  " <<  zoom2*ptr->second.y << std::endl;
		}		
	}
	else 
	{
		std::cerr << "Unable to open the file matchings."; 
	}

	file.close();


	
	// Write all the keypoints (row, col, scale, orientation, desciptor (128 integers)) to 
	// the file argv[6] (so that the users can match the keypoints with their own matching algorithm if they wish to)
	// keypoints in the 1st image
	std::ofstream file_key1(argv[6]);
	if (file_key1.is_open())
	{
		// Follow the same convention of David Lowe: 
		// the first line contains the number of keypoints and the length of the desciptors (128)
		file_key1 << num_keys1 << "  " << VecLength << "  " << std::endl;
		for (int tt = 0; tt < (int) keys1.size(); tt++)
		{
			for (int rr = 0; rr < (int) keys1[tt].size(); rr++)
			{
				keypointslist::iterator ptr = keys1[tt][rr].begin();
				for(int i=0; i < (int) keys1[tt][rr].size(); i++, ptr++)	
				{
					file_key1 << zoom1*ptr->x << "  " << zoom1*ptr->y << "  " << zoom1*ptr->scale << "  " << ptr->angle;
					
					for (int ii = 0; ii < (int) VecLength; ii++)
					{
						file_key1 << "  " << ptr->vec[ii];
					}
					
					file_key1 << std::endl;
				}
			}	
		}
	}
	else 
	{
		std::cerr << "Unable to open the file keys1."; 
	}

	file_key1.close();
	
	////// keypoints in the 2nd image
	std::ofstream file_key2(argv[7]);
	if (file_key2.is_open())
	{
		// Follow the same convention of David Lowe: 
		// the first line contains the number of keypoints and the length of the desciptors (128)
		file_key2 << num_keys2 << "  " << VecLength << "  " << std::endl;
		for (int tt = 0; tt < (int) keys2.size(); tt++)
		{
			for (int rr = 0; rr < (int) keys2[tt].size(); rr++)
			{
				keypointslist::iterator ptr = keys2[tt][rr].begin();
				for(int i=0; i < (int) keys2[tt][rr].size(); i++, ptr++)	
				{
					file_key2 << zoom2*ptr->x << "  " << zoom2*ptr->y << "  " << zoom2*ptr->scale << "  " << ptr->angle;
					
					for (int ii = 0; ii < (int) VecLength; ii++)
					{
						file_key2 << "  " << ptr->vec[ii];
					}					
					file_key2 << std::endl;
				}
			}	
		}
	}
	else 
	{
		std::cerr << "Unable to open the file keys2."; 
	}
	file_key2.close();
	
    return 0;
}
