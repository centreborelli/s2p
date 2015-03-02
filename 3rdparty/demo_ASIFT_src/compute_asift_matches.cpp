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
//*------------------------ compute_asift_matches-- -------------------------*/
// Match the ASIFT keypoints. 
// 
// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
// 
// Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image 
//            Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009. 
// Reference: ASIFT online demo (You can try ASIFT with your own images online.) 
//			  http://www.ipol.im/pub/algo/my_affine_sift/
/*---------------------------------------------------------------------------*/


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "compute_asift_matches.h"
#include "libMatch/match.h"
#include "orsa.h"


#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/* Remove the repetitive matches that appear in different simulations and retain only one */
void unique_match1(matchingslist &seg_in, matchingslist &seg_out, vector< vector <float> > &Minfoall_in, vector< vector <float> > &Minfoall_out)
{
  int i_in, i_out;
  float x1_in, x2_in, y1_in, y2_in, x1_out, x2_out, y1_out, y2_out;
  int flag_unique;
  float d1, d2;
  int Th2 = 2;
  
  seg_out.push_back(seg_in[0]);
  Minfoall_out.push_back(Minfoall_in[0]);	
	
  /* For other matches */
  if ( seg_in.size() > 1 )	
    {
      /* check if a match is unique. if yes, copy */
		matchingslist::iterator ptr_in = seg_in.begin();
		for ( i_in = 1; i_in < (int) seg_in.size(); i_in++, ptr_in++ )
			{				
			  x1_in = ptr_in->first.x;
			  y1_in = ptr_in->first.y;
			  x2_in = ptr_in->second.x; 
			  y2_in = ptr_in->second.y;	

			  flag_unique = 1;

			  matchingslist::iterator ptr_out = seg_out.begin(); 
			  for ( i_out = 0; i_out < (int) seg_out.size(); i_out++, ptr_out++ )
				{
					x1_out = ptr_out->first.x;
					y1_out = ptr_out->first.y;
					x2_out = ptr_out->second.x; 
					y2_out = ptr_out->second.y;	
					
					d1 = (x1_in - x1_out)*(x1_in - x1_out) + (y1_in - y1_out)*(y1_in - y1_out);
					d2 = (x2_in - x2_out)*(x2_in - x2_out) + (y2_in - y2_out)*(y2_in - y2_out);

	
					if ( ( d1 <= Th2) && ( d2 <= Th2) )
					 {
					   flag_unique = 0;
					   continue;
					 }	       
				}

			  if ( flag_unique == 1 )
				{									
					seg_out.push_back(seg_in[i_in]);
					Minfoall_out.push_back(Minfoall_in[i_in]);						
				}	  
		}
    }
}

/* Remove the ALL one-to-multiple matches. */
void clean_match1(matchingslist &seg_in, matchingslist &seg_out, vector< vector <float> > &Minfoall_in, vector< vector <float> > &Minfoall_out)
{
  int i1, i2;
  float x1_in, x2_in, y1_in, y2_in, x1_out, x2_out, y1_out, y2_out;

  // Guoshen Yu, 2010.09.22, Windows version
 // int flag_unique[seg_in.size()];
  int tmp_size = seg_in.size();
  int *flag_unique = new int[tmp_size];

  int sum_flag=0;
  float d1, d2;
  int Th1 = 1;
  int Th2 = 4;

 for ( i1 = 0; i1 < (int) seg_in.size(); i1++ )
    {
      flag_unique[i1] = 1;
    }

  /* Set the flag of redundant matches to 0. */
	matchingslist::iterator ptr_in = seg_in.begin();
	for ( i1 = 0; i1 < (int) seg_in.size() - 1; i1++, ptr_in++ )
    {					
		x1_in = ptr_in->first.x;
		y1_in = ptr_in->first.y;
		x2_in = ptr_in->second.x; 
		y2_in = ptr_in->second.y;

		matchingslist::iterator ptr_out = ptr_in+1;
		for ( i2 = i1 + 1; i2 < (int) seg_in.size(); i2++, ptr_out++ )
		{
		  x1_out = ptr_out->first.x;
		  y1_out = ptr_out->first.y;
		  x2_out = ptr_out->second.x; 
		  y2_out = ptr_out->second.y;	

		  d1 = (x1_in - x1_out)*(x1_in - x1_out) + (y1_in - y1_out)*(y1_in - y1_out);
		  d2 = (x2_in - x2_out)*(x2_in - x2_out) + (y2_in - y2_out)*(y2_in - y2_out);

		  /* If redundant, set flags of both elements to 0.*/
		  if ( ( ( d1 <= Th1) && ( d2 > Th2) ) || ( ( d1 > Th2) && ( d2 <= Th1) ) )
			{
			  flag_unique[i1] = 0;
			  flag_unique[i2] = 0;				  
			}
		}		
    }

	for ( i1 = 0; i1 < (int) seg_in.size(); i1++ )
    {
      sum_flag += flag_unique[i1];
    }
  
	/* Copy the matches that are not redundant */
	if ( sum_flag > 0 )
    {
	  for (  i1 = 0; i1 < (int) seg_in.size(); i1++ )
		{
		  if ( flag_unique[i1] == 1 )
			{				
				seg_out.push_back(seg_in[i1]);
				Minfoall_out.push_back(Minfoall_in[i1]);   
			}
		}
    }
  else
    {
      printf("Warning: all matches are redundant and are thus removed! This step of match cleaning is short circuited. (Normally this should not happen...)\n");
    }

	// Guoshen Yu, 2010.09.22, Windows version
	delete [] flag_unique;
}


/* Remove the ALL multiple-to-one matches */
void clean_match2(matchingslist &seg_in, matchingslist &seg_out, vector< vector <float> > &Minfoall_in, vector< vector <float> > &Minfoall_out)
{
	int i1, i2;
	float x1_in, x2_in, y1_in, y2_in, x1_out, x2_out, y1_out, y2_out;

  // Guoshen Yu, 2010.09.22, Windows version
  // int flag_unique[seg_in.size()];
  int tmp_size = seg_in.size();
  int *flag_unique = new int[tmp_size];

	int sum_flag=0;
	float d1, d2;
	int Th1 = 1;
	int Th2 = 4;

	for ( i1 = 0; i1 < (int) seg_in.size(); i1++ )
    {
      flag_unique[i1] = 1;
    }

	/* Set the flag of redundant matches to 0. */
	matchingslist::iterator ptr_in = seg_in.begin();
	for ( i1 = 0; i1 < (int) seg_in.size() - 1; i1++, ptr_in++ )
	{																				
		x1_in = ptr_in->first.x;
		y1_in = ptr_in->first.y;
		x2_in = ptr_in->second.x; 
		y2_in = ptr_in->second.y;	
															
		matchingslist::iterator ptr_out = ptr_in+1;
		for ( i2 = i1 + 1; i2 < (int) seg_in.size(); i2++, ptr_out++ )
		{
			
		  x1_out = ptr_out->first.x;
		  y1_out = ptr_out->first.y;
		  x2_out = ptr_out->second.x; 
		  y2_out = ptr_out->second.y;	

		  d1 = (x1_in - x1_out)*(x1_in - x1_out) + (y1_in - y1_out)*(y1_in - y1_out);
		  d2 = (x2_in - x2_out)*(x2_in - x2_out) + (y2_in - y2_out)*(y2_in - y2_out);


		  /* If redundant, set flags of both elements to 0.*/	
		  if ( ( d1 > Th2) && ( d2 <= Th1) ) 
			{
			  flag_unique[i1] = 0;
			  flag_unique[i2] = 0;				  
			}
		}
	}

  for ( i1 = 0; i1 < (int) seg_in.size(); i1++ )
  {
	  sum_flag += flag_unique[i1];
  }
  
  /* Copy the matches that are not redundant */
  if ( sum_flag > 0 )
  {
	  for (  i1 = 0; i1 < (int) seg_in.size(); i1++ )
	  {
		  if ( flag_unique[i1] == 1 )
		  {			 
				seg_out.push_back(seg_in[i1]);
				Minfoall_out.push_back(Minfoall_in[i1]);    
		  }
	  }
  }
  else
  {
      printf("Warning: all matches are redundant and are thus removed! This step of match cleaning is short circuited. (Normally this should not happen...)\n");
  }

  // Guoshen Yu, 2010.09.22, Windows version
  delete [] flag_unique;
}


// Normalize the coordinates of the matched points by compensating the simulate affine transformations
void compensate_affine_coor(matching &matching1, int w1, int h1, int w2, int h2, float t1, float t2, float Rtheta, float t_im2_1, float t_im2_2, float Rtheta2)
{
	float x_ori, y_ori;	
	float x_ori2, y_ori2, x_tmp, y_tmp;
	float x1, y1, x2, y2;
		
	Rtheta = Rtheta*PI/180;
	
	if ( Rtheta <= PI/2 )
    {
		x_ori = 0;
		y_ori = w1 * sin(Rtheta) / t1;
    }
	else
    {
		x_ori = -w1 * cos(Rtheta) / t2;
		y_ori = ( w1 * sin(Rtheta) + h1 * sin(Rtheta-PI/2) ) / t1;
    }
	
	Rtheta2 = Rtheta2*PI/180;
	
	if ( Rtheta2 <= PI/2 )
    {
		x_ori2 = 0;
		y_ori2 = w2 * sin(Rtheta2) / t_im2_1;
    }
	else
    {
		x_ori2 = -w2 * cos(Rtheta2) / t_im2_2;
		y_ori2 = ( w2 * sin(Rtheta2) + h2 * sin(Rtheta2-PI/2) ) / t_im2_1;
    }
	
	float sin_Rtheta = sin(Rtheta);
	float cos_Rtheta = cos(Rtheta);
	float sin_Rtheta2 = sin(Rtheta2);
	float cos_Rtheta2 = cos(Rtheta2);
	
	x1 = matching1.first.x;
	y1 = matching1.first.y;	
	x2 = matching1.second.x;
	y2 = matching1.second.y;
	
	/* project the coordinates of im1 to original image before tilt-rotation transform */
	/* Get the coordinates with respect to the 'origin' of the original image before transform */
	x1 = x1 - x_ori;
	y1 = y1 - y_ori;
	/* Invert tilt */
	x1 = x1 * t2;
	y1 = y1 * t1;
	/* Invert rotation (Note that the y direction (vertical) is inverse to the usual concention. Hence Rtheta instead of -Rtheta to inverse the rotation.) */
	x_tmp = cos_Rtheta*x1 - sin_Rtheta*y1;
	y_tmp = sin_Rtheta*x1 + cos_Rtheta*y1;
	x1 = x_tmp;
	y1 = y_tmp;
	
	/* Coordinate projection on image2  */
	
	/* Get the coordinates with respect to the 'origin' of the original image before transform */
	x2 = x2 - x_ori2;
	y2 = y2 - y_ori2;
	/* Invert tilt */
	x2 = x2 * t_im2_2;
	y2 = y2 * t_im2_1;
	/* Invert rotation (Note that the y direction (vertical) is inverse to the usual concention. Hence Rtheta instead of -Rtheta to inverse the rotation.) */
	x_tmp = cos_Rtheta2*x2 - sin_Rtheta2*y2;
	y_tmp = sin_Rtheta2*x2 + cos_Rtheta2*y2;
	x2 = x_tmp;
	y2 = y_tmp;
	
	matching1.first.x = x1;
	matching1.first.y = y1;	
	matching1.second.x = x2;
	matching1.second.y = y2;
}

int compute_asift_matches(int num_of_tilts1, int num_of_tilts2, int w1, int h1, int w2, int h2, int verb, vector< vector< keypointslist > >& keys1, vector< vector< keypointslist > >& keys2, matchingslist &matchings, siftPar &siftparameters)
// Match the ASIFT keypoints. 
// Input:
// num_of_tilts1, num_of_tilts2: number of tilts that have been simulated on the two images. (They can be different.)
// w1, h1, w2, h2: widht/height of image1/image2.
// verb: 1/0 --> show/don not show verbose messages. (1 for debugging) 
// keys1, keys2: ASIFT keypoints of image1/image2. (They should be calculated with compute_asift_keypoints.)
// matchings (output): the coordinates (col1, row1, col2, row2) of all the matching points. 
//
// Output: the number of matching points.
{  	
	float t_min, t_k, t;
	int num_tilt1, num_tilt2, tt, num_rot_t2, num_rot1, rr;
	int cc;

	int tt2, rr2, num_rot1_2;
	float t_im2;

	/* It stores the coordinates of ALL matches points of ALL affine simulations  */
	vector< vector <float> > Minfoall;

	int Tmin = 8;					
	float nfa_max = -2;

	num_rot_t2 = 10;

	t_min = 1;
	t_k = sqrt(2.);

	num_tilt1 = num_of_tilts1;
	num_tilt2 = num_of_tilts2;

	if ( ( num_tilt1 < 1 ) || ( num_tilt2 < 1 ) )
	{
		printf("Number of tilts num_tilt should be equal or larger than 1. \n");
		exit(-1);	
	}

	
	/* Initialize the vector structure for the matching points */	
	std::vector< vector< vector < vector < matchingslist > > > > matchings_vec(num_tilt1);	
	std::vector< vector< vector< vector< vector< vector <float> > > > > > Minfoall_vec(num_tilt1);	
	for (tt = 1; tt <= num_tilt1; tt++)
	{		
		t = t_min * pow(t_k, tt-1);
		if ( t == 1 )
		{			
			num_rot1 = 1;
		}
		else
		{
			num_rot1 = round(num_rot_t2*t/2);        
			if ( num_rot1%2 == 1 )
			{
				num_rot1 = num_rot1 + 1;
			}
			num_rot1 = num_rot1 / 2;
		}     
		
		matchings_vec[tt-1].resize(num_rot1);
		Minfoall_vec[tt-1].resize(num_rot1);
		
		for  ( rr = 1; rr <= num_rot1; rr++ ) 
		{   
			
			matchings_vec[tt-1][rr-1].resize(num_tilt2);
			Minfoall_vec[tt-1][rr-1].resize(num_tilt2);
			
			for (tt2 = 1; tt2 <= num_tilt2; tt2++)				
			{			
				t_im2 = t_min * pow(t_k, tt2-1);
				if ( t_im2 == 1 )
				{
					num_rot1_2 = 1;
				}
				else
				{
					num_rot1_2 = round(num_rot_t2*t_im2/2);        
					if ( num_rot1_2%2 == 1 )
					{
						num_rot1_2 = num_rot1_2 + 1;
					}
					num_rot1_2 = num_rot1_2 / 2;
				}
				
				matchings_vec[tt-1][rr-1][tt2-1].resize(num_rot1_2);		
				Minfoall_vec[tt-1][rr-1][tt2-1].resize(num_rot1_2);				
			}
		}
	}
	
	
	///*
	// * setup the tilt and rotation parameters
	// * for all the loops, this vector will hold
	// * the following parameters: 
	// * tt, num_rot1, rr, tt2, num_rot1_2, rr2
	// */
	//vector<int> tilt_rot;
	///* loop on tilts for image 1 */ 
	//for (int tt = 1; tt <= num_tilt1; tt++)
	//{
	//    float t = t_min * pow(t_k, tt-1);
	//    int num_rot1;
	//    /* if tilt t = 1, do not simulate rotation. */
	//    if ( 1 == tt )
	//		num_rot1 = 1;
	//    else
	//    {
	//		/* number of rotations to simulate */
	//		num_rot1 = round(num_rot_t2 * t / 2);        
	//		if ( num_rot1%2 == 1 )
	//			num_rot1 = num_rot1 + 1;
	//		num_rot1 = num_rot1 / 2;
	//    }
	//    /* loop on rotations for image 1 */ 
	//    for  (int rr = 1; rr <= num_rot1; rr++ ) 
	//    {   
	//		/* loop on tilts for image 2 */
	//		for (int tt2 = 1; tt2 <= num_tilt2; tt2++)
	//		{
	//			float t_im2 = t_min * pow(t_k, tt2-1);
	//			int num_rot1_2;
	//			if ( tt2 == 1 )
	//				num_rot1_2 = 1;
	//			else
	//			{
	//				num_rot1_2 = round(num_rot_t2 * t_im2 / 2);        
	//				if ( num_rot1_2%2 == 1 )
	//					num_rot1_2 = num_rot1_2 + 1;
	//				num_rot1_2 = num_rot1_2 / 2;
	//			}
	//			/* loop on rotations for image 2 */ 
	//			for  (int rr2 = 1; rr2 <= num_rot1_2; rr2++ ) 
	//			{
	//				tilt_rot.push_back(tt);
	//				tilt_rot.push_back(num_rot1);
	//				tilt_rot.push_back(rr);
	//				tilt_rot.push_back(tt2);
	//				tilt_rot.push_back(num_rot1_2);
	//				tilt_rot.push_back(rr2);
	//			}
	//		}
	//    }
	//}
		
	/* Calculate the number of simulations */
#ifdef _OPENMP
	omp_set_nested(1);
#endif 
	// loop on tilts for image 1. 
#pragma omp parallel for private(tt)
	for (int tt = 1; tt <= num_tilt1; tt++)
	{
		
		float t = t_min * pow(t_k, tt-1);
		
		/* Attention: the t1, t2 do not follow the same convention as in compute_asift_keypoints */   
		float t1 = t; 
		float t2 = 1;
		
		int num_rot1;
		
		// If tilt t = 1, do not simulate rotation. 	
		if ( tt == 1 )
		{
			num_rot1 = 1;
		}
		else
		{
			// The number of rotations to simulate under the current tilt.   
			num_rot1 = round(num_rot_t2*t/2);        
			if ( num_rot1%2 == 1 )
			{
				num_rot1 = num_rot1 + 1;
			}
			num_rot1 = num_rot1 / 2;
		}
		
		
		float delta_theta = PI/num_rot1;
		
		// Loop on rotations for image 1. 
#pragma omp parallel for private(rr)
		for  ( int rr = 1; rr <= num_rot1; rr++ ) 
		{   
			float theta = delta_theta * (rr-1);
			theta = theta * 180 / PI;
			
			/* Read the keypoints of image 1 */	       
			keypointslist keypoints1 = keys1[tt-1][rr-1];
			
			// loop on tilts for image 2.
#pragma omp parallel for private(tt2)
			for (int tt2 = 1; tt2 <= num_tilt2; tt2++)				
			{
				float t_im2 = t_min * pow(t_k, tt2-1);
				
				/* Attention: the t1, t2 do not follow the same convention as in asift_v1.c */
				float t_im2_1 = t_im2; 
				float t_im2_2 = 1;
				
				int num_rot1_2;
				
				if ( tt2 == 1 )
				{
					num_rot1_2 = 1;
				}
				else
				{
					num_rot1_2 = round(num_rot_t2*t_im2/2);        
					if ( num_rot1_2%2 == 1 )
					{
						num_rot1_2 = num_rot1_2 + 1;
					}
					num_rot1_2 = num_rot1_2 / 2;
				}
				
				float delta_theta2 = PI/num_rot1_2;
				
#pragma omp parallel for private(rr2)
				// Loop on rotations for image 2. 
				for  ( int rr2 = 1; rr2 <= num_rot1_2; rr2++ ) 
				{   
					float theta2 = delta_theta2 * (rr2-1);
					theta2 = theta2 * 180 / PI;
					
					/* Read the keypoints of image2. */	       
					keypointslist keypoints2 = keys2[tt2-1][rr2-1];
					
					
					// Match the keypoints of image1 and image2.
					matchingslist matchings1;					
					compute_sift_matches(keypoints1,keypoints2,matchings1,siftparameters);		       
					
					if ( verb ) 
					{
						printf("t1=%.2f, theta1=%.2f, num keys1 = %d, t2=%.2f, theta2=%.2f, num keys2 = %d, num matches=%d\n", t, theta, (int) keypoints1.size(), t_im2, theta2, (int) keypoints2.size(), (int) matchings1.size());
					}
					
					/* Store the matches */
					if ( matchings1.size() > 0 )
					{
						matchings_vec[tt-1][rr-1][tt2-1][rr2-1] = matchingslist(matchings1.size());
            Minfoall_vec[tt-1][rr-1][tt2-1][rr2-1].resize(matchings1.size());
						
						for ( int cc = 0; cc < (int) matchings1.size(); cc++ )		   
						{		      														
							///// In the coordinates the affine transformations have been normalized already in compute_asift_keypoints. So no need to normalize here. 
							// Normalize the coordinates of the matched points by compensating the simulate affine transformations
							//	compensate_affine_coor(matchings1[cc], w1, h1, w2, h2, t1, t2, theta, t_im2_1, t_im2_2, theta2);
							
							matchings_vec[tt-1][rr-1][tt2-1][rr2-1][cc] = matchings1[cc];
							
							vector<float> Minfo_1match(6);
							Minfo_1match[0] = t1;
							Minfo_1match[1] = t2;
							Minfo_1match[2] = theta;
							Minfo_1match[3] = t_im2_1;
							Minfo_1match[4] = t_im2_2;
							Minfo_1match[5] = theta2;									
							Minfoall_vec[tt-1][rr-1][tt2-1][rr2-1][cc] = Minfo_1match;
						} 												
					}
				}
			}
		}
	}
			
	// Move the matches to a 1D vector
	for (tt = 1; tt <= num_tilt1; tt++)
	{				
		t = t_min * pow(t_k, tt-1);

		if ( t == 1 )
		{			
			num_rot1 = 1;
		}
		else
		{
			num_rot1 = round(num_rot_t2*t/2);        
			if ( num_rot1%2 == 1 )
			{
				num_rot1 = num_rot1 + 1;
			}
			num_rot1 = num_rot1 / 2;
		}     
					
		for  ( rr = 1; rr <= num_rot1; rr++ ) 
		{   			
			for (tt2 = 1; tt2 <= num_tilt2; tt2++)				
			{			
				t_im2 = t_min * pow(t_k, tt2-1);
				if ( t_im2 == 1 )
				{
					num_rot1_2 = 1;
				}
				else
				{
					num_rot1_2 = round(num_rot_t2*t_im2/2);        
					if ( num_rot1_2%2 == 1 )
					{
						num_rot1_2 = num_rot1_2 + 1;
					}
					num_rot1_2 = num_rot1_2 / 2;
				}
		
				for  ( rr2 = 1; rr2 <= num_rot1_2; rr2++ ) 
				{ 
					for ( cc=0; cc < (int) matchings_vec[tt-1][rr-1][tt2-1][rr2-1].size(); cc++ )
					{
						matchings.push_back(matchings_vec[tt-1][rr-1][tt2-1][rr2-1][cc]);
						Minfoall.push_back(Minfoall_vec[tt-1][rr-1][tt2-1][rr2-1][cc]);
					}
				}
			}
		}
	}

	if ( verb )
	{
	  printf("The number of matches is %d \n", (int) matchings.size());
	}

	
	if ( matchings.size() > 0 )
	{
	  /* Remove the repetitive matches that appear in different simulations and retain only one. */		
		// Since tilts are simuated on both image 1 and image 2, it is normal to have repetitive matches. 
		matchingslist matchings_unique;	
		vector< vector<float> > Minfoall_unique;			
		unique_match1(matchings, matchings_unique, Minfoall, Minfoall_unique);      		
		matchings = matchings_unique;
		Minfoall = Minfoall_unique; 
		
		if ( verb )    
			{
			  printf("The number of unique matches is %d \n", (int) matchings.size());
			}
	  
		// There often appear to be some one-to-multiple/multiple-to-one matches (one point in image 1 matches with many points in image 2/vice versa). 
		// This is an artifact of SIFT on interpolated images, as the interpolation tends to create some auto-similar structures (steps for example). 
		// These matches need to be removed.		
		  /* Separating the removal of  multiple-to-one and one-to-multiple in two steps:
		 - first remove multiple-to-one
		 - then remove one-to-multiple
		 This allows to avoid removing some good matches: multiple-to-one matches is much more frequent than one-to-multiple. Sometimes some of the feature points in image 1 that take part in "multiple-to-one" bad matches have also correct matches in image 2. The modified scheme avoid removing these good matches. */
	 
		// Remove to multiple-to-one matches	
		matchings_unique.clear();
		Minfoall_unique.clear();		
		clean_match2(matchings, matchings_unique, Minfoall, Minfoall_unique);    		
		matchings = matchings_unique;
		Minfoall = Minfoall_unique;
		
		// Remove to one-to-multiple matches	
		matchings_unique.clear();
		Minfoall_unique.clear();		
		clean_match1(matchings, matchings_unique, Minfoall, Minfoall_unique);    		
		matchings = matchings_unique;
		Minfoall = Minfoall_unique;

				
		if ( verb )
		{
			printf("The number of final matches is %d \n", (int) matchings.size());
		}
		
		// If enough matches to do epipolar filtering
		if ( (int) matchings.size() >= Tmin )
		{			
			//////// Use ORSA to filter out the incorrect matches. 
			// store the coordinates of the matching points
			vector<Match> match_coor; 
			for ( cc = 0; cc < (int) matchings.size(); cc++ )	
			{
				Match match1_coor;
				match1_coor.x1 = matchings[cc].first.x;
				match1_coor.y1 = matchings[cc].first.y;
				match1_coor.x2 = matchings[cc].second.x;
				match1_coor.y2 = matchings[cc].second.y;
				
				match_coor.push_back(match1_coor);							
			}
			
			std::vector<float> index;
			// Guoshen Yu, 2010.09.23
			// index.clear();
			
			int t_value_orsa=10000;
			int verb_value_orsa=0;
			int n_flag_value_orsa=0;
			int mode_value_orsa=2;
			int stop_value_orsa=0;
			
			// epipolar filtering with the Moisan-Stival ORSA algorithm. 
//			float nfa = orsa(w1, h1, match_coor, index, t_value_orsa, verb_value_orsa, n_flag_value_orsa, mode_value_orsa, stop_value_orsa);
			float nfa = orsa((w1+w2)/2, (h1+h2)/2, match_coor, index, t_value_orsa, verb_value_orsa, n_flag_value_orsa, mode_value_orsa, stop_value_orsa);


			// if the matching is significant, register the good matches
			if ( nfa < nfa_max )
			{
				// extract meaningful matches
				matchings_unique.clear();
				Minfoall_unique.clear();
				for ( cc = 0; cc < (int) index.size(); cc++ )				
				{
					matchings_unique.push_back(matchings[(int)index[cc]]);
					Minfoall_unique.push_back(Minfoall[(int)index[cc]]);
				}
				matchings = matchings_unique;
				Minfoall = Minfoall_unique;     
				
				cout << "The two images match! " << matchings.size() << " matchings are identified. log(nfa)=" << nfa << "." << endl;
			}
			else 
			{
				matchings.clear();
				Minfoall.clear();
				cout << "The two images do not match. The matching is not significant: log(nfa)=" << nfa << "." << endl;
			}
		}
		else 
		{
			matchings.clear();
			Minfoall.clear();
			cout << "The two images do not match. Not enough matches to do epipolar filtering." << endl;		
		}
	}
	else 
	{
		cout << "The two images do not match.\n" << endl;
	}
	
	return matchings.size();

}

