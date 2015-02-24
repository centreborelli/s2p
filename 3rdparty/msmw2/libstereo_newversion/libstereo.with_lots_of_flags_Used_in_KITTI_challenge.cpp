#include "libstereo.h"
#include "smartparameter.h"
#include "cubic.h"

// Obscure parameters
// by default old IJCV setting
SMART_PARAMETER_INT(DONT_USE_SUBPIX_AT_LOW_SCALES,0)
SMART_PARAMETER_INT(USE_MODIFIED_MINDIST,0)
SMART_PARAMETER_INT(USE_OLD_SELFSIM,1)
SMART_PARAMETER_INT(USE_OLD_SELFSIM_ORDER,1)
SMART_PARAMETER_INT(DISABLE_MINDIFFDILATATION_AT_LOWER_SCALES,0)
SMART_PARAMETER_INT(USE_OLD_MERGE_LR,1)
SMART_PARAMETER_INT(USE_LR_CONSISTENT_SELECTION_OF_WINDOW_ORIENTATION,0)
SMART_PARAMETER_INT(DONT_REMOVE_ISOLATED_INNER,0)
 


   int global_window_index=0;
   int global_scale_index=0;

namespace libIIPStable {
	
    
    
    
    
    
    
    
    
    ///////////////////////
    ///////////////////////  Window related
    ///////////////////////
    
    
	//! Build Gaussian kernel taking alpha at boundary before normalization to have integral 1
	void  stereo_build_gaussian_kernel(flimage &prolate, float alpha, float beta)
	{
		int pwidth = prolate.w();
		int spwidth = (prolate.w()-1)/2;
		int spheight = (prolate.h()-1)/2;
		
		float sigmax2 = -(float) spwidth * (float) spwidth / (2.0f* log(alpha));
		float sigmay2 = -(float) spheight * (float) spheight / (2.0f* log(beta));
		
		for (int i=-spwidth; i <= spwidth; i++)
			for (int j=-spheight; j <= spheight; j++)
			{
				
				float x2 = (float)(i*i);
				float y2 = (float)(j*j);
                
                int l= (spheight + j) * pwidth + spwidth + i;
				prolate[l] = expf(-(0.5f * x2 / sigmax2) - (0.5f * y2 / sigmay2) );
			}
		
	}
    
    
    
    
    
    
    ///////////////////////
    ///////////////////////  Correlation
    ///////////////////////
    
    
    //! Pixelian correlation
    //! Takes a correlation window as input inside strPar
    //! oselfdist: best self match and olrfdist: subpixel reciprocity computed if flags activated
	void stereo_pixel_chain_one_direction(cflimage &input, cflimage &input2, flimage &imask, flimage &Dmin, flimage &Dmax, flimage &odisp, flimage &omask, flimage &odist, flimage &odist2, flimage &oselfdist, flimage &olrfdist, flimage &ocdif, strParameters &strPar)
	{
		
        
		//! Variable: prolate
		flimage corrwindow = strPar.prolate;
		corrwindow.normalizeL1();
		
        
		//! Variable: boundary
		int spwidth = (corrwindow.w() - 1) / 2;
		int spheight = (corrwindow.h() - 1) / 2;
		int boundary =  2*MAX(spwidth, spheight) + 1;
		
		
        
        
		//! Compute translation of second image for subpixel estimation
		//! It does not compute anything if inPrecisions = 1 but
        //! creates images2[0] = input2
        cflimage *images2 = new cflimage[strPar.inPrecisions];
		float step = 1.0 / (float) strPar.inPrecisions;
		images2[0] = input2;
		for (int ii=1; ii < strPar.inPrecisions; ii++)
		{
			images2[ii] = input2;
			images2[ii].fftShear(0.0, iipHorizontal, (float) ii * step, 0, 1);
			
		}
        
        
        //! Compute translation of first image for self similarity and left-right consistency
        cflimage *images1 = new cflimage[strPar.inPrecisions];
        images1[0] = input;
        if (strPar.flagSelfSimilarity || strPar.flagSubRecip)
        {
            for (int ii=1; ii < strPar.inPrecisions; ii++)
            {
                images1[ii] = input;
                images1[ii].fftShear(0.0, iipHorizontal, (float) ii * step, 0, 1);
                
            }
        }
        
        
        //! Compute mean of left and translated right image if idistype == L2M
		cflimage mean;
		cflimage *mean2 = new cflimage[strPar.inPrecisions];
		if (strPar.itypeDist == L2M)
		{
			if (strPar.flagListKernels) {
				mean= input.patchListMean(corrwindow);
				//mean= input.patchMean(corrwindow);
				for (int ii=0; ii < strPar.inPrecisions; ii++)  mean2[ii] = images2[ii].patchListMean(corrwindow);
				//for (int ii=0; ii < strPar.inPrecisions; ii++)  mean2[ii] = images2[ii].patchMean(corrwindow);
			} else {
				mean= input.patchMean(corrwindow);
				for (int ii=0; ii < strPar.inPrecisions; ii++)  mean2[ii] = images2[ii].patchMean(corrwindow);
			}
		}
		
        
        
        //! Compute mean of left translated images if idistype == L2M and flagSelfSimilarity activated
		cflimage *mean1 = new cflimage[strPar.inPrecisions];
		if (strPar.itypeDist == L2M && (strPar.flagSelfSimilarity || strPar.flagSubRecip))
		{
			if (strPar.flagListKernels) {
				for (int ii=0; ii < strPar.inPrecisions; ii++)  mean1[ii] = images1[ii].patchListMean(corrwindow);			
				//for (int ii=0; ii < strPar.inPrecisions; ii++)  mean1[ii] = images1[ii].patchMean(corrwindow);			
         }
			else
				for (int ii=0; ii < strPar.inPrecisions; ii++)  mean1[ii] = images1[ii].patchMean(corrwindow);

		}
		
        
        
        
        
		// Initialization
		omask = -1.0f;
		odist = fLarge;
		oselfdist = fLarge;
      cflimage oselfdisp = oselfdist;
      oselfdisp = 0;
        
        
        //! Begin for each pixel
#pragma omp parallel shared( images2, mean, mean2, mean1, images1,corrwindow)
		{
			
			
#pragma omp for schedule(dynamic) nowait
            for(int ipy = boundary ;  ipy < input.h() - boundary ; ipy++)
            {
                
                //! Compute correlation only at points with imask>0
                for(int ipx = boundary ;  ipx < input.w() - boundary; ipx++)
                    if (imask[ipy*imask.w() + ipx] > 0.0f)
                    {
                        
                        
                        
                        
                        
                        int l = ipy*input.w() + ipx;
                        
                        
                        //! correlation with first image if flagSelfSimilarity used
                        float fSelfBestDistance = fLarge;
                        float fSelfBestDisp = fLarge;
                        if (strPar.flagSelfSimilarity)
                        {
                            
                            //! range adapted to window size
         // The self-similarity range is centered at the reference pixel
         // This permits to spot self similar structures even if the search range is shifted
                           int r = fmax(ceil((Dmax[l]-Dmin[l])/2.0),3.0);
                           int imin = MAX(boundary, ipx - (int) r );
                           int imax = MIN(input.w() - boundary, ipx + (int) r );

                            if ( USE_OLD_SELFSIM() ) {
         // Use the range of the right image for the self-similarity test
         // This has a clear problem the image may not be self-similar in the considered range.
         // (this is the one used in the paper)
                               imin = MAX(boundary, ipx + (int) Dmin[l]);
                               imax = MIN(input.w() - boundary, ipx + (int) Dmax[l]);
                            }
                            
                            
                            //! compute distance
                            int iik=0;
                            for(int ci = imin; ci <= imax; ci++)
                                for(int ii = 0; ii < strPar.inPrecisions; ii++, iik++)
                                    if (abs(ci - ipx) > 1)
                                    {
                                        
                                        
                                        float fCurrentDistance = 0.0f;
                                        
                                        if (strPar.itypeDist == L2)
                                        {
                                            if (strPar.flagListKernels) 
											{
                                                fCurrentDistance = distancePatchListWL2(input, images1[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow);
											} 
											else 
											{
                                            if (strPar.flagWinWeighted)
                                                fCurrentDistance = distancePatchWL2(input, images1[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow);
                                            else
                                                fCurrentDistance = distancePatchL2(input, images1[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow.w(), corrwindow.h());
											}                                       
                                        }
                                        else if (strPar.itypeDist == L2M)
                                        {
											if (strPar.flagListKernels) 
											{
												fCurrentDistance = distancePatchListWL2M(input, images1[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow, mean, mean1[ii]);
											} 
											else 
											{
												if (strPar.flagWinWeighted)
													fCurrentDistance = distancePatchWL2M(input, images1[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow, mean, mean1[ii]);
												else
													fCurrentDistance = distancePatchL2M(input, images1[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow.w(), corrwindow.h(), mean, mean1[ii]);
											}


                                            
                                        } else {printf("... please enter a correct distance type"); exit(-1);}
                                        
                                        if (fCurrentDistance < fSelfBestDistance)  {
                                           fSelfBestDistance = fCurrentDistance;
                                           fSelfBestDisp = (float) ci + (float) ii * step -ipx;   
                                        }
                                    }
                            
                            oselfdist[l] = fSelfBestDistance;
                            oselfdisp[l] = fSelfBestDisp;
                            
                        }
                        
                        
                        
                        //! correlation with second image
                        float fBestDistance = fLarge;
                        float fSelectedPixel = 0.0;
                        float fPixelDistance = fLarge;
                        {
                            
                            //! range adapted to window size
                            int imin = MAX(boundary, ipx + (int) Dmin[l]);
                            int imax = MIN(input2.w() - boundary, ipx + (int) Dmax[l]);
                            
                            
                            //printf("sc x: %d   y: %d  w: %d  h: %d  imin: %d  imax: %d   (spw: %d  sph: %d)\n ", ipx, ipy, input.w(), input.h(), imin, imax, spwidth, spheight);
                            
                            
                            //! distance function taking into account translations
                            int ilength = strPar.inPrecisions * (imax-imin + 1);
                            if (imax<imin) 
                              ilength = strPar.inPrecisions * (imin-imax + 1);
                            //float *fldistance=new float[ilength];
                            float fldistance[ilength];
                            for(int ii=0; ii < ilength; ii++) fldistance[ii] = fLarge;
                            
                            
                            //! compute distance
                            int iik=0;
                            for(int ci = imin; ci <= imax; ci++)
                                for(int ii = 0; ii < strPar.inPrecisions; ii++, iik++)
                                {
                                    
                                    
                                    
                                    float fCurrentDistance = 0.0f;
                                    
                                    if (strPar.itypeDist == L2)
                                    {
                                        if (strPar.flagListKernels) 
                                        {
                                            fCurrentDistance = distancePatchListWL2(input, images2[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow);
                                        } 
                                        else 
                                        {                                        
                                            if (strPar.flagWinWeighted)
                                                fCurrentDistance = distancePatchWL2(input,images2[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow);
                                            else
                                                fCurrentDistance = distancePatchL2(input,images2[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow.w(), corrwindow.h());
                                        }
                                    }
                                    else if (strPar.itypeDist == L2M)
                                    {
                                        if (strPar.flagListKernels) 
                                        {
                                            fCurrentDistance = distancePatchListWL2M(input, images2[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow, mean, mean2[ii]);
                                        } 
                                        else 
                                        {     
                                            if (strPar.flagWinWeighted)
                                                fCurrentDistance = distancePatchWL2M(input,images2[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow, mean, mean2[ii]);
                                            else
                                                fCurrentDistance = distancePatchL2M(input,images2[ii],ipx-spwidth,ipy-spheight,ci-spwidth,ipy-spheight,corrwindow.w(), corrwindow.h(), mean, mean2[ii]);
                                        }
                                        
                                    } else {printf("... please enter a correct distance type"); exit(-1);}
                                    
                                    
                                    fldistance[iik] = fCurrentDistance;
                                    
                                    if (fCurrentDistance < fBestDistance)
                                    {
                                        fBestDistance = fCurrentDistance;
                                        fSelectedPixel = (float) ci + (float) ii * step;
                                        fPixelDistance =  distanceL2(input,images2[ii],ipx ,ipy , ci ,ipy);
                                    }
                                    
                                    
                                }


                            //! refine the minimum by a factor 4 using cubic interpolation
                            for (int ii=1;ii<iik-2;ii++) {
                               float tdist, tx;
                               float increment = 1.0 / strPar.cubicRefinementFactor;
                               for(tx=increment; tx<1; tx+=increment) {
                                  tdist = cubicInterpolate ( fldistance+(ii-1), tx);
                                  if(fBestDistance > tdist ) {
                                     fSelectedPixel = (float) imin + (float) ii*step + tx*step;
                                     fBestDistance  = tdist;
                                  }
                               }
                               // we can go to the analytic minimum, but is is not stable.
//                               if(fBestDistance > cubicMinimum(  fldistance+(ii-1), &tdist,&tx)) {
////                                  if(tdist>=0) {
//                                    fSelectedPixel = (float) imin + (float) ii*step + tx*step;
//                                    fBestDistance  = tdist;
////                                  }
//                               }
                            }

                            
                            //! check for first and second minima of distance
                            {
                                float dif0, dif1;
                                stereo_find_local_minima(fldistance, ilength, dif1, dif0);
                                odist[l] = dif0;
                                odist2[l] = dif1;
                                
                            }
                            
                            
                            // check if minima at extremum of range
                            //if ( fSelectedPixel == (float) imin || fSelectedPixel == (float) imax)
                            //{omask[l] = 0.0f;}
                            
                            
                           //delete[] fldistance;
                        }
                        
                        
                        
                        //! left-right with subpixel precision
                        if (strPar.flagSubRecip)
                        {
                            
                            
                            
                            int iSelectedPixel = (int) floorf(fSelectedPixel);
                            int iiSel = rint((fSelectedPixel - (float) iSelectedPixel) / step );
                            
                            
                            //! range adapted to window size
                            int imin = MAX(boundary, iSelectedPixel + (int) Dmin[l]);
                            int imax = MIN(input.w() - boundary, iSelectedPixel + (int) Dmax[l]);
                            
                            
                            //! compute distance
                            int iik=0;
                            float fLRBestDistance = fLarge;
                            float fLRSelectedPixel;
                            for(int ci = imin; ci <= imax; ci++)
                                for(int ii = 0; ii < strPar.inPrecisions; ii++, iik++)
                                {
                                    
                                    
                                    float fCurrentDistance = 0.0f;
                                    
                                    if (strPar.itypeDist == L2)
                                    {
                                        if (strPar.flagListKernels) 
                                        {
                                            fCurrentDistance = distancePatchListWL2(images1[ii], images2[iiSel], ci-spwidth, ipy-spheight, iSelectedPixel - spwidth ,ipy-spheight,corrwindow);
                                        } 
                                        else 
                                        {         
                                        if (strPar.flagWinWeighted)
                                            fCurrentDistance = distancePatchWL2(images1[ii], images2[iiSel],  ci-spwidth, ipy-spheight, iSelectedPixel - spwidth ,ipy-spheight,corrwindow);
                                        else
                                            fCurrentDistance = distancePatchL2(images1[ii], images2[iiSel], ci-spwidth, ipy-spheight, iSelectedPixel - spwidth ,ipy-spheight,corrwindow.w(), corrwindow.h());
                                        }
                                    }
                                    else if (strPar.itypeDist == L2M)
                                    {
                                        if (strPar.flagListKernels) 
                                        {
                                            fCurrentDistance = distancePatchListWL2M(images1[ii], images2[iiSel], ci-spwidth, ipy-spheight, iSelectedPixel - spwidth ,ipy-spheight,corrwindow,  mean1[ii], mean2[iiSel]);
                                        } 
                                        else 
                                        {  
                                        if (strPar.flagWinWeighted)
                                            fCurrentDistance = distancePatchWL2M(images1[ii], images2[iiSel],  ci-spwidth, ipy-spheight, iSelectedPixel - spheight ,ipy-spheight,corrwindow,  mean1[ii], mean2[iiSel]);
                                        else
                                            fCurrentDistance = distancePatchL2M(images1[ii], images2[iiSel],ci-spwidth,ipy-spheight, iSelectedPixel - spheight ,ipy-spheight,corrwindow.w(), corrwindow.h() , mean1[ii], mean2[iiSel]);
                                        }
                                        
                                    } else {printf("... please enter a correct distance type"); exit(-1);}
                                    
                                    if (fCurrentDistance < fLRBestDistance) { fLRBestDistance = fCurrentDistance;  fLRSelectedPixel = (float) ci + (float) ii * step; }
                                    
                                }
                            
                            
                            olrfdist[l] = fabsf( (float) ipx - fLRSelectedPixel);
                            
                            
                        }
                        
                        
                        
                        //! if pixel not invalidated
                        odisp[l] = fSelectedPixel - (float) ipx;
                        omask[l] = 1.0f;
                        odist[l] = fBestDistance;
                        ocdif[l] = fPixelDistance;
                        
                    }
            }
            
        }

            char str[1024];
            sprintf(str, "ssimilarity_disp%d_%d.tif",global_scale_index,global_window_index);
            oselfdisp.save(str);
            sprintf(str, "ssimilarity_dist%d_%d.tif",global_scale_index,global_window_index);
            oselfdist.save(str);
            sprintf(str, "ssimilarity_pdisp%d_%d.tif",global_scale_index,global_window_index);
            odisp.save(str);
            sprintf(str, "ssimilarity_pdist%d_%d.tif",global_scale_index,global_window_index);
            odist.save(str);

		 delete[] mean1;
         delete[] mean2;
	//! free translated images
         delete [] images2;
         delete [] images1;
        
	}
    
    
    
    
    
    
    
    ///////////////////////
    ///////////////////////  STROBE EFFECT
    ///////////////////////
    
	void stereo_find_local_minima(float *v, int n, float &dif1, float &dif0)
	{
		
		float dist0 = fLarge;
		float dist1 = fLarge;
		int nminima = 0;
		
		for (int i=1; i < n-1; i++)
			if (v[i-1] >= v[i] && v[i] <= v[i+1])
			{
				if (v[i] < dist0)
				{
					dist1 = dist0;
					dist0 = v[i];
					
				} else if (v[i] < dist1)
				{
					dist1 = v[i];
				}
				
				nminima++;
			}
		
		
		if (nminima > 1)
		{
			dif0 =  dist0;
			dif1 =  dist1;
		}
		else if (nminima == 1)	// Single minima
		{
			dif1 = fLarge;
			dif0 = dist0;
			
		} else	// Absolut minima is at boundary
		{
			dif1 = fLarge;
			dif0 = fLarge;
		}
		
		
	}
	
    
    ////////////////////////////
    ////////////////////////////
    //////////////////////////// MultiScale
    ////////////////////////////
    ////////////////////////////
    
        void combine_scales ( flimage &out_hr, flimage &mask_hr, float factor, flimage &out_low, flimage &mask_low , flimage &out, flimage &mask, float threshold)  {
        
           int w = out_hr.w();
           int h = out_hr.h();
           int lw = out_low.w();
           int lh = out_low.h();
           for(int y=0;y<h;y++) {
              for(int x=0;x<w;x++) {
                 if( mask_hr[x+w*y] == 0 ) {   // if the pixel is rejected at the current scale
                    int lx=x/factor;              // corresponding piels at the previous scale
                    int ly=y/factor;
                    float minv = INFINITY;
                    float maxv = -INFINITY;
                    float average = 0;
                    int count = 0;
                    // scan the 4-neighborhood of the pixel in the previous scale
                    for (int j=0;j<factor;j++){
                       for (int i=0;i<factor;i++){
                          if ( lx+i>=0 && ly+j>=0 && lx+i<lw  && ly+j<lh)  {
                             int pos = lx+i+lw*(ly+j);
                             if (mask_low[pos] == 0) {
                                maxv = INFINITY;      // invalidate the point
                             } else {
                                maxv = fmax(maxv,out_low[pos]);
                                minv = fmin(minv,out_low[pos]);
                                average += out_low[pos];
                                count ++;
                             }
                          }
                       } 
                    } 
                    if (count>0) average = average/count;
                    // if the value is stable then copy it to the current scale
                    if ( (maxv - minv)*factor < threshold && count >= 4 ) {
                       out[x+w*y] = average*factor;
                       mask[x+w*y] = 1;
                    }

                 }
              }
           }
        
        
        }

    void stereo_pixel_multiscale_chain_recursive(cflimage &input, cflimage &input2, int &in, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &out, flimage &out2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar, flimage scale_output[], flimage scale_mask[], flimage scale_output2[], flimage scale_mask2[], flimage scale_win[])
    {
        
        
        // initial scale is 0
        int currentScale = in;
        
        
        if (in < strPar.nScales)
        {
            
            // subsample by factor 2: input, input2, dmin and dmax
            float fGaussianConvol = 0.8f;
            cflimage sinput  = input.subSample(2, fGaussianConvol,1);
            cflimage sinput2 = input2.subSample(2, fGaussianConvol,1);
            
            flimage sDmin = Dmin.subSample(2, fGaussianConvol,1);
            flimage sDmax = Dmax.subSample(2, fGaussianConvol,1);
            
            flimage siDmin = iDmin.subSample(2, fGaussianConvol,1);
            flimage siDmax = iDmax.subSample(2, fGaussianConvol,1);
            
            
            // multiply dmin and dmax by 0.5
            sDmin *= 0.5f;   sDmax *= 0.5f;
            siDmin *= 0.5f;   siDmax *= 0.5f;
            
            
            // compute multiscale
            in++;
            
            flimage somask(sinput.w(), sinput.h());
            flimage somask2(sinput2.w(), sinput2.h());
            flimage sout(sinput.w(), sinput.h());
            flimage sout2(sinput2.w(), sinput2.h());
            
            flimage sodist(sinput.w(), sinput.h());
            flimage sodistr(sinput2.w(), sinput2.h());
            
            
            
            stereo_pixel_multiscale_chain_recursive(sinput, sinput2, in, sDmin, sDmax, siDmin, siDmax, sout, sout2, sodist, sodistr, somask, somask2, strPar, scale_output, scale_mask, scale_output2, scale_mask2, scale_win);
            
            
            
            // sDmin and sDmax could be updated here instead of inside stereo_pixel_chain
            flimage aDmin = sDmin.upSampleSplines(2.0f, 0);  aDmin *= 2.0f; aDmin -= 2.0f;
            flimage aDmax = sDmax.upSampleSplines(2.0f, 0);  aDmax *= 2.0f; aDmax += 2.0f;
            
            flimage aiDmin = siDmin.upSampleSplines(2.0f, 0);  aiDmin *= 2.0f; aiDmin -= 2.0f;
            flimage aiDmax = siDmax.upSampleSplines(2.0f, 0);  aiDmax *= 2.0f; aiDmax += 2.0f;
            
            
            // possible incoherence since upsampled range can have different size than the original one because of
            // quantization to integer dimensions. Pad the image with zeros
            //Dmin = aDmin.padding(input.w(), input.h(), 0.0f);
            //Dmax = aDmax.padding(input.w(), input.h(), 0.0f);
            //iDmin = aiDmin.padding(input2.w(), input2.h(), 0.0f);
            //iDmax = aiDmax.padding(input2.w(), input2.h(), 0.0f);

            flimage tDmin = aDmin.padding(input.w(), input.h(), 0.0f);
            flimage tDmax = aDmax.padding(input.w(), input.h(), 0.0f);

            flimage tiDmin = aiDmin.padding(input2.w(), input2.h(), 0.0f);
            flimage tiDmax = aiDmax.padding(input2.w(), input2.h(), 0.0f);

            // the updated range may not exceed the input one
            for(int i=0;i<input.wh();i++){
               Dmin[i] = fmax(Dmin[i],tDmin[i]);
               Dmax[i] = fmin(Dmax[i],tDmax[i]);
            }
            for(int i=0;i<input2.wh();i++){
               iDmin[i] = fmax(iDmin[i],tiDmin[i]);
               iDmax[i] = fmin(iDmax[i],tiDmax[i]);
            }

        }
        
        
        
        // build parameter structure for current scale
        strParameters csPar;
        set_strParameters_for_current_scale(strPar, csPar, currentScale);



        // change the precision for the scale for combine_scales
        if(strPar.flagCombineLastScale && currentScale==2) csPar.inPrecisions = csPar.inPrecisions*2;


        // for lower scales it is always better to have a subpixel estimate (and it is very cheap)
        if (DONT_USE_SUBPIX_AT_LOW_SCALES()==0)  {
           if (currentScale>1 && csPar.inPrecisions<4 ) csPar.inPrecisions = 4;
        }
        else { 
           printf("USING THE OBSCURE PARAMETER: DONT_USE_SUBPIX_AT_LOW_SCALES=%d ", DONT_USE_SUBPIX_AT_LOW_SCALES());
        }
        
   global_scale_index = currentScale;
        // apply stereo chain: dmin, dmax, idmin and idmax are updated inside stereo_pixel_chain
        if(strPar.flagMultiWin) 
           stereo_pixel_chain_multi_window(input,input2,Dmin,Dmax, iDmin, iDmax, out, out2, odist, odistr, omask, omask2, csPar);
        else
           stereo_pixel_chain(input,input2,Dmin,Dmax, iDmin, iDmax, out, out2, odist, odistr, omask, omask2, csPar);

        scale_output[currentScale] = out;
        scale_output2[currentScale] = out2;
        scale_mask[currentScale] = omask;
        scale_mask2[currentScale] = omask2;
        // scale_win
        
        // combine latest scale with the previous scale
        if (strPar.flagCombineLastScale && currentScale ==1 && strPar.nScales > 1 ) {
           printf("combining %d and %d\n", currentScale, currentScale+1);
           combine_scales ( scale_output[currentScale], scale_mask[currentScale], 2.0, scale_output[currentScale+1], scale_mask[currentScale+1], out, omask, 1);
           combine_scales ( scale_output2[currentScale], scale_mask2[currentScale], 2.0, scale_output2[currentScale+1], scale_mask2[currentScale+1], out2, omask2, 1);


		      //! Check left right consistency
		      //! reciprocity flag
            // THIS WORSEN THE RESULTS
		      if(1){
		      	
		      	libIIPStable::flimage oBin(input.w(), input.h());     oBin = 255.0f;
		      	libIIPStable::flimage oTMP(input.w(), input.h());     oTMP = fLarge;
		      	
		      	libIIPStable::flimage oBin2(input2.w(), input2.h());     oBin2 = 255.0f;
		      	libIIPStable::flimage oTMP2(input2.w(), input2.h());     oTMP2 = fLarge;
		      	
		      	//stereo_check_pixelian_reciprocity(out, omask, out2, omask2, oBin, strPar.valueRecip, oTMP);
		      	//stereo_check_pixelian_reciprocity(out2, omask2, out, omask, oBin2, strPar.valueRecip, oTMP2);
		      	stereo_check_pixelian_reciprocity_remove_only_if_inconsistent(out, omask, out2, omask2, oBin, strPar.valueRecip, oTMP);
		      	stereo_check_pixelian_reciprocity_remove_only_if_inconsistent(out2, omask2, out, omask, oBin2, strPar.valueRecip, oTMP2);
		      	
		      	
		      	for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i]*10;      // TAG THESE POINTS
		      	for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i]*10;
		      	
		      }
        }



    }
    
    
    
    
    
    void stereo_pixel_multiscale_chain(cflimage &input, cflimage &input2, int &in, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &out, flimage &out2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar) {


      // make an array of images and masks for the intermediate scales
      libIIPStable::flimage *scale_output;
      scale_output = new libIIPStable::flimage[strPar.nScales+1];
      for (int i=0;i<strPar.nScales+1;i++) scale_output[i].create(1, 1);

      libIIPStable::flimage *scale_mask;
      scale_mask = new libIIPStable::flimage[strPar.nScales+1];
      for (int i=0;i<strPar.nScales+1;i++) scale_mask[i].create(1, 1);

      libIIPStable::flimage *scale_output2;
      scale_output2 = new libIIPStable::flimage[strPar.nScales+1];
      for (int i=0;i<strPar.nScales+1;i++) scale_output2[i].create(1, 1);

      libIIPStable::flimage *scale_mask2;
      scale_mask2 = new libIIPStable::flimage[strPar.nScales+1];
      for (int i=0;i<strPar.nScales+1;i++) scale_mask2[i].create(1, 1);

      libIIPStable::flimage *scale_win;
      scale_win = new libIIPStable::flimage[strPar.nScales+1];
      for (int i=0;i<strPar.nScales+1;i++) scale_win[i].create(1, 1);
       
      stereo_pixel_multiscale_chain_recursive(input, input2, in, Dmin, Dmax, iDmin, iDmax, out, out2, odist, odistr, omask, omask2, strPar, scale_output, scale_mask, scale_output2, scale_mask2, scale_win);

      delete[] scale_output;
      delete[] scale_mask;
      delete[] scale_output2;
      delete[] scale_mask2;
      delete[] scale_win;
    }
    
    
    
    void set_strParameters_for_current_scale(strParameters & strIn, strParameters & strOut, int iiScale)
    {
        
        
        
        //! Type of distance and window size
        strOut.itypeDist = strIn.itypeDist;
        strOut.flagWinWeighted = strIn.flagWinWeighted;
        strOut.prolate = strIn.prolate;
        
        
        //! Number of precisions used for single scale pixelian estimation
		strOut.inPrecisions = strIn.inPrecisions;
		strOut.cubicRefinementFactor = strIn.cubicRefinementFactor;
		
        
        //! Number of scales
        strOut.nScales = strIn.nScales;
        strOut.currentScale = iiScale;
        float scaleAux = powf(2.0, (float) iiScale - 1);
        
        //! Range min and max
        strOut.dmin0 = strIn.dmin0 / scaleAux;
        strOut.idmin0 = strIn.idmin0 / scaleAux;
        
        strOut.dmax0 = strIn.dmax0 / scaleAux;
        strOut.idmax0 = strIn.idmax0 / scaleAux;
        
        
        strOut.flagSubRecip = strIn.flagSubRecip;
        strOut.valueSubRecip = strIn.valueSubRecip;
       
        
        strOut.flagSelfSimilarity = strIn.flagSelfSimilarity;
        strOut.valueSelfSimilarity = strIn.valueSelfSimilarity;
        
        
        strOut.flagRemoveIsolated = strIn.flagRemoveIsolated;
        strOut.valueRemoveIsolated = strIn.valueRemoveIsolated;
        
        // reduce the size of the grain at lower scales 
        strOut.flagGrainFilter = strIn.flagGrainFilter;
        strOut.valueGrainFilter = strIn.valueGrainFilter/scaleAux;

        strOut.flagCostFilter = strIn.flagCostFilter;
        strOut.valueCostFilter = strIn.valueCostFilter/scaleAux;
        
        strOut.flagHorInfo = strIn.flagHorInfo;
        strOut.valueHorInfo = strIn.valueHorInfo;
        
        
        //! Noise standard deviation
        float scaleNoiseAux = powf(2.0, (float) iiScale - 1);
        strOut.fSigma = strIn.fSigma / scaleNoiseAux;
        
        
        
        strOut.flagRecip   =  strIn.flagRecip;
        strOut.valueRecip   =  strIn.valueRecip;
        strOut.flagRecipDual   =  strIn.flagRecipDual;
        strOut.flagiRecip  =  strIn.flagiRecip;
        strOut.flagMultiWin =  strIn.flagMultiWin;
        strOut.flagListKernels= strIn.flagListKernels;
		
        
        
        strOut.flagVariance = strIn.flagVariance;
        strOut.valueVariance = strIn.valueVariance;
        
        strOut.flagMinDist =  strIn.flagMinDist;
        strOut.valueMinDist =  strIn.valueMinDist;
        strOut.flagMinDistDilatation = strIn.flagMinDistDilatation;
        strOut.flagRegressionMinDist = strIn.flagRegressionMinDist;
        
        strOut.flagCombineLastScale=  strIn.flagCombineLastScale;
        
    }
    
    
    
	
	// multi window inside multiscale chain:
	// this function implements multi window for a given scale
	void stereo_pixel_chain_multi_window(cflimage &input, cflimage &input2, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &odisp, flimage &odisp2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar)
	{
		
		
		// GENERATE THE KERNELS 
		// THESE SHOULD BE LIST KERNELS  -- 

      int win= strPar.prolate.w(); //window size
      // number of window orientations
      int nprol = 5; //  by default 5 orientations
      nprol=strPar.flagListKernels; // 5 or 9

        
      // the maximum number of window orientations we will ever consider is 9
      libIIPStable::flimage *prolate;
      prolate = new libIIPStable::flimage[10];

      // this is the case of the old non-list windows which are simpler to maintain
      // lists are only coded for 5x5 windows 
      // however lits a slightly faster 
      if( win!=5 || strPar.flagListKernels == 0)  
//      if(1) // always run without lists (easyer life)
      {
         nprol = fmax(strPar.flagListKernels,1); // nprol can't be 0
         strPar.flagListKernels = 0;   // disable the list kernels


         prolate[0].create(win, win); prolate[0] = 1.0f; prolate[0].normalizeL1();

         if (win == 3)       //! 3x3
         {
            prolate[1].create(9,1); prolate[1] = 1.0; prolate[1].normalizeL1();
            prolate[2].create(1,9); prolate[2] = 1.0; prolate[2].normalizeL1();

//            stereo_diagonal_flat_window(9, 0, 0, prolate[3]);  prolate[3].normalizeL1();
//            stereo_diagonal_flat_window(9, 0, 1, prolate[4]);  prolate[4].normalizeL1();
            nprol=fmin(nprol,3);



         }else if (win == 5)   //! 5x5
         {
            prolate[1].create(9,3); prolate[1] = 1.0; prolate[1].normalizeL1();
            prolate[2].create(3,9); prolate[2] = 1.0; prolate[2].normalizeL1();

            stereo_diagonal_flat_window(9, 1, 0, prolate[3]);  prolate[3].normalizeL1();
            stereo_diagonal_flat_window(9, 1, 1, prolate[4]);  prolate[4].normalizeL1();

            // crazy experiment
            //        prolate[5].create(1,25); prolate[5] = 1.0; prolate[5].normalizeL1();
            //        prolate[6].create(25,1); prolate[6] = 1.0; prolate[6].normalizeL1();

            // 25.5 degrees windows
            int w5[81] = { 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           1 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           1 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           0 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,1 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 };
            int w6[81] = { 0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 };
            int w7[81] = { 0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 };
            int w8[81] = { 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,
                           0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,1 ,
                           0 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0 ,
                           1 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,
                           1 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                           0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 };
            prolate[5].create(9,9); for(int i=0;i<81;i++) prolate[5][i] = w5[i]; prolate[5].normalizeL1();
            prolate[6].create(9,9); for(int i=0;i<81;i++) prolate[6][i] = w6[i]; prolate[6].normalizeL1();
            prolate[7].create(9,9); for(int i=0;i<81;i++) prolate[7][i] = w7[i]; prolate[7].normalizeL1();
            prolate[8].create(9,9); for(int i=0;i<81;i++) prolate[8][i] = w8[i]; prolate[8].normalizeL1();
            nprol=fmin(nprol,9);

         } else if (win == 7) // remove
         {

            prolate[1].create(11,5); prolate[1] = 1.0; prolate[1].normalizeL1();
            prolate[2].create(5,11); prolate[2] = 1.0; prolate[2].normalizeL1();

            stereo_diagonal_flat_window(11, 2, 0, prolate[3]);  prolate[3].normalizeL1();
            stereo_diagonal_flat_window(11, 2, 1, prolate[4]);  prolate[4].normalizeL1();
            nprol=fmin(nprol,5);


            //// crazy experimento
            //prolate[5].create(3,17); prolate[5] = 1.0; prolate[5].normalizeL1();
            //prolate[6].create(17,3); prolate[6] = 1.0; prolate[6].normalizeL1();

            //stereo_diagonal_flat_window(17, 1, 0, prolate[7]);  prolate[7].normalizeL1();
            //stereo_diagonal_flat_window(17, 1, 1, prolate[8]);  prolate[8].normalizeL1();


         } else if (win == 9)
         {

            prolate[1].create(15,5); prolate[1] = 1.0; prolate[1].normalizeL1();
            prolate[2].create(5,15); prolate[2] = 1.0; prolate[2].normalizeL1();

            stereo_diagonal_flat_window(15, 2, 0, prolate[3]);  prolate[3].normalizeL1();
            stereo_diagonal_flat_window(15, 2, 1, prolate[4]);  prolate[4].normalizeL1();
            nprol=fmin(nprol,5);


         } else if (win == 13)
         {

            prolate[1].create(21,7); prolate[1] = 1.0; prolate[1].normalizeL1();
            prolate[2].create(7,21); prolate[2] = 1.0; prolate[2].normalizeL1();

            stereo_diagonal_flat_window(21, 3, 0, prolate[3]);  prolate[3].normalizeL1();
            stereo_diagonal_flat_window(21, 3, 1, prolate[4]);  prolate[4].normalizeL1();
            nprol=fmin(nprol,5);


         } else if (win == 17)
         {

            prolate[1].create(29,9); prolate[1] = 1.0; prolate[1].normalizeL1();
            prolate[2].create(9,29); prolate[2] = 1.0; prolate[2].normalizeL1();

            stereo_diagonal_flat_window(29, 4, 0, prolate[3]);  prolate[3].normalizeL1();
            stereo_diagonal_flat_window(29, 4, 1, prolate[4]);  prolate[4].normalizeL1();
            nprol=fmin(nprol,5);

         } else if (win == 21)
         {

            prolate[1].create(33,13); prolate[1] = 1.0; prolate[1].normalizeL1();
            prolate[2].create(13,33); prolate[2] = 1.0; prolate[2].normalizeL1();

            stereo_diagonal_flat_window(33, 6, 0, prolate[3]);  prolate[3].normalizeL1();
            stereo_diagonal_flat_window(33, 6, 1, prolate[4]);  prolate[4].normalizeL1();
            nprol=fmin(nprol,5);

         } else
         {
            printf("error :: please specify a correct window size (3,5,7,9,13,17,21) \n"); abort();
         }

         // make sure the lists are disabled
         for (int f=0; f<nprol;f++) {
            prolate[f].list_len = 0;
            prolate[f].offx     = NULL;
            prolate[f].offy     = NULL;
            prolate[f].offval   = NULL;
         }


      }
      else
      {

         //square
         float w5_0[][25] = {
            {-2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2},
            {-2, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2},
            {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04} 
         };

         // 90 degrees
         float w5_1[][27] = {
            {-4, -3, -2, -1, 0, 1, 2, 3, 4, -4, -3, -2, -1, 0, 1, 2, 3, 4, -4, -3, -2, -1, 0, 1, 2, 3, 4}, 
            { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            { 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037}
         };
         float w5_2[][27] = {
            {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1}, 
            {-4, -4, -4, -3, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4},
            {0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037, 0.037037} 
         };

         // 45 degrees
         float w5_3[][25] = {
            {-4, -3, -4, -3, -2, -3, -2, -1, -2, -1, 0, -1, 0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4}, 
            {-4, -4, -3, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4}, 
            {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04}
         };
         float w5_4[][25] = {
            {3, 4, 2, 3, 4, 1, 2, 3, 0, 1, 2, -1, 0, 1, -2, -1, 0, -3, -2, -1, -4, -3, -2, -4, -3},
            {-4, -4, -3, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
            {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04}
         };

         // 22.5 degrees
         float w5_5[][25] = {
            {-4, -3, -2, -4, -3, -2, -1, 0, 1, -3, -2, -1, 0, 1, 2, 3, -1, 0, 1, 2, 3, 4, 2, 3, 4},
            {-2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2}, 
            {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04} 
         };
         float w5_6[][25] = {
            {-2, -1, -2, -1, 0, -2, -1, 0, -1, 0, 1, -1, 0, 1, -1, 0, 1, 0, 1, 2, 0, 1, 2, 1, 2}, 
            {-4, -4, -3, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4}, 
            {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04}, 
         };
         float w5_7[][25] = {
            {1, 2, 0, 1, 2, 0, 1, 2, -1, 0, 1, -1, 0, 1, -1, 0, 1, -2, -1, 0, -2, -1, 0, -2, -1}, 
            {-4, -4, -3, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4}, 
            {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04}, 
         };
         float w5_8[][25] = {
            {2, 3, 4, -1, 0, 1, 2, 3, 4, -3, -2, -1, 0, 1, 2, 3, -4, -3, -2, -1, 0, 1, -4, -3, -2}, 
            {-2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2}, 
            {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04}, 
         };


         int   w5_lens[] = {25,27,27, 25,25,25, 25,25,25};

         // for lanted windows with ~25^2 pixels the support is larger than 5
         win=9;


         // this is loading the windows as lists
         float* w5[9];
         w5[0]=(float*)w5_0; w5[1]=(float*)w5_1; w5[2]=(float*)w5_2;
         w5[3]=(float*)w5_3; w5[4]=(float*)w5_4; w5[5]=(float*)w5_5;
         w5[6]=(float*)w5_6; w5[7]=(float*)w5_7; w5[8]=(float*)w5_8;

         for(int f=0;f<nprol;f++){
            prolate[f].create(win,win);
            prolate[f] = 0;
            prolate[f].list_len = w5_lens[f];
            prolate[f].offx     = new int[prolate[f].list_len];
            prolate[f].offy     = new int[prolate[f].list_len];
            prolate[f].offval   = new float[prolate[f].list_len];
            for (int i=0;i<prolate[f].list_len;i++)  {
               int spd_w = win/2;
               int spd_h = win/2;
               prolate[f].offx[i]   = (int) w5[f][0*w5_lens[f] + i]+spd_w;
               prolate[f].offy[i]   = (int) w5[f][1*w5_lens[f] + i]+spd_h;
               prolate[f].offval[i] =       w5[f][2*w5_lens[f] + i];
               // also copy the kernel as a patch
               prolate[f]( prolate[f].offx[i] , 
                     prolate[f].offy[i] ) = w5[f][2*w5_lens[f] + i];
            }
            char stri[1024];
            sprintf(stri, "p%d.tif", f);
            prolate[f].save(stri);
         }

      }


      // Memory for disparity and mask of selected pixel taking left image as reference
      libIIPStable::flimage *tmpodisp = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *tmpodisp2 = new  libIIPStable::flimage[nprol];

      libIIPStable::flimage *tmpomask = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *tmpomask2 = new  libIIPStable::flimage[nprol];

      libIIPStable::flimage *tmpodist = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *tmpodist2 = new  libIIPStable::flimage[nprol];


      for (int ii=0; ii < nprol; ii++)
      {

         tmpodisp[ii].create(input.w(), input.h());         tmpodisp[ii] = 0.0f;
         tmpodisp2[ii].create(input2.w(), input2.h());      tmpodisp2[ii] = 0.0f;

         tmpomask[ii].create(input.w(), input.h());         tmpomask[ii] = 0.0f;
         tmpomask2[ii].create(input2.w(), input2.h());      tmpomask2[ii] = 0.0f;

         tmpodist[ii].create(input.w(), input.h());         tmpodist[ii] = fLarge;
         tmpodist2[ii].create(input2.w(), input2.h());      tmpodist2[ii] = fLarge;

      }


      /*libIIPStable::flimage tmpDmin(input.w(), input.h());     tmpDmin = Dmin; //Dmin =  fLarge;
        libIIPStable::flimage tmpDmax(input.w(), input.h());     tmpDmax = Dmax; //Dmax = -fLarge;
        libIIPStable::flimage tmpiDmin(input2.w(), input2.h());    tmpiDmin = iDmin; //iDmin =  fLarge;
        libIIPStable::flimage tmpiDmax(input2.w(), input2.h());    tmpiDmax = iDmax; //iDmax = -fLarge;*/

      for (int ii=0; ii < nprol; ii++)
      {
         //! change the prolate
         strPar.prolate = prolate[ii];
         printf("multiwindow each scale --> w: %d h: %d l: %d s: %d \n", strPar.prolate.w(), strPar.prolate.h(), strPar.prolate.list_len, strPar.currentScale);
         strPar.flagWinWeighted = 1;

         global_window_index=ii;
			
			libIIPStable::flimage ttmpDmin(input.w(), input.h());     ttmpDmin = Dmin;//tmpDmin;
			libIIPStable::flimage ttmpDmax(input.w(), input.h());     ttmpDmax = Dmax;//tmpDmax;
			libIIPStable::flimage ttmpiDmin(input2.w(), input2.h());    ttmpiDmin = iDmin;//tmpiDmin;
			libIIPStable::flimage ttmpiDmax(input2.w(), input2.h());    ttmpiDmax = iDmax; //tmpiDmax;
			
			stereo_pixel_chain(input, input2, ttmpDmin, ttmpDmax, ttmpiDmin, ttmpiDmax, tmpodisp[ii], tmpodisp2[ii], tmpodist[ii], tmpodist2[ii], tmpomask[ii], tmpomask2[ii], strPar);
			//        stereo_pixel_chain(input, input2, Dmin, Dmax, iDmin, iDmax, tmpodisp[ii], tmpodisp2[ii], tmpodist[ii], tmpodist2[ii], tmpomask[ii], tmpomask2[ii], strPar);
			
			/*
			 for (int iidx=0; iidx < Dmin.w()*Dmin.h() ; iidx++) {
			 Dmin[iidx] =fmin(Dmin[iidx], ttmpDmin[iidx]);
			 Dmax[iidx] =fmax(Dmax[iidx], ttmpDmax[iidx]);
			 }
			 for (int iidx=0; iidx < iDmin.w()*iDmin.h() ; iidx++) {
			 iDmin[iidx] =fmin(iDmin[iidx], ttmpiDmin[iidx]);
			 iDmax[iidx] =fmax(iDmax[iidx], ttmpiDmax[iidx]);
			 }
			 */
			//           char stri[1024];
			//           sprintf(stri, "Dmin%d.tif", ii);
			//           Dmin.save(stri);
		}
		
		
		
		//! Computing choice and additional criteria
		odist = fLarge;
		odisp = 0.0f;
		omask = 0.0f;
		
		odistr = fLarge;
		odisp2 = 0.0f;
		omask2 = 0.0f;
		
		libIIPStable::flimage oChoice(input.w(), input.h());      oChoice = -1.0f;

      // needed for costFilter
      libIIPStable::flimage costFilter_mindist(odist);
      libIIPStable::flimage costFilter_mindistr(odistr);

      //! Select minimum distance
      for (int ii=0; ii < nprol; ii++)
      {


         //// ADVANCED SELECTION OF MINIMUM DISTANCE CONSISTENT WITH LEFT AND RIGHT MAPS
         if(USE_LR_CONSISTENT_SELECTION_OF_WINDOW_ORIENTATION() )  {

            for(int jj=0; jj < tmpodisp[ii].wh(); jj++) {
               float candidatedist = tmpodist[ii][jj];
               float candidatedisp = tmpodisp[ii][jj];
               float candidatemask = tmpomask[ii][jj];

               if ( candidatemask > 0.0f ) {
                  int jjRev = jj +  (int) rintf(tmpodisp[ii][jj]);
                  float candidatedist2 = tmpodist2[ii][jjRev];
                  float candidatedisp2 = tmpodisp2[ii][jjRev];
                  float candidatemask2 = tmpomask2[ii][jjRev];

                  if ( candidatemask2 > 0.0f && (candidatedist + candidatedist2)/2.0 < odist[jj] ) {
                     odist[jj] = (candidatedist+candidatedist2)/2.0;
                     odisp[jj] = candidatedisp;
                     omask[jj] = candidatemask;
                     oChoice[jj] = (float) ii;

                     odistr[jjRev] = (candidatedist+candidatedist2)/2.0;
                     odisp2[jjRev] = candidatedisp2;
                     omask2[jjRev] = candidatemask2;

                     costFilter_mindist[jj] = (candidatedist+candidatedist2)/2.0;
                  } else {
                     if(omask[jj] <=0 && candidatedist < costFilter_mindist[jj])
                        costFilter_mindist[jj] = candidatedist;
                  }

               } else { 
                  if(omask[jj] <=0 && candidatedist < costFilter_mindist[jj])
                     costFilter_mindist[jj] = candidatedist;
               }
            }


            for(int jj=0; jj < tmpodisp2[ii].wh(); jj++) {
               float candidatedist2 = tmpodist2[ii][jj];
               float candidatedisp2 = tmpodisp2[ii][jj];
               float candidatemask2 = tmpomask2[ii][jj];

               if (  candidatemask2 > 0.0f ) {
                  int jjRev = jj +  (int) rintf(tmpodisp2[ii][jj]);
                  float candidatedist = tmpodist[ii][jjRev];
                  float candidatedisp = tmpodisp[ii][jjRev];
                  float candidatemask = tmpomask[ii][jjRev];

                  if ( candidatemask > 0.0f && (candidatedist2 + candidatedist)/2.0 < odistr[jj]) {
                     odistr[jj] = (candidatedist2 + candidatedist)/2.0;
                     odisp2[jj] = candidatedisp2;
                     omask2[jj] = candidatemask2;

                     odist[jjRev] = (candidatedist2 + candidatedist)/2.0;
                     odisp[jjRev] = candidatedisp;
                     omask[jjRev] = candidatemask;
                     oChoice[jjRev] = (float) ii;

                     costFilter_mindistr[jj] = (candidatedist+candidatedist2)/2.0;
                  } else {
                     if(omask2[jj] <=0 && candidatedist2 < costFilter_mindistr[jj])
                        costFilter_mindistr[jj] = candidatedist2;
                  }
               } else {
                  if(omask2[jj] <=0 && candidatedist2 < costFilter_mindistr[jj])
                     costFilter_mindistr[jj] = candidatedist2;
               }
            }



         } else {
            // JUST TAKE THE MINIMUM DISTANCE FOR EACH PIXEL
            for(int jj=0; jj < tmpodisp[ii].wh(); jj++)
               if (  tmpomask[ii][jj] > 0.0f  && tmpodist[ii][jj] < odist[jj])
               { // do not update dist and disp unless the pixel is unmasked
                  odist[jj] = tmpodist[ii][jj];
                  odisp[jj] = tmpodisp[ii][jj];
                  omask[jj] = tmpomask[ii][jj];
                  oChoice[jj] = (float) ii;
                  costFilter_mindist[jj] = tmpodist[ii][jj];
               } else if(omask[jj] <=0 && tmpodist[ii][jj] < costFilter_mindist[jj]){
                  costFilter_mindist[jj] = tmpodist[ii][jj];
               }


            for(int jj=0; jj < tmpodisp2[ii].wh(); jj++)
               if (  tmpomask2[ii][jj] > 0.0f  && tmpodist2[ii][jj] < odistr[jj])
               {
                  odistr[jj] = tmpodist2[ii][jj];
                  odisp2[jj] = tmpodisp2[ii][jj];
                  omask2[jj] = tmpomask2[ii][jj];
                  //oChoice[jj] = (float) ii;
                  costFilter_mindistr[jj] = tmpodist2[ii][jj];
               } else if(omask2[jj] <=0 && tmpodist2[ii][jj] < costFilter_mindistr[jj]){
                  costFilter_mindistr[jj] = tmpodist2[ii][jj];
               }


         }


      }



		//! Check left right consistency
		//! reciprocity flag
		if (strPar.flagRecip && strPar.valueRecip >= 0.0f)
		{
			
			libIIPStable::flimage oBin(input.w(), input.h());     oBin = 255.0f;
			libIIPStable::flimage oTMP(input.w(), input.h());     oTMP = fLarge;
			
			libIIPStable::flimage oBin2(input2.w(), input2.h());     oBin2 = 255.0f;
			libIIPStable::flimage oTMP2(input2.w(), input2.h());     oTMP2 = fLarge;
			
		   // it's the second time we're applying LR in this pipeline, to prevent 
         // further eliminations we remove only matches that are explicitly inconsistent	
         if (USE_OLD_MERGE_LR()) {
			stereo_check_pixelian_reciprocity(odisp, omask, odisp2, omask2, oBin, strPar.valueRecip, oTMP);
			stereo_check_pixelian_reciprocity(odisp2, omask2, odisp, omask, oBin2, strPar.valueRecip, oTMP2);
         } else {
         stereo_check_pixelian_reciprocity_remove_only_if_inconsistent(odisp, omask, odisp2, omask2, oBin, strPar.valueRecip, oTMP);
         stereo_check_pixelian_reciprocity_remove_only_if_inconsistent(odisp2, omask2, odisp, omask, oBin2, strPar.valueRecip, oTMP2);
         }
			
			
			for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
			for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
			
		}
		
		
		

		//! Remove isolated points
		if (strPar.flagRemoveIsolated)
		{
			
			float valueRemoveIsolated;
			
			
			libIIPStable::flimage oBin(input.w(), input.h());     oBin = 255.0f;
			libIIPStable::flimage oBin2(input2.w(), input2.h());     oBin2 = 255.0f;
			
			
         // THIS ONE IS WITH SQUARED WINDOWS BY DESIGN! WITH SIZE DOUBLE OF THE WINDOW
			int pwin = strPar.prolate.w()/2;
         pwin = strPar.prolate.w();
			stereo_remove_isolated_points(omask, pwin, oBin, strPar.valueRemoveIsolated);
			stereo_remove_isolated_points(omask2, pwin, oBin2, strPar.valueRemoveIsolated);
			
			for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
			for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
			
		}
		

        // Perform Grain Filter
        if (strPar.flagGrainFilter)
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            int minarea = strPar.valueGrainFilter;

            stereo_grain_filter(omask, minarea, oBin);
            stereo_grain_filter(omask2, minarea, oBin2);
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }


        // Filter: Cost histogram analysis 
        if (strPar.flagCostFilter)
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            // default 0.05 
            float discard_threshold = strPar.valueCostFilter; 

            stereo_cost_filter(costFilter_mindist, omask, discard_threshold, oBin);
            stereo_cost_filter(costFilter_mindistr, omask2, discard_threshold, oBin2);
            
            oBin.save("debug_cost_histogram0.tif");
            omask.save("debug_cost_histogram1.tif");

            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];

            omask.save("debug_cost_histogram2.tif");
            
        }





		// TODO:  DEBUG REMOVE THIS
		char stri[1024];
		 sprintf(stri, "win%d.tif", strPar.currentScale);
		oChoice.save(stri);
		 sprintf(stri, "disp%d.tif", strPar.currentScale);
      odisp.save(stri);
		 sprintf(stri, "cost%d.tif", strPar.currentScale);
      odist.save(stri);
		 sprintf(stri, "mask%d.tif", strPar.currentScale);
      omask.save(stri);
		
		
		
      // Update min and max
      // MODIFICATION: this processing resembles more what is written in the description
      // compute it for each window and then combine into a signle Dmin/Dmax
		libIIPStable::flimage newDmin(input.w(), input.h());     newDmin =  INFINITY;//tmpDmin;
		libIIPStable::flimage newDmax(input.w(), input.h());     newDmax =  -INFINITY;//tmpDmax;
		libIIPStable::flimage newiDmin(input2.w(), input2.h());  newiDmin = INFINITY;//tmpiDmin;
		libIIPStable::flimage newiDmax(input2.w(), input2.h());  newiDmax = -INFINITY; //tmpiDmax;
		for (int ii=0; ii < nprol; ii++)
		{

			libIIPStable::flimage ttmpDmin(input.w(), input.h());     ttmpDmin = Dmin;//tmpDmin;
			libIIPStable::flimage ttmpDmax(input.w(), input.h());     ttmpDmax = Dmax;//tmpDmax;
			libIIPStable::flimage ttmpiDmin(input2.w(), input2.h());    ttmpiDmin = iDmin;//tmpiDmin;
			libIIPStable::flimage ttmpiDmax(input2.w(), input2.h());    ttmpiDmax = iDmax; //tmpiDmax;

		   update_dmin_dmax_window(ttmpDmin,ttmpDmax,odisp,omask, strPar.dmin0, strPar.dmax0, prolate[ii]);
		   update_dmin_dmax_window(ttmpiDmin,ttmpiDmax,odisp2,omask2, strPar.idmin0, strPar.idmax0, prolate[ii]);

			for (int i=0; i < Dmin.wh(); i++)  {
            newDmin[i] = fmin (newDmin[i] , ttmpDmin[i]);
            newDmax[i] = fmax (newDmax[i] , ttmpDmax[i]);
         }
			for (int i=0; i < iDmin.wh(); i++) {
            newiDmin[i] = fmin (newiDmin[i] , ttmpiDmin[i]);
            newiDmax[i] = fmax (newiDmax[i] , ttmpiDmax[i]);
         }
			
      }
      Dmin  = newDmin;
      Dmax  = newDmax;
      iDmin = newiDmin;
      iDmax = newiDmax;
		
		
//    OLD
//		// Update min and max
//		update_dmin_dmax(Dmin,Dmax,odisp,omask, strPar.dmin0, strPar.dmax0, strPar.prolate.w(), strPar.prolate.h());
//		update_dmin_dmax(iDmin,iDmax,odisp2,omask2, strPar.idmin0, strPar.idmax0, strPar.prolate.w(), strPar.prolate.h());
		
		
		delete[] tmpodisp;
		delete[] tmpodisp2;
		
		delete[] tmpomask;
		delete[] tmpomask2;
		
		delete[] tmpodist;
		delete[] tmpodist2;
		
      delete[] prolate;

	} 




    void update_dmin_dmax_window(flimage &Dmin, flimage &Dmax, flimage &out, flimage &omask, float minim, float maxim, flimage &prolate)
    {
        
        
        // using min and max on neighborhood
        int spwidth =  (prolate.w()-1)/2; //spwidth++;
        int spheight = (prolate.h()-1)/2; //spheight++;
        
        
        //int imin = (int) rintf((float) (2*spwidth+1) * (float) (2*spheight + 1) * 0.5f);
        
        
        for (int jj=spheight; jj < omask.h() - spheight; jj++)
            for (int ii=spwidth; ii < omask.w() - spwidth; ii++)
                
                
                if (omask[jj*omask.w() + ii] > 0.0f)
                {
                    
                    
                    float fmax = -fLarge;
                    float fmin = fLarge;
                    
                    int iCount = 0;
                    for (int ss=-spheight; ss <= spheight; ss++)
                        for (int rr=-spwidth; rr <= spwidth; rr++)
                        {
                            int irr = ii + rr;
                            int jss = jj + ss;
                            int iindex = jss*out.w()+irr;
                            
                            int prol_i = rr + spwidth; //prolate.w()/2;
                            int prol_j = ss + spheight;//prolate.h()/2;

                            assert(jss>=0);
                            assert(jss<out.h());
                            assert(irr>=0);
                            assert(irr<out.w());
                            assert(prol_i>=0); assert(prol_j>=0);
                            assert(prol_i<prolate.w()); assert(prol_j<prolate.h());
                            
                            if (prolate(prol_i, prol_j) > 0.0f && omask[iindex] > 0.0f)
                            {
                                if (out[iindex] < fmin) fmin = out[iindex];
                                
                                if (out[iindex] > fmax) fmax = out[iindex];
                                
                                iCount++;
                            }
                            
                            
                        }
                    
                    
                    
                    // if (iCount > imin)
                    //{
                    
                    
                    Dmin[jj*Dmin.w()+ii] = fmin;
                    Dmax[jj*Dmax.w()+ii] = fmax;
                    
                    //}
                    
                    
                    
                } else {
                    
                    Dmin[jj*omask.w() + ii] = minim;
                    Dmax[jj*omask.w() + ii] = maxim;
                    
                }
        
        
        
    }
    
    
    
    
    







	
    
    
    void stereo_pixel_chain(cflimage &input, cflimage &input2, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &odisp, flimage &odisp2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar)
    {
        

        
        //! Memory for disparity and mask of selected pixel taking left image as reference
        odist = fLarge;
        libIIPStable::flimage odist2(input.w(), input.h());    odist2 = fLarge;
        libIIPStable::flimage oselfdist(input.w(), input.h());    oselfdist = fLarge;
        libIIPStable::flimage olrdist(input.w(), input.h());    olrdist = fLarge;
        libIIPStable::flimage ocdist(input.w(), input.h());    ocdist = fLarge;
        
        {
            libIIPStable::flimage imask(input.w(), input.h()); imask=255.0f;
            stereo_pixel_chain_one_direction( input,  input2, imask,  Dmin, Dmax, odisp,  omask,  odist, odist2, oselfdist, olrdist, ocdist, strPar);
        }
        
        
        //! Memory for disparity and mask of selected pixel taking right image as reference
        odistr = fLarge;
        libIIPStable::flimage odistr2(input2.w(), input2.h());    odistr2 = fLarge;
        libIIPStable::flimage oselfdistr(input2.w(), input2.h());   oselfdistr = fLarge;
        libIIPStable::flimage olrdistr(input2.w(), input2.h());    olrdistr = fLarge;
        libIIPStable::flimage ocdistr(input2.w(), input2.h());    ocdistr = fLarge;
        {
            libIIPStable::flimage imask(input2.w(), input2.h()); imask=255.0f;
            stereo_pixel_chain_one_direction( input2,  input, imask,  iDmin, iDmax, odisp2,  omask2,  odistr, odistr2, oselfdistr, olrdistr, ocdistr, strPar);
        }
        
        
        
        //! Save original masks
        libIIPStable::flimage omask0 =omask;
        libIIPStable::flimage omask20 =omask2;
        
        
        //! Perform Reciprocal checking if flag activated
        if (strPar.flagRecip)
        {
            
            float fThreshold = strPar.valueRecip;  // to assure pixelian precision
            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
            
            flimage oTMP(input.w(), input.h()); 	oTMP = fLarge;
            flimage oTMP2(input2.w(), input2.h()); 	oTMP2 = fLarge;
            
            stereo_check_pixelian_reciprocity(odisp, omask0, odisp2, omask20, oBin, fThreshold, oTMP);
            stereo_check_pixelian_reciprocity(odisp2, omask20, odisp, omask0, oBin2, fThreshold, oTMP2);
            
            //oBin.save("st_reciprocity.pmf");
            
            char str[1024];
            sprintf(str, "st_reciprocity%d_%d.tif",global_scale_index, global_window_index);
            oBin.save(str);
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        
        
        
        
        //! Perform Reciprocal checking if flag activated
        if (strPar.flagRecipDual)
        {
            
            float fThreshold = 1.25;  // to assure pixelian precision
            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
            
            stereo_check_pixelian_reciprocity_dual(odisp, omask0, odisp2, omask20, oBin, fThreshold);
            stereo_check_pixelian_reciprocity_dual(odisp2, omask20, odisp, omask0, oBin2, fThreshold);
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        
        
        
        //! Perform Reciprocal checking if flag activated
        if (strPar.flagiRecip)
        {
            
            float fThreshold = 1.25;  // to assure pixelian precision
            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
            
            stereo_check_inverse_pixelian_reciprocity(odisp, omask0, odisp2, omask20, oBin, fThreshold);
            stereo_check_inverse_pixelian_reciprocity(odisp2, omask20, odisp, omask0, oBin2, fThreshold);
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        
        
        //! Perform Variance test
        if (strPar.flagVariance)
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
            
            flimage oTMP(input.w(), input.h()); 	oTMP = fLarge;
            flimage oTMP2(input2.w(), input2.h()); 	oTMP2 = fLarge;
            
            
            stereo_check_horizontal_variance( input, oBin, strPar.valueVariance, strPar.fSigma, 1 ,strPar.prolate.w(), strPar.prolate.h(), oTMP);
            
            stereo_check_horizontal_variance( input2, oBin2, strPar.valueVariance, strPar.fSigma, 1 , strPar.prolate.w(), strPar.prolate.h(), oTMP2);
            
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        
        
        
        //! Perform Derivative test
        if (strPar.flagHorInfo)
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
            
            stereo_check_integral_derivatives(input, oBin, strPar.valueHorInfo, strPar.fSigma, strPar.prolate);
            stereo_check_integral_derivatives(input2, oBin2, strPar.valueHorInfo, strPar.fSigma, strPar.prolate);
            // dilate the rejected points by the size of the correlation window : invalidates the whole patch
            dilate_mask_by_window_support(oBin, strPar.prolate) ;
            dilate_mask_by_window_support(oBin2, strPar.prolate) ;
            
            
            //oBin.save("st_horinfo.pmf");
            
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        
        
        
        
        
        
        //! Perform SelfSimilarity
        if(! USE_OLD_SELFSIM_ORDER()) 
        if (strPar.flagSelfSimilarity)
        {
            
            float fTransFactor = 1.0 / (float) (2 * strPar.inPrecisions);
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oTMP(input.w(), input.h()); 	oTMP = 255.0f;
            
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            
            stereo_check_strobe_and_self_simililarity_effect(input, odist, oselfdist, omask0, oBin, 1, strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP,strPar.flagListKernels);
            stereo_check_strobe_and_self_simililarity_effect(input2, odistr, oselfdistr, omask20, oBin2, 1 , strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP, strPar.flagListKernels);
            
            
            char str[1024];
            sprintf(str, "st_ssimilarity%d_%d.tif",global_scale_index, global_window_index);
            oBin.save(str);

            //oBin.save("st_ssimilarity.tif");

            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
            
        }
        




        //! Perform MinDist
        if (strPar.flagMinDist)
        {
            
            
            //float fThreshold = 1.0 / (float) strPar.inPrecisions;
            float fThreshold = strPar.valueMinDist;
            
            
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            flimage oTMP(input.w(), input.h()); 	oTMP = fLarge;
            flimage oTMP2(input2.w(), input2.h()); 	oTMP2 = fLarge;
            
            //libIIPStable::flimage prolate;
            //prolate.create(3,3);
            //prolate=1;
        
            //stereo_check_min_dist( odisp, odist,  omask0,  oBin, fThreshold, strPar.prolate, oTMP);
            //stereo_check_min_dist( odisp2, odistr,  omask20,  oBin2, fThreshold, strPar.prolate, oTMP2);
            // CONSIDER THE REJECTIONS OF THE CURRENT omask 
            if (strPar.flagRegressionMinDist) {
               stereo_check_regression_min_dist( odisp, odist,  omask,  oBin, fThreshold, strPar.prolate, oTMP);
               stereo_check_regression_min_dist( odisp2, odistr,  omask2,  oBin2, fThreshold, strPar.prolate, oTMP2);
            } else {
               stereo_check_min_dist( odisp, odist,  omask,  oBin, fThreshold, strPar.prolate, oTMP);
               stereo_check_min_dist( odisp2, odistr,  omask2,  oBin2, fThreshold, strPar.prolate, oTMP2);
            }
            
            
            char str[1024];
            sprintf(str, "st_mindist%d_%d.tif",global_scale_index, global_window_index);
            oBin.save(str);
            //! oBin.save("st_mindist.pmf");
            
            //! Dilate of 1 pixel
            if (strPar.flagMinDistDilatation && ( ( DISABLE_MINDIFFDILATATION_AT_LOWER_SCALES() <= 0 ) || strPar.currentScale==1) )
            {
                flimage itmp = oBin;
                flimage itmp2 = oBin2;
                float fRadius = 1.0;
                fiPatchMin(itmp.v(), oBin.v(), fRadius, itmp.w(), itmp.h());
                fiPatchMin(itmp2.v(), oBin2.v(), fRadius, itmp2.w(), itmp2.h());
            }
            
            
            //oBin.save("st_mindist.pmf");

            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        
        
        //! Perform SelfSimilarity
        if(USE_OLD_SELFSIM_ORDER()) 
        if (strPar.flagSelfSimilarity)
        {
            
            float fTransFactor = 1.0 / (float) (2 * strPar.inPrecisions);
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oTMP(input.w(), input.h()); 	oTMP = 255.0f;
            
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            
            stereo_check_strobe_and_self_simililarity_effect(input, odist, oselfdist, omask0, oBin, 1, strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP,strPar.flagListKernels);
            stereo_check_strobe_and_self_simililarity_effect(input2, odistr, oselfdistr, omask20, oBin2, 1 , strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP, strPar.flagListKernels);
            
            
            char str[1024];
            sprintf(str, "st_ssimilarity%d_%d.tif",global_scale_index, global_window_index);
            oBin.save(str);

            //oBin.save("st_ssimilarity.tif");

            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
            
        }
        





        
        
        // Perform Remove Isolated
        if (strPar.flagRemoveIsolated && ! DONT_REMOVE_ISOLATED_INNER())
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            
            
            // USE THE CURRENT WINDOW
            //int win = (MAX(strPar.prolate.w(), strPar.prolate.h()) - 1) / 2;
            //stereo_remove_isolated_points(omask, win, oBin, strPar.valueRemoveIsolated);
            //stereo_remove_isolated_points(omask2, win, oBin2, strPar.valueRemoveIsolated);
            stereo_remove_isolated_points_window(omask, strPar.prolate, oBin, strPar.valueRemoveIsolated);
            stereo_remove_isolated_points_window(omask2, strPar.prolate, oBin2, strPar.valueRemoveIsolated);
            
            char str[1024];
            sprintf(str, "st_isolated%d_%d.tif",global_scale_index, global_window_index);
            oBin.save(str);
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }




        // Perform Grain Filter
        if (strPar.flagGrainFilter && ! DONT_REMOVE_ISOLATED_INNER())
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            
            int minarea = strPar.valueGrainFilter;
//            printf("minarea %d\n", minarea);

            stereo_grain_filter(omask, minarea, oBin);
            stereo_grain_filter(omask2, minarea, oBin2);

            char str[1024];
            sprintf(str, "st_grain%d_%d.tif",global_scale_index, global_window_index);
            oBin.save(str);
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        
       
        update_dmin_dmax(Dmin,Dmax,odisp,omask, strPar.dmin0, strPar.dmax0, strPar.prolate.w(), strPar.prolate.h());
        update_dmin_dmax(iDmin,iDmax,odisp2,omask2, strPar.idmin0, strPar.idmax0, strPar.prolate.w(), strPar.prolate.h());
        
        
    }
    
    
    
    
    void update_dmin_dmax(flimage &Dmin, flimage &Dmax, flimage &out, flimage &omask, float minim, float maxim, int pwidth, int pheight)
    {
        
        
        // using min and max on neighborhood
        int spwidth = (pwidth-1)/2; spwidth++;
        int spheight = (pheight-1)/2; spheight++;
        
        
        //int imin = (int) rintf((float) (2*spwidth+1) * (float) (2*spheight + 1) * 0.5f);
        
        
        for (int jj=spheight; jj < omask.h() - spheight-1; jj++)
            for (int ii=spwidth; ii < omask.w() - spwidth-1; ii++)
                
                
                if (omask[jj*omask.w() + ii] > 0.0f)
                {
                    
                    
                    float fmax = -fLarge;
                    float fmin = fLarge;
                    
                    int iCount = 0;
                    for (int ss=-spheight; ss <= spheight; ss++)
                        for (int rr=-spwidth; rr <= spwidth; rr++)
                        {
                            int irr = ii + rr;
                            int jss = jj + ss;
                            int iindex = jss*out.w()+irr;
                            
                            
                            if (omask[iindex] > 0.0f)
                            {
                                if (out[iindex] < fmin) fmin = out[iindex];
                                
                                if (out[iindex] > fmax) fmax = out[iindex];
                                
                                iCount++;
                            }
                            
                            
                        }
                    
                    
                    
                    // if (iCount > imin)
                    //{
                    
                    
                    Dmin[jj*Dmin.w()+ii] = fmin;
                    Dmax[jj*Dmax.w()+ii] = fmax;
                    
                    //}
                    
                    
                    
                } else {
                    
                    Dmin[jj*omask.w() + ii] = minim;
                    Dmax[jj*omask.w() + ii] = maxim;
                    
                }
        
        
        
    }
    
    
    
    
    
    
    
    
    
    //************** Not here yet ****************//
    
    
    
    
    ///////////////////////
    ///////////////////////   Compute Window
    ///////////////////////
	
    
    
    
    void stereo_diagonal_flat_window(int win, int d, int inverse, flimage &kernel)
    {
        
        
        kernel.create(win, win); kernel = 0.0f;
        
        
        if (!inverse)
        {
            
            for(int ii=0; ii < win; ii++)
            {
                
                for(int jj=-d; jj<=d; jj++)
                    if (ii+jj >= 0 && ii+jj < win)  kernel[ ii * win + ii + jj] = 1.0f;
                
            }
            
        } else
        {
            
            for(int ii=0; ii < win; ii++)
            {
                
                for(int jj=-d; jj<=d; jj++)
                    if (ii+jj >= 0 && ii+jj < win)  kernel[ (win - ii - 1) * win + ii + jj] = 1.0f;
                
            }
            
            
        }
        
        
        //! Count number of ones
        // int iSum = 0.0f; for(int ii=0; ii < kernel.wh(); ii++) if (kernel[ii] == 1.0f) iSum++; printf("Number of ones: %d \n", iSum);
        
        kernel.normalizeL1();
        
        
        
    }
    
    
    
    
    
	
	
    
    
	//! Applies Rafa strategy to avoid micro adhesion
	//! Only gradient larger than fRafaThresh is taken into account
	void stereo_adapt_window_rafa(int ipx,int ipy, flimage &corrwindow, flimage &grad, float fRafaThresh)
    {
        
        //! half window sizes
        int spheight = (corrwindow.h()-1)/2;
        int spwidth = (corrwindow.w()-1)/2;
        
        //! Left top corner
        int iqx = ipx - spwidth;
        int iqy = ipy - spheight;
        
        fRafaThresh *= fRafaThresh;
        
        int ll=0;
        for(int s=0; s < corrwindow.h(); s++)
            for(int r=0; r < corrwindow.w(); r++, ll++)
            {
                
                float gradient = MAX( fRafaThresh,  grad[iqx + r +  (iqy + s) * grad.w()]);
                if (corrwindow[ll] > 0.0f) corrwindow[ll] /= gradient;
            }
    }
    
    
    
    
    
    
    
    //! Cancels weight of pixel no having the same pixelian depth. Points without available disparity are also canceled
    void stereo_adapt_window_to_depth(int ipx, int ipy, flimage &corrwindow, flimage &idisparity, flimage &imask, float fprec)
    {
        
        
        //! half window sizes
        int spwidth = (corrwindow.w() - 1) / 2;
        int spheight = (corrwindow.w() - 1) / 2;
        
        
        //! Left top corner
        int iqx = ipx - spwidth;
        int iqy = ipy - spheight;
        
        //! Central pixel index
        int lp = ipy * idisparity.w() + ipx;
        
        
        int l=0;
        for(int s=0; s < corrwindow.h(); s++)
            for(int r=0; r < corrwindow.w(); r++, l++)
            {
                
                int lq = (iqy + s) * idisparity.w() + (iqx + r);
                
                if (imask[lq] > 0.0)
                {
                    float fDif = fabsf(idisparity[lp] - idisparity[lq]);
                    
                    if (fDif > fprec)
                        corrwindow[l] = 0.0f;
                    
                } else
                    corrwindow[l] = 0.0f;
                
            }
    }
    
    
    
    
    ///////////////////////
    ///////////////////////   Close segments
    ///////////////////////
    void stereo_close_segments(flimage &imask, flimage &omask, int iLength)
	{
        
        omask=imask;
        for (int jj=0; jj < imask.h(); jj++)
            for (int ii=0; ii < imask.w()-1; )
                if (imask[jj*imask.w()+ii] <= 0.0f && imask[jj*imask.w()+ii+1] > 0.0f)
                {
                    
                    int mm = 0;
                    int kk= ii + 1;
                    while (kk <  imask.w() && imask[jj*imask.w() + kk] > 0.0f) {kk++; mm++; }
                    
                    if (kk == imask.w() || mm > iLength) ii = kk;
                    else
                    {
                        for(int ll = ii+1; ll < kk; ll++) omask[jj*imask.w()+ll] = 0.0f;
                        ii = kk;
                    }
                    
                } else ii++;
        
        
    }
    
    
    
    
    void stereo_detect_occlusions(flimage &idisparity, flimage &imask, flimage &omask)
	{
        
        omask=0.0f;
        for (int jj=0; jj < imask.h(); jj++)
            for (int ii=0; ii < imask.w()-1; )
                if (imask[jj*imask.w()+ii] > 0.0f && imask[jj*imask.w()+ii+1] <= 0.0f)
                {
                    
                    int mm = 0;
                    int kk= ii + 1;
                    while (kk <  imask.w() && imask[jj*imask.w() + kk] <= 0.0f) {kk++; mm++; }
                    
                    if (kk == imask.w()) ii = kk;
                    else
                    {
                        int ivalue = abs(ii + rintf(idisparity[jj*imask.w()+ii]) -  kk - rintf(idisparity[jj*imask.w()+kk]));
                        
                        for(int ll = ii+1; ll < kk; ll++) omask[jj*imask.w()+ll] = (float) mm - (float) ivalue;
                        
                        ii = kk;
                    }
                    
                } else ii++;
        
    }
    
    
    
    
    
    
    
    void stereo_remove_isolated_points_window(flimage &imask, flimage &prolate, flimage &omask, float fPercent)
    {
        
        // initialization
        omask = imask;
        
        
        // Variable: boundary
		  int pwidth = prolate.w();
		  int pheight = prolate.h();
		  
		  int spheight = (pheight-1)/2;
		  int spwidth = (pwidth-1)/2;
        int boundary = MAX(spwidth, spheight) +  1;

        int numpix = 0;
        for (int i = 0; i < pwidth*pheight; i++) 
           if(prolate[i] > 0.0f) numpix ++;
        
        
        int mincounter = (int) rintf( (float) numpix * fPercent);
        
        
        for (int y = boundary; y < imask.h()-boundary; y++)
            for (int x = boundary; x < imask.w()-boundary; x++)
                if (imask[y * imask.w() + x] > 0.0f)
                {
                    
                    
                    int icounter = 0;
                    int l0 = y * imask.w() + x;
                    
                    for (int s=-spheight; s<=spheight; s++)
                        for (int r=-spwidth; r<=spwidth; r++)
                        {
                            int l = (y+s) * imask.w() + x + r;
                            
                            if (imask[l] > 0.0f && prolate[(s+spheight)*pwidth + r+spwidth] > 0.0f)		icounter++;
                        }
                    
                    
                    if (icounter <= mincounter)		omask[l0] = 0.0f;
                    
                    
                    
                }
        
    }
    
    void stereo_remove_isolated_points(flimage &imask, int iradius, flimage &omask, float fPercent)
    {
        
        // initialization
        omask = imask;
        
        
        // Variable: boundary
        int boundary =  iradius +  3;
        
        
        int mincounter = (int) rintf( (float) (2*iradius+1) * (float) (2*iradius + 1) * fPercent);
        
        
        for (int y = boundary; y < imask.h()-boundary; y++)
            for (int x = boundary; x < imask.w()-boundary; x++)
                if (imask[y * imask.w() + x] > 0.0f)
                {
                    
                    
                    int icounter = 0;
                    int l0 = y * imask.w() + x;
                    
                    for (int r=-iradius; r<=iradius; r++)
                        for (int s=-iradius; s<=iradius; s++)
                        {
                            int l = (y+s) * imask.w() + x + r;
                            
                            if (imask[l] > 0.0f)		icounter++;
                        }
                    
                    
                    if (icounter <= mincounter)		omask[l0] = 0.0f;
                    
                    
                    
                }
        
    }
    
    ///////////////////////
    ///////////////////////   Cost histogram analysis filter
    ///////////////////////

    /* auxiliar for the qsort */
    #include <stdlib.h>
    int cmp_qsort (const void * a, const void * b)
    { 
       if ( *((float*)a) > *((float*)b) ) return 1; 
       if ( *((float*)a) < *((float*)b) ) return -1;
       return 0;
    }

    float stereo_cost_filter(flimage &cost, flimage &imask, float discard_statistic , flimage &omask)
    {

       // initialization
       omask = imask;
       int wh = cost.wh();

       float *tmp = (float*) calloc( wh, sizeof(*tmp));
       float *tmpRej = (float*) calloc( wh, sizeof(*tmp));

       int NN=0;
       int NNRej=0;
       for(int i = 0 ;  i< wh ; i++){
          if (cost[i] < fLarge) {
             if (imask[i] > 0.0f) {
                tmp[NN] = cost[i];
                NN++;
             }else {
                tmpRej[NNRej]=cost[i];
                NNRej++;
             }
          }
       }



       // sort
       qsort( tmp, NN , sizeof(float) , *cmp_qsort );
       qsort( tmpRej, NNRej , sizeof(float) , *cmp_qsort );

       

       // estimate the threshold from the accepted points
       // remove the noisy tails from the data
       // and fit a gaussian to it
       // then compute the threshold as mu + 5std
       //
       int robustNN = NN*(1-discard_statistic);
       float mu=0;
       for(int i = 0 ;  i< robustNN ; i++) mu+=tmp[i];
       mu/=robustNN;
       float var=0;
       for(int i = 0 ;  i< robustNN ; i++) var+=(tmp[i] - mu)*(tmp[i] - mu);
       var/=robustNN;
//       printf("%d %d : %f %f\n", NN, robustNN, mu, var);

       float computed_threshold = mu + 5*sqrt(var);

       // classify the pixels
       for(int i = 0 ;  i< wh ; i++) {
          omask[i]  = (cost[i] < computed_threshold); 
       }


//       // estimate the threshold as the balance point 
//       //   P(err>t | accept) = P(err<t | reject)
//       // means that the probability of being accepted with 
//       // errors larger than t equals the probability of 
//       // being rejected while having errors smaller than t 
//       //
//       // the empiric probabilities are trained from the
//       // costs of accepted and rejected pixels, after removing the tails
//       
//       int robustNN = NN*(1-discard_statistic);
//       int robustNNRej = NNRej*(1-discard_statistic);
//       int i=0;
//       int j=NNRej-1;
//       float x1 = tmp[i];
//       float x2 = tmpRej[j];
//       while (x2 > x1) {
//          i++;
//          x1 = tmp[i];
//          x2 = tmpRej[j];
//
//          float p1 = ((float)i+1)/robustNN;
//          float p2 = ((float)(NNRej-j))/robustNNRej;
//          while (p2 < p1) {
//            j--;
//            p2 = ((float)(NNRej-j))/robustNNRej;
//          }
//       }
//       float computed_threshold = x1;
//       printf("balance threshold %f\n", computed_threshold);
//
//       // classify the pixels
//       for(int i = 0 ;  i< wh ; i++) {
//          omask[i]  = (cost[i] < computed_threshold); 
//       }






       free(tmp);
       free(tmpRej);





       // generate histograms for gnuplot
       //
       //
       //
       int nbins=10000;
       float histo[nbins];
       float histoRej[nbins];

       float vmax=-INFINITY, vmin=+INFINITY;
       for(int i = 0 ;  i< wh ; i++){
          if (cost[i] < fLarge) {
             vmax = fmax(vmax, cost[i]);
             vmin = fmin(vmin, cost[i]);
          }
       }

       NN=NNRej=0;
       for (int i=0;i<nbins;i++) histo[i]=histoRej[i]=0;
       for(int i = 0 ;  i< wh ; i++) {
          if (cost[i] < fLarge) {
             if(imask[i] > 0) {
                int bin = (cost[i]-vmin)/(vmax-vmin)*(nbins-1);
                histo[bin]++;
                NN++;
             }else{
                int bin = (cost[i]-vmin)/(vmax-vmin)*(nbins-1);
                histoRej[bin]++;
                NNRej++;
             }
          }
       }

       FILE *fd;
       fd = fopen("histo.gnuplot", "w");
       fprintf(fd, "b=0\ncummulative_sum2(x)=(b=b+x,b)\n");
       fprintf(fd, "a=0\ncummulative_sum(x)=(a=a+x,a)\nset logscale x\n");
       fprintf(fd, "set xrange [%f:%f]\n", vmin, vmax);
       //fprintf(fd, "plot '-' using (1-cummulative_sum($1)) w steps, '-' using (cummulative_sum2($1)) w steps\n");
       fprintf(fd, "plot '-' using (1-($1)) w steps, '-' using ($1) w steps\n");

       float a=0;
       fprintf(fd,"%f\n", 0.);
       for (int i=0;i<nbins;i++){ a+=(histo[i]/((float)NN)); fprintf(fd, "%f\n", a); }
       fprintf(fd, "EOF\n");
       a=0;
       fprintf(fd,"%f\n", 0.);
       for (int i=0;i<nbins;i++){ a+=(histoRej[i]/((float)NNRej)); fprintf(fd, "%f\n", a); }

       fclose(fd);
       printf("visualize the densities calling> gnuplot --persist histo.gnuplot\n");



       return computed_threshold;

    }

    ///////////////////////
    ///////////////////////   Grain Filter
    ///////////////////////
    
extern "C" {
#include "../libmw/fgrain.h"
}

    void stereo_grain_filter(flimage &imask, int grainarea, flimage &omask)
    {
        
        // initialization
        omask = imask;
        
        fgrain(grainarea, imask.v(), imask.w(), imask.h(), omask.v());
        
    }
    
    
    
    
    ///////////////////////
    ///////////////////////   Reciprocity
    ///////////////////////
    
    
	void stereo_check_pixelian_reciprocity(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold, flimage &oDif)
	{
		
		omask = imask;
        oDif=fLarge;
        
		for(int ipy = 0 ;  ipy < input.h() ; ipy++)
			for(int ipx = 0 ;  ipx < input.w() ; ipx++)
				if (imask[ipy*imask.w() + ipx] > 0.0f)
				{
					
					int l = ipy* input.w() + ipx;
					
					// pixel in input2 for ipx
					int ipxRev = ipx +  (int) rintf(input[l]);
					
					if (ipxRev >= 0 && ipxRev < input2.w() && imask2[ipy * imask2.w() + ipxRev] > 0.0f)
					{
						
						// pixel in input for ipxRev
						int ipxRev2 = ipxRev + (int) rintf(input2[ipy * input2.w() + ipxRev]);
						
						float aux = fabsf((float) (ipx - ipxRev2));
						
						if (aux <= fThreshold) omask[l] = 1;
                        else omask[l] = 0.0;
						
                        oDif[l] = aux;
                        
					} else
                    {
						
						omask[l] = 0.0;
                        
                    }
					
					
					
				}
		
	}
    
	
// 'remove_only_if_inconsistent' : if the pixel on the second image is masked don't remove it in the first
// this is to prevent LR from removing more matches than necessary in a second application of the filter
	void stereo_check_pixelian_reciprocity_remove_only_if_inconsistent(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold, flimage &oDif)
	{
		
		omask = imask;
        oDif=fLarge;
        
		for(int ipy = 0 ;  ipy < input.h() ; ipy++)
			for(int ipx = 0 ;  ipx < input.w() ; ipx++)
				if (imask[ipy*imask.w() + ipx] > 0.0f)
				{
					
					int l = ipy* input.w() + ipx;
					
					// pixel in input2 for ipx
					int ipxRev = ipx +  (int) rintf(input[l]);
					
					if (ipxRev >= 0 && ipxRev < input2.w() && imask2[ipy * imask2.w() + ipxRev] > 0.0f)
					{
						
						// pixel in input for ipxRev
						int ipxRev2 = ipxRev + (int) rintf(input2[ipy * input2.w() + ipxRev]);
						
						float aux = fabsf((float) (ipx - ipxRev2));
						
						if (aux <= fThreshold) omask[l] = 1;
                        else omask[l] = 0.0;
						
                        oDif[l] = aux;
                        
					} else
                    {
						
						//omask[l] = 0.0;
                        
                    }
					
					
					
				}
		
	}
    
    
// 'dual' : When inconsistent remove the pixel in both images left and right
    void stereo_check_pixelian_reciprocity_dual(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold)
	{
		
		omask = imask;
		for(int ipy = 0 ;  ipy < input.h() ; ipy++)
			for(int ipx = 0 ;  ipx < input.w() ; ipx++)
				if (imask[ipy*imask.w() + ipx] > 0.0f && omask[ipy*imask.w() + ipx] > 0.0f)
				{
					
					int l = ipy* input.w() + ipx;
					
					// pixel in input2 for ipx
					int ipxRev = ipx +  (int) rintf(input[l]);
					
					if (ipxRev >= 0 && ipxRev < input2.w() && imask2[ipy * imask2.w() + ipxRev] > 0.0f)
					{
						
						// pixel in input for ipxRev
						int ipxRev2 = ipxRev + (int) rintf(input2[ipy * input2.w() + ipxRev]);
                        int l2 = ipy* input.w() + ipxRev2;
                        
						float aux = fabsf((float) (ipx - ipxRev2));
						
						if (aux > fThreshold) { omask[l] = 0.0; omask[l2] = 0.0;}
						
					} else
                    {
						
						omask[l] = 0.0;
                        
                    }
					
					
					
				}
		
	}
    
    
    
    
    
    
    void stereo_check_inverse_pixelian_reciprocity(flimage &input, flimage &imask, flimage &input2, flimage &imask2, flimage &omask, float fThreshold)
	{
		
		omask = imask;
		for(int ipy = 0 ;  ipy < input2.h() ; ipy++)
			for(int ipx = 0 ;  ipx < input2.w() ; ipx++)
				if (imask2[ipy*imask2.w() + ipx] > 0.0f)
				{
					
					int l = ipy* input2.w() + ipx;
					
					// pixel in input2 for ipx
					int ipxRev = ipx +  (int) rintf(input2[l]);
					
					if (ipxRev >= 0 && ipxRev < input.w() && imask[ipy * imask.w() + ipxRev] > 0.0f)
					{
						
						// pixel in input for ipxRev
						int ipxRev2 = ipxRev + (int) rintf(input[ipy * input.w() + ipxRev]);
						
                        int l2 = ipy * input.w() + ipxRev;
                        
                        
						float aux = fabsf((float) (ipx - ipxRev2));
						
						if (aux > fThreshold) omask[l2] = 0.0;
						
					} else
                    {
						
						omask[l] = 0.0;
					}
					
					
					
				}
		
	}
    
    
    
    
    
    ///////////////////////
    ///////////////////////   FLAT
    ///////////////////////
    
    
    //! Compute integral derivative
	float  stereo_compute_integral_derivative(int ipx, int ipy,  flimage &corrwindow,  cflimage &grad)
	{
		
		//! half window size
		int pwidth = corrwindow.w();
		int pheight = corrwindow.h();
		
		int spheight = (pheight-1)/2;
		int spwidth = (pwidth-1)/2;
		
		
		float e1 = 0.0, e2 = 0.0;
		
		for (int cc=0; cc < grad.c(); cc++)
        {
            
            
            //! for each line
            for(int s=-spheight ; s <= spheight; s++)
            {
                
                int i1 = (spheight+s)*pwidth + spwidth;	//! index in window of first pixel in line
                int i2 = cc * grad.wh() + (ipy + s)*grad.w() + ipx;		//! index in image of first pixel in line
                
                for(int r=-spwidth;r <= spwidth; r++)
                {
                    
                    int ii = i1 + r;
                    float gradient = grad[i2 + r];
                    
                    e1 += corrwindow[ii] * gradient;	//! += phi u'^2
                    e2 += corrwindow[ii];				//! += phi
                }
            }
            
        }
		
        
		
		return e1 / e2;
	}
	
    
    void dilate_mask_by_window_support(flimage &mask, flimage &corrwindow) 
   {

        //! Variable: boundary
        int spwidth = corrwindow.w() / 2;
        int spheight = corrwindow.h() / 2;
        int boundary = MAX(spwidth, spheight) +  1;

        flimage tmpimage(mask.w(), mask.h()); 

        // dilate the rejected points by the size of the correlation window
        tmpimage = mask; // copy the current mask before applying the dilation
        for (int y = boundary; y < mask.h()-boundary; y++)
        {
            for (int x = boundary; x < mask.w()-boundary; x++)
            {
               if ( tmpimage[ y * mask.w() + x] == 0.0 ) 
               {
         		//! half window size
         		int pwidth = corrwindow.w();
         		int pheight = corrwindow.h();
         		
         		int spheight = (pheight-1)/2;
         		int spwidth = (pwidth-1)/2;

               //! for each line
               for(int s=-spheight ; s <= spheight; s++)
                   for(int r=-spwidth;r <= spwidth; r++)
                   {
                        int iwin = (spheight+s)*pwidth + spwidth + r;
                        int iimg = (s+y)*mask.w() + r + x;
                        if(corrwindow[iwin]>0) 
                           mask[iimg] = 0;
                   }
               }

                
            }
         } 

   }
    
    void stereo_check_integral_derivatives(cflimage &input, flimage &omask, float fThreshold, float fSigma, flimage &corrwindow)
    {
        
        
        //! initialization
        omask = 0.0;
        
        
        //! Variable: boundary
        int spwidth = corrwindow.w() / 2;
        int spheight = corrwindow.h() / 2;
        int boundary = MAX(spwidth, spheight) +  1;
        
        
        //! compute square image gradient
        cflimage xgrad(input.w(), input.h(), input.c());
        {
            cflimage iiptmp = input.xgradient('f');
            for(int i=0; i < xgrad.whc(); i++)  xgrad[i] = iiptmp[i] * iiptmp[i] / 2.0f;   // since iitmp[i] / sqrt(2) has std sigma
        }

        
        
        fThreshold = fThreshold * fSigma * fSigma;
        
        
        
        flimage tmpimage(input.w(), input.h()); tmpimage = 0.0f;
        
        for (int y = boundary; y < input.h()-boundary; y++) 
        {
            for (int x = boundary; x < input.w()-boundary; x++)
            {
                
                float fprec =  stereo_compute_integral_derivative(x, y, corrwindow, xgrad);
                
                tmpimage[y * input.w() + x] = fprec;
                
                if (fprec > fThreshold)  omask[ y * input.w() + x] = 1;
                else    omask[ y * input.w() + x] = 0.0;
            }
         }

        
    }
    
    
    
    
    
    
    
    
    
    void stereo_compute_horizontal_variance(cflimage & input, cflimage & variance, int pwidth, int pheight)
    {
        
        if (!variance.isSameSize(input))
        {
            variance.erase();
            variance.create(input.w(), input.h(), input.c());
        }
        variance = 0.0f;
        
        
        
        //! Parameters
        int wRadius = pwidth / 2;
        int hRadius = pheight / 2;
        
        for (int iC = 0; iC < input.c(); iC++)
        {
            
            float *fpI = input.v(iC);
            float *fpV = variance.v(iC);
            
            
            for(int x = wRadius; x < input.w()-wRadius; x++)
                for(int y = hRadius; y < input.h()-hRadius; y++)
                {
                    
                    float var=0.0f;
                    for(int j=-hRadius;j<=hRadius;j++)
                    {
                        
                        
                        float mean = 0.0f;
                        float mean2 = 0.0f;
                        
                        for(int i=-wRadius;i<=wRadius;i++)
                        {
                            mean += fpI[ (y+j) * input.w() + x + i ];
                            mean2 += fpI[ (y+j) * input.w() + x + i ] *  fpI[ (y+j) * input.w() + x + i ] ;
                        }
                        
                        mean2 = mean2 / (float) (2*wRadius+1);
                        mean = mean / (float) (2*wRadius+1);
                        
                        var +=  mean2 - mean*mean;
                    }
                    
                    var /= (float)  (2*hRadius+1);
                    fpV[ y * input.w() + x] = var;
                    
                    
                }
            
        }
        
        
    }
    
    
    
    
    
	void stereo_check_horizontal_variance(cflimage &input, flimage &mask, float fThreshold, float fSigma, int flagWholePatch, int pwidth, int pheight, flimage &oVAR)
	{
		
        
		//! Variable: boundary
		int spwidth = (pwidth - 1) / 2;
		int spheight = (pheight - 1) / 2;
		int boundary =  spwidth +  3;
        
        
        //! Variance
		//! New variance per lines
        cflimage variance(input.w(), input.h(), input.c());
        stereo_compute_horizontal_variance(input, variance, pwidth, pheight);
        
        
        //variance.save("variance.pmf");
		
        
        
		float fSigmaThr = fThreshold * fSigma * fSigma;
        
		oVAR=0.0f;
		mask=0.0f;
		for(int ipy = boundary ;  ipy < input.h() - boundary ; ipy++)
			for(int ipx = boundary ;  ipx < input.w() - boundary; ipx++)
			{
				
				//! Multichannel variance
				float fVar = 0.0f;
				for	(int i = 0 ; i < input.c(); i++)
				{
					fVar += variance[i * input.wh() + ipy * input.w() + ipx];
				}
				
				fVar /= (float) input.c();
				
				oVAR[ipy * input.w() + ipx ] = fVar;
                
                //! Test Variance and option for invalidating the whole patch
				if (fVar  <  fSigmaThr)
				{
					
                    if (flagWholePatch)
                    {
                        for (int j=-spwidth; j <= spwidth; j++)
                            for (int i=-spwidth; i <= spwidth; i++)
                            {
                                mask[(ipy + j) * input.w() + ipx + i] = 0.0;
                            }
                    } else
                        mask[ipy * input.w() + ipx ] = 0.0;
                    
				} else
                {
                    mask[ipy * input.w() + ipx ] = 1.0f;
                }
				
				
			}
		
	}
	
	
    
    
    
    
     
    float stereo_compute_pixel_patch_list_distance(cflimage &input, cflimage &input2, int ipx, int ipy, int iqx, int iqy, cflimage &mean, cflimage &mean2, flimage &corrwindow, int itypeDist)
    {
        
        float fDistance;
        if ( itypeDist == L2)
            fDistance = distancePatchListWL2(input,input2,ipx, ipy, iqx, iqy, corrwindow);

        else if (itypeDist == L2M)
            fDistance = distancePatchListWL2M(input, input2,ipx, ipy, iqx, iqy, corrwindow, mean, mean2);

        else {printf("... please enter a correct distance type"); exit(-1);}
        
        
        return fDistance;
        
    }
    
    
    
	
	
	
    
    float stereo_compute_pixel_patch_distance(cflimage &input, cflimage &input2, int ipx, int ipy, int iqx, int iqy, cflimage &mean, cflimage &mean2, flimage &corrwindow, int itypeDist)
    {
        
        float fDistance;
        if ( itypeDist == L2)
            fDistance = distancePatchWL2(input,input2,ipx, ipy, iqx, iqy, corrwindow);
        
        else if (itypeDist == L2M)
            fDistance = distancePatchWL2M(input, input2,ipx, ipy, iqx, iqy, corrwindow, mean, mean2);
        
        else {printf("... please enter a correct distance type"); exit(-1);}
        
        
        return fDistance;
        
    }
    
    
    
    
    
    
    void stereo_check_strobe_and_self_simililarity_effect(cflimage &input, flimage &odist1, flimage &odist2, flimage &imask, flimage &omask,  int iflagTrans, int itypeDist, float fTrans, float fSigma, float fMult, flimage &prolate, flimage &oDif, int flagListKernels=0)
    {
        
        
        
        //! Compute translations and means
        cflimage translated1,  translated2;
        cflimage mean, meant1,  meant2;
        if (iflagTrans)
        {
            
            translated1 = input;
            translated1.fftShear(0.0, iipHorizontal, fTrans, 0, 1 );
            
            translated2 = input;
            translated2.fftShear(0.0, iipHorizontal, -fTrans, 0, 1 );
            
            
            if (itypeDist == L2M)
            {
                if (flagListKernels) {
                    mean = input.patchListMean(prolate);
                    meant1 = translated1.patchListMean(prolate);
                    meant2 = translated2.patchListMean(prolate);    
                } else {
                    mean = input.patchMean(prolate);
                    meant1 = translated1.patchMean(prolate);
                    meant2 = translated2.patchMean(prolate);                    
                }
            }
            
        }
        
        
        
        
        
        //! Variable: Noise variance
        float fSigma2 = fSigma  *  fSigma;
        
        
        //! Variable: boundary
        int spwidth = (prolate.w() - 1) / 2;
        int spheight = (prolate.h() - 1) / 2;
        int boundary =  2 * MAX(spwidth, spheight) +  1;
        
        omask=0.0f;
        oDif = fLarge;
        
        for(int ipy = boundary ;  ipy < input.h() - boundary ; ipy++)
            for(int ipx = boundary ;  ipx < input.w() - boundary; ipx++)
                if (imask[ipy*imask.w() + ipx] > 0.0f)
                {
                    
                    int l = ipy * input.w() + ipx;
                    
                    float fDistTrans;
                    if (iflagTrans)
                    {
                        
                        float fDistance1 = 0.0f;
                        if (flagListKernels)
                            fDistance1 = stereo_compute_pixel_patch_list_distance(input,translated1,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant1,  prolate,  itypeDist);
                        else
                            fDistance1 = stereo_compute_pixel_patch_distance(input,translated1,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant1,  prolate,  itypeDist);
                        
                        float fDistance2 = 0.0f;
                        if (flagListKernels)
                            fDistance1 = stereo_compute_pixel_patch_list_distance(input,translated2,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant2,  prolate,  itypeDist);
                        else
                            fDistance2 = stereo_compute_pixel_patch_distance(input,translated2,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant2,  prolate,  itypeDist);
                        
                        fDistTrans =   MAX(fDistance1, fDistance2);
                        
                    } else
                        fDistTrans = 0.0f;
                    
                    
                    float fDifference = odist2[l] - odist1[l];
                    oDif[l] = fDifference - fDistTrans;
                    if (fDifference - fDistTrans > 2.0f * fMult * fSigma2) omask[l] = 1.0f;
                    
                    
                }
        
    }
    
    
    
    
    
    void stereo_check_distance(cflimage &input, flimage &idist1, flimage &imask, flimage &omask,  int iflagTrans, int itypeDist, float fTrans, float fSigma, float fMult, flimage &prolate, flimage &oDif, int flagListKernels=0)
    {
        
        
        
        //! Compute translations and means
        cflimage translated1,  translated2;
        cflimage mean, meant1,  meant2;
        if (iflagTrans)
        {
            
            translated1 = input;
            translated1.fftShear(0.0, iipHorizontal, fTrans, 0, 1 );
            
            translated2 = input;
            translated2.fftShear(0.0, iipHorizontal, -fTrans, 0, 1 );
            
            
            if (itypeDist == L2M)
            {
                if (flagListKernels) {
                    mean = input.patchListMean(prolate);
                    meant1 = translated1.patchListMean(prolate);
                    meant2 = translated2.patchListMean(prolate);    
 //                   mean = input.patchMean(prolate);
 //                   meant1 = translated1.patchMean(prolate);
 //                   meant2 = translated2.patchMean(prolate);    
                } else {
                    mean = input.patchMean(prolate);
                    meant1 = translated1.patchMean(prolate);
                    meant2 = translated2.patchMean(prolate);                    
                }
            }
            
        }
        
        
        
        
        
        //! Variable: Noise variance
        float fSigma2 = fSigma  *  fSigma;
        
        
        //! Variable: boundary
        int spwidth = (prolate.w() - 1) / 2;
        int spheight = (prolate.h() - 1) / 2;
        int boundary = 2* MAX(spwidth, spheight) +  1;
        
        
        omask=0.0f;
        oDif = fLarge;
        
        for(int ipy = boundary ;  ipy < input.h() - boundary ; ipy++)
            for(int ipx = boundary ;  ipx < input.w() - boundary; ipx++)
                if (imask[ipy*imask.w() + ipx] > 0.0f)
                {
                    
                    int l = ipy * input.w() + ipx;
                    
                    float fDistTrans;
                    if (iflagTrans)
                    {
                        float fDistance1 = 0.0f;
                        if (flagListKernels)
                            fDistance1 = stereo_compute_pixel_patch_list_distance(input,translated1,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant1,  prolate,  itypeDist);
                        else
                            fDistance1 = stereo_compute_pixel_patch_distance(input,translated1,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant1,  prolate,  itypeDist);
                        
                        float fDistance2 = 0.0f;
                        if (flagListKernels)
                            fDistance1 = stereo_compute_pixel_patch_list_distance(input,translated2,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant2,  prolate,  itypeDist);
                        else
                            fDistance2 = stereo_compute_pixel_patch_distance(input,translated2,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight, mean, meant2,  prolate,  itypeDist);
                        
                        fDistTrans =   MAX(fDistance1, fDistance2);
                        
                    } else
                        fDistTrans = 0.0f;
                    
                    
                    float fDifference = idist1[l];
                    oDif[l] = fDifference - fDistTrans;
                    
                    if (fDifference - fDistTrans < 2.0f * fMult * fSigma2) omask[l] = 1.0f;
                    
                    
                }
        
    }
    
    
    
    
    
#include "ransacFIX.h"
    
    
    ///////////////////////
    ///////////////////////  MINDIST FILTER
    ///////////////////////
    
    void stereo_check_regression_min_dist(flimage &idisparity, flimage &idistance, flimage &imask, flimage &omask, float fThreshold, flimage &prolate, flimage &oDif)
	{
		
		// initialization
		omask = 0.0;
		
		
		// Variable: boundary
        int pwidth = prolate.w();
        int pheight = prolate.h();
        
		int spwidth = pwidth / 2;
		int spheight = pheight / 2;
		int boundary = MAX(spwidth, spheight) +  1;
		
		
		for (int y = boundary; y < idisparity.h()-boundary; y++)
			for (int x = boundary; x < idisparity.w()-boundary; x++)
				if (imask[ y * idisparity.w() + x] > 0.0f)
				{
					
					
					int l0 = y * idisparity.w() + x;
               float fMinDist = fLarge;
					float sdisparity = idisparity[l0];
					float sdist      = 0.0001;
               int sidx=0;
               int sidy=0;

					int l2=0;
					for (int r=-spheight; r<=spheight; r++)
						for (int s=-spwidth; s<=spwidth; s++, l2++)
						{

							int l = (y + r) * idisparity.w() + x + s;
							if (imask[l] > 0.0f && prolate[l2] > 0.0f)   
							{

								if (idistance[l] < fMinDist)
								{
									fMinDist = idistance[l];
									sdisparity = idisparity[l];
                           sdist = fmax(0.0001,hypot(r,s));
                           sidx = s;         // store the position of the min
                           sidy = r;
								}
							}
							
						}

					float dif = fabsf(sdisparity - idisparity[l0]);


					oDif[l0] = dif;
                    
               // If the mindiff is within the threshold, then keep it
               // This is reasonable because we know that the points
               // accepted by mindiff are correct, the problem
               // is that it also discards many correct points
					if (dif <= fThreshold) 
               {
						omask[l0] = 1.0f;
               } 
               else 
               {
                  // If the center doesn't pass the mindiff text, 
                  // maybe it's because is on a slanted plane. 
                  // Let's check if the regression can explain it, 
                  // otherwise reject 
						omask[l0] = 0.0f;

                  // store the points of the window in a structure for ransac
                  float points[pwidth*pheight][3];
                  int FIXED_inliers[pwidth*pheight];
                  int num_FIXED =0;
                  int pidx=0;
                  int idx_center = -1;
                  int idx_min = -1;
					   int l2=0;
   					for (int r=-spheight; r<=spheight; r++)
   						for (int s=-spwidth; s<=spwidth; s++, l2++)
   						{
   							int l = (y + r) * idisparity.w() + x + s;
   							if (imask[l] > 0.0f && prolate[l2] > 0.0f)
   							{
                           points[pidx][0] = s;
                           points[pidx][1] = r;
                           points[pidx][2] = idisparity[l];
                           // fix the minimum
                           if((s==sidx && r==sidy)) {
//                              FIXED_inliers[num_FIXED] = pidx;
//                              num_FIXED++;
                              idx_min = pidx;
                           }
                           // fix the center
                           if((s==0 && r==0)) {
                              FIXED_inliers[num_FIXED] = pidx;
                              num_FIXED++;
                              idx_center = pidx;
                           }
                           pidx++;
   							}
   							
   						}
                  
                  
                  if(pidx>3) {
                     float model[3];
                     int numbestcons=0;
                     float err;

                     // the best fitting model with at least pix/2 inliers. the center point is fixed num_FIXED=1
                     int ransac_success = ransacFIX(3, pidx, points[0], 3, 100, fThreshold, pidx/2, model, &numbestcons, NULL, &err, num_FIXED, FIXED_inliers);


                     float dif = INFINITY;
                     if ( ransac_success ) dif = fabs (test_model(3, points[idx_min] , model));
      					
      					oDif[l0] = dif;
      
      					if (dif < fThreshold) {
      						omask[l0] = 3.0f;
                     }
      					else {
      						omask[l0] = 0.0f;
                     }
   
                  }
   
                  }
                    
					
				}
		
	}
    
	
    void stereo_check_min_dist(flimage &idisparity, flimage &idistance, flimage &imask, flimage &omask, float fThreshold, flimage &prolate, flimage &oDif)
	{
		
		
		// initialization
		omask = 0.0;
		
		
		// Variable: boundary
        int pwidth = prolate.w();
        int pheight = prolate.h();
        
		int spwidth = pwidth / 2;
		int spheight = pheight / 2;
		int boundary = MAX(spwidth, spheight) +  1;
		
		
		for (int y = boundary; y < idisparity.h()-boundary; y++)
			for (int x = boundary; x < idisparity.w()-boundary; x++)
				if (imask[ y * idisparity.w() + x] > 0.0f)
				{
					
					
					int l0 = y * idisparity.w() + x;
                    float fMinDist = fLarge;
					float sdisparity = idisparity[l0];
					float sdist      = 0.0001;
					
					int l2=0;
					for (int r=-spheight; r<=spheight; r++)
						for (int s=-spwidth; s<=spwidth; s++, l2++)
						{
							int l = (y + r) * idisparity.w() + x + s;
							if (imask[l] > 0.0f && prolate[l2] > 0.0f)
							{
								if (idistance[l] < fMinDist)
								{
									fMinDist = idistance[l];
									sdisparity = idisparity[l];
                           sdist = fmax(0.0001,hypot(r,s));
								}
							}
							
						}
					
					float dif = fabsf(sdisparity - idisparity[l0]);

               // The modified mindist imposes a constraint on the maximum slant 
               // of the window. It si done by dividing to the distance at the center.
               if( USE_MODIFIED_MINDIST() ) 
                  dif = fabsf(sdisparity - idisparity[l0])/sdist;

					oDif[l0] = dif;
                    
					if (dif <= fThreshold)
						omask[l0] = 1.0f;
					
					else
						omask[l0] = 0.0f;
					
				}
		
	}
	
    
    
    
	
	
	
	
	
    
	void stereo_apply_mask(flimage &flDisparity, flimage &flBin, cflimage &out, float flMin = fLarge, float flMax = -fLarge)
	{
		
		if (flMin == fLarge && flMax == -fLarge)
		for (int i=0; i < flDisparity.wh(); i++)
		{
			if (flBin[i] > 0.0f && flDisparity[i] > flMax) flMax = flDisparity[i];
			if (flBin[i] > 0.0f && flDisparity[i] < flMin) flMin = flDisparity[i];
		}
		
		
		for (int i=0; i < flDisparity.wh(); i++)
		{
			
			if (flBin[i] > 0.0f)
			{
				float value = 255.0f * (flDisparity[i] - flMin) / (flMax - flMin);
				
				out[i] = value;
				out[flDisparity.wh() + i ]= value;
				out[2 * flDisparity.wh() + i ]= value;
				
			} else {
				
				out[i] = 128.0f;
				out[flDisparity.wh() + i ]= 0.0f ;
				out[2 * flDisparity.wh() + i ]= 0.0f ;
				
			}
			
		}
	}
	
	
    
    
	
    ///////////////////////
    /////////////////////// Distance computation
    ///////////////////////
	
	
    void stereo_subpixel_computation(cflimage &input, cflimage &input2, flimage &idisparity, flimage &imask, flimage &out,flimage &omask, flimage &odist, strSubPixelPar &strPar)
    {
        
        
        //! Build prolate
        flimage prolate = strPar.prolate;
        
        
        
        //! compute square image gradient
        flimage  gxgrad(input.w(), input.h());
        if (strPar.sflagRafa)
        {
            cflimage xgrad(input.w(), input.h(), input.c());
            
            cflimage iiptmp = input.xgradient('f');
            for(int i=0; i < xgrad.whc(); i++)  xgrad[i] = iiptmp[i] * iiptmp[i];
            
            gxgrad = xgrad.getGray();
            
        }
        
        
        
        //! Compute image translations of right image by using the shear function
        //! For a fixed number of scales n: we increase precision as 1, 1/2, 1/4, ... , 1/2^n
        //! Precision for examples 1/4 can be achieved by using  2^{n-2} * 1 / 2^n.  So we only need to translate image of factor 1 / 2^n
        int nimages = (int) powf(2.0f, strPar.fnScales);
        float compstep = 1.0 / (float) nimages;
        float *steps = new float[strPar.fnScales+1];
        cflimage *images2 = new cflimage[nimages+1];
        {
            
            steps[0] = 1.0;
            for(int i=1; i <= strPar.fnScales; i++) steps[i] = steps[i-1] / 2.0f;
            
            
            images2[0] = input2;
            for(int i=1; i <= nimages; i++)
            {
                float step = (float) i * steps[strPar.fnScales];
                images2[i] = input2;
                images2[i].fftShear(0.0, iipHorizontal, step, 0, 1);
            }
            
        }
        
        
        
        
        //! Compute translations of left image if right-left option activated
        cflimage *images1 = NULL;
        /*if (strPar.sflagRecip)
         {
         images1 = new cflimage[nimages+1];
         images1[0] = input;
         
         for(int i=1; i <= nimages; i++)
         {
         
         float step = (float) i * steps[strPar.snscales];
         images1[i] = input;
         images1[i].fftShear(0.0, iipHorizontal, step, 0, 1 SimFlag);
         }
         }
         */
        
        
        
        //! Processing
        odist = 0.0f;
        omask = 0.0f;
        out=0.0f;
        
        
        //! Compute precision of current window
        //if (strPar.sflagPrecision)   strPar.sPrecision = 1.0f;
        //if (strPar.sflagRecip)   strPar.sReciprocity = 1.0f;
        
        
        
        //! Boundary
        int spwidth  = (prolate.w()-1)/2;
        int spheight = (prolate.h()-1)/2;
        int	boundary = prolate.w();
        
        
        //! Begin for each pixel with available pixel disparity
#pragma omp parallel shared(images1, images2, gxgrad, prolate)
        {
            
            
#pragma omp for schedule(dynamic) nowait
            for(int ipy = boundary ;  ipy < input.h() - boundary ; ipy++)
            {
                for(int ipx = boundary ;  ipx < input.w() - boundary; ipx++)
                    if (imask[ipy*imask.w()+ipx] > 0.0f)
                    {
                        
                        //! index of (ipx,ipy)
                        int l=ipy*input.w()+ipx;
                        
                        //! Correlation window
                        flimage corrwindow = prolate;
                        
                        
                        //! Initial subpixelian scale, depending in precision of initial disparity
                        int initscale =  strPar.inScales;
                        
                        
                        //! corresponding point in second image (updated during the subpixelian refinement)
                        //! initially truncate to initscale (1, 1/2, 1/4, 1/8, 1/16, ....)
                        //! initscale=0     truncate to 1
                        //! initscale=2^i   truncate to 1/2^i
                        
                        float fValueTMP = fabsf(idisparity[l]);
                        float fAUX = 0.0f;
                        for (int i=0; i <= initscale; i++)
                        {
                            float fAUX0 = floorf(fValueTMP / steps[i]) * steps[i];
                            fAUX += fAUX0;
                            fValueTMP = fValueTMP - fAUX0;
                        }
                        
                        if (idisparity[l] < 0.0f) fAUX = -fAUX;
                        float iselected = fAUX + (float) ipx;
                        
                        
                        
                        //! Compute adaptive window DEPTH
                        corrwindow = prolate;
                        
                        float sdPrec = 1.5f * steps[strPar.inScales];
                        if (strPar.sflagSameDepth) stereo_adapt_window_to_depth(ipx, ipy, corrwindow, idisparity, imask, sdPrec);
                        corrwindow.normalizeL1();
                        
                        
                        //! Compute adaptive window RAFA
                        if (strPar.sflagRafa) stereo_adapt_window_rafa( ipx, ipy, corrwindow, gxgrad, strPar.sValueRafa * strPar.sfSigma);
                        corrwindow.normalizeL1();
                        
                        
                        //! Compute Bernard precision formula of current window
                        //if (strPar.sflagPrecision)
                        //{
                        //    strPar.sPrecision[l] = 2.0f * strPar.sfSigma *  strPar.sfSigma * stereo_compute_precision(ipx, ipy, corrwindow, xgrad);
                        //}
                        
                        
                        
                        //! Correlate with second image
                        int pindex=0;	//! index of translated image giving the best match
                        int ppoint=0;	//! x component of pixel giving the best match in translated image
                        if (iselected > boundary && iselected < input2.w() - boundary)
                        {
                            // printf("(%d, %d) : idisp: %f  fAux: %f  iselected: %f \n", ipx, ipy, idisparity[l], fAUX, iselected);
                            //getchar();
                            
                            
                            float second_best2 = 1000000000.0f;
                            //! we check  -subi/2^i, ..., -1/2^i ,0, 1/2^i, ... , -subi/2^i centered at best pixel at previous scale
                            int subi=3;
                            
                            for(int i = initscale+1; i <= strPar.fnScales; i++)
                            {
                                
                                
                                float step = steps[i];					//! step = 1 / 2^i, at first pass step = 1/2
                                float selected0 = iselected;			//! at first i, equals pixelian disparity
                                
                                for(int k = -subi ; k <= subi; k++)
                                {
                                    
                                    
                                    //! subpixel coordinates of current point
                                    float k0 = selected0 + (float) k * step;
                                    
                                    
                                    int point = (int) floor(k0);
                                    int index =  (int) ((k0 - (float) point) /  compstep);
                                    
                                    
                                    float distance;
                                    distance =   distancePatchWL2(input,images2[index],ipx-spwidth,ipy-spheight,point-spwidth,ipy-spheight,corrwindow);
                                    distance /= (float) input.c();
                                    
                                    if (distance < second_best2) {second_best2 = distance ; iselected = k0; ppoint=point; pindex=index;}
                                    
                                }
                            }
                            
                            
                            
                            /*
                             //! Check reciprocity if selected
                             if (strPar.sflagRecip)			//! Begin of reciprocity flag activated
                             {
                             
                             //! pixel coordinate of selected point in right image
                             float irselected = rintf(iselected);
                             
                             
                             //! selected point correlated in right - left?
                             int leftcorrespondingCorrelated = (int) strPar.imaskr[ipy * input2.w() + (int) irselected];
                             
                             //! pixel coordinate of match in left image
                             irselected = irselected + strPar.idispr[ipy * input2.w() + (int) irselected];
                             
                             
                             
                             //! if pixel in right image correlated and correspoding in left image agrees with boundary conditions
                             if ( leftcorrespondingCorrelated && irselected > boundary && irselected < input.w() - boundary)
                             {
                             
                             
                             float second_best1 = 1000000000.0f;
                             for(int i = initscale; i <= strPar.snscales; i++)
                             {
                             
                             
                             float step = steps[i];
                             float selected0 = irselected;
                             
                             for(int k = -subi ; k <= subi; k++)
                             {
                             
                             float k0 = selected0 + (float) k * step;
                             
                             int point = (int) floor(k0);
                             int index =  (int) ((k0 - (float) point) /  compstep);
                             
                             float distance;
                             
                             
                             distance = distancePatchWL2(images1[index],images2[pindex],point-spwidth,ipy-spheight,ppoint-spwidth,ipy-spheight,corrwindow);
                             
                             distance /= (float) input.c();
                             
                             if (distance < second_best1) {second_best1 = distance ; irselected = k0;}
                             
                             }
                             
                             
                             }
                             
                             strPar.sReciprocity[l] = fabsf((float) ipx - irselected);
                             
                             
                             }  //! END of pixel correlated and no boundary problems
                             
                             
                             
                             }  //! END of reciprocity flag activated
                             */
                            
                            
                            omask[l] = 255.0;
                            out[l] = iselected - (float) ipx;
                            
                            //float dif1 = distancePatchWL2(input,image1dist1,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight,corrwindow) / (float) input.c();
                            //float dif2 = distancePatchWL2(input,image1dist2,ipx-spwidth,ipy-spheight,ipx-spwidth,ipy-spheight,corrwindow) / (float) input.c();
                            //odist[l] =  sqrtf( MAX(second_best2 - MAX(dif1, dif2),0.0f) / (strPar.sfSigma * strPar.sfSigma) ) ;
                            
                            //!odist[l] =  sqrtf( second_best2  / (strPar.sfSigma * strPar.sfSigma) ) ;
                            
                        }
                        
                        
                        
                    }
                
            } 	//! End for each pixel with available pixel disparity
            
        }
      delete[] images2;
      delete[] steps;
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
	
	
}





