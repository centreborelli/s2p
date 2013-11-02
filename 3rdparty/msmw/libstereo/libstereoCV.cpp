#include "libstereo.h"
#include "smartparameter.h"
#include "cubic.h"
#include "math.h"

// Obscure parameters
SMART_PARAMETER_INT(DONT_USE_SUBPIX_AT_LOW_SCALES,0)
SMART_PARAMETER_INT(USE_MODIFIED_MINDIST,0)


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
                        if (strPar.flagSelfSimilarity)
                        {
                            
                            //! range adapted to window size
                            int imin = MAX(boundary, ipx + (int) Dmin[l]);
                            int imax = MIN(input.w() - boundary, ipx + (int) Dmax[l]);
                            
                            
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
                                        
                                        if (fCurrentDistance < fSelfBestDistance)  fSelfBestDistance = fCurrentDistance;
                                    }
                            
                            oselfdist[l] = fSelfBestDistance;
                            
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
	
    
    
    void stereo_pixel_chain_multiwinCV_simpler(cflimage &input, cflimage &input2, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, 
          flimage &odisp, flimage &odispr, flimage &odist, flimage &odistr, flimage &omask, flimage &omaskr, strParameters &strPar); // << attention these are also outputs
    
    
    ////////////////////////////
    ////////////////////////////
    //////////////////////////// MultiScale
    ////////////////////////////
    ////////////////////////////
    
    void stereo_pixel_multiscale_chain(cflimage &input, cflimage &input2, int &in, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage &out, flimage &out2, flimage &odist, flimage &odistr, flimage &omask, flimage &omask2, strParameters &strPar)
    {
        
        
        // initial scale is 0
        int currentScale = in;
        
        
        if (in < strPar.nScales)
        {
            
            // subsample by factor 2: input, input2, dmin and dmax
            float fGaussianConvol = 1.2f;
            cflimage sinput  = input.subSample(2, fGaussianConvol);
            cflimage sinput2 = input2.subSample(2, fGaussianConvol);
            
            flimage sDmin = Dmin.subSample(2, fGaussianConvol);
            flimage sDmax = Dmax.subSample(2, fGaussianConvol);
            
            flimage siDmin = iDmin.subSample(2, fGaussianConvol);
            flimage siDmax = iDmax.subSample(2, fGaussianConvol);
            
            
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
            
            
            
            stereo_pixel_multiscale_chain(sinput, sinput2, in, sDmin, sDmax, siDmin, siDmax, sout, sout2, sodist, sodistr, somask, somask2, strPar);
            
            
            
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
        
        // for lower scales it is always better to have a subpixel estimate (and it is very cheap)
        if (DONT_USE_SUBPIX_AT_LOW_SCALES()==0)  {
           if (currentScale>1 && csPar.inPrecisions<4 ) csPar.inPrecisions = 4;
        }
        else { 
           printf("USING THE OBSCURE PARAMETER: DONT_USE_SUBPIX_AT_LOW_SCALES=%d ", DONT_USE_SUBPIX_AT_LOW_SCALES());
        }
        
        // apply stereo chain: dmin, dmax, idmin and idmax are updated inside stereo_pixel_chain
        if(strPar.flagMultiWin) 
           stereo_pixel_chain_multiwinCV_simpler(input,input2,Dmin,Dmax, iDmin, iDmax, out, out2, odist, odistr, omask, omask2, csPar);
           //stereo_pixel_chain_multi_window(input,input2,Dmin,Dmax, iDmin, iDmax, out, out2, odist, odistr, omask, omask2, csPar);
        else
           stereo_pixel_chain(input,input2,Dmin,Dmax, iDmin, iDmax, out, out2, odist, odistr, omask, omask2, csPar);
        
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
        
        
    }
    
    
    
	




    void update_dmin_dmax_window(flimage &Dmin, flimage &Dmax, flimage &out, flimage &omask, float minim, float maxim, flimage &prolate)
    {
        
        
        // using min and max on neighborhood
        int spwidth =  (prolate.w()-1)/2; //spwidth++;
        int spheight = (prolate.h()-1)/2; //spheight++;
        
        
        //int imin = (int) rintf((float) (2*spwidth+1) * (float) (2*spheight + 1) * 0.5f);
        
#pragma omp parallel for
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
    
    


    void fiComputeIntegralImageINPLACE(float *d, int width, int height)
    {
       // recursivity in the rows
       for(int jj=0; jj < height; jj++)
          for(int ii=1; ii < width; ii++)
             d[width*jj+ii] = d[width*jj+ii-1] + d[width*jj+ii];

       // recursivity in the columns
       for(int ii=0; ii < width; ii++)
          for(int jj=1; jj < height; jj++)
             d[width*jj+ii] = d[width*(jj-1)+ii] + d[width*jj+ii];
    }





    void apply_filters_chain(cflimage &input, cflimage &input2, 
          flimage &odisp, flimage &odisp2, flimage &odist, flimage &odistr, flimage &oselfdist , flimage &oselfdistr, 
          flimage &omask, flimage &omask2, strParameters &strPar)
    {

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
            stereo_check_min_dist( odisp, odist,  omask,  oBin, fThreshold, strPar.prolate, oTMP);
            stereo_check_min_dist( odisp2, odistr,  omask2,  oBin2, fThreshold, strPar.prolate, oTMP2);
            
            
            //! oBin.save("st_mindist.pmf");
            
            //! Dilate of 1 pixel
            if (strPar.flagMinDistDilatation)
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
        if (strPar.flagSelfSimilarity)
        {
            
            float fTransFactor = 1.0 / (float) (2 * strPar.inPrecisions);
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oTMP(input.w(), input.h()); 	oTMP = 255.0f;
            
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            
            stereo_check_strobe_and_self_simililarity_effect(input, odist, oselfdist, omask0, oBin, 1, strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP,strPar.flagListKernels);
            stereo_check_strobe_and_self_simililarity_effect(input2, odistr, oselfdistr, omask20, oBin2, 1 , strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP, strPar.flagListKernels);
            
            
            //oBin.save("st_ssimilarity.pmf");

            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
            
        }
        
        
        
        // Perform Remove Isolated
        if (strPar.flagRemoveIsolated)
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            
            
            // USE THE CURRENT WINDOW
            //int win = (MAX(strPar.prolate.w(), strPar.prolate.h()) - 1) / 2;
            //stereo_remove_isolated_points(omask, win, oBin, strPar.valueRemoveIsolated);
            //stereo_remove_isolated_points(omask2, win, oBin2, strPar.valueRemoveIsolated);
            stereo_remove_isolated_points_window(omask, strPar.prolate, oBin, strPar.valueRemoveIsolated);
            stereo_remove_isolated_points_window(omask2, strPar.prolate, oBin2, strPar.valueRemoveIsolated);
            
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }




        // Perform Grain Filter
        if (strPar.flagGrainFilter)
        {
            
            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
            
            
            int minarea = strPar.valueGrainFilter;
//            printf("minarea %d\n", minarea);

            stereo_grain_filter(omask, minarea, oBin);
            stereo_grain_filter(omask2, minarea, oBin2);
            
            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
            
        }
        
        



    }



    static int good_modulus(int n, int p)
    {
       assert(p);
       if (p < 0) return good_modulus(n, -p);

       int r;
       if (n >= 0)
          r = n % p;
       else
       {
          r = p - (-n) % p;
          if (r == p)
             r = 0;
       }
       assert(r >= 0);
       assert(r < p);
       return r;
    }

void gen_raw_costvolume_tranche_simpler(cflimage *images1, cflimage *images2, int isteps,int y_base , int y_end, cflimage &cv, int dmin, int dmax) 
{
         cv = 0;
         int w1 = images1[0].w();
         int w2 = images2[0].w();
         int wh1 = images1[0].wh();
         int wh2 = images2[0].wh();
         int cv_w = cv.w();
         int cv_h = cv.h();
         int imch = images1[0].c();

         assert(cv_w == w1);
         assert(w1 == w2);

#pragma omp parallel for
         for (int d=dmin;d<=dmax;d++)
         for (int s2=0;s2<isteps;s2++) 
         { 
            int c = (d - dmin)*isteps + s2;

            if(c >= cv.c() || c < 0 ) { printf("MIERDA!\n"); continue; }
            float *curr_cv = cv.v(c);
            const float *im1 = images1[0].v();
            const float *im2 = images2[s2].v();

            //    compute cv(x,d)
            for (int y=y_base;y<y_end;y++)
            { 
               if(y-y_base >= cv_h) { printf("PUTA MIERDA!\n"); continue; }

               // x range adjusted to d
               for (int x=0; x<w1; x++) 
               {   
                  if(x+d<0 || x+d>=w2) continue;

                  double diff = 0;
                  for (int ch=0; ch<imch; ch++) {
                     double t = im1[ x +     y*w1 + wh1 * ch ] -
                                im2[ x + d + y*w2 + wh2 * ch ];
                     diff += t*t;  // pixel differences
                  }
                  curr_cv[ x + (y - y_base)*cv_w ] = diff;
               }

            }
         }
         // At this point cv contains the raw cost volume

}













void fiComputeIntegralImageINPLACEstride(float *d, int width, int height, int stridex, int stridey)
{
   // recursivity in the rows
   for(int jj=0; jj < height; jj++)
      for(int ii=1*stridex; ii < width; ii++)
         d[width*jj+ii] = d[width*jj+ii-1*stridex] + d[width*jj+ii];

   // recursivity in the columns
   for(int ii=0; ii < width; ii++)
      for(int jj=1*stridey; jj < height; jj++)
         d[width*jj+ii] = d[width*(jj-1*stridey)+ii] + d[width*jj+ii];
}

void fiComputeIntegralImageINPLACEstride_double(double *d, int width, int height, int stridex, int stridey)
{
   // recursivity in the rows
   for(int jj=0; jj < height; jj++)
      for(int ii=1*stridex; ii < width; ii++)
         d[width*jj+ii] = d[width*jj+ii-1*stridex] + d[width*jj+ii];

   // recursivity in the columns
   for(int ii=0; ii < width; ii++)
      for(int jj=1*stridey; jj < height; jj++)
         d[width*jj+ii] = d[width*(jj-1*stridey)+ii] + d[width*jj+ii];
}

void fiComputeIntegralImageINPLACE_double(double *d, int width, int height)
{
   // recursivity in the rows
   for(int jj=0; jj < height; jj++)
      for(int ii=1; ii < width; ii++)
         d[width*jj+ii] = d[width*jj+ii-1] + d[width*jj+ii];

   // recursivity in the columns
   for(int ii=0; ii < width; ii++)
      for(int jj=1; jj < height; jj++)
         d[width*jj+ii] = d[width*(jj-1)+ii] + d[width*jj+ii];
}

/**
 * This function only filters all the cost volume channels with a ww x wh window
 * The strides correspond to the subpixel scaling
 * */
void filter_costvolume_integralimage_stride(cflimage &cv, cflimage &ccv, int winsize[2], int strides[2]) 
         {
         //    filter cv(x,.) with ((3-1)*isteps+1)x3 windows->ccv(x,d)
         
         ccv=fLarge;

         const int ww = winsize[0];
         const int wh = winsize[1];
         const int stridex = strides[0];
         const int stridey = strides[1];

         int cv_h = cv.h();
         int cv_w = cv.w();
         int cv_wh = cv.wh();
#pragma omp parallel for
         for (int c=0;c<cv.c();c++)
         {

            float *curr_cv = cv.v(c);
            float *curr_ccv = ccv.v(c);

            double *tmp_cv = (double*) malloc(cv_wh*sizeof(double));
            for (int i=0;i<cv_wh;i++) tmp_cv[i] = curr_cv[i];

//            fiComputeIntegralImageINPLACEstride(curr_cv, cv_w, cv_h, stridex, stridey);
            fiComputeIntegralImageINPLACEstride_double(tmp_cv, cv_w, cv_h, stridex, stridey);


            // subwindow size
            const int halfww = ww/2;
            const int halfwh = wh/2;

            for (int y=stridey*halfwh;y<cv_h-stridey*halfwh;y++) 
               for (int x=stridex*halfww;x<cv_w-stridex*halfww;x++) {

                  // test if d is within Dmin Dmax
//                  int idx  = (int)(x*step) +     (y+y_base)*in1.w();
//                  int iidx = (int)(x*step) + d + (y+y_base)*in2.w();
                  //if ( (y+y_base) < in1.h() &&  
                  //     ( (fd <= Dmax[idx]  && fd >= Dmin[idx] ) ||
                  //       (-fd <= iDmax[iidx] && -fd >= iDmin[iidx] ) ))
                  if ( 1 )
                  {
                     int xp = x+stridex*halfww; 
                     int xm = x-stridex*(halfww+1); 
                     int yp = y+stridey*halfwh;
                     int ym = y-stridey*(halfwh+1);

                     double retval = tmp_cv[xp + yp*cv_w];
                     if(xm>=0)
                        retval -= tmp_cv[xm + yp*cv_w];
                     if(ym>=0)
                        retval -= tmp_cv[xp + ym*cv_w];
                     if(xm>=0 && ym>=0)
                        retval += tmp_cv[xm + ym*cv_w];

                     curr_ccv[x+y*cv_w] = retval;   // FIXME  error indexing when isteps != 1
                  }
               }
            free(tmp_cv);
         }
         // At this point ccv contains the 3x3 aggregated cost volume


}



/**
 * This function only filters all the cost volume channels with a ww x wh window
 * */
void filter_costvolume(cflimage &cv, cflimage &ccv, int winsize[2]) 
         {
         //    filter cv(x,.) with ((3-1)*isteps+1)x3 windows->ccv(x,d)
         
         ccv=fLarge;

         const int ww = winsize[0];
         const int wh = winsize[1];

         int cv_h = cv.h();
         int cv_w = cv.w();
#pragma omp parallel for
         for (int c=0;c<cv.c();c++)
         {

            float *curr_cv = cv.v(c);
            float *curr_ccv = ccv.v(c);

            // subwindow size
            const int halfww = ww/2;
            const int halfwh = wh/2;

            for (int y=halfwh;y<cv_h-halfwh;y++) 
               for (int x=halfww;x<cv_w-halfww;x++) {
                  double retval= 0;

                  for(int j=-halfwh;j<=halfwh;j++)
                  for(int i=-halfww;i<=halfww;i++) {
                     retval += curr_cv[x+i + (y+j)*cv_w];
                  }
                  curr_ccv[x + y*cv_w] = retval;
               }
         }
         // At this point ccv contains the 3x3 aggregated cost volume

}

/**
 * This function only filters all the cost volume channels with a ww x wh window
 * */
void filter_costvolume_integralimage(cflimage &cv, cflimage &ccv, int winsize[2]) 
         {
         //    filter cv(x,.) with ((3-1)*isteps+1)x3 windows->ccv(x,d)
         
         ccv=fLarge;

         const int ww = winsize[0];
         const int wh = winsize[1];

         int cv_h = cv.h();
         int cv_w = cv.w();
         int cv_wh = cv.wh();
#pragma omp parallel for
         for (int c=0;c<cv.c();c++)
         {

            float *curr_cv = cv.v(c);
            float *curr_ccv = ccv.v(c);

            double *tmp_cv = (double*) malloc(cv_wh*sizeof(double));
            for (int i=0;i<cv_wh;i++) tmp_cv[i] = curr_cv[i];

//            fiComputeIntegralImageINPLACEstride(curr_cv, cv_w, cv_h, stridex, stridey);
            fiComputeIntegralImageINPLACE_double(tmp_cv, cv_w, cv_h);


            // subwindow size
            const int halfww = ww/2;
            const int halfwh = wh/2;

            for (int y=halfwh;y<cv_h-halfwh;y++) 
               for (int x=halfww;x<cv_w-halfww;x++) {

                  // test if d is within Dmin Dmax
//                  int idx  = (int)(x*step) +     (y+y_base)*in1.w();
//                  int iidx = (int)(x*step) + d + (y+y_base)*in2.w();
                  //if ( (y+y_base) < in1.h() &&  
                  //     ( (fd <= Dmax[idx]  && fd >= Dmin[idx] ) ||
                  //       (-fd <= iDmax[iidx] && -fd >= iDmin[iidx] ) ))
                  if ( 1 )
                  {
                     int xp = x+halfww; 
                     int xm = x-(halfww+1); 
                     int yp = y+halfwh;
                     int ym = y-(halfwh+1);

                     double retval = tmp_cv[xp + yp*cv_w];
                     if(xm>=0)
                        retval -= tmp_cv[xm + yp*cv_w];
                     if(ym>=0)
                        retval -= tmp_cv[xp + ym*cv_w];
                     if(xm>=0 && ym>=0)
                        retval += tmp_cv[xm + ym*cv_w];

                     curr_ccv[x+y*cv_w] = retval;   // FIXME  error indexing when isteps != 1
                  }
               }
            free(tmp_cv);
         }
         // At this point ccv contains the 3x3 aggregated cost volume


}


void cleanup_ccv_simpler( cflimage &ccv, int subwindowsize[2], int isteps, int y_base, int y_end, int dmin, int dmax){
///////////////// CLEANUP ccv
         int cv_c = ccv.c();
         int cv_w = ccv.w();
         int cv_h = ccv.h();

#pragma omp parallel for
         for (int c=0;c<cv_c;c++) 
         {
            float *curr_ccv = ccv.v(c);

            float fd = dmin + ((float)c)/isteps;

            // subwindow size
            const int ww = subwindowsize[0] / 2;
            const int wh = subwindowsize[1] / 2;

            for (int y=0;y<cv_h;y++)  
            for (int x=0;x<cv_w;x++) {

                  //// estimate equivalent disparity for this setup of d,s1,s2
                  //float fd = d*step;
                  //// test if d is within Dmin Dmax
                  //int idx  = (int)(x*step) +     (y+y_base)*in1.w();
                  //int iidx = (int)(x*step) + d + (y+y_base)*in2.w();
                  //     ( (fd <= Dmax[idx]  && fd >= Dmin[idx] ) ||
                  //       (-fd <= iDmax[iidx] && -fd >= iDmin[iidx] ) ))
                  
                  if ( (y+y_base) >= y_end ||
                       y-wh<0 || y+wh>=cv_h || 
                       x-ww<0 || x+ww>=cv_w ||
                       x+fd-ww < 0  ||                //TODO: this line was here but seems to be wrong ..  x+fd-(ww+1) < 0  ||
                       x+fd+ww >= cv_w
                     ) 
                     curr_ccv[x+y*cv_w] = fLarge;   

            }
         }
}


/** 
 * ZSSDc is computed in constant time for any patch size by using integral images 
 * The different windows are computed by combining the costs of smaller subwindows.
 * Let SSD(P,Q) denote the SSD cost for two blocks P and Q (for color images the
 * costs of all channels are addedi), mu_i(P) is the mean value of the block for the channel i, 
 * and B the size (in pixels) of the block.
 * Then : 
 *    ZSSDc (P,Q) = SSD(P,Q) - B \sum_i |mu_i(P) - mu_i(Q)|^2
 *
 * And ZSSD for the union of two regions can be computed as:
 *    ZSSDc (P \cup P', Q \cup Q') = SSD(P,Q) + SSD(P',Q') - B/2 \sum_i |(mu_i(P) - mu_i(Q)) + (mu_i(P') - mu_i(Q'))|^2
 * */
void compute_winners_simplerZSSDc(cflimage &ccv, cflimage *means1, cflimage *means2, int subwindowsize[2],int isteps, int y_base, int y_end, int dmin, int dmax, flimage* tmpodisp, flimage* tmpodispr, flimage* tmpodist, flimage* tmpodistr, flimage* tmpomask, flimage* tmpomaskr, int nprol, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, int use_ZSSD_instead_of_SSD , flimage *prolates) 
{
   int w1 = tmpodisp[0].w();
   int w2 = tmpodispr[0].w();
   int wh1 = tmpodisp[0].wh();
   int wh2 = tmpodispr[0].wh();
   int cv_w = ccv.w();
   int cv_h = ccv.h();
   int cv_wh = ccv.wh();
   float step = 1.0/isteps;
   int nch = (means1[0]).c();

   const int swin = subwindowsize[0];
   const int win =swin* 3;
   const int halfwin = win / 2;



   // COMPUTE THE WINNER FOR EACH PIXEL 
   //    for w in window_configurations (max 9) 
   for (int c=0;c<ccv.c();c++) 
   {
      // estimate equivalent disparity for this setup of d,s1,s2
      float fd = c*step+dmin;
      int s2 = good_modulus(c,isteps);
      int d  = dmin + (c - s2)/isteps;

      float *curr_ccv = ccv.v(c);
      //float *curr_cv  = cv.v(c);


      float *mu1 = means1[0].v();
      float *mu2 = means2[s2].v();

      cflimage meandifferences = cflimage(cv_w,cv_h,nch);
      meandifferences = 0;
#pragma omp parallel for
      for (int y=y_base;y<y_end;y++)
         for (int x=0;x<cv_w;x++) 
         {
            //if(x+d>=halfwin && x+d<cv_w-halfwin && x>=halfwin && x+halfwin<cv_w ) {
            if(x+d>=0 && x+d<cv_w) {
               for(int t=0;t<nch ;t++) {
                  float tmp = mu1[x     + y*w1 + wh1*t ] 
                     - mu2[x + d + y*w2 + wh2*t ];
                  meandifferences[x + (y-y_base)*cv_w + cv_wh*t ] = tmp;
               }
            }
         }



#pragma omp parallel for
         for (int y=y_base+halfwin;y<y_end-halfwin;y++)
            for (int x=halfwin;x<cv_w-halfwin;x++) 
               if(x+d>=halfwin && x+d<cv_w-halfwin && x>=halfwin && x+halfwin<cv_w )   // near the boundaries the SSD aren't computed while averages are averages are, this may produce costs lower than zero. To avoid this we only compute if the support window is completely contained in the image.
               {

                  int idx  = x    + y * w1;
                  int iidx = x+fd + y * w2;

                  double costs_win[9];

                  // compute: SSD(P,Q) + SSD(P',Q')
                  for (int i=0;i<nprol;i++) {
                     costs_win[i] = 0;
                     flimage *curr_prolate = &prolates[i];
                     for (int p=0;p<curr_prolate->list_len;p++) 
                     {
                        int dx = curr_prolate->offx[p];
                        int dy = curr_prolate->offy[p];
                        costs_win[i] += curr_ccv[ (x+dx) + (y+dy - y_base)*cv_w ];
                     }
                  }

                  // compute: - |B|/2 \sum_i |(mu_i(P) - mu_i(Q)) + (mu_i(P') - mu_i(Q'))|^2
                  // where |B| is the number of pixels in the block
                  // and 2 is the nuber of blocsk being aggregated
                  if(use_ZSSD_instead_of_SSD) {
                     double blocksize = swin*swin;
                     for(int t=0;t<nch ;t++) 
                     {


                        for (int i=0;i<nprol;i++) {
                           flimage *curr_prolate = &prolates[i];
                           double tmp = 0;
                           for (int p=0;p<curr_prolate->list_len;p++) 
                           {
                              int dx = curr_prolate->offx[p];
                              int dy = curr_prolate->offy[p];
                              tmp +=  meandifferences[ (x+ dx) + (y+dy- y_base)*cv_w + cv_wh*t ];
                           }
                           costs_win[i] -= tmp*tmp*blocksize/curr_prolate->list_len;
                        }

                     }
                  }


                  // normalize costs
                  double nn = 1.0/((double)(swin*swin*nch));
                  for(int i=0;i<nprol;i++){
                     costs_win[i]*=nn;
                     costs_win[i]/=prolates[i].list_len;
                     // max with 0
                     //   costs_win[i] = fmax(0, costs_win[i]);
                  }



                  //       update odisp[d] ocost[d] odispr[d] ocostr[d]
                  for (int i=0;i<nprol;i++){
                     //#pragma omp critical   <<<< only needed if for c is parallel
                     if(tmpodist[i][idx] > costs_win[i]) 
                        if (fd <= Dmax[idx]  && fd >= Dmin[idx] ) {
                           tmpomask[i][idx] = 1;
                           tmpodist[i][idx] = costs_win[i];
                           tmpodisp[i][idx] = fd;
                        } else {
//                           tmpomask[i][idx] = 2;
//                           tmpodist[i][idx] = costs_win[i];
//                           tmpodisp[i][idx] = fd;
                        }


                     if(s2 == 0 && x+fd >= 0 && x+fd < cv_w)  {
                        // This inequality is complicated to explain.
                        // When the image is textureless the left disparity map keeps the lowest disparity
                        // and the right map the opposite. To avoid this corner case, we allow the right 
                        // map to be take the lowest value so that these regions are removed 
                        //#pragma omp critical   <<<< only needed if for c is parallel
                        if(tmpodistr[i][iidx] >= costs_win[i]) 
                           if (-fd <=iDmax[iidx] && -fd >= iDmin[iidx] )  {
                              tmpomaskr[i][iidx] = 1;
                              tmpodistr[i][iidx] = costs_win[i];
                              tmpodispr[i][iidx] = -fd;
                           } else {
//                              tmpomaskr[i][iidx] = 2;
//                              tmpodistr[i][iidx] = costs_win[i];
//                              tmpodispr[i][iidx] = -fd;
                           }
                     }
                  }

               }

         }
   }





























void compute_subpixel_match_multiple_windows_one_direction( cflimage *images1, cflimage *images2, cflimage *means1, cflimage *means2, int dmin, int dmax, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, flimage *tmpodisp, flimage *tmpodispr, flimage *tmpodist, flimage *tmpodistr, flimage *tmpomask, flimage *tmpomaskr, strParameters &strPar, flimage *prolates)
{

      // number of window orientations
      const int nprol = fmin(9,strPar.flagListKernels);

      int isteps = strPar.inPrecisions;

      int win= strPar.prolate.w(); //window size

      int subwindowsize[2] = {3,3};
      if (win==5) {
         subwindowsize[0] = subwindowsize[1] = 3;
      }else if (win==9) {
         subwindowsize[0] = subwindowsize[1] = 5;
      }else if (win==13) {
         subwindowsize[0] = subwindowsize[1] = 7;
      }else if (win==17) {
         subwindowsize[0] = subwindowsize[1] = 9;
      }else if (win==21) {
         subwindowsize[0] = subwindowsize[1] = 11;
      }

      int support_window = subwindowsize[0] * 3;

      printf("windows=%d support=%d subwindow=%d drange=(%d %d)\t lines:", nprol, support_window, subwindowsize[0], dmin, dmax);

      // check if the range is feasible
      if(dmax<dmin) {
         printf(" nothing to do\n");
         return;
      }

      // costvolume buffer parameters 
      const int cv_h = support_window*6;
      const int cv_w = images1[0].w();
      const int cv_c = (dmax-dmin+1)*isteps; 
      const int cv_wh=cv_w*cv_h;

      const int image_height = images1[0].h();
   
      // allocate costvolume buffer
      libIIPStable::cflimage cv = cflimage(cv_w,cv_h,cv_c);
      libIIPStable::cflimage ccv = cflimage(cv_w,cv_h,cv_c);


      // process by bands // win  is the overlap needed between bands
      // TODO REVIEW UPPER LIMIT
      for (int y_base=0;y_base<image_height;y_base+=cv_h-support_window) 
      {

         int y_end = fmin(y_base+cv_h,image_height);

         printf("[%d %d] ", y_base, y_end); fflush(stdout);

         // POPULATE the raw cost volume cv for all subpixel shifts: d in [dmin..dmax] (subpixel)
         cv = 0;
         gen_raw_costvolume_tranche_simpler(images1, images2, isteps, y_base , y_end, cv, dmin, dmax);
         
         filter_costvolume_integralimage(cv, ccv, subwindowsize);
         //filter_costvolume(cv, ccv, subwindowsize);

         // Remove ccv values that are incorrect because the get out of the image
         cleanup_ccv_simpler( ccv, subwindowsize, isteps, y_base, y_end, dmin, dmax);

         compute_winners_simplerZSSDc(ccv, means1, means2, subwindowsize,isteps, y_base, y_end, dmin, dmax, tmpodisp, tmpodispr, tmpodist, tmpodistr, tmpomask, tmpomaskr, nprol, Dmin, Dmax, iDmin, iDmax, strPar.itypeDist == L2M, prolates);

      }
      printf("\n");

}
      




    void stereo_pixel_chain_multiwinCV_simpler(cflimage &input, cflimage &input2, flimage &Dmin, flimage &Dmax, flimage &iDmin, flimage &iDmax, 
          flimage &odisp, flimage &odispr, flimage &odist, flimage &odistr, flimage &omask, flimage &omaskr, strParameters &strPar) // << attention these are also outputs
   {
      const int dmin = strPar.dmin0;
      const int dmax = strPar.dmax0;

      // number of window orientations
      const int nprol=fmin(9,strPar.flagListKernels);


      int win= strPar.prolate.w(); //window size

      int subwindowsize[2] = {3,3};
      if (win==5) {
         subwindowsize[0] = subwindowsize[1] = 3;
      }else if (win==9) {
         subwindowsize[0] = subwindowsize[1] = 5;
      }else if (win==13) {
         subwindowsize[0] = subwindowsize[1] = 7;
      }else if (win==17) {
         subwindowsize[0] = subwindowsize[1] = 9;
      }else if (win==21) {
         subwindowsize[0] = subwindowsize[1] = 11;
      }

      int support_window = subwindowsize[0] * 3;
      int half_support_window = support_window/2;

      printf("%d %d %d\n", nprol, support_window, subwindowsize[0]);

      // this is used for update Dmin Dmax
      flimage prolate; prolate.create(support_window, support_window);  prolate = 1.0f;  prolate.normalizeL1();



      // the maximum number of window orientations we will ever consider is 9
      libIIPStable::flimage *prolates = new libIIPStable::flimage[9];
      for (int i=0;i<9;i++) {
         prolates[i].create(support_window, support_window); prolates[i] = 0;
      }



      ////// GENERATE THE PROLATE WINDOWS AND THE LISTS NEEDED FOR EVALUATING THEM
      //////
      //////
      {

         int wrange = (support_window/2 - (subwindowsize[0]/2)) ;

         // generate the lists 
         for (int i=0;i<nprol;i++) {
            prolates[i].create(support_window, support_window); 
            prolates[i] = 0;
            prolates[i].offx     = new int[5];
            prolates[i].offy     = new int[5];
            prolates[i].offval   = new float[5];
            prolates[i].list_len = 0;

            // special case : square, centered window with 4,5 or 9 elements
            if(i==0) {
               int DEBUG_CASE_USE_ONLY_ONE_CENTERED_SUBWINDOW=0;           // TODO REMOVE
               if(DEBUG_CASE_USE_ONLY_ONE_CENTERED_SUBWINDOW) {
                  int ox0[5] = {0};
                  int oy0[5] = {0};
                  int ov0[5] = {1};

                  prolates[i].list_len = 1;        // JUST THE CENTRAL SUBWINDOW! 

                  for(int t=0;t<prolates[i].list_len;t++){
                     prolates[i].offx[t]   = ox0[t]*wrange/2;
                     prolates[i].offy[t]   = oy0[t]*wrange/2;
                     prolates[i].offval[t] = ov0[t];
                  }
               }else{
                  int ox0[9] = {-1, 1,-1, 1, 0,-1,1, 0,0};
                  int oy0[9] = {-1,-1, 1, 1, 0, 0,0,-1,1};
                  int ov0[9] = { 1, 1, 1, 1, 1, 1,1, 1,1};

                  prolates[i].list_len = 4; // 4 yields a FLATTER window rather than 5 or 9

                  for(int t=0;t<prolates[i].list_len;t++){
                     prolates[i].offx[t]   = ox0[t]*wrange/2;
                     prolates[i].offy[t]   = oy0[t]*wrange/2;
                     prolates[i].offval[t] = ov0[t];
                  }
               }
            } else { // all the other cases 3 elements

               prolates[i].list_len = 3;

               float fangle = 0;
               int transpose=0;
               if(i==1 || i==2) fangle = 0;
               if(i==3) fangle = M_PI/4;
               if(i==4) fangle = -M_PI/4;
               if(i==5 || i==6) fangle = M_PI/8;
               if(i==7 || i==8) fangle = -M_PI/8;

               if(i==2 || i==3 || i==6 || i==7) transpose =1;

               int t=0;
               for(int j=-wrange;j<=wrange;j+=wrange) { 
                  int ctrx = round(j*cosf(fangle));
                  int ctry = round(j*sinf(fangle));
                  prolates[i].offx[t]   = transpose?ctry:ctrx;
                  prolates[i].offy[t]   = transpose?ctrx:ctry;
                  prolates[i].offval[t] = 1;
                  t++;
               }
            }
         }


         // generate the prolates
         for(int f=0; f<nprol; f++)   {
            for(int t=0;t<prolates[f].list_len;t++) { 
               int ctrx = half_support_window + prolates[f].offx[t];
               int ctry = half_support_window + prolates[f].offy[t];

               int ww = subwindowsize[0]/2;
               for(int y=-ww;y<=ww;y++)
                  for(int x=-ww;x<=ww;x++)
                     prolates[f][ctrx+x + (ctry+y)*support_window] += 1;
            }
            prolates[f].normalizeL1();
         }



         for(int f=0; f<nprol; f++) {
            char stri[1024];
            sprintf(stri, "prolate%d.tif", f);
            prolates[f].save(stri);
         }

      }






      





      // initialize disp and costs vectors for each window configuration : odisp[0-5] ocost[0-5]   odispr[0-5] ocostr[0-5] 
      libIIPStable::flimage *tmpodisp = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *tmpodispr = new  libIIPStable::flimage[nprol];

      libIIPStable::flimage *tmpomask = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *tmpomaskr = new  libIIPStable::flimage[nprol];

      libIIPStable::flimage *tmpodist = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *tmpodistr = new  libIIPStable::flimage[nprol];


      for (int ii=0; ii < nprol; ii++)
      {

         tmpodisp[ii].create(input.w(), input.h());         tmpodisp[ii] = 0.0f;
         tmpodispr[ii].create(input2.w(), input2.h());      tmpodispr[ii] = 0.0f;

         tmpomask[ii].create(input.w(), input.h());         tmpomask[ii] = 0.0f;
         tmpomaskr[ii].create(input2.w(), input2.h());      tmpomaskr[ii] = 0.0f;

         tmpodist[ii].create(input.w(), input.h());         tmpodist[ii] = fLarge;
         tmpodistr[ii].create(input2.w(), input2.h());      tmpodistr[ii] = fLarge;

      }






		//! Compute translation of second image for subpixel estimation
		//! It does not compute anything if inPrecisions = 1 but
        //! creates images2[0] = input2
      int isteps = strPar.inPrecisions;
		float step = 1.0 / (float) isteps;

      cflimage *images2 = new cflimage[isteps];
		images2[0] = input2;
		for (int ii=1; ii < isteps; ii++)
		{
			images2[ii] = input2;
			images2[ii].fftShear(0.0, iipHorizontal, (float) ii * step, 0, 1);
		}
        
      //! Compute translation of first image for subpixel estimation
      cflimage *images1 = new cflimage[isteps];
      images1[0] = input;
      for (int ii=1; ii < isteps; ii++)
      {
          images1[ii] = input;
          images1[ii].fftShear(0.0, iipHorizontal, (float) ii * step, 0, 1);
      }

      // precompute means for ZSSD
      cflimage *means1 = new cflimage[isteps];
      cflimage *means2 = new cflimage[isteps];
      for (int ii=0; ii < isteps; ii++)
      {
         flimage ker; ker.create(subwindowsize[0], subwindowsize[1]);  ker= 1.0f;  ker.normalizeL1();
         means1[ii] = images1[ii].patchMean(ker);
         means2[ii] = images2[ii].patchMean(ker);
      }





      // MATCH WITH ALL THE WINDOWS

      compute_subpixel_match_multiple_windows_one_direction( images1, images2, means1, means2, dmin, dmax, Dmin, Dmax, iDmin, iDmax, tmpodisp, tmpodispr, tmpodist, tmpodistr, tmpomask, tmpomaskr, strPar, prolates);


      if(isteps>1)  {
         for (int ii=0;ii<nprol;ii++) { tmpodistr[ii] = fLarge; tmpomaskr[ii]=0; }
         compute_subpixel_match_multiple_windows_one_direction( images2, images1, means2, means1, -dmax, -dmin, iDmin, iDmax, Dmin, Dmax, tmpodispr, tmpodisp, tmpodistr, tmpodist, tmpomaskr, tmpomask, strPar, prolates);
      }



      





      //// SELF SIMILARITY

      libIIPStable::flimage *tmposelfdist  = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *tmposelfdistr = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *ttrash        = new  libIIPStable::flimage[nprol];

      for (int ii=0; ii < nprol; ii++)
      {

         tmposelfdist[ii].create(input.w(), input.h());         tmposelfdist[ii] = fLarge;
         tmposelfdistr[ii].create(input2.w(), input2.h());      tmposelfdistr[ii] = fLarge;

      }


      if(strPar.flagSelfSimilarity)
      {

      libIIPStable::flimage *ttmpodisp = new  libIIPStable::flimage[nprol];
      libIIPStable::flimage *ttmpodispr = new  libIIPStable::flimage[nprol];

      for (int ii=0; ii < nprol; ii++)
      {

         ttmpodisp[ii].create(input.w(), input.h());         ttmpodisp[ii] = 0.0f;
         ttmpodispr[ii].create(input2.w(), input2.h());      ttmpodispr[ii] = 0.0f;

         ttrash[ii].create(input2.w(), input2.h());      ttrash[ii] = fLarge;

      }




      int old_selfsim = 0;

      if (old_selfsim)
      {
         // Use the range of the right image for the self-similarity test
         // This has a clear problem the image may not be self-similar in the considered range.
         // (this is the one used in the paper)

         if(dmax > 0)
            compute_subpixel_match_multiple_windows_one_direction( images1, images1, means1, means1, 
                  fmax(dmin,2) , dmax, Dmin, Dmax, iDmin, iDmax, 
                  ttmpodisp, ttmpodispr, tmposelfdist, ttrash, ttrash, ttrash, strPar, prolates);
         if(dmin < 0)
            compute_subpixel_match_multiple_windows_one_direction( images1, images1, means1, means1, 
                  dmin , fmin(dmax,-2), Dmin, Dmax, iDmin, iDmax, 
                  ttmpodisp, ttmpodispr, tmposelfdist, ttrash, ttrash, ttrash, strPar, prolates);


         if(-dmin > 0)
            compute_subpixel_match_multiple_windows_one_direction( images2, images2, means2, means2, 
                  fmax(-dmax,2) , -dmin, Dmin, Dmax, iDmin, iDmax, 
                  ttmpodispr, ttmpodisp, tmposelfdistr, ttrash, ttrash, ttrash, strPar, prolates);
         if(-dmax < 0)
            compute_subpixel_match_multiple_windows_one_direction( images2, images2, means2, means2, 
                  -dmax, fmin(-dmin,-2), iDmin, iDmax, Dmin, Dmax, 
                  ttmpodispr, ttmpodisp, tmposelfdistr, ttrash, ttrash, ttrash, strPar, prolates);

      } 
      else 
      { 
         // The self-similarity range is centered at the reference pixel
         // This permits to spot self similar structures even if the search range is shifted

         flimage sDmin = Dmin;
         flimage sDmax = Dmax;
         for (int i=0;i<sDmin.wh();i++){
            int r = fmax(ceil((Dmax[i]-Dmin[i])/2.0),3.0);
            sDmin[i] = -r;
            sDmax[i] = r;
         }

         compute_subpixel_match_multiple_windows_one_direction( images1, images1, means1, means1, 
               2 , dmax , sDmin, sDmax, sDmin, sDmax, 
               ttmpodisp, ttmpodispr, tmposelfdist, ttrash, ttrash, ttrash, strPar, prolates);
         compute_subpixel_match_multiple_windows_one_direction( images1, images1, means1, means1, 
               dmin, -2 , sDmin, sDmax, sDmin, sDmax, 
               ttmpodisp, ttmpodispr, tmposelfdist, ttrash, ttrash, ttrash, strPar, prolates);



         sDmin = iDmin;
         sDmax = iDmax;
         for (int i=0;i<sDmin.wh();i++){
            int r = fmax(ceil((iDmax[i]-iDmin[i])/2.0),3.0);
            sDmin[i] = -r;
            sDmax[i] = r;
         }

         compute_subpixel_match_multiple_windows_one_direction( images2, images2, means2, means2,  
               2 , dmax, sDmin, sDmax, sDmin, sDmax, 
               ttmpodispr, ttmpodisp, tmposelfdistr, ttrash, ttrash, ttrash, strPar, prolates);
         compute_subpixel_match_multiple_windows_one_direction( images2, images2, means2, means2, 
               dmin, -2, sDmin, sDmax, sDmin, sDmax, 
               ttmpodispr, ttmpodisp, tmposelfdistr, ttrash, ttrash, ttrash, strPar, prolates);
      }


      }




      // APPLYING THE FILTERS TO ALL THE ALL THE MAPS

      for (int i=0;i<nprol;i++){
         strPar.prolate = prolates[i];
         strPar.flagListKernels = 0;

         apply_filters_chain(input, input2, 
            tmpodisp[i], tmpodispr[i], tmpodist[i], tmpodistr[i], tmposelfdist[i], tmposelfdistr[i], tmpomask[i], tmpomaskr[i], strPar);
      }





      // COMBINING THE MAPS FOR DIFFERENT WINDOW ORIENTATIONS


		//! Computing choice and additional criteria
		odist = fLarge;
		odisp = 0.0f;
		omask = 0.0f;
		
		odistr = fLarge;
		odispr = 0.0f;
		omaskr = 0.0f;
		
		libIIPStable::flimage oChoice(input.w(), input.h());      oChoice = -1.0f;

		
		//! Select minimum distance
		for (int ii=0; ii < nprol; ii++)
		{
			
         for(int jj=0; jj < tmpodisp[ii].wh(); jj++) {
            if ( tmpodist[ii][jj] < odist[jj]) {
               odist[jj] = tmpodist[ii][jj];
               if (  tmpomask[ii][jj] > 0.0f )
               {
                  odisp[jj] = tmpodisp[ii][jj];
                  omask[jj] = tmpomask[ii][jj];
                  oChoice[jj] = (float) ii;
               }
            }
         }
			
			
         for(int jj=0; jj < tmpodispr[ii].wh(); jj++) {
            if ( tmpodistr[ii][jj] < odistr[jj]) {
               odistr[jj] = tmpodistr[ii][jj];
               if (  tmpomaskr[ii][jj] > 0.0f )
               {
                  odispr[jj] = tmpodispr[ii][jj];
                  omaskr[jj] = tmpomaskr[ii][jj];
                  //oChoice[jj] = (float) ii;
               }
            }
         }
			
		}

		
      odist.save("debug_dists_left.tif");
      odistr.save("debug_dists_right.tif");
		


		//! Check left right consistency
		//! reciprocity flag
		if (strPar.flagRecip && strPar.valueRecip >= 0.0f)
		{
			
			libIIPStable::flimage oBin(input.w(), input.h());     oBin = 255.0f;
			libIIPStable::flimage oTMP(input.w(), input.h());     oTMP = fLarge;
			
			libIIPStable::flimage oBin2(input2.w(), input2.h());     oBin2 = 255.0f;
			libIIPStable::flimage oTMP2(input2.w(), input2.h());     oTMP2 = fLarge;
			
			
			stereo_check_pixelian_reciprocity(odisp, omask, odispr, omaskr, oBin, strPar.valueRecip, oTMP);
			stereo_check_pixelian_reciprocity(odispr, omaskr, odisp, omask, oBin2, strPar.valueRecip, oTMP2);
			
			
			for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
			for (int i=0; i < omaskr.wh(); i++) if (oBin2[i] <= 0) omaskr[i] = oBin2[i];
			
		}
		
		
		

		//! Remove isolated points
		if (strPar.flagRemoveIsolated)
		{
			
			float valueRemoveIsolated;
			
			
			libIIPStable::flimage oBin(input.w(), input.h());     oBin = 255.0f;
			libIIPStable::flimage oBin2(input2.w(), input2.h());     oBin2 = 255.0f;
			
			
         // THIS ONE IS WITH SQUARED WINDOWS BY DESIGN! WITH SIZE DOUBLE OF THE WINDOW
			int pwin = strPar.prolate.w();
			stereo_remove_isolated_points(omask, pwin, oBin, strPar.valueRemoveIsolated);
			stereo_remove_isolated_points(omaskr, pwin, oBin2, strPar.valueRemoveIsolated);
			
			for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
			for (int i=0; i < omaskr.wh(); i++) if (oBin2[i] <= 0) omaskr[i] = oBin2[i];
			
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
		
		
		//// Update min and max
      // STILL WORSE THAN THE MODIFICATION
		//update_dmin_dmax(Dmin,Dmax,odisp,omask, strPar.dmin0, strPar.dmax0, support_window, support_window);
		//update_dmin_dmax(iDmin,iDmax,odispr,omaskr, strPar.idmin0, strPar.idmax0, support_window, support_window);



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

		   update_dmin_dmax_window(ttmpDmin,ttmpDmax,odisp,omask, strPar.dmin0, strPar.dmax0, prolates[ii]);
		   update_dmin_dmax_window(ttmpiDmin,ttmpiDmax,odispr,omaskr, strPar.idmin0, strPar.idmax0, prolates[ii]);

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





      
      // delete stuff
		delete[] tmpodisp;
		delete[] tmpodispr;
		
		delete[] tmpomask;
		delete[] tmpomaskr;
		
		delete[] tmpodist;
		delete[] tmpodistr;

      delete [] images2;
      delete [] images1;

      delete[] means2;
      delete[] means1;
		
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
        
        
         apply_filters_chain(input, input2, 
            odisp, odisp2, odist, odistr, oselfdist , oselfdistr, omask, omask2, strPar);

       
        update_dmin_dmax(Dmin,Dmax,odisp,omask, strPar.dmin0, strPar.dmax0, strPar.prolate.w(), strPar.prolate.h());
        update_dmin_dmax(iDmin,iDmax,odisp2,omask2, strPar.idmin0, strPar.idmax0, strPar.prolate.w(), strPar.prolate.h());

        
         return;



//        //! Save original masks
//        libIIPStable::flimage omask0 =omask;
//        libIIPStable::flimage omask20 =omask2;
//        
//        
//        //! Perform Reciprocal checking if flag activated
//        if (strPar.flagRecip)
//        {
//            
//            float fThreshold = strPar.valueRecip;  // to assure pixelian precision
//            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
//            
//            flimage oTMP(input.w(), input.h()); 	oTMP = fLarge;
//            flimage oTMP2(input2.w(), input2.h()); 	oTMP2 = fLarge;
//            
//            stereo_check_pixelian_reciprocity(odisp, omask0, odisp2, omask20, oBin, fThreshold, oTMP);
//            stereo_check_pixelian_reciprocity(odisp2, omask20, odisp, omask0, oBin2, fThreshold, oTMP2);
//            
//            //oBin.save("st_reciprocity.pmf");
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//        }
//        
//        
//        
//        
//        
//        //! Perform Reciprocal checking if flag activated
//        if (strPar.flagRecipDual)
//        {
//            
//            float fThreshold = 1.25;  // to assure pixelian precision
//            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
//            
//            stereo_check_pixelian_reciprocity_dual(odisp, omask0, odisp2, omask20, oBin, fThreshold);
//            stereo_check_pixelian_reciprocity_dual(odisp2, omask20, odisp, omask0, oBin2, fThreshold);
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//        }
//        
//        
//        
//        
//        //! Perform Reciprocal checking if flag activated
//        if (strPar.flagiRecip)
//        {
//            
//            float fThreshold = 1.25;  // to assure pixelian precision
//            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
//            
//            stereo_check_inverse_pixelian_reciprocity(odisp, omask0, odisp2, omask20, oBin, fThreshold);
//            stereo_check_inverse_pixelian_reciprocity(odisp2, omask20, odisp, omask0, oBin2, fThreshold);
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//        }
//        
//        
//        
//        //! Perform Variance test
//        if (strPar.flagVariance)
//        {
//            
//            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
//            
//            flimage oTMP(input.w(), input.h()); 	oTMP = fLarge;
//            flimage oTMP2(input2.w(), input2.h()); 	oTMP2 = fLarge;
//            
//            
//            stereo_check_horizontal_variance( input, oBin, strPar.valueVariance, strPar.fSigma, 1 ,strPar.prolate.w(), strPar.prolate.h(), oTMP);
//            
//            stereo_check_horizontal_variance( input2, oBin2, strPar.valueVariance, strPar.fSigma, 1 , strPar.prolate.w(), strPar.prolate.h(), oTMP2);
//            
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//        }
//        
//        
//        
//        
//        //! Perform Derivative test
//        if (strPar.flagHorInfo)
//        {
//            
//            flimage oBin(input.w(), input.h()); 	oBin = 255.0f;
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 255.0f;
//            
//            stereo_check_integral_derivatives(input, oBin, strPar.valueHorInfo, strPar.fSigma, strPar.prolate);
//            stereo_check_integral_derivatives(input2, oBin2, strPar.valueHorInfo, strPar.fSigma, strPar.prolate);
//            // dilate the rejected points by the size of the correlation window : invalidates the whole patch
//            dilate_mask_by_window_support(oBin, strPar.prolate) ;
//            dilate_mask_by_window_support(oBin2, strPar.prolate) ;
//            
//            
//            //oBin.save("st_horinfo.pmf");
//            
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//        }
//        
//        
//        
//        
//        
//        //! Perform MinDist
//        if (strPar.flagMinDist)
//        {
//            
//            
//            //float fThreshold = 1.0 / (float) strPar.inPrecisions;
//            float fThreshold = strPar.valueMinDist;
//            
//            
//            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
//            
//            flimage oTMP(input.w(), input.h()); 	oTMP = fLarge;
//            flimage oTMP2(input2.w(), input2.h()); 	oTMP2 = fLarge;
//            
//            //libIIPStable::flimage prolate;
//            //prolate.create(3,3);
//            //prolate=1;
//        
//            //stereo_check_min_dist( odisp, odist,  omask0,  oBin, fThreshold, strPar.prolate, oTMP);
//            //stereo_check_min_dist( odisp2, odistr,  omask20,  oBin2, fThreshold, strPar.prolate, oTMP2);
//            // CONSIDER THE REJECTIONS OF THE CURRENT omask 
//            stereo_check_min_dist( odisp, odist,  omask,  oBin, fThreshold, strPar.prolate, oTMP);
//            stereo_check_min_dist( odisp2, odistr,  omask2,  oBin2, fThreshold, strPar.prolate, oTMP2);
//            
//            
//            //! oBin.save("st_mindist.pmf");
//            
//            //! Dilate of 1 pixel
//            if (strPar.flagMinDistDilatation)
//            {
//                flimage itmp = oBin;
//                flimage itmp2 = oBin2;
//                float fRadius = 1.0;
//                fiPatchMin(itmp.v(), oBin.v(), fRadius, itmp.w(), itmp.h());
//                fiPatchMin(itmp2.v(), oBin2.v(), fRadius, itmp2.w(), itmp2.h());
//            }
//            
//            
//            //oBin.save("st_mindist.pmf");
//
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//        }
//        
//        
//        
//        //! Perform SelfSimilarity
//        if (strPar.flagSelfSimilarity)
//        {
//            
//            float fTransFactor = 1.0 / (float) (2 * strPar.inPrecisions);
//            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
//            flimage oTMP(input.w(), input.h()); 	oTMP = 255.0f;
//            
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
//            
//            
//            stereo_check_strobe_and_self_simililarity_effect(input, odist, oselfdist, omask0, oBin, 1, strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP,strPar.flagListKernels);
//            stereo_check_strobe_and_self_simililarity_effect(input2, odistr, oselfdistr, omask20, oBin2, 1 , strPar.itypeDist, fTransFactor, strPar.fSigma,  strPar.valueSelfSimilarity,  strPar.prolate, oTMP, strPar.flagListKernels);
//            
//            
//            //oBin.save("st_ssimilarity.pmf");
//
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//            
//        }
//        
//        
//        
//        // Perform Remove Isolated
//        if (strPar.flagRemoveIsolated)
//        {
//            
//            flimage oBin(input.w(), input.h()); 	oBin = 0.0f;
//            flimage oBin2(input2.w(), input2.h()); 	oBin2 = 0.0f;
//            
//            
//            
//            // USE THE CURRENT WINDOW
//            //int win = (MAX(strPar.prolate.w(), strPar.prolate.h()) - 1) / 2;
//            //stereo_remove_isolated_points(omask, win, oBin, strPar.valueRemoveIsolated);
//            //stereo_remove_isolated_points(omask2, win, oBin2, strPar.valueRemoveIsolated);
//            stereo_remove_isolated_points_window(omask, strPar.prolate, oBin, strPar.valueRemoveIsolated);
//            stereo_remove_isolated_points_window(omask2, strPar.prolate, oBin2, strPar.valueRemoveIsolated);
//            
//            
//            for (int i=0; i < omask.wh(); i++) if (oBin[i] <= 0) omask[i] = oBin[i];
//            for (int i=0; i < omask2.wh(); i++) if (oBin2[i] <= 0) omask2[i] = oBin2[i];
//            
//        }
//        
//       
//        update_dmin_dmax(Dmin,Dmax,odisp,omask, strPar.dmin0, strPar.dmax0, strPar.prolate.w(), strPar.prolate.h());
//        update_dmin_dmax(iDmin,iDmax,odisp2,omask2, strPar.idmin0, strPar.idmax0, strPar.prolate.w(), strPar.prolate.h());
        
        
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
    
    
    
    
    
    
    
    ///////////////////////
    ///////////////////////  MINDIST FILTER
    ///////////////////////
    
    
	
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





