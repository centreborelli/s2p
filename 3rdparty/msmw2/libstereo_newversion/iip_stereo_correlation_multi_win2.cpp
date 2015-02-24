#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../libraryStable/libraryBasic.h"

#include "libstereo.h"

using namespace std;




int main(int argc, char **argv)
{
	
	
	vector <OptStruct *> options;
    
    
    OptStruct om = {"m:", 0, NULL, NULL,  "minimum displacement (can be an image)"};  options.push_back(&om);
	OptStruct oM = {"M:", 0, NULL, NULL,  "maximum displacement (can be an image)"};  options.push_back(&oM);
    
    //! Option for distance, scale and precision
    OptStruct oSeg = {"S:", 0, NULL, NULL, "optional input segmentation (NOT WORKING YET)"}; options.push_back(&oSeg);
	OptStruct oK = {"k:", 0, NULL, NULL, "input correlation window"}; options.push_back(&oK);
	
    OptStruct oX = {"x:", 0, NULL, NULL, "window x size"}; options.push_back(&oX);
	OptStruct oY = {"y:", 0, NULL, NULL, "window y size"}; options.push_back(&oY);
	OptStruct oW = {"w:", 0, "0", NULL, "flag for weighting window"}; options.push_back(&oW);
	OptStruct oWL = {"W:", 0, "5", NULL, "flag for using windows as lists (5x5 windows only)\nthe number (!=0) indicates how many of the orientations should be considered "}; options.push_back(&oWL);
    OptStruct oI = {"i:", 0, "1", NULL, "type of distance"}; options.push_back(&oI);
	OptStruct oP = {"p:", 0, "1", NULL, "number of precisions for single scale"}; options.push_back(&oP);
	OptStruct oPrefine = {"P:", 0, "1", NULL, "factor of disparity refinement by cubic interpolation"}; options.push_back(&oPrefine);
	OptStruct oN = {"n:", 0, "3", NULL, "number of scales"}; options.push_back(&oN);
    
    //! Noise std
    OptStruct oF = {"f:", 0, NULL, NULL, "noise standard deviation"}; options.push_back(&oF);
    
    //! Options for checking LR-RL
    OptStruct oR = {"r:", 0, NULL, NULL, "reciprocity flag and value"}; options.push_back(&oR);
    OptStruct oG = {"g:", 0, NULL, NULL, "subpixel reciprocity flag"}; options.push_back(&oG);
    OptStruct oRD = {"R:", 0, NULL, NULL, "reciprocity dual flag"}; options.push_back(&oRD);
    OptStruct oIR = {"l:", 0, NULL, NULL, "inverse reciprocity flag"}; options.push_back(&oIR);
    
    //! Other checking
    OptStruct od = {"d:", 0, NULL, NULL, "mindist flag and value"}; options.push_back(&od);
    OptStruct oT = {"t:", 0, NULL, NULL, "mindist dilatation flag"}; options.push_back(&oT);
    OptStruct oD = {"D:", 0, NULL, NULL, "REGRESSION MINDISS flag"}; options.push_back(&oD);
    
    OptStruct oS = {"s:", 0, NULL, NULL, "self similarity flag and value"}; options.push_back(&oS);
    OptStruct oB = {"b:", 0, NULL, NULL, "integral of derivatives"}; options.push_back(&oB);
    OptStruct oV = {"v:", 0, NULL, NULL, "variance flag and value"}; options.push_back(&oV);
    
    OptStruct oO = {"o:", 0, "0", NULL, "remove isolated flag"}; options.push_back(&oO);
    OptStruct oOg = {"O:", 0, "0", NULL, "remove isolated grain (# pixels)"}; options.push_back(&oOg);
    OptStruct oC = {"C:", 0, "-1", NULL, "filter using the cost, train removing a fraction of the accepted points (ej. 0.05)"}; options.push_back(&oC);
	
	OptStruct oA = {"a:", 0, NULL, NULL, "use laplacian of the image instead of the image itself"}; options.push_back(&oA);
    
   //! postprocessing option : combine
	OptStruct oc = {"c:", 0, NULL, NULL, "combine last scale with the previous one to densify the result"}; options.push_back(&oc);
    
	vector<ParStruct2 *> parameters;
	ParStruct2 pinput = {"image1", NULL, "intput left image",0}; parameters.push_back(&pinput);
	ParStruct2 pinput2 = {"image2", NULL, "input right image",0}; parameters.push_back(&pinput2);
	ParStruct2 poutput = {"output", NULL, "output disparity",0}; parameters.push_back(&poutput);
	ParStruct2 pmask = {"mask", NULL, "outptut mask",0}; parameters.push_back(&pmask);
	ParStruct2 poutput2 = {"outputr", NULL, "right output disparity",1}; parameters.push_back(&poutput2);
	ParStruct2 pmask2 = {"maskr", NULL, "right outptut mask",1}; parameters.push_back(&pmask2);
    
	
	
	if (!parsecmdline2((char*) "iip_stereo_multi_scale_pixel_computation",(char*) "Pixel precision correlation", argc, argv, options, parameters))
		return 0;
	
	
	
	//! Parameters: Input Images
	libIIPStable::cflimage input;
	input.load(pinput.value);
	
	libIIPStable::cflimage input2;
	input2.load(pinput2.value);
	

   //! Set to 0 all nans and inf in the input images 
   for (int i=0;i<input.whc();i++) input[i] = isfinite(input[i]) ? input[i] : 0;
   for (int i=0;i<input2.whc();i++) input2[i] = isfinite(input2[i]) ? input2[i] : 0;
	
	
    //! Parameters: structure
	libIIPStable::strParameters strPar;
	
    
    
    //! Correlation window
    if (!oK.flag)
    {
        
        
        if (!oX.flag || !oY.flag) { printf("error :: please specify window size\n"); return -1;}
        strPar.flagWinWeighted= atoi(oW.value);

        // LISTS!! 
        strPar.flagListKernels = atoi(oWL.value);

        // BY DEFAULT USE MULTI WINDOW
        strPar.flagMultiWin = 1; 
        
        libIIPStable::flimage prolate(atoi(oX.value), atoi(oY.value));
        prolate=1.0f;
        if (strPar.flagWinWeighted)
        {
            stereo_build_gaussian_kernel(prolate, 0.05f, 0.05f);
        }
        prolate.normalizeL1();
        strPar.prolate = prolate;
        if (STEREOVERBOSE) {printf("info :: saving used prolate\n "); prolate.save("stereo_prolate.pmf");}
        
        
    } else
    {
        libIIPStable::flimage prolate;
        prolate.load(oK.value);
        prolate.normalizeL1();
        strPar.flagWinWeighted = 1;
        strPar.flagListKernels = 0; // no lists if the kernel is provided by the user
        strPar.flagMultiWin = 0; // no multi window if the kernel is provided
        strPar.prolate = prolate;
    }
    
    
    
    //! type of distance
    strPar.itypeDist= atoi(oI.value);

    //! number of precisions for single scale
    strPar.inPrecisions = atoi(oP.value);
	
    //! factor of disparity refinement by cubic interpolation
    strPar.cubicRefinementFactor = atoi(oPrefine.value);
	
    //! number of scales
    strPar.nScales = atoi(oN.value);
	
    
    //! reciprocity flag
    strPar.flagRecip = oR.flag;
    if (oR.flag)
    {
        if (atof(oR.value) < 0.0f) strPar.flagRecip = 0;
        else strPar.valueRecip = atof(oR.value);
    }
    
    
    //! reciprocity flag
    strPar.flagRecipDual = oRD.flag;
    if (oRD.flag && atof(oRD.value) < 0.0f) strPar.flagRecipDual = 0;
    
    
    //! inverse reciprocity flag
    strPar.flagiRecip = oIR.flag;
    if (oIR.flag && atof(oIR.value) < 0.0f) strPar.flagiRecip = 0;
    
    
    //! flag subrecip flag
    strPar.flagSubRecip = oG.flag;
    if (oG.flag)
    {
        if (atof(oG.value) < 0.0) strPar.flagSubRecip = 0;
        else strPar.valueSubRecip = atof(oG.value);
    }
    
    
    
    //! Remove isolated
    strPar.flagRemoveIsolated = oO.flag;
    if (strPar.flagRemoveIsolated)
    {
        if (atof(oO.value) <= 0.0f) strPar.flagRemoveIsolated = 0;
        else    strPar.valueRemoveIsolated = atof(oO.value);
    }
    
    //! Remove isolated grain filter
    strPar.flagGrainFilter = oOg.flag;
    if (strPar.flagGrainFilter)
    {
        if (atof(oOg.value) <= 0.0f) strPar.flagGrainFilter = 0;
        else    strPar.valueGrainFilter = atof(oOg.value);
    }
    
    //! Filter based on the cost 
    strPar.flagCostFilter = oC.flag;
    if (strPar.flagCostFilter)
    {
        if (atof(oC.value) < 0.0f) strPar.flagCostFilter = 0;
        else    strPar.valueCostFilter = atof(oC.value);
    }
    
    
    //! mindist flag
    strPar.flagMinDist = od.flag;
    if (od.flag)
    {
        if (atof(od.value) < 0.0f) strPar.flagMinDist=0;
        else strPar.valueMinDist = atof(od.value);
	 }

    //! regression mindist flag
    strPar.flagRegressionMinDist = oD.flag;
    if (strPar.flagRegressionMinDist)
    {
        if (atof(oD.value) <= 0.0f) strPar.flagRegressionMinDist=0;
        else strPar.flagRegressionMinDist=1;
	 }
    
    
    strPar.flagMinDistDilatation = oT.flag;
    if (oT.flag && atof(oT.value) < 0.0f) strPar.flagMinDistDilatation = 0;
    

    
    //! self similarity flag
    strPar.flagSelfSimilarity = oS.flag;
    if (oS.flag)
    {
        if (atof(oS.value) < 0.0f) strPar.flagSelfSimilarity = 0;
        else strPar.valueSelfSimilarity = atof(oS.value);
    }
    
    
    //! variance flag
    strPar.flagVariance = oV.flag;
    if (oV.flag)
    {
        if (atof(oV.value) < 0.0f) strPar.flagVariance = 0;
        else strPar.valueVariance = atof(oV.value);
    }
    
    
    //! integral derivatives flag
    strPar.flagHorInfo = oB.flag;
    if (oB.flag)
    {
        if (atof(oB.value) < 0.0f) strPar.flagHorInfo = 0;
        else strPar.valueHorInfo = atof(oB.value);
    }
    
    
    
    //! noise standard deviation
    if (oF.flag)
        strPar.fSigma = atof(oF.value);
    else if (oV.flag || oS.flag || oT.flag)
    {
        printf("error :: please specify noise std");
        return -1;
    }
    
    
    
    //! Compute laplacian of the image
    if (oA.flag && atof(oA.value) >= 0.0f)
    {
    
        //! build kernel
        int kwidth = 7;
        int skwidth = kwidth / 2;
        libIIPStable::flimage kernel(kwidth,kwidth); kernel=1.0f;
        stereo_build_gaussian_kernel(kernel, 0.05f, 0.05f); kernel.normalizeL1();
        
        kernel[skwidth * kwidth + skwidth] -= 1.0;
        
        
        //! normalize to variance equal 1
        float fSum = 0.0f;
        for (int ii=0; ii < kernel.wh(); ii++) fSum += kernel[ii] * kernel[ii];
        for (int ii=0; ii < kernel.wh(); ii++) kernel[ii] /= sqrtf(fSum);
      

        
        input = input.convolve(kernel);
        input2 = input2.convolve(kernel);
        
        
    }

    //! Flag to allow combining restul of last nth scale with the previous one 
    //  This densifies the result
    strPar.flagCombineLastScale = oc.flag;
    if (strPar.flagCombineLastScale)
    {
        if (atof(oc.value) <= 0.0f) strPar.flagCombineLastScale=0;
        else strPar.flagCombineLastScale=atoi(oc.value);
	 }
    
    
	// Parameters: dmin, dmax images
	// Range of search for each pixel of the first image in the epipolar line of the second image (x+dmin,x+dmax)
	// Range is controlled to be inside the image in the correlation process
   // GF MODIFICATION : reads images if needed
	libIIPStable::flimage Dmin(input.w(),input.h());
	libIIPStable::flimage Dmax(input.w(),input.h());
	{
		
		if (!om.flag)   Dmin = - (float) input2.w();
		else // Dmin = atof(om.value);
      {
         // if it is a number read it
         char *end;
         float val = strtof(om.value, &end);
         if(end!=om.value) Dmin = val;
         else { 
            // otherwise it could be a file, then read it!
            Dmin.load(om.value);
            if (Dmin.w() != input.w() || Dmin.h() != input.h() ) {
               printf("error :: the Dmin/Dmax images must have the same size as the input");
               return -1;
            }
         }
      }
		
		if (!oM.flag)   Dmax = (float) input2.w();
		else // Dmax = atof(oM.value);
      {
         // if it is a number read it
         char *end;
         float val = strtof(oM.value, &end);
         if(end!=oM.value) Dmax = val;
         else { 
            // otherwise it could be a file, then read it!
            Dmax.load(oM.value);
            if (Dmax.w() != input.w() || Dmax.h() != input.h() ) {
               printf("error :: the Dmin/Dmax images must have the same size as the input");
               return -1;
            }
         }
      }
        
      strPar.dmin0 = Dmin.min();
		strPar.dmax0 = Dmax.max();
		
      for(int i=0;i<input.wh();i++){
         if( isnan(Dmin[i] ) ) Dmin[i] = Dmin.min();
         if( isnan(Dmax[i] ) ) Dmax[i] = Dmax.max();
      }

    }
	

    //! load the segmentation
    if (oSeg.flag )
    {
        strPar.segments.load(oSeg.value); 
	 }else {
        strPar.segments= NULL;
    }
    
    
    
	// Parameters: idmin, idmax images
	// Range of search for each pixel of the second image in the epipolar line of the first image
	// Computed as the inverse range of the first image
	// Range is controlled to be inside the image in the correlation function
	// Used only if reciprocity flag is used and correlation is performed in both directions
	libIIPStable::flimage iDmin(input2.w(),input2.h());
	libIIPStable::flimage iDmax(input2.w(),input2.h());
   {

		if (!oM.flag)
			iDmin = - (float) input.w();
		else
			iDmin = - Dmax.max(); //atof(oM.value);
		
		
		if (!om.flag)
			iDmax = (float) input.w();
		else
			iDmax= - Dmin.min(); //atof(om.value);
		
		
		strPar.idmin0 = iDmin.min();
		strPar.idmax0 = iDmax.max();
		
   }
    
    
    
    // Memory for disparity and mask of selected pixel taking left image as reference
	int iiScale = 1;
    libIIPStable::flimage odisp(input.w(), input.h());     odisp = 0.0f;
    libIIPStable::flimage odist(input.w(), input.h());     odist = fLarge;
	libIIPStable::flimage omask(input.w(), input.h());     omask = 0.0f;
	
    libIIPStable::flimage odisp2(input2.w(), input2.h());     odisp = 0.0f;
    libIIPStable::flimage odistr(input2.w(), input2.h());     odistr = fLarge;
	libIIPStable::flimage omask2(input2.w(), input2.h());     omask = 0.0f;

   // SET MULTI WINDOW
   strPar.flagMultiWin = 1;

	
    stereo_pixel_multiscale_chain(input, input2, iiScale, Dmin, Dmax, iDmin, iDmax, odisp, odisp2, odist, odistr, omask, omask2, strPar);
    
    
    
    
	odisp.save(poutput.value);
	omask.save(pmask.value);
   if (pmask2.value) {
	   odisp2.save(poutput2.value);
	   omask2.save(pmask2.value);
   }
    
	
}



