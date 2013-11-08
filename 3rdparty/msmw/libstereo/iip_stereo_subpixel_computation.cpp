#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../libraryStable/libraryBasic.h"

#include "libstereo.h"

using namespace std;



int main(int argc, char **argv)
{
	

	vector <OptStruct *> options;
	OptStruct oX = {"x:", 0, "9", NULL, "window size"}; options.push_back(&oX);
    OptStruct oW = {"w:", 0, "0", NULL, "flag for weighting window"}; options.push_back(&oW);
    OptStruct oI = {"i:", 0, "0", NULL, "initial scale"}; options.push_back(&oI);
    OptStruct oF = {"f:", 0, "6", NULL, "final scale"}; options.push_back(&oF);
    OptStruct oD = {"d:",  0, "0", NULL, "flag for same depth window"}; options.push_back(&oD);
    OptStruct oR = {"r:", 0, NULL, NULL, "adaptive window by Rafa"}; options.push_back(&oR);
    OptStruct oS = {"s:", 0, "1", NULL, "noise standard deviation"}; options.push_back(&oS);

	
    /*
	OptStruct oP = {"p:", 0, NULL, NULL, "image of precision"}; options.push_back(&oP);
	OptStruct oC = {"c:", 0, NULL, NULL, "consistency flag"}; options.push_back(&oC);
	OptStruct oM = {"m:", 0, NULL, NULL, "right left mask used with consistency flag"}; options.push_back(&oM);
	OptStruct oI = {"i:", 0, NULL, NULL, "right left disparity used with consistency flag"}; options.push_back(&oI);
	*/
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image1", NULL, "intput left image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image2", NULL, "input right image"}; parameters.push_back(&pinput2);
	
	ParStruct pidisp = {"idisparity", NULL, "input disparity"}; parameters.push_back(&pidisp);
	ParStruct pimask = {"imask", NULL, "input mask"}; parameters.push_back(&pimask);

	ParStruct podisp = {"odisp", NULL, "output disparity"}; parameters.push_back(&podisp);
	//ParStruct podist = {"odist", NULL, "output disparity"}; parameters.push_back(&podist);
	ParStruct pmask = {"omask", NULL, "outptut mask"}; parameters.push_back(&pmask);
	//ParStruct pparameters = {"parameters", NULL, "input parameters file"}; parameters.push_back(&pparameters);
	
	
	
	if (!parsecmdline("iip_stereo_subpixel_computation","SubPixel precision correlation", argc, argv, options, parameters))
		return 0;

	

	//! Parameters: Input Images
	libIIPStable::cflimage input;	
	input.load(pinput.value);

	libIIPStable::cflimage input2;
	input2.load(pinput2.value);
	
	libIIPStable::flimage idisparity;
	idisparity.load(pidisp.value);

	
	
	//! Parameters: Initial mask
	libIIPStable::flimage imask;	
	imask.load(pimask.value);

	
	//! Process left-right direction
	libIIPStable::flimage odisp(input.w(), input.h());     odisp = 0.0f;
	libIIPStable::flimage odist(input.w(), input.h());     odist = 0.0f;
	libIIPStable::flimage omask(input.w(), input.h());     omask = 0.0f;
	
	

	//! Parameters: structure
	libIIPStable::strSubPixelPar strPar;
    
    strPar.sflagWinWeighted = atoi(oW.value);
    libIIPStable::flimage prolate(atoi(oX.value), atoi(oX.value));
    prolate=1.0f;
    if (strPar.sflagWinWeighted)
    {
        stereo_build_gaussian_kernel(prolate, 0.05f, 0.05f);
    }
    prolate.normalizeL1();
    strPar.prolate = prolate;
    
    


	strPar.inScales = atoi(oI.value);
    strPar.fnScales = atoi(oF.value);
	strPar.sflagSameDepth = atoi(oD.value);
    strPar.sfSigma = atof(oS.value);
    
    
    if (oR.flag)
	{
		strPar.sflagRafa = 1;
		strPar.sValueRafa = atof(oR.value);
        
        if (!oS.flag)
        {
            printf("... -r option needs of -s option to be activated");
            return -1;
        }
	}
    
    
    
    /*


	strPar.sflagPrecision = oP.flag;
	if (oP.flag) strPar.sPrecision.create(input.w(), input.h());
	
	
	strPar.sflagRecip = oC.flag;
	if (oC.flag)
	{
		
		
		if (!oM.flag) {printf("error :: iip_stereo_subpixel_computation :: -m option needed if -c activated\n"); exit(-1);}
		strPar.imaskr.load(oM.value);
		
		
		if (!oI.flag) {printf("error :: iip_stereo_subpixel_computation :: -i option needed if -c activated\n"); exit(-1);}
		strPar.idispr.load(oI.value);
		
		strPar.sReciprocity.create(input.w(), input.h());
	}
	*/
	
	
	
	libIIPStable::stereo_subpixel_computation(input, input2, idisparity, imask, odisp, omask, odist, strPar);
	

	//! Process: save outputs	
	odisp.save(podisp.value);
	//odist.save(podist.value);
	omask.save(pmask.value);

    
	/*
	if (oP.flag)
	{
		strPar.sPrecision.save(oP.value);
	}
	
	
	if (oC.flag)
	{
		strPar.sReciprocity.save(oC.value);
	}
	*/
	
	
}

