/**
 * @file main.cpp
 *
 * @brief Main function to compute the MSMW algorithm.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <stdio.h>
#include <stdlib.h>


//! Local includes
#include "Utilities/Parameters.h"
#include "Utilities/Time.h"
#include "LibImages/LibImages.h"
#include "LibMSMW/LibMSMW.h"

using namespace std;

// ./msmw -il a.png -ir b.png -dl pixell.tif -kl maskl.tif -dr pixelr.tif -kr maskr.tif -n 4 -W 5 -x 9 -y 9 -m -87 -M 63

int main(int argc, char **argv) {

    //! Retrieve the arguments
    Parameters params;
    if (params.checkArgs(argc, argv) != EXIT_SUCCESS)
        return EXIT_FAILURE;

    //! Time
    Time time;

    //! Read left and right input images
    Image imL, imR;
    imL.read(params.inpLeft());
    imR.read(params.inpRight());
    time.get_time("Read images");

    //! Set to 0 all NaNs and Inf in the input images
    imL.removeNaN();
    imR.removeNaN();
    time.get_time("Remove NaN");

    //! Output images
    Image imDispL, imDispR, imDistL, imDistR, imMaskL, imMaskR;

    //! Initialize MSMW
    MSMW msmw;
    msmw.init(imL, imR, params);
    time.get_time("Initialize msmw object");

    //! Launch the stereo multiscale chain
    msmw.run(imDispL, imDispR, imDistL, imDistR, imMaskL, imMaskR);
    time.get_time("StereoMultiscaleChain");

    //! Save images
    imDispL.write(params.outDispL());
    imMaskL.write(params.outMaskL());
    if (params.outDispR() != NULL)
        imDispR.write(params.outDispR());
    if (params.outMaskR() != NULL)
        imMaskR.write(params.outMaskR());
    time.get_time("Write images");

    return EXIT_SUCCESS;
}
