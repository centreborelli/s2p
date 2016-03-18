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
#include <omp.h>
#include <fftw3.h> // to remove


//! Local includes
#include "Utilities/Parameters.h"
#include "Utilities/Time.h"
#include "LibImages/LibImages.h"
#include "LibMSMW/LibMSMW.h"
#include "Utilities/Memory.h" // to remove

using namespace std;

// ./MSMW -il /home/ibal/Images/input/stereo/a.png -ir /home/ibal/Images/input/stereo/b.png -dl /home/ibal/Images/output/S2P/pixell.tif -kl maskl.tif -dr pixelr.tif -kr maskr.tif -n 4 -W 5 -x 9 -y 9 -m -8.759999999999999432e+01 -M 6.359999999999999432e+01

int main(int argc, char **argv) {

  //! Retrieve the arguments
  Parameters params;
  if (params.checkArgs(argc, argv) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Time
  Time time;

  //! Read left and right input images
  Image imL, imR;
  imL.read(params.inpLeft());
  imR.read(params.inpRight());
  time.getTime("Read images");

  //! Set to 0 all NaNs and Inf in the input images
  imL.removeNaN();
  imR.removeNaN();
  time.getTime("Remove NaN");

  //! Output images
  Image imDispL, imDispR, imDistL, imDistR, imMaskL, imMaskR;

  //! Initialize MSMW
  MSMW msmw;
  msmw.init(imL, imR, params);
  time.getTime("Initialize msmw object");

  //! Launch the stereo multiscale chain
  msmw.run(imDispL, imDispR, imDistL, imDistR, imMaskL, imMaskR);
  time.getTime("StereoMultiscaleChain");

  //! Save images
  imDispL.write(params.outDispL());
  if (params.outDispR() != NULL) {
    imDispR.write(params.outDispR());
  }
  imMaskL.write(params.outMaskL());
  if (params.outMaskR() != NULL) {
    imMaskR.write(params.outMaskR());
  }
  time.getTime("Write images");

  return EXIT_SUCCESS;
}



