/**
 * @file main.cpp
 *
 * @brief Main function for Landscape Evolution Modelling.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

//! Global includes
#include <ctime>


//! Local includes
#include "LibImages/LibImages.h"
#include "Utilities/Time.h"
#include "Utilities/Utilities.h"
#include "Utilities/Parameters.h"
#include "LibHomography/Homography.h"


using namespace std;


/**
 * @brief Main function.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/
int main(
  int argc,
  char** argv) {

  //! Retrieve the arguments
  Parameters params;
  if (params.checkArgs(argc, argv) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Initialization of time
  Time time;

  //! Read the input image
  Image imI;
  imI.read(params.inpName());
  if (params.verbose()) time.getTime("Read image");

  //! Read the homography
  double mat[9];
  readHomography(params.inpHomo(), mat);
  if (params.verbose()) time.getTime("Read homography");

  //! Call the mapping function
  Image imO;
  runHomography(imI, mat, imO, params);
  if (params.verbose()) time.getTime("Apply the homography");

  //! Write the image
  imO.write(params.outName());
  if (params.verbose()) time.getTime("Write image");

  //! Exit the main function
  return EXIT_SUCCESS;
}
