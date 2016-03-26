#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED


//! Global includes
#include <cstring>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <xmmintrin.h>
#include <x86intrin.h>
#include <vector>


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief Return the optimal cut of the lines of an image according to
 *        the number of threads available.
 *
 * @param i_height: the height of the image;
 * @param o_heights: will contain the lines position (must be already allocated)
 * @param p_nbThreads: the number of threads available.
 **/
void initializeHeights(
  const size_t i_height,
  size_t* o_heights,
  const size_t p_nbThreads);


/**
 * @brief Test if x is not NaN.
 **/
inline bool isNumber(
  const float i_x) {
    return (i_x == i_x);
}


/**
 * @brief Read an homography: [a b c;
                               d e f;
                               g h i]
 **/
void readHomography(
  const char* i_fileName,
  double o_mat[9]);



#endif // UTILITIES_H_INCLUDED
