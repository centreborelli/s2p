#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED


//! Global includes


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
 * @brief Compute the x modulus y.
 **/
float mod(
  const float i_x,
  const float i_y);


/**
 * @brief Get the bins from the orientation.
 **/
int ori2bin(
  const float i_ori,
  const int i_nbBins);


/**
 * @brief Get the orientation from the bins.
 **/
float bin2ori(
  const float i_bin,
  const int i_nbBins);


/**
 * @brief Iterative box filter of w 3 bins.
 *
 * @param p_nbBins: number of bins;
 * @param p_nbIter: number of iterations;
 * @param io_hist: histogram to smooth.
 **/
void smoothCircularHistogram(
  const int p_nbBins,
  const int p_nbIter,
  float* io_hist);


#endif // UTILITIES_H_INCLUDED
