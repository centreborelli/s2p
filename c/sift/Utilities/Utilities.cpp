/**
 * @file Utilites.cpp
 *
 * @brief Miscellaneous functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <png.h>
#include <tiffio.h>


//! Local includes
#include "Utilities.h"
#include "Memory.h"


using namespace std;


//! Return the optimal cut of the lines of an image according to the number of
//! threads available.
void initializeHeights(
  const size_t i_height,
  size_t* o_heights,
  const size_t p_nbThreads) {

  //! Get the initial number of lines for each thread to deal with.
  const size_t nb = p_nbThreads;
  const size_t h = ceil(float(i_height) / float(nb));

  //! Initialize the number of lines
  size_t* lines = (size_t*) malloc(nb * sizeof(size_t));
  for (size_t n = 0; n < nb; n++) {
    lines[n] = h;
  }

  //! Remove one by one the lines that are not needed
  int nCur = nb - 1;
  int k = int(h * nb) - int(i_height);
  while (k > 0) {
    lines[nCur]--;
    if (nCur > 0) {
      nCur--;
    }
    else {
      nCur = nb - 1;
    }
    k--;
  }

  //! Get the lines position
  o_heights[0] = 0;
  for (size_t n = 0; n < nb; n++) {
    o_heights[n + 1] = o_heights[n] + lines[n];
  }

  //! Release memory
  free(lines);
}


//! Compute the x modulus y.
float mod(
  const float i_x,
  const float i_y) {

  float z = i_x;
  int n;
  if(z < 0){
      n = (int)((-z)/i_y)+1;
      z += n*i_y;
  }
  n = (int)(z/i_y);
  z -= n*i_y;
  return z;
}


//! Get the bins from the orientation.
int ori2bin(
  const float i_ori,
  const int i_nbBins) {

  return (int) ((i_ori < 0 ? i_ori + 2 * M_PI : i_ori) / (2 * M_PI) * i_nbBins + 0.5) % i_nbBins;
}


//! Get the orientation from the bins.
float bin2ori(
  const float i_bin,
  const int i_nbBins) {

  const float ori = (i_bin + 0.5) * 2 * M_PI / float(i_nbBins);
  return ori > M_PI ? ori - 2 * M_PI : ori;
}


//! Iterative box filter of w 3 bins.
void smoothCircularHistogram(
  const int p_nbBins,
  const int p_nbIter,
  float* io_hist) {

  //! Initialization
  float* tmp = (float*) memalloc(16, (p_nbBins + 2) * sizeof(float));

  //! Convolution with box filters
  for (int i = 0; i < p_nbIter; i++) {

    //! Copy the current histogram
    for (int n = 0; n < p_nbBins; n++) {
      tmp[n + 1] = io_hist[n];
    }

    //! Circular border conditions
    tmp[0] = io_hist[p_nbBins - 1];
    tmp[p_nbBins + 1] = io_hist[0];

    //! Apply the convolution
    for (int n = 0; n < p_nbBins; n++) {
      io_hist[n] = (tmp[n] + tmp[n + 1] + tmp[n + 2]) / 3.f;
    }
  }

  //! Release memory
  memfree(tmp);
}










