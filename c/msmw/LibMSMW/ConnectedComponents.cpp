/**
 * @file ConnectedComponents.cpp
 *
 * @brief Small functions to deal with connected components.
 *
 * @author Original: Toni Buades, Gabriele Facciolo, Enric Meinhardt-Llopis
 * @author Modified: Marc Lebrun <marc.lebrun.ik@gmail.com>.
 **/


//! Global includes
#include <vector>


//! Local includes
#include "ConnectedComponents.h"


//! Connected components of positive pixels of the image rep (8 connected)
void updatePositiveConnectedComponent(
  int *io_rep,
  const int p_width,
  const int p_height) {

  for (int i = 0; i < p_height - 1; i++) {
    int* iR0 = io_rep +  i      * p_width;
    int* iR1 = io_rep + (i + 1) * p_width;

    for (int j = 0; j < p_width; j++) {
      if (iR0[j] >= 0) {
        if (iR0[j + 1] >= 0) {
          joinCC(io_rep, i * p_width + j, i * p_width + j + 1);
        }

        if (iR1[j] >= 0) {
          joinCC(io_rep, i * p_width + j, (i + 1) * p_width + j);
        }

        if (iR1[j + 1] >= 0) {
          joinCC(io_rep, i * p_width + j, (i + 1) * p_width + j + 1);
        }

        if (j > 0 && iR1[j - 1] >= 0) {
          joinCC(io_rep, i * p_width + j, (i + 1) * p_width + j - 1);
        }
      }
    }
  }

  for (int k = 0; k < p_width * p_height; k++) {
    if (io_rep[k] >= 0) {
      io_rep[k] = findCC(io_rep, k);
    }
  }
}


//! Find index such as i_array[index] = index and update i_array[p_value].
int findCC(
  int *io_array,
  const int p_value) {

  //! Recursive call
  if (p_value != io_array[p_value]) {
    io_array[p_value] = findCC(io_array, io_array[p_value]);
  }

  //! Return the value
  return io_array[p_value];
}


//! Join connected components.
void joinCC(
  int *io_array,
  const int p_a,
  const int p_b) {

  //! Find positions
  const int a = findCC(io_array, p_a);
  const int b = findCC(io_array, p_b);

  //! Make link if needed
  if (a != b) {
    if (a < b) { // arbitrary choice
      io_array[b] = a;
    }
    else {
      io_array[a] = b;
    }
  }
}


//! Initialize index of an image.
void initRep(
  const Image& i_im,
  int* o_rep) {

  //! For convenience
  const size_t w = i_im.width();

  for (size_t i = 0, k = 0; i < i_im.height(); i++) {
    const float* iI = i_im.getPtr(0, i);
    int* oR = o_rep + i * w;

    for (size_t j = 0; j < w; j++, k++) {
      oR[j] = iI[j] > 0 ? (int) k : -1;
    }
  }
}


//! Remove small regions of an image based on connected components.
void removeRegions(
  Image& io_im,
  int* io_rep,
  const size_t p_minArea) {

  //! For convenience
  const size_t h = io_im.height();
  const size_t w = io_im.width ();

  //! Count connected components
  int* ccIdx = (int*) malloc(w * h * sizeof(int));
  int nbCC = 0;
  for (int k = 0; k < int(w * h); k++) {
    if (io_rep[k] == k) {
      ccIdx[nbCC++] = k;
    }
  }

  //! Store the cc_pixels
  std::vector< std::vector< struct pix_new > > cc_pixels(w * h);

  for (size_t i = 0; i < h; i++) {
    const float* iI = io_im.getPtr(0, i);
    const int* iR = io_rep + i * w;

    for (size_t j = 0; j < w; j++) {
      const int cc = iR[j];

      //! If it's a region
      if (cc >= 0) {
        struct pix_new pixel;
        pixel.i = i;
        pixel.j = j;
        pixel.v = iI[j];
        cc_pixels[cc].push_back(pixel);
      }
    }
  }

  //! For each connected component
  for (int n = 0; n < nbCC; n++) {
    const int cc = ccIdx[n];
    std::vector<pix_new>* pPix = &cc_pixels[cc];

    //! Remove small regions
    if (pPix->size() < p_minArea) {
      for (size_t k = 0; k < pPix->size(); k++) {
        io_im.getPtr(0, (*pPix)[k].i)[(*pPix)[k].j] = 0.f;
      }
    }
  }

  //! Release memory
  free(ccIdx);
}
