/**
 * @file UtilitiesMSMW.cpp
 *
 * @brief Functions for distance computation and best value selection in SSE/AVX
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <stdlib.h>


//! Local includes
#include "UtilitiesMSMW.h"


//! Compute the distance between the current patch in i_im1 and four patches in
//! i_im2, where four images are stored columns aside.
__m128 computePatchDistance4(
  const Image &i_im1,
  const Image &i_im2,
  const int ipx,
  const int ipy,
  const int iqx,
  const int iqy,
  const Image& i_kernel,
  const bool p_normalize) {

  //! Initialization
  __m128 dist = _mm_setzero_ps();
  const int kh = i_kernel.height();
  const int kw = i_kernel.width ();

  //! Compute the distance over all channels
  for (size_t c = 0; c < i_im1.channels(); c++) {

    //! Compute the mean of the current patch
    __m128 m1 = _mm_setzero_ps(), m2 = _mm_setzero_ps();
    if (p_normalize) {
      for (int i = 0; i < kh; i++) {
        const float* iI1 = i_im1.getPtr(c, ipy + i) + ipx;
        const float* iI2 = i_im2.getPtr(c, iqy + i) + iqx * 4;
        const float* iK  = i_kernel.getPtr(0, i);

        for (int j = 0; j < kw; j++) {
          m1 += _mm_set1_ps(iI1[j] * iK[j]);
          m2 += _mm_loadu_ps(iI2 + 4 * j) * _mm_set1_ps(iK[j]);
        }
      }
    }

    //! Loop over the patch
    for (int i = 0; i < kh; i++) {
      const float* iI1 = i_im1.getPtr(c, ipy + i) + ipx;
      const float* iI2 = i_im2.getPtr(c, iqy + i) + iqx * 4;
      const float* iK  = i_kernel.getPtr(0, i);

      for (int j = 0; j < kw; j++) {
        // Here using _mm_broadcast_ss should be faster than using
        // _mm256_set1_ps because the pointer is already available.
        const __m128 xVal = _mm_broadcast_ss(iI1 + j) - m1
                          - _mm_loadu_ps(iI2 + 4 * j) + m2;
        dist += _mm_set1_ps(iK[j]) * xVal * xVal;
      }
    }
  }

  //! Return the normalized result
  return dist / _mm_set1_ps(float(i_im1.channels()));
}


//! Compute the distance between two current patches in i_im1 and 2 x four
//! patches in i_im2, where four images are stored columns aside.
__m256 computePatchDistance8(
  const Image &i_im1,
  const Image &i_im2,
  const int ipx,
  const int ipy,
  const int iqx,
  const int iqy,
  const Image& i_kernel,
  const bool p_normalize) {

  //! Initialization
  __m256 dist = _mm256_setzero_ps();
  const int kh = i_kernel.height();
  const int kw = i_kernel.width ();

  //! Compute the distance over all channels
  for (size_t c = 0; c < i_im1.channels(); c++) {

    //! Compute the mean of the current patch
    __m256 m1 = _mm256_setzero_ps(), m2 = _mm256_setzero_ps();
    if (p_normalize) {
      for (int i = 0; i < kh; i++) {
        const float* iI1 = i_im1.getPtr(c, ipy + i) + ipx;
        const float* iI2 = i_im2.getPtr(c, iqy + i) + iqx * 4;
        const float* iK  = i_kernel.getPtr(0, i);

        for (int j = 0; j < kw; j++) {
          m1 += _mm256_broadcast_ss(iI1 + j) * _mm256_set1_ps(iK[j]);
          m2 += _mm256_loadu_ps(iI2 + 4 * j) * _mm256_set1_ps(iK[j]);
        }
      }
    }

    //! Compute the distance over the patch
    for (int i = 0; i < kh; i++) {
      const float* iI1 = i_im1.getPtr(c, ipy + i) + ipx;
      const float* iI2 = i_im2.getPtr(c, iqy + i) + iqx * 4;
      const float* iK  = i_kernel.getPtr(0, i);

      for (int j = 0; j < kw; j++) {
        // Here using _mm256_broadcast_ss should be faster than using
        // _mm256_set1_ps because the pointer is already available.
        const __m256 xVal = _mm256_broadcast_ss(iI1 + j) - m1
                          - _mm256_loadu_ps(iI2 + 4 * j) + m2;
        dist += _mm256_set1_ps(iK[j]) * xVal * xVal;
      }
    }
  }

  //! Return the normalized result
  return dist / _mm256_set1_ps(float(i_im1.channels()));
}


//! Select the smallest value in i_dist.
void getBest4(
  const __m128 i_dist,
  const __m128 i_pos,
  float& io_dist,
  float& o_pos) {

  //! Get the values from the __m256 data
  float dist[4];
  _mm_storeu_ps(dist, i_dist);
  float pos [4];
  _mm_storeu_ps(pos , i_pos );

  //! Select the best one
  if (dist[0] < io_dist) {
    io_dist = dist[0];
    o_pos   = pos [0];
  }
  if (dist[1] < io_dist) {
    io_dist = dist[1];
    o_pos   = pos [1];
  }
  if (dist[2] < io_dist) {
    io_dist = dist[2];
    o_pos   = pos [2];
  }
  if (dist[3] < io_dist) {
    io_dist = dist[3];
    o_pos   = pos [3];
  }
}


//! Select the smallest value in i_dist.
void getBest8(
  const __m256 i_dist,
  const __m256 i_pos,
  float& io_dist,
  float& o_pos) {

  //! Get the values from the __m256 data
  float dist[8];
  _mm256_storeu_ps(dist, i_dist);
  float pos [8];
  _mm256_storeu_ps(pos , i_pos );

  //! Select the best one
  for (size_t n = 0; n < 8; n++) {
    if (dist[n] < io_dist) {
      io_dist = dist[n];
      o_pos   = pos [n];
    }
  }
}


//! Initialize a set of windows with different size and orientations.
void fillWindows(
  Image* o_windows,
  const Parameters& p_params,
  size_t& io_nb) {

  //! For convenience
  const int win = p_params.window()->width();

  //! At least one window
  io_nb = std::max((int) p_params.orientation(), 1);

  //! First window is size win x win, full of 1
  o_windows[0].init(win, win);
  o_windows[0] = 1.f;
  o_windows[0].normalizeL1();

  //! 3x3 window
  if (win == 3) {
    o_windows[1].init(9, 1);
    o_windows[1] = 1.f;
    o_windows[1].normalizeL1();

    o_windows[2].init(1, 9);
    o_windows[2] = 1.f;
    o_windows[2].normalizeL1();

    //! Update the number of orientation
    io_nb = std::min(io_nb, (size_t) 3);
  }

  //! 5x5
  else if (win == 5) {

    //! First two are uniforms
    o_windows[1].init(9, 3);
    o_windows[1] = 1.f;
    o_windows[1].normalizeL1();

    o_windows[2].init(3, 9);
    o_windows[2] = 1.f;
    o_windows[2].normalizeL1();

    //! Second two are diagonals
    fillDiagonal(9, 1, false, o_windows[3]);
    fillDiagonal(9, 1, true , o_windows[4]);

    //! 25.5 degrees windows
    int w5[81] = {0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  1 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  1 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  0 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,1 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 };
    int w6[81] = {0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 };
    int w7[81] = {0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 };
    int w8[81] = {0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,
                  0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,1 ,
                  0 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0 ,
                  1 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,
                  1 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
                  0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 };

    //! Initialize the last four windows
    o_windows[5].init(9, 9);
    o_windows[6].init(9, 9);
    o_windows[7].init(9, 9);
    o_windows[8].init(9, 9);

    //! Fill them
    for (size_t i = 0; i < 9; i++) {
      for (size_t j = 0; j < 9; j++) {
        o_windows[5].getPtr(0, i)[j] = w5[i * 9 + j];
        o_windows[6].getPtr(0, i)[j] = w6[i * 9 + j];
        o_windows[7].getPtr(0, i)[j] = w7[i * 9 + j];
        o_windows[8].getPtr(0, i)[j] = w8[i * 9 + j];
      }
    }

    //! Normalization
    o_windows[5].normalizeL1();
    o_windows[6].normalizeL1();
    o_windows[7].normalizeL1();
    o_windows[8].normalizeL1();

    //! Update the number of windows
    io_nb = std::min(io_nb, (size_t) 9);
  }

  //! 7x7
  else if (win == 7 || win == 9 || win == 13 || win == 17 || win == 21) {

    //! For convenience
    const size_t n1 = win ==  7 ? 11 :
                      win ==  9 ? 15 :
                      win == 13 ? 21 :
                      win == 17 ? 29 : 33;
    const size_t n2 = win ==  7 ?  5 :
                      win ==  9 ?  5 :
                      win == 13 ?  7 :
                      win == 17 ?  9 : 13;
    const size_t n3 = win ==  7 ?  2 :
                      win ==  9 ?  2 :
                      win == 13 ?  3 :
                      win == 17 ?  4 : 6;

    //! First two are uniforms
    o_windows[1].init(n1, n2);
    o_windows[1] = 1.f;
    o_windows[1].normalizeL1();

    o_windows[2].init(n2, n1);
    o_windows[2] = 1.f;
    o_windows[2].normalizeL1();

    //! Second two are diagonals
    fillDiagonal(n1, n3, false, o_windows[3]);
    o_windows[3].normalizeL1();
    fillDiagonal(n1, n3, true , o_windows[4]);
    o_windows[4].normalizeL1();

    //! Update the number of windows
    io_nb = std::min(io_nb, (size_t) 5);
  }

  //! This case shouldn't exist
  else {
    std::cout << "fillProlate: please specify a correct window size between {3,\
    5, 7, 9, 13, 17, 21}.\n Abort." << std::endl;
    exit(EXIT_FAILURE);
  }
}


//! Compute 4 translations of an image: 0, 0.25, 0.5, 0.75.
void getTranslated(
  const Image& i_im,
  Image& o_im) {

  //! For convenience
  const size_t chnls = i_im.channels();
  const size_t h     = i_im.height();
  const size_t w     = i_im.width();

  //! Compute translation for subpixel estimation
  Image im1(i_im), im2(i_im), im3(i_im);
  im1.translate(0.25);
  im2.translate(0.50);
  im3.translate(0.75);

  //! Initialize the size of the output image
  o_im.init(4 * w, h, chnls);

  //! Put all the images together
  for (size_t c = 0; c < chnls; c++) {
    for (size_t i = 0; i < h; i++) {
      const float* iI0 = i_im.getPtr(c, i);
      const float* iI1 = im1 .getPtr(c, i);
      const float* iI2 = im2 .getPtr(c, i);
      const float* iI3 = im3 .getPtr(c, i);
      float* oI = o_im.getPtr(c, i);

      for (size_t j = 0; j < w; j++) {
        oI[4 * j + 0] = iI0[j];
        oI[4 * j + 1] = iI1[j];
        oI[4 * j + 2] = iI2[j];
        oI[4 * j + 3] = iI3[j];
      }
    }
  }
}


//! Fill a window with a diagonal.
void fillDiagonal(
  const size_t i_size,
  const int p_d,
  const bool p_inverse,
  Image &o_window) {

  //! Initialize the kernel
  o_window.init(i_size, i_size);
  o_window = 0.f;

  //! Normal diagonal
  if (!p_inverse) {
    for (int i = 0; i < (int) i_size; i++) {
      float* oK = o_window.getPtr(0, i);

      for (int j = -p_d; j <= p_d; j++) {
        if (i + j >= 0 && i + j < (int) i_size) {
          oK[i + j] = 1.0f;
        }
      }
    }
  }

  //! Inverse diagonal
  else {
    for (int i = 0; i < (int) i_size; i++) {
      float* oK = o_window.getPtr(0, i_size - i - 1);

      for (int j = -p_d; j <= p_d; j++) {
        if (i + j >= 0 && i + j < (int) i_size) {
          oK[i + j] = 1.0f;
        }
      }
    }
  }

  //! Normalize the kernel
  o_window.normalizeL1();
}
