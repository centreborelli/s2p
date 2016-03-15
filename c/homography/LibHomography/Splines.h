#ifndef SPLINES_H_INCLUDED
#define SPLINES_H_INCLUDED


//! Global includes
#include <cstdlib>
#include <cmath>
#include <xmmintrin.h>
#include <x86intrin.h>


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief Prepare the image for cardinal spline interpolation.
 */
void prepareSpline(
  Image& io_im);


/**
 * @brief Spline interpolation.
 **/
void interpolateSpline(
  const Image& i_im,
  Image& o_im,
  const float p_x,
  const float p_y,
  const size_t p_i,
  const size_t p_j);


/**
 * @brief Init the forward recursion for spline application.
 *        normal version.
 **/
float initForward(
  const float* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z);


/**
 * @brief Init the forward recursion for spline application.
 *        SSE version.
 **/
__m128 initForward(
  const __m128* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z);


/**
 * @brief Init the backward recursion for spline application.
 *        normal version.
 **/
inline float initBackward(
  const float* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z);


/**
 * @brief Init the backward recursion for spline application.
 *        SSE version.
 **/
inline __m128 initBackward(
  const __m128* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z);


/**
 * @brief Apply the 1D spline interpolation.
 *        normal version.
 **/
void applySpline(
  float* io_vec,
  const size_t p_step,
  const size_t p_size,
  const double* z,
  const size_t p_nPoles);


/**
 * @brief Apply the 1D spline interpolation.
 *        SSE version.
 **/
void applySpline(
  __m128* io_vec,
  const size_t p_step,
  const size_t p_size,
  const double* z,
  const size_t p_nPoles);


/**
 * @brief Initialize the coefficients of the spline of order 5.
 **/
void initSpline(
  float io_w[6],
  const float p_val);


/**
 * @brief For convenience return x^5.
 **/
inline float pow5(
  const float x) {

  const float x2 = x * x;
  return x2 * x2 * x;
}


/**
 * @brief Symetrization of the index of an image out of domain.
 **/
inline int symi(
  const int p_x,
  const int p_w) {

  //! Left / top
  if (p_x < 0) {
    return -p_x - 1;
  }

  //! Right / bottom
  if (p_x >= p_w) {
    return 2 * p_w - p_x - 1;
  }

  //! Center
  return p_x;
}


#endif // SPLINES_H_INCLUDED
