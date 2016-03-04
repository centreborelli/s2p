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
 *        Template is for both normal and SSE version.
 **/
template <typename T>
T initForward(
  const T* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z) {

  //! Initialization
  double zk  = p_z;
  double iz  = 1.0 / p_z;
  double z2k = pow(p_z, double(p_size - 1));
  T sum      = i_vec[0] + float(z2k) * i_vec[p_step * (p_size - 1)];
  z2k = z2k * z2k * iz;

  //! Loop over the pixels
  for (size_t k = 1; k < p_size - 1; k++) {
    sum += float(zk + z2k) * i_vec[p_step * k];
    zk  *= p_z;
    z2k *= iz;
  }

  //! Get the result
  return sum / float(1.0 - zk * zk);
}


/**
 * @brief Init the backward recursion for spline application.
 *        Template is for both normal and SSE version.
 **/
template <typename T>
inline T initBackward(
  const T* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z) {

  return float(p_z / (p_z * p_z - 1.0)) * (float(p_z) *
    i_vec[p_step * (p_size - 2)] + i_vec[p_step * (p_size - 1)]);
}


/**
 * @brief Apply the 1D spline interpolation.
 *        Template is for both normal and SSE version.
 **/
template <typename T>
void applySpline(
  T* io_vec,
  const size_t p_step,
  const size_t p_size,
  const double* z,
  const size_t p_nPoles) {

  //! Normalization
  const float lambda = 1.430575  * (1.0 + 1.0 / 0.430575 ) *
                       1.0430963 * (1.0 + 1.0 / 0.0430963);

  //! Initialization
  for (size_t k = 0; k < p_size; k++) {
    io_vec[k * p_step] *= lambda;
  }

  //! Loop on poles
  for (size_t n = 0; n < p_nPoles; n++) {
    const float zn = z[n];

    //! Forward recursion
    io_vec[0] = initForward(io_vec, p_step, p_size, zn);
    T sum = io_vec[0];
    for (size_t k = 1; k < p_size; k++) {
      sum = zn * sum + io_vec[k * p_step];
      io_vec[k * p_step] = sum;
    }

    //! Backward recursion
    io_vec[(p_size - 1) * p_step] = initBackward(io_vec, p_step, p_size, zn);
    sum = io_vec[(p_size - 1) * p_step];
    for (int k = int(p_size) - 2; k >= 0; k--) {
      sum = zn * (sum - io_vec[k * p_step]);
      io_vec[k * p_step] = sum;
    }
  }
}


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
