/**
 * @file Splines.cpp
 *
 * @brief Function to prepare and apply splines of size 5 to an image.
 *
 * @author Pascal Monasse (original version),
 * @author Gabriele Facciolo (modified version),
 * @author Marc Lebrun (final version) <marc.lebrun.ik@gmail.com>
 **/


//! Global includes


//! Local includes
#include "Splines.h"
#include "../LibImages/LibImages.h"
#include "../Utilities/Memory.h"
#include "../Utilities/Utilities.h"


using namespace std;


//! Prepare image for cardinal spline interpolation.
void prepareSpline(
  Image& io_im) {

  //! For convenience
  const size_t chnls  = io_im.channels();
  const size_t width  = io_im.width();
  const size_t height = io_im.height();

  //! Replace NaN with 0
  for (size_t c = 0; c < chnls; c++) {
    for (size_t i = 0; i < height; i++) {
      float* iI = io_im.getPtr(c, i);

      for (size_t j = 0; j < width; j++) {
        if (!isNumber(iI[j])) {
          iI[j] = 0.f;
        }
      }
    }
  }

  //! Initialize poles of associated z-filter
  const double z[2] = {-0.430575, -0.0430963};
  const size_t nPoles = 2;

  //! For each channels
  for (size_t c = 0; c < chnls; c++) {

    //! Allocation for SSE version
    __m128* vecI = (__m128*) memalloc(16, width * sizeof(__m128));
    size_t i = 0;

    //! Apply the interpolation over the lines - SSE version
    for (; i < height - 4; i += 4) {
      float* iI0 = io_im.getPtr(c, i + 0);
      float* iI1 = io_im.getPtr(c, i + 1);
      float* iI2 = io_im.getPtr(c, i + 2);
      float* iI3 = io_im.getPtr(c, i + 3);

      //! Copy data
      for (size_t j = 0; j < width; j++) {
        vecI[j] = _mm_set_ps(iI0[j], iI1[j], iI2[j], iI3[j]);
      }

      //! Apply the spline
      applySpline(vecI, 1, width, z, nPoles);

      //! Copy data
      for (size_t j = 0; j < width; j++) {
        float value[4];
        _mm_storeu_ps(value, vecI[j]);
        iI0[j] = value[3];
        iI1[j] = value[2];
        iI2[j] = value[1];
        iI3[j] = value[0];
      }
    }

    //! Release memory
    memfree(vecI);

    //! Apply the interpolation over the lines - normal version
    for (; i < height; i++) {
      applySpline(io_im.getPtr(c, i), 1, width, z, nPoles);
    }

    //! Allocation for SSE version
    __m128* vecJ = (__m128*) memalloc(16, height * sizeof(__m128));
    size_t j = 0;

    //! Apply the interpolation of the columns - SSE version
    for (; j < width - 4; j += 4) {

      //! Copy the data
      for (size_t i = 0; i < height; i++) {
        vecJ[i] = _mm_loadu_ps(io_im.getPtr(c, i) + j);
      }

      //! Apply the spline
      applySpline(vecJ, 1, height, z, nPoles);

      //! Copy the data
      for (size_t i = 0; i < height; i++) {
        _mm_storeu_ps(io_im.getPtr(c, i) + j, vecJ[i]);
      }
    }

    //! Release memory
    memfree(vecJ);

    //! Apply the interpolation of the columns - normal version
    for (; j < width; j++) {
      applySpline(io_im.getPtr(c, 0) + j, width, height, z, nPoles);
    }
  }
}


//! Spline interpolation.
void interpolateSpline(
  const Image& i_im,
  Image& o_im,
  const float p_x,
  const float p_y,
  const size_t p_i,
  const size_t p_j) {

  //! For convenience
  const int w = i_im.width();
  const int h = i_im.height();
  const size_t chnls = i_im.channels();
  const float NaN = sqrtf(-1.f);

  //! Check if the pixel is inside
  if (p_x < 0.f || p_x > float(w) || p_y < 0.f || p_y > float(h)) {
    for (size_t c = 0; c < chnls; c++) {

      //! Return NaN
      o_im.getPtr(c, p_i)[p_j] = NaN;
    }
    return;
  }

  //! Position initialization
  const float x  = p_x - 0.5f;
  const float y  = p_y - 0.5f;
  const int xi   = x < 0 ? -1 : int(x);
  const int yi   = y < 0 ? -1 : int(y);
  const float ux = x - float(xi);
  const float uy = y - float(yi);

  //! Precomputation of the spline of order 5
  float cx[6], cy[6];
  initSpline(cx, ux);
  initSpline(cy, uy);
  const __m128 xx = _mm_set_ps(cx[2], cx[3], cx[4], cx[5]);

  //! This test saves computational time
  if (xi >= 2 && xi < w - 3 && yi >= 2 && yi < h - 3) {

    for (size_t c = 0; c < chnls; c++) {

      const float* iI0 = i_im.getPtr(c, yi - 2);
      const float* iI1 = i_im.getPtr(c, yi - 1);
      const float* iI2 = i_im.getPtr(c, yi    );
      const float* iI3 = i_im.getPtr(c, yi + 1);
      const float* iI4 = i_im.getPtr(c, yi + 2);
      const float* iI5 = i_im.getPtr(c, yi + 3);

      const __m128 xVal = xx * (_mm_set1_ps(cy[5]) * _mm_loadu_ps(iI0 + xi - 2) +
                                _mm_set1_ps(cy[4]) * _mm_loadu_ps(iI1 + xi - 2) +
                                _mm_set1_ps(cy[3]) * _mm_loadu_ps(iI2 + xi - 2) +
                                _mm_set1_ps(cy[2]) * _mm_loadu_ps(iI3 + xi - 2) +
                                _mm_set1_ps(cy[1]) * _mm_loadu_ps(iI4 + xi - 2) +
                                _mm_set1_ps(cy[0]) * _mm_loadu_ps(iI5 + xi - 2));

      const float value = cy[5] * (iI0[xi + 2] * cx[1] + iI0[xi + 3] * cx[0]) +
                          cy[4] * (iI1[xi + 2] * cx[1] + iI1[xi + 3] * cx[0]) +
                          cy[3] * (iI2[xi + 2] * cx[1] + iI2[xi + 3] * cx[0]) +
                          cy[2] * (iI3[xi + 2] * cx[1] + iI3[xi + 3] * cx[0]) +
                          cy[1] * (iI4[xi + 2] * cx[1] + iI4[xi + 3] * cx[0]) +
                          cy[0] * (iI5[xi + 2] * cx[1] + iI5[xi + 3] * cx[0]);

      float tmp[4];
      _mm_storeu_ps(tmp, xVal);
      o_im.getPtr(c, p_i)[p_j] = value + tmp[0] + tmp[1] + tmp[2] + tmp[3];
    }
  }
  else {
    for (size_t c = 0; c < chnls; c++) {
      float value = 0.f;

      for (int di = -2; di <= 3; di++) {
        const float* iI = i_im.getPtr(c, yi + di <  0 ? -yi - di - 1 :
                                         yi + di >= h ? 2 * h - yi - di - 1 : yi + di);
        value += cy[3 - di] * (iI[symi(xi - 2, w)] * cx[5] +
                               iI[symi(xi - 1, w)] * cx[4] +
                               iI[symi(xi    , w)] * cx[3] +
                               iI[symi(xi + 1, w)] * cx[2] +
                               iI[symi(xi + 2, w)] * cx[1] +
                               iI[symi(xi + 3, w)] * cx[0]);
      }
      o_im.getPtr(c, p_i)[p_j] = value;
    }
  }
}


//! Initialize the coefficients of the spline of order 5.
void initSpline(
  float io_w[6],
  const float p_val) {

  const float ak[6] = {1. / 120., -0.05, 0.125, -1. / 6., 0.125, -0.05};

  io_w[0] = pow5(p_val) * ak[0];
  io_w[1] = pow5(p_val) * ak[1] + pow5(p_val + 1) * ak[0];
  io_w[2] = pow5(p_val) * ak[2] + pow5(p_val + 1) * ak[1] + pow5(p_val + 2) * ak[0];
  io_w[3] = pow5(p_val) * ak[3] + pow5(p_val + 1) * ak[2] + pow5(p_val + 2) * ak[1] + pow5(p_val + 3) * ak[0];
  io_w[4] = pow5(p_val) * ak[4] + pow5(p_val + 1) * ak[3] + pow5(p_val + 2) * ak[2] + pow5(p_val + 3) * ak[1] + pow5(p_val + 4) * ak[0];
  io_w[5] = pow5(p_val) * ak[5] + pow5(p_val + 1) * ak[4] + pow5(p_val + 2) * ak[3] + pow5(p_val + 3) * ak[2] + pow5(p_val + 4) * ak[1] + pow5(p_val + 5) * ak[0];
}


/**
 * @brief Apply the 1D spline interpolation.
 *        normal version.
 **/
void applySpline(
  float* io_vec,
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
    float sum = io_vec[0];
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
 * @brief Apply the 1D spline interpolation.
 *        SSE version.
 **/
void applySpline(
  __m128* io_vec,
  const size_t p_step,
  const size_t p_size,
  const double* z,
  const size_t p_nPoles) {

  //! Normalization
  const __m128 lambda = _mm_set1_ps(1.430575  * (1.0 + 1.0 / 0.430575 ) *
                       1.0430963 * (1.0 + 1.0 / 0.0430963));

  //! Initialization
  for (size_t k = 0; k < p_size; k++) {
    io_vec[k * p_step] *= lambda;
  }

  //! Loop on poles
  for (size_t n = 0; n < p_nPoles; n++) {
    const __m128 zn = _mm_set1_ps((float) z[n]);

    //! Forward recursion
    io_vec[0] = initForward(io_vec, p_step, p_size, z[n]);
    __m128 sum = io_vec[0];
    for (size_t k = 1; k < p_size; k++) {
      sum = zn * sum + io_vec[k * p_step];
      io_vec[k * p_step] = sum;
    }

    //! Backward recursion
    io_vec[(p_size - 1) * p_step] = initBackward(io_vec, p_step, p_size, z[n]);
    sum = io_vec[(p_size - 1) * p_step];
    for (int k = int(p_size) - 2; k >= 0; k--) {
      sum = zn * (sum - io_vec[k * p_step]);
      io_vec[k * p_step] = sum;
    }
  }
}


/**
 * @brief Init the forward recursion for spline application.
 *        normal version.
 **/
float initForward(
  const float* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z) {

  //! Initialization
  double zk  = p_z;
  double iz  = 1.0 / p_z;
  double z2k = pow(p_z, double(p_size - 1));
  float sum      = i_vec[0] + float(z2k) * i_vec[p_step * (p_size - 1)];
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
 * @brief Init the forward recursion for spline application.
 *        and SSE version.
 **/
__m128 initForward(
  const __m128* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z) {

  //! Initialization
  double zk  = p_z;
  double iz  = 1.0 / p_z;
  double z2k = pow(p_z, double(p_size - 1));
  __m128 sum      = i_vec[0] + _mm_set1_ps((float) z2k) * i_vec[p_step * (p_size - 1)];
  z2k = z2k * z2k * iz;

  //! Loop over the pixels
  for (size_t k = 1; k < p_size - 1; k++) {
    sum += _mm_set1_ps(float(zk + z2k)) * i_vec[p_step * k];
    zk  *= p_z;
    z2k *= iz;
  }

  //! Get the result
  return sum / _mm_set1_ps(float(1.0 - zk * zk));
}


/**
 * @brief Init the backward recursion for spline application.
 *        normal version.
 **/
inline float initBackward(
  const float* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z) {

  return float(p_z / (p_z * p_z - 1.0)) * (float(p_z) *
    i_vec[p_step * (p_size - 2)] + i_vec[p_step * (p_size - 1)]);
}
/**
 * @brief Init the backward recursion for spline application.
 *        SSE version.
 **/
inline __m128 initBackward(
  const __m128* i_vec,
  const size_t p_step,
  const size_t p_size,
  const double p_z) {

  return _mm_set1_ps(float(p_z / (p_z * p_z - 1.0))) * (_mm_set1_ps(float(p_z)) *
    i_vec[p_step * (p_size - 2)] + i_vec[p_step * (p_size - 1)]);
}
