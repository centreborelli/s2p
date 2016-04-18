#ifndef UTILITIESMSMW_H_INCLUDED
#define UTILITIESMSMW_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"
#include "../Utilities/Parameters.h"


/**
 * @brief Compute a cubic interpolation. WARNING: the inline is really important
 *        for speed purpose.
 *
 * @param i_array: ptr to the first of four values stored in an array;
 * @param p_x: coefficient for the cubic interpolation.
 **/
inline float cubicInterpolate(
  const float* i_array,
  const float p_x) {
  return i_array[1] + 0.5f * p_x * (i_array[2] - i_array[0] + p_x * (2.f *
    i_array[0] - 5.f * i_array[1] + 4.f * i_array[2] - i_array[3] + p_x * (3.f *
    (i_array[1] - i_array[2]) + i_array[3] - i_array[0])));
}


/**
 * @brief Compute four cubic interpolations at once (SSE version). WARNING: the
 *        inline is really important for speed purpose.
 *
 * @param i_array: ptr to the first of 4x4 values stored in an array;
 * @param p_x: coefficient for the cubic interpolation.
 **/
inline __m128 cubicInterpolate128(
  const float* i_array,
  const float p_x) {

  const __m128 a0 = _mm_loadu_ps(i_array + 0);
  const __m128 a1 = _mm_loadu_ps(i_array + 1);
  const __m128 a2 = _mm_loadu_ps(i_array + 2);
  const __m128 a3 = _mm_loadu_ps(i_array + 3);
  const __m128 xx = _mm_set1_ps(p_x);

  return a1 + _mm_set1_ps(0.5) * xx * (a2 - a0 + xx * (_mm_set1_ps(2.f) * a0 -
    _mm_set1_ps(5.f) * a1 + _mm_set1_ps(4.f) * a2 - a3 + xx * (_mm_set1_ps(3.f)
    * (a1 - a2) + a3 - a0)));
}


/**
 * @brief Compute eight cubic interpolations at once (AVX version). WARNING: the
 *        inline is really important for speed purpose.
 *
 * @param i_array: ptr to the first of 8x4 values stored in an array;
 * @param p_x: coefficient for the cubic interpolation.
 **/
inline __m256 cubicInterpolate256(
  const float* i_array,
  const float p_x) {

  const __m256 a0 = _mm256_loadu_ps(i_array + 0);
  const __m256 a1 = _mm256_loadu_ps(i_array + 1);
  const __m256 a2 = _mm256_loadu_ps(i_array + 2);
  const __m256 a3 = _mm256_loadu_ps(i_array + 3);
  const __m256 xx = _mm256_set1_ps(p_x);

  return a1 + _mm256_set1_ps(0.5) * xx * (a2 - a0 + xx * (_mm256_set1_ps(2.f) *
    a0 - _mm256_set1_ps(5.f) * a1 + _mm256_set1_ps(4.f) * a2 - a3 + xx * (
    _mm256_set1_ps(3.f) * (a1 - a2) + a3 - a0)));
}


/**
 * @brief Compute the distance between the current patch in i_im1 and four
 *        patches in i_im2, where four images are stored columns aside.
 *
 * @param i_im1: reference image;
 * @param i_im2: contains four images stored columns aside;
 * @param p_px, p_py: pixel coordinates for the reference image;
 * @param p_qx, p_qy: pixel coordinates for the four second image;
 * @param i_kernel: weighting window for the distance computation;
 * @param p_normalize: if true substract the mean of the patches before distance
 *                     computation.
 **/
__m128 computePatchDistance4(
  const Image &i_im1,
  const Image &i_im2,
  const int p_px,
  const int p_py,
  const int p_qx,
  const int p_qy,
  const Image& i_kernel,
  const bool p_normalize);


/**
 * @brief Compute the distance between two current patches in i_im1 and 2 x four
 *        patches in i_im2, where four images are stored columns aside.
 *
 * @param i_im1: reference image;
 * @param i_im2: contains four images stored columns aside;
 * @param p_px, p_py: pixel coordinates for the reference image;
 * @param p_qx, p_qy: pixel coordinates for the four second image;
 * @param i_kernel: weighting window for the distance computation;
 * @param p_normalize: if true substract the mean of the patches before distance
 *                     computation.
 **/
__m256 computePatchDistance8(
  const Image &i_im1,
  const Image &i_im2,
  const int p_px,
  const int p_py,
  const int p_qx,
  const int p_qy,
  const Image& i_kernel,
  const bool p_normalize);


/**
 * @brief Select the smallest value in i_dist.
 *
 * @param i_dist: contains four distance values;
 * @param i_pos: contains the four corresponding positions;
 * @param io_dist: will contain the smallest value of i_dist if the new distance
 *                 is smaller than the old one;
 * @param o_pos: update the position according to the distance.
 **/
void getBest4(
  const __m128 i_dist,
  const __m128 i_pos,
  float& io_dist,
  float& o_pos);


/**
 * @brief Select the smallest value in i_dist.
 *
 * @param i_dist: contains eight distance values;
 * @param i_pos: contains the eight corresponding positions;
 * @param io_dist: will contain the smallest value of i_dist if the new distance
 *                 is smaller than the old one;
 * @param o_pos: update the position according to the distance.
 **/
void getBest8(
  const __m256 i_dist,
  const __m256 i_pos,
  float& io_dist,
  float& o_pos);


/**
 * @brief Initialize a set of windows with different size and orientations.
 *
 * @param o_windows: will contain a set of window. Must be allocated.
 * @param p_params: see Parameters;
 * @param io_nb: will contain the number of window.
 **/
void fillWindows(
  Image* o_windows,
  const Parameters& p_params,
  size_t& io_nb);


/**
 * @brief Compute 4 translations of an image: 0, 0.25, 0.5, 0.75.
 *        Those 4 images will be stored into one, column by column.
 **/
void getTranslated(
  const Image& i_im,
  Image& o_im);


/**
 * @brief Fill a window with a diagonal.
 *
 * @param i_size: size of the window;
 * @param p_d: size of the diagonal;
 * @param p_inverse: fill a left or right diagonal accordingly;
 * @param o_window: will be filled and normalized.
 **/
void fillDiagonal(
  const size_t i_size,
  const int p_d,
  const bool p_inverse,
  Image &o_window);


#endif // UTILITIESMSMW_H_INCLUDED
