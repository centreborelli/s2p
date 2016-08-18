#ifndef LIBSSE_H_INCLUDED
#define LIBSSE_H_INCLUDED


//! Global includes
#ifdef __SSE__
#include <xmmintrin.h>
#endif // __SSE__
#ifdef __AVX__
#include <immintrin.h>
#endif // __AVX__

/**
 * @brief Arctangent of y / x.
 **/
#ifdef __SSE__
__m128 atan2_128(
  const __m128& y,
  const __m128& x);
#endif // __SSE__
#ifdef __AVX__
__m256 atan2_256(
  const __m256& y,
  const __m256& x);
#endif // __AVX__


/**
 * @brief Compute the x modulus y.
 **/
 #ifdef __SSE__
__m128 modulus_128(
  const __m128& x,
  const __m128& y);
#endif // __SSE__
#ifdef __AVX__
__m256 modulus_256(
  const __m256& x,
  const __m256& y);
#endif // __AVX__


/**
 * @brief Compute the hypotenuse of a right-angled triange.
 **/
 #ifdef __SSE__
inline __m128 hypot_128(
  const __m128& x,
  const __m128& y) {
  return _mm_sqrt_ps(x * x + y * y);
}
#endif // __SSE__
#ifdef __AVX__
inline __m256 hypot_256(
  const __m256& x,
  const __m256& y) {
  return _mm256_sqrt_ps(x * x + y * y);
}
#endif // __AVX__


/**
 * @brief Compute the exponential of a value.
 **/
 #ifdef __SSE__
__m128 exp_128(
  const __m128& x);
#endif // __SSE__
#ifdef __AVX__
__m256 exp_256(
  const __m256& x);
#endif // __AVX__


/**
 * @brief Compute the absolute value of a value.
 **/
 #ifdef __SSE__
inline __m128 abs_128(
  const __m128& x) {
  return _mm_max_ps(x, -x);
}
#endif // __SSE__
#ifdef __AVX__
inline __m256 abs_256(
  const __m256& x) {
  return _mm256_max_ps(x, -x);
}
#endif // __AVX__


/**
 * @brief Compute the orientation bin.
 **/
#ifdef __SSE__
__m128 ori_to_bin_128(
  const __m128& ori,
  const int nbins);
#endif // __SSE__
#ifdef __AVX__
__m256 ori_to_bin_256(
  const __m256& ori,
  const int nbins);
#endif // __AVX__


/**
 * @brief If else then SSE implementation.
 **/
#ifdef __SSE__
inline __m128 applyMask_ps(
  const __m128 i_mask,
  const __m128 i_then,
  const __m128 i_else){

  // return mask ? then : else
  return _mm_or_ps(_mm_and_ps(i_mask, i_then), _mm_andnot_ps(i_mask, i_else));
}
#endif // __SSE__


/**
 * @brief If else then AVX implementation.
 **/
#ifdef __AVX__
inline __m256 applyMask256_ps(
  const __m256 i_mask,
  const __m256 i_then,
  const __m256 i_else){

  // return mask ? then : else
  return _mm256_or_ps(_mm256_and_ps(i_mask, i_then),
                      _mm256_andnot_ps(i_mask, i_else));
}
#endif // __AVX__


#endif // LIBSSE_H_INCLUDED
