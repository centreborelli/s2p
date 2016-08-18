/**
 * @brief SSE functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

//! Global includes
#include <cmath>


//! Local includes
#include "LibSSE.h"
#include "../Utilities/Utilities.h"


//! Arctangent of y / x.
#ifdef __SSE__
__m128 atan2_128(
  const __m128& y,
  const __m128& x) {

  //! For convenience
  float a[4];
  float b[4];
  _mm_storeu_ps(a, x);
  _mm_storeu_ps(b, y);

  //! Compute the arc tangent
  a[0] = atan2(b[0], a[0]);
  a[1] = atan2(b[1], a[1]);
  a[2] = atan2(b[2], a[2]);
  a[3] = atan2(b[3], a[3]);

  //! Get the result
  return _mm_loadu_ps(a);
}
#endif // __SSE__
#ifdef __AVX__
__m256 atan2_256(
  const __m256& y,
  const __m256& x) {

  //! For convenience
  float a[8];
  float b[8];
  _mm256_storeu_ps(a, x);
  _mm256_storeu_ps(b, y);

  //! Compute the arc tangent
  a[0] = atan2(b[0], a[0]);
  a[1] = atan2(b[1], a[1]);
  a[2] = atan2(b[2], a[2]);
  a[3] = atan2(b[3], a[3]);
  a[4] = atan2(b[4], a[4]);
  a[5] = atan2(b[5], a[5]);
  a[6] = atan2(b[6], a[6]);
  a[7] = atan2(b[7], a[7]);

  //! Get the result
  return _mm256_loadu_ps(a);
}
#endif // __AVX__


//! Compute the x modulus y.
#ifdef __SSE__
__m128 modulus_128(
  const __m128& x,
  const __m128& y) {

  __m128 z = x;
  __m128 n = _mm_round_ps((-z) / y, _MM_FROUND_TO_ZERO) + _mm_set1_ps(1.f);
  __m128 mask = _mm_cmplt_ps(z, _mm_set1_ps(0.f));
  z = applyMask_ps(mask, z + n * y, z);

  n = _mm_round_ps(z / y, _MM_FROUND_TO_ZERO);
  return z - n * y;
}
#endif // __SSE__
#ifdef __AVX__
__m256 modulus_256(
  const __m256& x,
  const __m256& y) {

  __m256 z = x;
  __m256 n = _mm256_round_ps((-z) / y, _MM_FROUND_TO_ZERO) + _mm256_set1_ps(1.f);
  __m256 mask = _mm256_cmp_ps(z, _mm256_set1_ps(0.f), _CMP_LT_OS);
  z = applyMask256_ps(mask, z + n * y, z);

  n = _mm256_round_ps(z / y, _MM_FROUND_TO_ZERO);
  return z - n * y;
}
#endif // __AVX__


//! Compute the exponential of a value.
#ifdef __SSE__
__m128 exp_128(
  const __m128& x) {

  //! Clip the value
  __m128 y = _mm_max_ps(_mm_min_ps(x, _mm_set1_ps(88.3762626647949f)),
                                      _mm_set1_ps(-88.3762626647949f));

  //! Express exp(x) as exp(g + n * log(2))
  __m128 fx = y * _mm_set1_ps(1.44269504088896341) + _mm_set1_ps(0.5f);

  //! Floor
  const __m128 tmp = _mm_round_ps(fx, _MM_FROUND_TO_ZERO);

  //! If greater, substract 1
  const __m128 mask = _mm_and_ps(_mm_cmpgt_ps(tmp, fx), _mm_set1_ps(1.f));
  fx = tmp - mask;

  y -= fx * _mm_set1_ps(0.693359375 - 2.12194440e-4);
  const __m128 z = y * y;


  const __m128 t = (((((_mm_set1_ps(1.9875691500E-4)  * y +
                        _mm_set1_ps(1.3981999507E-3)) * y +
                        _mm_set1_ps(8.3334519073E-3)) * y +
                        _mm_set1_ps(4.1665795894E-2)) * y +
                        _mm_set1_ps(1.6666665459E-1)) * y +
                        _mm_set1_ps(5.0000001201E-1)) * z + y + _mm_set1_ps(1.f);

  //! Build 2^n
  const __m128i emm0 = _mm_add_epi32(_mm_cvttps_epi32(fx), _mm_set1_epi32(0x7f));

  //! Return the result
  return t * _mm_castsi128_ps(_mm_slli_epi32(emm0, 23));
}
#endif // __SSE__
#ifdef __AVX__
__m256 exp_256(
  const __m256& x) {

  //! Clip the value
  __m256 y = _mm256_max_ps(_mm256_min_ps(x, _mm256_set1_ps(88.3762626647949f)),
                                            _mm256_set1_ps(-88.3762626647949f));

  //! Express exp(x) as exp(g + n * log(2))
  __m256 fx = y * _mm256_set1_ps(1.44269504088896341) + _mm256_set1_ps(0.5f);

  //! Floor
  const __m256 tmp = _mm256_round_ps(fx, _MM_FROUND_TO_ZERO);

  //! If greater, substract 1
  const __m256 mask = _mm256_and_ps(_mm256_cmp_ps(tmp, fx, _CMP_GT_OS),
                                    _mm256_set1_ps(1.f));
  fx = tmp - mask;

  y -= fx * _mm256_set1_ps(0.693359375 - 2.12194440e-4);
  const __m256 z = y * y;


  const __m256 t = (((((_mm256_set1_ps(1.9875691500E-4)  * y +
                        _mm256_set1_ps(1.3981999507E-3)) * y +
                        _mm256_set1_ps(8.3334519073E-3)) * y +
                        _mm256_set1_ps(4.1665795894E-2)) * y +
                        _mm256_set1_ps(1.6666665459E-1)) * y +
                        _mm256_set1_ps(5.0000001201E-1)) * z + y +
                        _mm256_set1_ps(1.f);

  //! Build 2^n (split it into two SSE array, since AVX2 equivalent functions
  //! aren't available.
  const __m128i emm0 = _mm_add_epi32(_mm_cvttps_epi32(_mm256_castps256_ps128(fx)), _mm_set1_epi32(0x7f));
  const __m128i emm1 = _mm_add_epi32(_mm_cvttps_epi32(_mm256_extractf128_ps(fx, 1)), _mm_set1_epi32(0x7f));

  fx = _mm256_castps128_ps256(_mm_castsi128_ps(_mm_slli_epi32(emm0, 23)));
  fx = _mm256_insertf128_ps(fx, _mm_castsi128_ps(_mm_slli_epi32(emm1, 23)), 1);

  //! Return the result
  return t * fx;
}
#endif // __AVX__


//! Compute the orientation bin.
#ifdef __SSE__
__m128 ori_to_bin_128(
  const __m128& ori,
  const int nbins) {

  //! For convenience
  const __m128 x2PI = _mm_set1_ps(2 * M_PI);
  const __m128 xbins = _mm_set1_ps(nbins);

  //! Get it positive
  const __m128 mask = _mm_cmplt_ps(ori, _mm_setzero_ps());

  //! Get the value
  const __m128 val = _mm_round_ps(applyMask_ps(mask, ori + x2PI, ori)
    / x2PI * xbins + _mm_set1_ps(0.5f), _MM_FROUND_TO_ZERO);

  //! Return the modulo of it
  return val - xbins * _mm_round_ps(val / xbins, _MM_FROUND_TO_ZERO);
}
#endif // __SSE__
#ifdef __AVX__
__m256 ori_to_bin_256(
  const __m256& ori,
  const int nbins) {

  //! For convenience
  const __m256 x2PI  = _mm256_set1_ps(2 * M_PI);
  const __m256 xbins = _mm256_set1_ps(nbins);

  //! Get it positive
  const __m256 mask = _mm256_cmp_ps(ori, _mm256_setzero_ps(), _CMP_LT_OS);

  //! Get the value
  const __m256 val = _mm256_round_ps(applyMask256_ps(mask, ori + x2PI, ori)
    / x2PI * xbins + _mm256_set1_ps(0.5f), _MM_FROUND_TO_ZERO);

  //! Return the modulo of it
  return val - xbins * _mm256_round_ps(val / xbins, _MM_FROUND_TO_ZERO);
}
#endif // __SSE__














