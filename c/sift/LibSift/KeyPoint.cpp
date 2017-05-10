/**
 * @file KeyPoint.cpp
 *
 * @brief Class to handle the key points.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <iostream>
#include <math.h>

#include <cassert>

//! Local includes
#include "KeyPoint.h"
#include "../Utilities/Memory.h"
#include "../Utilities/Utilities.h"
#include "../LibSSE/LibSSE.h"


using namespace std;


//! Default constructor
KeyPoint::KeyPoint() :

  m_x(0.f),
  m_y(0.f),
  m_sigma(0.f),
  m_theta(0.f),

  m_o(0),
  m_s(0),
  m_i(0),
  m_j(0),

  m_val(0.f),
  m_edgeResp(0.f),
  m_nbHist(0),
  m_nbOri(0),
  m_descriptors(NULL),

  m_nbBins(0),
  m_oriHist(NULL) {
}


//! Copy constructor.
KeyPoint::KeyPoint(
  const KeyPoint& i_keyPoint) :

  m_x(i_keyPoint.m_x),
  m_y(i_keyPoint.m_y),
  m_sigma(i_keyPoint.m_sigma),
  m_theta(i_keyPoint.m_theta),

  m_o(i_keyPoint.m_o),
  m_s(i_keyPoint.m_s),
  m_i(i_keyPoint.m_i),
  m_j(i_keyPoint.m_j),

  m_val(i_keyPoint.m_val),
  m_edgeResp(i_keyPoint.m_edgeResp),
  m_nbHist(i_keyPoint.m_nbHist),
  m_nbOri(i_keyPoint.m_nbOri),
  m_descriptors((float*) memalloc(16, m_nbOri * (m_nbHist + 1) * (m_nbHist + 1) * sizeof(float))),

  m_nbBins(i_keyPoint.m_nbBins),
  m_oriHist((float*) memalloc(16, m_nbBins * sizeof(float))) {

  //! Copy the data
  for (size_t n = 0; n < m_nbOri * m_nbHist * m_nbHist; n++) {
    m_descriptors[n] = i_keyPoint.m_descriptors[n];
  }

  for (size_t n = 0; n < m_nbBins; n++) {
    m_oriHist[n] = i_keyPoint.m_oriHist[n];
  }
}


//! Overload constructor.
KeyPoint::KeyPoint(
  const size_t i_nbOri,
  const size_t i_nbHist,
  const size_t i_nbBins) :

  m_x(0.f),
  m_y(0.f),
  m_sigma(0.f),
  m_theta(0.f),

  m_o(0),
  m_s(0),
  m_i(0),
  m_j(0),

  m_val(0.f),
  m_edgeResp(0.f),
  m_nbHist(i_nbHist),
  m_nbOri(i_nbOri),
  m_descriptors((float*) memalloc(16, m_nbOri * m_nbHist * m_nbHist * sizeof(float))),

  m_nbBins(i_nbBins),
  m_oriHist((float*) memalloc(16, m_nbBins * sizeof(float))) {

  //! Fill the data of zeros
  for (size_t n = 0; n < m_nbOri * m_nbHist * m_nbHist; n++) {
    m_descriptors[n] = 0.f;
  }

  for (size_t n = 0; n < m_nbBins; n++) {
    m_oriHist[n] = 0.f;
  }
}


//! Operator overload.
KeyPoint& KeyPoint::operator=(
  const KeyPoint& i_keyPoint) {

  if (&i_keyPoint == this) {
    return *this;
  }

  releaseMemory();

  new (this) KeyPoint(i_keyPoint);
  return *this;
}


//! Default destructor
KeyPoint::~KeyPoint() {
  releaseMemory();
}


//! Extract principal orientations from gradient orientation histogram, and
//! update the histogram.
int KeyPoint::extractPrincipalOrientations(
  const int p_nbBins,
  const float p_threshold,
  float* o_principalOrientations) {

  //! Number of principal orientations ( the return value).
  int nbOri = 0;

  //! Smooth histogram : 6 iterated box filters
  smoothCircularHistogram(p_nbBins, 6, m_oriHist);

  //! What is the value of the global maximum
  float maxValue = m_oriHist[0];
  for (int n = 1; n < p_nbBins; n++) {
    maxValue = std::max(maxValue, m_oriHist[n]);
  }

  //! Search for local extrema in the histogram
  for (int i = 0; i < p_nbBins; i++) {
    const float hc = m_oriHist[i];
    const float hp = m_oriHist[(i - 1 + p_nbBins) % p_nbBins];
    const float hn = m_oriHist[(i + 1           ) % p_nbBins];

    if (hc > p_threshold * maxValue && hc > hp && hc > hn) {

        //! Quadratic interpolation of the position of each local maximum
        const float di = (hp - hn) / (2 * (hp + hn - 2 * hc));

        //! Add to vector of principal orientations (expressed in [0, 2pi])
        o_principalOrientations[nbOri++] = bin2ori(float(i) + di, p_nbBins);
    }
  }

  //! Return the number of principal orientations
  return nbOri;
}


//! Extract feature vector and update the current keypoint.
void KeyPoint::extractFeatureVector(
  const ScaleSpace* i_sx,
  const ScaleSpace* i_sy,
  const Parameters& p_params) {

  //! For convenience
  const int h        = i_sx->getOctave(m_o)->getHeight();
  const int w        = i_sx->getOctave(m_o)->getWidth ();
  const float delta  = i_sx->getOctave(m_o)->getDelta();
  const float xk     = m_x / delta;
  const float yk     = m_y / delta;
  const float sk     = m_sigma / delta;
  const int    nHist = p_params.nbHist();
  const size_t nOri  = p_params.nbOri();
  const float lambda = p_params.lambdaDes();
  const float cosT   = cos(-m_theta);
  const float sinT   = sin(-m_theta);
  const float t2     = 1.f / (2.f * (lambda * sk) * (lambda * sk));

#ifdef __AVX__
  const __m256 zTheta  = _mm256_set1_ps(m_theta);
  const __m256 z2PI    = _mm256_set1_ps(2 * M_PI);
  const __m256 zt2     = _mm256_set1_ps(t2);
  const __m256 zNorm1  = _mm256_set1_ps(1.f / (2 * lambda * sk / nHist));
  const __m256 zNorm2  = _mm256_set1_ps(nOri / (2 * M_PI));
  const __m256 zHist   = _mm256_set1_ps((nHist - 1.0) / 2.0);
  const __m256 z0      = _mm256_setzero_ps();
  const __m256 z1      = _mm256_set1_ps(1.f);
#endif // __AVX__
#ifdef __SSE__
  const __m128 xTheta  = _mm_set1_ps(m_theta);
  const __m128 x2PI    = _mm_set1_ps(2 * M_PI);
  const __m128 xt2     = _mm_set1_ps(t2);
  const __m128 norm1   = _mm_set1_ps(1.f / (2 * lambda * sk / nHist));
  const __m128 norm2   = _mm_set1_ps(nOri / (2 * M_PI));
  const __m128 dHist   = _mm_set1_ps((nHist - 1.0) / 2.0);
  const __m128 x0      = _mm_setzero_ps();
  const __m128 x1      = _mm_set1_ps(1.f);
#endif // __SSE__

  //! Initialize descriptors tab
  float* iD = m_descriptors;
  for (size_t n = 0; n < nHist * nHist * nOri; n++) {
    iD[n] = 0.f;
  }

  //! Contributing pixels are inside a patch [siMin;siMax]X[sjMin;sjMax] of
  //! width 2*lambda_descr*sigma_key*(nhist+1)/nhist
  const float R  = (1 + 1 / float(nHist)) * lambda * sk;
  const float Rp = M_SQRT2 * R;
  const size_t siMin = std::max(0, int(xk - Rp + 0.5));
  const size_t sjMin = std::max(0, int(yk - Rp + 0.5));
  const size_t siMax = std::min(int(xk + Rp + 0.5), h - 1);
  const size_t sjMax = std::min(int(yk + Rp + 0.5), w - 1);

  //! To speed up
  const size_t N = (siMax - siMin + 1) * (sjMax - sjMin + 1);
  float* gx = (float*) memalloc(16, N * sizeof(float));
  float* gy = (float*) memalloc(16, N * sizeof(float));
  float* rx = (float*) memalloc(16, N * sizeof(float));
  float* ry = (float*) memalloc(16, N * sizeof(float));
  size_t nbPix = 0;

  /// For each pixel inside the patch.
  for (size_t si = siMin; si < siMax; si++) {
    const float* iX = i_sx->getOctave(m_o)->getPtr(m_s, si);
    const float* iY = i_sy->getOctave(m_o)->getPtr(m_s, si);

    for (size_t sj = sjMin; sj < sjMax; sj++) {

      //! Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
      const float X = cosT * (float(si) - xk) - sinT * (float(sj) - yk);
      const float Y = sinT * (float(si) - xk) + cosT * (float(sj) - yk);

      //! Does this sample fall inside the descriptor area ?
      if (fabs(X) < R && fabs(Y) < R) {
        gx[nbPix] = iX[sj];
        gy[nbPix] = iY[sj];
        rx[nbPix] = X;
        ry[nbPix] = Y;
        nbPix++;
      }
    }
  }

  float* arrI = (float*) memalloc(16, nbPix * sizeof(float));
  float* arrJ = (float*) memalloc(16, nbPix * sizeof(float));
  float* arrG = (float*) memalloc(16, nbPix * sizeof(float));
  float* arrM = (float*) memalloc(16, nbPix * sizeof(float));
  float* arr0 = (float*) memalloc(16, nbPix * sizeof(float));
  float* arr1 = (float*) memalloc(16, nbPix * sizeof(float));
  float* arr2 = (float*) memalloc(16, nbPix * sizeof(float));
  float* arr3 = (float*) memalloc(16, nbPix * sizeof(float));

  //! Loop over all the valid pixels
  size_t n = 0;

#ifdef __AVX__
  //! AVX version
  for (; n < nbPix - 8; n += 8) {

    //! Compute the gradient orientation (theta) on keypoint referential.
    const __m256 dx = _mm256_loadu_ps(gx + n);
    const __m256 dy = _mm256_loadu_ps(gy + n);
    const __m256 ori = modulus_256(atan2_256(dy, dx) - zTheta, z2PI);
    const __m256 X = _mm256_loadu_ps(rx + n);
    const __m256 Y = _mm256_loadu_ps(ry + n);

    //! Compute the gradient magnitude and apply a Gaussian weighing to give
    //! less emphasis to distant sample
    const __m256 M = hypot_256(dx, dy) * exp_256(-(X * X + Y * Y) * zt2);

    //! Bin indices, Compute the (tri)linear weightings ...
    const __m256 alpha = X * zNorm1 + zHist;
    const __m256 beta  = Y * zNorm1 + zHist;
    const __m256 gamma = ori * zNorm2;

    //! ...and add contributions to respective bins in different histograms.
    //! a loop with 1 or two elements
    const __m256 i0 = _mm256_floor_ps(alpha);
    const __m256 j0 = _mm256_floor_ps(beta );

    //! Masks
    const __m256 m0 = _mm256_cmp_ps(i0, z0, _CMP_LT_OS);
    const __m256 m1 = _mm256_cmp_ps(_mm256_set1_ps(nHist - 2), i0, _CMP_LT_OS);
    const __m256 m2 = _mm256_cmp_ps(j0, z0, _CMP_LT_OS);
    const __m256 m3 = _mm256_cmp_ps(_mm256_set1_ps(nHist - 2), j0, _CMP_LT_OS);

    //! Store the intermediate values
    _mm256_storeu_ps(arrI + n, i0);
    _mm256_storeu_ps(arrJ + n, j0);
    _mm256_storeu_ps(arrG + n, gamma);
    _mm256_storeu_ps(arrM + n, M);
    _mm256_storeu_ps(arr0 + n, applyMask256_ps(m0, z0, z1 - abs_256(i0 - alpha)));
    _mm256_storeu_ps(arr1 + n, applyMask256_ps(m1, z0, z1 - abs_256(i0 + z1 - alpha)));
    _mm256_storeu_ps(arr2 + n, applyMask256_ps(m2, z0, z1 - abs_256(j0 - beta)));
    _mm256_storeu_ps(arr3 + n, applyMask256_ps(m3, z0, z1 - abs_256(j0 + z1 - beta)));
  }
#endif // __AVX__

#ifdef __SSE__
  //! SSE version
  for (; n < nbPix - 4; n += 4) {

    //! Compute the gradient orientation (theta) on keypoint referential.
    const __m128 dx = _mm_loadu_ps(gx + n);
    const __m128 dy = _mm_loadu_ps(gy + n);
    const __m128 ori = modulus_128(atan2_128(dy, dx) - xTheta, x2PI);
    const __m128 X = _mm_loadu_ps(rx + n);
    const __m128 Y = _mm_loadu_ps(ry + n);

    //! Compute the gradient magnitude and apply a Gaussian weighing to give
    //! less emphasis to distant sample
    const __m128 M = hypot_128(dx, dy) * exp_128(-(X * X + Y * Y) * xt2);

    //! Bin indices, Compute the (tri)linear weightings ...
    const __m128 alpha = X * norm1 + dHist;
    const __m128 beta  = Y * norm1 + dHist;
    const __m128 gamma = ori * norm2;

    //! ...and add contributions to respective bins in different histograms.
    //! a loop with 1 or two elements
    const __m128 i0 = _mm_floor_ps(alpha);
    const __m128 j0 = _mm_floor_ps(beta );

    //! Masks
    const __m128 m0 = _mm_cmplt_ps(i0, x0);
    const __m128 m1 = _mm_cmplt_ps(_mm_set1_ps(nHist - 2), i0);
    const __m128 m2 = _mm_cmplt_ps(j0, x0);
    const __m128 m3 = _mm_cmplt_ps(_mm_set1_ps(nHist - 2), j0);

    //! Store the intermediate values
    _mm_storeu_ps(arrI + n, i0);
    _mm_storeu_ps(arrJ + n, j0);
    _mm_storeu_ps(arrG + n, gamma);
    _mm_storeu_ps(arrM + n, M);
    _mm_storeu_ps(arr0 + n, applyMask_ps(m0, x0, x1 - abs_128(i0 - alpha)));
    _mm_storeu_ps(arr1 + n, applyMask_ps(m1, x0, x1 - abs_128(i0 + x1 - alpha)));
    _mm_storeu_ps(arr2 + n, applyMask_ps(m2, x0, x1 - abs_128(j0 - beta)));
    _mm_storeu_ps(arr3 + n, applyMask_ps(m3, x0, x1 - abs_128(j0 + x1 - beta)));
  }
#endif // __SSE__

  //! Normal version
  for (; n < nbPix; n++) {

    //! Compute the gradient orientation (theta) on keypoint referential.
    const float dx = gx[n];
    const float dy = gy[n];
    const float ori = mod(atan2(dy, dx) - m_theta, 2 * M_PI);
    const float X = rx[n];
    const float Y = ry[n];

    //! Compute the gradient magnitude and apply a Gaussian weighing to give
    //! less emphasis to distant sample
    const float M = hypot(dx, dy) * exp(-(X * X + Y * Y) * t2);

    //! Bin indices, Compute the (tri)linear weightings ...
    const float alpha = X / (2 * lambda * sk / nHist) + (nHist - 1.0) / 2.0;
    const float beta  = Y / (2 * lambda * sk / nHist) + (nHist - 1.0) / 2.0;
    const float gamma = ori / (2 * M_PI) * nOri;

    //! ...and add contributions to respective bins in different histograms.
    //! a loop with 1 or two elements
    const int i0 = floor(alpha);
    const int j0 = floor(beta);

    //! Store the precomputed values
    arrI[n] = i0;
    arrJ[n] = j0;
    arr0[n] = i0 < 0 ? 0.f : 1.0 - fabs(float(i0) - alpha);
    arr1[n] = i0 > nHist - 2 ? 0.f : 1.0 - fabs(float(i0 + 1) - alpha);
    arr2[n] = j0 < 0 ? 0.f : 1.0 - fabs(float(j0) - beta);
    arr3[n] = j0 > nHist - 2 ? 0.f : 1.0 - fabs(float(j0 + 1) - beta);
    arrG[n] = gamma;
    arrM[n] = M;
  }

  for (size_t n = 0; n < nbPix; n++) {
    const float gamma = arrG[n];
    const float M = arrM[n];
    const float gl    = 1.0 - (gamma - floor(gamma));
    const float gr    = 1.0 - (floor(gamma) + 1 - gamma);
    const int i0 = arrI[n];
    const int j0 = arrJ[n];
    const int ia = std::max(0, i0);
    const int ib = std::min(ia + 1, nHist - 1);
    const int ja = std::max(0, j0);
    const int jb = std::min(ja + 1, nHist - 1);
    const size_t bl = (size_t) ( int(gamma)      % nOri);
    const size_t br = (size_t) ((int(gamma) + 1) % nOri);
    const float v0 = arr0[n];
    const float v1 = arr1[n];
    const float v2 = arr2[n];
    const float v3 = arr3[n];
    iD[(ia * nHist + ja) * nOri + bl] += gl * v0 * v2 * M;
    iD[(ia * nHist + ja) * nOri + br] += gr * v0 * v2 * M;

    iD[(ia * nHist + jb) * nOri + bl] += gl * v0 * v3 * M;
    iD[(ia * nHist + jb) * nOri + br] += gr * v0 * v3 * M;
    
    iD[(ib * nHist + ja) * nOri + bl] += gl * v1 * v2 * M;
    iD[(ib * nHist + ja) * nOri + br] += gr * v1 * v2 * M;

    iD[(ib * nHist + jb) * nOri + bl] += gl * v1 * v3 * M;
    iD[(ib * nHist + jb) * nOri + br] += gr * v1 * v3 * M;
  }

  //! Release memory
  memfree(gx);
  memfree(gy);
  memfree(rx);
  memfree(ry);
  memfree(arrI);
  memfree(arrJ);
  memfree(arrG);
  memfree(arrM);
  memfree(arr0);
  memfree(arr1);
  memfree(arr2);
  memfree(arr3);
}


//! Threshold and quantize the descriptors
void KeyPoint::thresholdAndQuantizeFeatureVector(
  const size_t p_n) {

  //! Threshold parameter
  const float threshold = 0.2f;

  //! Normalize
  float l2norm = 0;
  for (size_t i = 0; i < p_n; i++){
    l2norm += m_descriptors[i] * m_descriptors[i];
  }
  l2norm = sqrt(l2norm);

  //! Threshold bins
  const float thresh = l2norm * threshold;
  l2norm = 0.f;
  for (size_t i = 0; i < p_n; i++){
    const float d = std::min(m_descriptors[i], thresh);
    m_descriptors[i] = d;
    l2norm += d * d;
  }

  //! Quantization
  const float tQ = 512.f / sqrt(l2norm);
  for (size_t i = 0; i < p_n; i++){
    m_descriptors[i] = std::min((int)(m_descriptors[i] * tQ), 255);
  }
}


//! Release memory.
void KeyPoint::releaseMemory() {
  memfree(m_descriptors);
  memfree(m_oriHist);
}









