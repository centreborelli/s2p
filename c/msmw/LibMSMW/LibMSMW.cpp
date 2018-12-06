/**
 * @brief Launch the MSMW algorithm.
 *
 * @author Original: Toni Buades, Gabriele Facciolo, Enric Meinhardt-Llopis
 * @author Modified: Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
//#include <omp.h>


//! Local includes
#include "LibMSMW.h"
#include "UtilitiesMSMW.h"
#include "ConnectedComponents.h"
#include "../Utilities/Memory.h"
#include "../Utilities/Utilities.h"


//! Namespace
using namespace std;


//! Default constructor
MSMW::MSMW() :

  //! Images
  m_im1(NULL),
  m_im2(NULL),
  m_imDmin1(NULL),
  m_imDmax1(NULL),
  m_imDmin2(NULL),
  m_imDmax2(NULL),
  m_imDisp1(NULL),
  m_imDisp2(NULL),
  m_imDist1(NULL),
  m_imDist2(NULL),
  m_imMask1(NULL),
  m_imMask2(NULL),

  //! Convenience images: precomputed translation
  m_imTrans1(NULL),
  m_imTrans2(NULL),
  m_imStrobe11(NULL),
  m_imStrobe12(NULL),
  m_imStrobe21(NULL),
  m_imStrobe22(NULL),

  //! Windows
  m_windows(NULL),

  //! Size
  m_channels(0),
  m_height  (0),
  m_width   (0),

  //! Parameters
  m_params(NULL),
  m_currentScale(1),
  m_nbWindows(0) {

}


//! Copy constructor.
MSMW::MSMW(
  const MSMW& i_msmw) :

  //! Images
  m_im1(new Image(*i_msmw.m_im1)),
  m_im2(new Image(*i_msmw.m_im2)),
  m_imDmin1(new Image(*i_msmw.m_imDmin1)),
  m_imDmax1(new Image(*i_msmw.m_imDmax1)),
  m_imDmin2(new Image(*i_msmw.m_imDmin2)),
  m_imDmax2(new Image(*i_msmw.m_imDmax2)),
  m_imDisp1(new Image(*i_msmw.m_imDisp1)),
  m_imDisp2(new Image(*i_msmw.m_imDisp2)),
  m_imDist1(new Image(*i_msmw.m_imDist1)),
  m_imDist2(new Image(*i_msmw.m_imDist2)),
  m_imMask1(new Image(*i_msmw.m_imMask1)),
  m_imMask2(new Image(*i_msmw.m_imMask2)),

  //! Convenience images: precomputed translation
  m_imTrans1(new Image(*i_msmw.m_imTrans1)),
  m_imTrans2(new Image(*i_msmw.m_imTrans2)),
  m_imStrobe11(new Image(*i_msmw.m_imStrobe11)),
  m_imStrobe12(new Image(*i_msmw.m_imStrobe12)),
  m_imStrobe21(new Image(*i_msmw.m_imStrobe21)),
  m_imStrobe22(new Image(*i_msmw.m_imStrobe22)),

  //! Windows
  m_windows(new Image[i_msmw.m_nbWindows]),

  //! Size
  m_channels(i_msmw.m_channels),
  m_height  (i_msmw.m_height  ),
  m_width   (i_msmw.m_width   ),

  //! Parameters
  m_params(new Parameters(*i_msmw.m_params)),
  m_currentScale(i_msmw.m_currentScale),
  m_nbWindows(i_msmw.m_nbWindows) {

  //! Copy the windows
  for (size_t n = 0; n < m_nbWindows; n++) {
    m_windows[n] = i_msmw.m_windows[n];
  }
}


//! Default destructor
MSMW::~MSMW() {

  //! Release the memory
  releaseMemory();
}


//! Operator overload.
MSMW& MSMW::operator=(
  const MSMW& i_msmw) {

  if (&i_msmw == this) {
    return *this;
  }

  releaseMemory();

  new (this) MSMW(i_msmw);
  return *this;
}


//! Initialize the images and the parameters.
void MSMW::init(
  const Image& i_im1,
  const Image& i_im2,
  const Parameters& p_params) {

  //! Release memory
  this->releaseMemory();

  //! Parameters initialization
  m_params = new Parameters(p_params);
  m_nbWindows = m_params->orientation();

  //! Size initialization
  m_channels = i_im1.channels();
  m_height   = i_im1.height();
  m_width    = i_im1.width();

  //! Check that the size are correct
  if (i_im2.channels() != m_channels || i_im2.height() != m_height ||
      i_im2.width() != m_width) {
    cout << "im1 and im2 must have the same size. Abort.";
    exit(EXIT_FAILURE);
  }

  //! Image initialization
  m_im1     = new Image(i_im1);
  m_im2     = new Image(i_im2);

  //! Size initialization
  m_imDmin1 = new Image(m_width, m_height);
  m_imDmax1 = new Image(m_width, m_height);
  m_imDmin2 = new Image(m_width, m_height);
  m_imDmax2 = new Image(m_width, m_height);

  m_imDisp1 = new Image(m_width, m_height);
  m_imDisp2 = new Image(m_width, m_height);
  m_imDist1 = new Image(m_width, m_height);
  m_imDist2 = new Image(m_width, m_height);
  m_imMask1 = new Image(m_width, m_height);
  m_imMask2 = new Image(m_width, m_height);

  //! Value initialization
  *m_imDmin1 =  m_params->minDisp();
  *m_imDmax1 =  m_params->maxDisp();
  *m_imDmin2 = -m_params->maxDisp();
  *m_imDmax2 = -m_params->minDisp();

  *m_imDist1 = INFINITY;
  *m_imDist2 = INFINITY;
  *m_imDisp1 = 0.f;
  *m_imDisp2 = 0.f;
  *m_imMask1 = 0.f;
  *m_imMask2 = 0.f;

  //! Windows initialization
  m_windows = new Image[m_nbWindows];
  fillWindows(m_windows, *m_params, m_nbWindows);

  //! Translated images initialization
  //! Get the 4 translated images (inPrecision = 4) for stereo one direction
  m_imTrans1 = new Image();
  getTranslated(*m_im1, *m_imTrans1);
  m_imTrans2 = new Image();
  getTranslated(*m_im2, *m_imTrans2);

  //! Get the 2 translated images for strobe effect detection
  m_imStrobe11 = new Image(*m_im1);
  m_imStrobe11->translate(+0.125);
  m_imStrobe12 = new Image(*m_im1);
  m_imStrobe12->translate(-0.125);
  m_imStrobe21 = new Image(*m_im2);
  m_imStrobe21->translate(+0.125);
  m_imStrobe22 = new Image(*m_im2);
  m_imStrobe22->translate(-0.125);
}


//! Run the MSMW estimator.
void MSMW::run(
  Image& o_imDisp1,
  Image& o_imDisp2,
  Image& o_imDist1,
  Image& o_imDist2,
  Image& o_imMask1,
  Image& o_imMask2) {

  //! Check that the initialization have been made
  if (m_params == NULL) {
    cout << "MSMW: the initialization must be done before calling run. Abort\n";
    exit(EXIT_FAILURE);
  }

  //! Launch the algorithm
  this->runMultiscaleChainRecursive();

  //! Get the results
  o_imDisp1 = *m_imDisp1;
  o_imDisp2 = *m_imDisp2;
  o_imDist1 = *m_imDist1;
  o_imDist2 = *m_imDist2;
  o_imMask1 = *m_imMask1;
  o_imMask2 = *m_imMask2;
}


//! Recursive call to the MSMW multiscale algorithm.
void MSMW::runMultiscaleChainRecursive() {

  //! Initial scale is 0
  const int currentScale = m_currentScale;

  if (m_currentScale < m_params->nbScales()) {

    //! Subsample by factor 2
    Image imSub1, imSub2;
    m_im1->subSample(imSub1);
    m_im2->subSample(imSub2);

    //! Create new class
    MSMW msmw;
    msmw.init(imSub1, imSub2, *m_params);

    //! Subsample the min / max images
    m_imDmin1->subSample(*msmw.m_imDmin1);
    m_imDmax1->subSample(*msmw.m_imDmax1);
    m_imDmin2->subSample(*msmw.m_imDmin2);
    m_imDmax2->subSample(*msmw.m_imDmax2);

    //! Multiply dmin / dmax by 0.5
    *msmw.m_imDmin1 *= 0.5;
    *msmw.m_imDmax1 *= 0.5;
    *msmw.m_imDmin2 *= 0.5;
    *msmw.m_imDmax2 *= 0.5;

    //! Update the current scale
    msmw.m_currentScale = m_currentScale + 1;

    //! Recursive call
    msmw.runMultiscaleChainRecursive();

    const size_t hs = msmw.m_imDmin1->height();
    const size_t ws = msmw.m_imDmin1->width ();

    //! Up-scale and update the images
    for (size_t i = 0; i < m_height; i++) {
      const size_t is = std::min((size_t) floor(i * 0.5f), hs - 1);
      const float* sMin1 = msmw.m_imDmin1->getPtr(0, is);
      const float* sMax1 = msmw.m_imDmax1->getPtr(0, is);
      const float* sMin2 = msmw.m_imDmin2->getPtr(0, is);
      const float* sMax2 = msmw.m_imDmax2->getPtr(0, is);

      float* oMin1 = m_imDmin1->getPtr(0, i);
      float* oMax1 = m_imDmax1->getPtr(0, i);
      float* oMin2 = m_imDmin2->getPtr(0, i);
      float* oMax2 = m_imDmax2->getPtr(0, i);

      for (size_t j = 0; j < m_width; j++) {
        const size_t js = std::min((size_t) floor(j * 0.5f), ws - 1);
        oMin1[j] = std::max(oMin1[j], 2.f * sMin1[js] - 2.f);
        oMin2[j] = std::max(oMin2[j], 2.f * sMin2[js] - 2.f);
        oMax1[j] = std::min(oMax1[j], 2.f * sMax1[js] + 2.f);
        oMax2[j] = std::min(oMax2[j], 2.f * sMax2[js] + 2.f);
      }
    }
  }

  //! Build parameter structure for current scale
  m_params->update(currentScale);

  //! Apply stereo chain: dmin, dmax, idmin and idmax are updated inside stereo_pixel_chain
  this->runChainMultiWindow();
}


//! This function implements multi-window for a current scale.
void MSMW::runChainMultiWindow() {

  //! Loop over the orientations
  for (size_t n = 0; n < m_nbWindows; n++) {

    //! Update the prolate according to the current orientation
    *m_params->window() = m_windows[n];

    //! Print some informations
    cout << "multiwindow each scale --> w:" << m_params->window()->width();
    cout << " h: " << m_params->window()->height() << " s: ";
    cout << m_currentScale << endl;

    //! Call the function
    this->runPixelChain();
  }

  //! Check left right consistency
  this->checkPixelianReciprocity();

  //! Remove isolated points
  this->removeIsolatedPoints();

  //! Perform Grain Filter
  this->performGrainFilter();

  //! Initialize Dmin and Dmax
  *m_imDmin1 = INFINITY; *m_imDmax1 = -INFINITY;
  *m_imDmin2 = INFINITY; *m_imDmax2 = -INFINITY;

  //! Update Dmin/Dmax for each window
  for (size_t n = 0; n < m_nbWindows; n++) {
    this->updateWindowBoundary(m_windows[n]);
  }
}


//! Run the chain for a given window.
void MSMW::runPixelChain() {

  //! Keep track of previous results
  Image imDisp1 = *m_imDisp1;
  Image imDisp2 = *m_imDisp2;
  Image imDist1 = *m_imDist1;
  Image imDist2 = *m_imDist2;
  Image imMask1 = *m_imMask1;
  Image imMask2 = *m_imMask2;

  //! Images initialization
  Image imSelf1(m_width, m_height, m_channels);
  Image imSelf2(m_width, m_height, m_channels);

  //! Memory for disparity and mask of selected pixel taking left image as reference
  this->runPixelChain1D(imSelf1, true);

  //! Memory for disparity and mask of selected pixel taking right image as reference
  this->runPixelChain1D(imSelf2, false);

  //! Perform MinDist: consider the rejections of the current mask
  this->checkMinDist();

  //! Perform SelfSimilarity
  this->checkSelf(imSelf1, true);
  this->checkSelf(imSelf2, false);

  //! Perform Reciprocal checking if flag activated
  this->checkPixelianReciprocity();

  //! Perform Grain Filter
  this->performGrainFilter();

  //! Select minimum distance
  for (size_t i = 0; i < m_height; i++) {
    const float* iM1 = imMask1.getPtr(0, i);
    const float* iM2 = imMask2.getPtr(0, i);
    const float* iT1 = imDist1.getPtr(0, i);
    const float* iT2 = imDist2.getPtr(0, i);
    const float* iP1 = imDisp1.getPtr(0, i);
    const float* iP2 = imDisp2.getPtr(0, i);

    float* oM1 = m_imMask1->getPtr(0, i);
    float* oM2 = m_imMask2->getPtr(0, i);
    float* oT1 = m_imDist1->getPtr(0, i);
    float* oT2 = m_imDist2->getPtr(0, i);
    float* oP1 = m_imDisp1->getPtr(0, i);
    float* oP2 = m_imDisp2->getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      if (!(oM1[j] > 0.f && oT1[j] < iT1[j])) {
        oT1[j] = iT1[j];
        oP1[j] = iP1[j];
        oM1[j] = iM1[j];
      }
      if (!(oM2[j] > 0.f && oT2[j] < iT2[j])) {
        oT2[j] = iT2[j];
        oP2[j] = iP2[j];
        oM2[j] = iM2[j];
      }
    }
  }
}


//! Pixelian correlation over one direction.
void MSMW::runPixelChain1D(
  Image &o_imSelf,
  const bool p_isFirst) {

  //! Parameters
  const float inPrecisions           = 4.f;
  const size_t cubicRefinementFactor = 32;
  const bool normalize               = m_params->dist();

  //! Variable: the correlation window must be normalized
  // is it necessary ? Yep. Strange.
  Image window = *m_params->window();
  window.normalizeL1();

  //! Variable: boundary
  const int wc = (window.width () - 1) / 2;
  const int hc = (window.height() - 1) / 2;
  const int boundary = 2 * std::max(wc, hc) + 1;

  //! Initialization
  (p_isFirst ? *m_imDist1 : *m_imDist2) = INFINITY;
  o_imSelf = INFINITY;
  const float step = 1.0 / inPrecisions;
  const float increment = 1.f / float(cubicRefinementFactor);

  //! For convenience
  const int wi = m_width;
  const int hi = m_height;
  const __m256 xPos256 = _mm256_set_ps(7.f * step, 6.f * step, 5.f * step,
    4.f * step, 3.f * step, 2.f * step, step, 0.f);
  const __m128 xPos128 = _mm_set_ps(3.f * step, 2.f * step, step, 0.f);
  const __m256 x2_256 = _mm256_set1_ps(2.f);
  const __m128 x2_128 = _mm_set1_ps(2.f);
  const __m256 xStep256 = _mm256_set1_ps(step);
  const __m128 xStep128 = _mm_set1_ps(step);

//! Begin for each pixel
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for schedule(dynamic) nowait
#endif // _OPENMP
  for (int i = boundary; i < hi - boundary ; i++) {
    const float* iMin = (p_isFirst ? m_imDmin1 : m_imDmin2)->getPtr(0, i);
    const float* iMax = (p_isFirst ? m_imDmax1 : m_imDmax2)->getPtr(0, i);
    float* oSD = o_imSelf.getPtr(0, i);
    float* oP  = (p_isFirst ? m_imDisp1 : m_imDisp2)->getPtr(0, i);
    float* oT  = (p_isFirst ? m_imDist1 : m_imDist2)->getPtr(0, i);

    //! Compute correlation only at points with imask>0
    for (int j = boundary; j < wi - boundary; j++) {

      //! Correlation with first image
      {
        //! range adapted to window size
        //! The self-similarity range is centered at the reference pixel
        //! This permits to spot self similar structures even if the search
        //! range is shifted
        const int r    = std::max(std::ceil((iMax[j] - iMin[j]) / 2.0), 3.0);
        const int imin = std::max(boundary, j - (int) r);
        const int imax = std::min(wi - boundary, j + (int) r );

        //! Compute distance : WARNING: only works if inPrecision = 4.
        float fSelfBestDistance = INFINITY;
        int k = imin;

        //! AVX version
        __m256 xBest256 = _mm256_set1_ps(INFINITY);
        for (; k <= imax - 2; k += 2) {

          //! Compute the distance over all channels for the 4 translated images
          const __m256 xDist = computePatchDistance8(p_isFirst ? *m_im1 : *m_im2
            , p_isFirst ? *m_imTrans1 : *m_imTrans2, j - wc, i - hc, k - wc,
            i - hc, window, normalize);

          //! Keep the best one
          const __m256 xDelta = _mm256_set1_ps(k - j) + xPos256;
          const __m256 xMask = _mm256_and_ps(_mm256_cmp_ps(xDist, xBest256, _CMP_LT_OS),
            _mm256_cmp_ps(_mm256_max_ps(xDelta, -xDelta), x2_256, _CMP_GE_OS));
          xBest256 = applyMask256_ps(xMask, xDist, xBest256);
        }

        //! Update the best self distance
        getBest8(xBest256, xBest256, fSelfBestDistance, fSelfBestDistance);

        //! SSE version
        __m128 xBest128 = _mm_set1_ps(INFINITY);
        for (; k <= imax; k++) {

          //! Compute the distance over all channels for the 4 translated images
          const __m128 xDist = computePatchDistance4(p_isFirst ? *m_im1 : *m_im2
            , p_isFirst ? *m_imTrans1 : *m_imTrans2, j - wc, i - hc, k - wc,
            i - hc, window, normalize);

          //! Keep the best one
          const __m128 xDelta = _mm_set1_ps(k - j) + xPos128;
          const __m128 xMask  = _mm_and_ps(_mm_cmplt_ps(xDist, xBest128),
            _mm_cmpge_ps(_mm_max_ps(xDelta, -xDelta), x2_128));
          xBest128 = applyMask_ps(xMask, xDist, xBest128);
        }

        //! Update the best self distance
        getBest4(xBest128, xBest128, fSelfBestDistance, fSelfBestDistance);

        //! Store the best distance
        oSD[j] = fSelfBestDistance;
      }

      //! Correlation with second image
      {
        //! Range adapted to window size
        const int imin = std::max(boundary, j + (int) iMin[j]);
        const int imax = std::min(wi - boundary, j + (int) iMax[j]);

        //! Distance function taking into account translations
        const int length = inPrecisions * (abs(imax - imin) + 1);
        float* distance = (float*) memalloc(16, length * sizeof(float));
        for (int l = 0; l < length; l++) {
          distance[l] = INFINITY;
        }

        //! Compute distance - Initialization
        float fBestDistance = INFINITY;
        float fSelectedPixel = 0.0;
        int k = imin;

        //! AVX version
        __m256 xBest256 = _mm256_set1_ps(INFINITY);
        __m256 xSel256  = _mm256_set1_ps(0.f);
        for (; k <= imax - 2; k += 2) {

          //! Compute the distance of two pixels for all 4 precision at once
          const __m256 xDist = computePatchDistance8(p_isFirst ? *m_im1 : *m_im2
            , p_isFirst ? *m_imTrans2 : *m_imTrans1, j - wc, i - hc, k - wc,
            i - hc, window, normalize);

          //! Store the distances
          _mm256_storeu_ps(distance + 4 * (k - imin), xDist);

          //! Keep the best one
          const __m256 xMask = _mm256_cmp_ps(xDist, xBest256, _CMP_LT_OS);
          xBest256 = applyMask256_ps(xMask, xDist, xBest256);
          xSel256  = applyMask256_ps(xMask, _mm256_set1_ps(k) + xPos256, xSel256);
        }

        //! Keep the best
        getBest8(xBest256, xSel256, fBestDistance, fSelectedPixel);

        //! SSE version
        __m128 xBest128 = _mm_set1_ps(INFINITY);
        __m128 xSel128  = _mm_set1_ps(0.f);
        for (; k <= imax; k++) {

          //! Compute the distance for all 4 precision at once
          const __m128 xDist = computePatchDistance4(p_isFirst ? *m_im1 : *m_im2
            , p_isFirst ? *m_imTrans2 : *m_imTrans1, j - wc, i - hc, k - wc,
            i - hc, window, normalize);

          //! Store the distances
          _mm_storeu_ps(distance + 4 * (k - imin), xDist);

          //! Keep the best one
          const __m128 xMask = _mm_cmplt_ps(xDist, xBest128);
          xBest128 = applyMask_ps(xMask, xDist, xBest128);
          xSel128  = applyMask_ps(xMask, _mm_set1_ps(k) + xPos128, xSel128);
        }

        //! Keep the best
        getBest4(xBest128, xSel128, fBestDistance, fSelectedPixel);

        //! Refine the minimum by a factor 4 using cubic interpolation
        int l = 1;

        //! AVX version
        xBest256 = _mm256_set1_ps(fBestDistance);
        xSel256  = _mm256_set1_ps(fSelectedPixel);
        for (; l < length - 2 - 8; l += 8) {
          const __m256 xl = _mm256_set_ps(l + 7, l + 6, l + 5, l + 4, l + 3,
            l + 2, l + 1, l + 0) * xStep256 + _mm256_set1_ps(imin); // imin + l * step

          for (float x = increment; x < 1; x += increment) {

            //! Compute the interpolated distance
            const __m256 xDist = cubicInterpolate256(distance + l - 1, x);

            //! Check if it's better than what we have
            const __m256 xMask = _mm256_cmp_ps(xBest256, xDist, _CMP_GT_OS);

            //! Get new best results if needed
            xBest256 = applyMask256_ps(xMask, xDist, xBest256);
            xSel256  = applyMask256_ps(xMask, xl + _mm256_set1_ps(x) * xStep256, xSel256);
          }
        }

        //! Select only one best result among 8
        getBest8(xBest256, xSel256, fBestDistance, fSelectedPixel);

        //! SSE version
        for (; l < length - 2 - 4; l += 4) {
          const __m128 xl = _mm_set_ps(l + 3, l + 2, l + 1, l + 0) * xStep128 +
            _mm_set1_ps(imin);

          for (float x = increment; x < 1; x += increment) {

            //! Compute the interpolated distance
            const __m128 xDist = cubicInterpolate128(distance + l - 1, x);

            //! Check if it's better than what we have
            const __m128 xMask = _mm_cmp_ps(xBest128, xDist, _CMP_GT_OS);

            //! Get new best results if needed
            xBest128 = applyMask_ps(xMask, xDist, xBest128);
            xSel128  = applyMask_ps(xMask, xl + _mm_set1_ps(x) * xStep128, xSel128);
          }
        }

        //! Select only one best result among 4
        getBest4(xBest128, xSel128, fBestDistance, fSelectedPixel);

        //! Normal version
        for (; l < length - 2; l++) {
          for (float x = increment; x < 1; x += increment) {

            //! Compute the interpolated distance
            const float dist = cubicInterpolate(distance + l - 1, x);

            //! Check if it's better than what we have
            if(fBestDistance > dist) {
              fSelectedPixel = (float) imin + (float(l) + x) * step;
              fBestDistance  = dist;
            }
          }
        }

        //! If pixel not invalidated
        oP[j] = fSelectedPixel - (float) j;
        oT[j] = fBestDistance;

        //! Release memory
        memfree(distance);
      }
    }
  }
#ifdef _OPENMP
}
#endif // _OPENMP
}


//! Check the pixelian reciprocity on both disparity images.
void MSMW::checkPixelianReciprocity() {

  //! For convenience
  const float threshold = m_params->valueRecip();

  //! Loop over the lines
  for (size_t i = 0; i < m_height; i++) {
    const float* iI1 = m_imDisp1->getPtr(0, i);
    const float* iI2 = m_imDisp2->getPtr(0, i);
    float* iM1 = m_imMask1->getPtr(0, i);
    float* iM2 = m_imMask2->getPtr(0, i);

    //! Keep track of both mask values for the current line
    float* m1 = (float*) memalloc(16, m_width * sizeof(float));
    float* m2 = (float*) memalloc(16, m_width * sizeof(float));
    for (size_t j = 0; j < m_width; j++) {
      m1[j] = iM1[j];
      m2[j] = iM2[j];
    }

    //! Loop over all the pixels of the current line
    for(int j = 0; j < (int) m_width; j++) {

      //! Update the first mask according to the second one
      if (m1[j] > 0.f) {
        const float dL = iI1[j];
        const int jr = j + int(rintf(dL));

        //! Update the mask
        if (jr >= 0 && jr < (int) m_width) {
          iM1[j] = m2[jr] > 0.f ? fabsf(dL + iI2[jr]) <= threshold : 0.f;
        }
      }

      //! Update the second mask according to the first one
      if (m2[j] > 0.f) {
        const float dL = iI2[j];
        const int jr = j + int(rintf(dL));

        //! Update the mask
        if (jr >= 0 && jr < (int) m_width) {
          iM2[j] = m1[jr] > 0.f ? fabsf(dL + iI1[jr]) <= threshold : 0.f;
        }
      }
    }

    //! Release memory
    memfree(m1);
    memfree(m2);
  }
}


//! Update both mask to remove isolated points.
void MSMW::removeIsolatedPoints() {

  //! For convenience
  const size_t N = m_params->window()->width();
  const float percent = m_params->valueRemoveIsolated();
  const int boundary = N + 3;
  const float Thresh = rintf(float((2 * N + 1) * (2 * N + 1)) * percent);

  //! Integral images.
  Image iim1(m_width, m_height), iim2(m_width, m_height);

  //! Pre-compute the integral image
  for (size_t i = 0; i < m_height - 1; i++) {
    const float* iM1  = m_imMask1->getPtr(0, i);
    const float* iM2  = m_imMask2->getPtr(0, i);
    const float* oC1 = iim1.getPtr(0, i);
    const float* oC2 = iim2.getPtr(0, i);
    float* oB1 = iim1.getPtr(0, i + 1);
    float* oB2 = iim2.getPtr(0, i + 1);
    float sum1 = 0.f;
    float sum2 = 0.f;

    //! Compute the integral images
    for (size_t j = 0; j < m_width - 1; j++) {
      sum1 += iM1[j] > 0.f;
      oB1[j + 1] = sum1 + oC1[j + 1];
      sum2 += iM2[j] > 0.f;
      oB2[j + 1] = sum2 + oC2[j + 1];
    }
  }

  //! Loop over the pixels
  for (size_t i = boundary; i < m_height - boundary; i++) {
    const float* iIT1 = iim1.getPtr(0, i - N - 1);
    const float* iIT2 = iim2.getPtr(0, i - N - 1);
    const float* iIB1 = iim1.getPtr(0, i + N);
    const float* iIB2 = iim2.getPtr(0, i + N);
    float* iM1 = m_imMask1->getPtr(0, i);
    float* iM2 = m_imMask2->getPtr(0, i);

    for (size_t j = boundary; j < m_width - boundary; j++) {
      const int jt = j - N - 1;
      const int jb = j + N;

      if (iM1[j] > 0.f && iIT1[jt] - iIT1[jb] - iIB1[jt] + iIB1[jb] <= Thresh) {
        iM1[j] = 0.f;
      }
      if (iM2[j] > 0.f && iIT2[jt] - iIT2[jb] - iIB2[jt] + iIB2[jb] <= Thresh) {
        iM2[j] = 0.f;
      }
    }
  }
}


//! Perform grain filter. Update both mask at once.
void MSMW::performGrainFilter() {

  //! Initialize the rep image to -1 or 1 for both mask
  int* rep1 = (int*) memalloc(16, m_width * m_height * sizeof(int));
  int* rep2 = (int*) memalloc(16, m_width * m_height * sizeof(int));
  initRep(*m_imMask1, rep1);
  initRep(*m_imMask2, rep2);

  //! Identify the connected components of positive values
  updatePositiveConnectedComponent(rep1, m_width, m_height);
  updatePositiveConnectedComponent(rep2, m_width, m_height);

  //! Remove small regions
  removeRegions(*m_imMask1, rep1, m_params->grainArea());
  removeRegions(*m_imMask2, rep2, m_params->grainArea());

  //! Release memory
  memfree(rep1);
  memfree(rep2);
}


//! Perform MinDist: consider the rejections of the current mask.
void MSMW::checkMinDist() {

  //! For convenience
  const float threshold = m_params->valueMinDist();
  const int wi = m_width;
  const int hi = m_height;

  //! Initialization
  *m_imMask1 = -1.f;
  *m_imMask2 = -1.f;

  //! Variable: boundary
  const int wc = m_params->window()->width ();
  const int hc = m_params->window()->height();

  //! For convenience
  const int wc2 = wc / 2;
  const int hc2 = hc / 2;
  // same boundary as in stereo_pixel_chain_one_direction
  const int boundary = 2 * std::max(wc2, hc2) + 1;

  //! SSE convenience
  const __m128 xm = _mm_set1_ps(-1.f);
  const __m128 xM = _mm_set1_ps( 1.f);
  const __m128 xThresh = _mm_set1_ps(threshold);

#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
  for (int i = boundary; i < hi - boundary; i++) {
    const float* iD1 = m_imDisp1->getPtr(0, i);
    const float* iD2 = m_imDisp2->getPtr(0, i);
    float* oM1 = m_imMask1->getPtr(0, i);
    float* oM2 = m_imMask2->getPtr(0, i);
    int j = boundary;

    //! SSE version
    for (; j < wi - boundary - 4; j += 4) {

      //! Initialization
      __m128 xMinDist1 = _mm_set1_ps(INFINITY);
      __m128 xMinDist2 = _mm_set1_ps(INFINITY);
      __m128 xDisp1    = _mm_loadu_ps(iD1 + j);
      __m128 xDisp2    = _mm_loadu_ps(iD2 + j);

      //! Looking window
      for (int p = -hc2; p <= hc2; p++) {
        const float* pDist1 = m_imDist1->getPtr(0, i + p) + j;
        const float* pDist2 = m_imDist2->getPtr(0, i + p) + j;
        const float* pDisp1 = m_imDisp1->getPtr(0, i + p) + j;
        const float* pDisp2 = m_imDisp2->getPtr(0, i + p) + j;
        const float* pProl  = m_params->window()->getPtr(0, p + hc2);

        for (int q = -wc2; q <= wc2; q++) {
          if (pProl[q + wc2] > 0.f) {
            const __m128 xDist1 = _mm_loadu_ps(pDist1 + q);
            const __m128 xDist2 = _mm_loadu_ps(pDist2 + q);
            const __m128 mask1  = _mm_cmplt_ps(xDist1, xMinDist1);
            const __m128 mask2  = _mm_cmplt_ps(xDist2, xMinDist2);

            //! Update values
            xMinDist1 = applyMask_ps(mask1, xDist1, xMinDist1);
            xMinDist2 = applyMask_ps(mask2, xDist2, xMinDist2);
            xDisp1    = applyMask_ps(mask1, _mm_loadu_ps(pDisp1 + q), xDisp1);
            xDisp2    = applyMask_ps(mask2, _mm_loadu_ps(pDisp2 + q), xDisp2);
          }
        }
      }

      //! Update the mask
      const __m128 xVal1 = xDisp1 - _mm_loadu_ps(iD1 + j);
      const __m128 xVal2 = xDisp2 - _mm_loadu_ps(iD2 + j);
      const __m128 mask1 = _mm_cmple_ps(_mm_max_ps(xVal1, -xVal1), xThresh);
      const __m128 mask2 = _mm_cmple_ps(_mm_max_ps(xVal2, -xVal2), xThresh);
      _mm_storeu_ps(oM1 + j, applyMask_ps(mask1, xM, xm));
      _mm_storeu_ps(oM2 + j, applyMask_ps(mask2, xM, xm));
    }

    //! Normal version
    for (; j < wi - boundary; j++) {

      //! Initialization
      float minDist1 = INFINITY;
      float minDist2 = INFINITY;
      float disp1    = iD1[j];
      float disp2    = iD2[j];

      //! Looking window
      for (int p = -hc2; p <= hc2; p++) {
        const float* pDist1 = m_imDist1->getPtr(0, i + p) + j;
        const float* pDist2 = m_imDist2->getPtr(0, i + p) + j;
        const float* pDisp1 = m_imDisp1->getPtr(0, i + p) + j;
        const float* pDisp2 = m_imDisp2->getPtr(0, i + p) + j;
        const float* pProl  = m_params->window()->getPtr(0, p + hc2);

        for (int q = -wc2; q <= wc2; q++) {
          if (pProl[q + wc2] > 0.f && pDist1[q] < minDist1) {
            minDist1 = pDist1[q];
            disp1    = pDisp1[q];
          }
          if (pProl[q + wc2] > 0.f && pDist2[q] < minDist2) {
            minDist2 = pDist2[q];
            disp2    = pDisp2[q];
          }
        }
      }

      //! Update the mask
      oM1[j] = fabsf(disp1 - iD1[j]) <= threshold ? 1.f : -1.f;
      oM2[j] = fabsf(disp2 - iD2[j]) <= threshold ? 1.f : -1.f;
    }
  }
}


//! Check strobe and self similarity effect.
void MSMW::checkSelf(
  const Image &i_imSelf,
  const bool p_isFirst) {

  //! Variable: boundary
  const int wp  = m_params->window()->width ();
  const int hp  = m_params->window()->height();
  const int wp2 = (wp - 1) / 2;
  const int hp2 = (hp - 1) / 2;
  const int boundary = 2 * std::max(wp2, hp2) + 1;

#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
  for(size_t i = boundary ;  i < m_height - boundary ; i++) {
    const float* iD1 = (p_isFirst ? m_imDist1 : m_imDist2)->getPtr(0, i);
    const float* iD2 = i_imSelf.getPtr(0, i);
    float* iM = (p_isFirst ? m_imMask1 : m_imMask2)->getPtr(0, i);

    for(size_t j = boundary ;  j < m_width - boundary; j++) {
      if (iM[j] > 0.f) {

        //! Initialization
        float d1 = 0.f, d2 = 0.f;
        __m128 xD1 = _mm_setzero_ps(), xD2 = _mm_setzero_ps();

        //! Compute the distance over all channels
        for (size_t c = 0; c < m_channels; c++) {
          for (int p = 0; p < hp; p++) {
            const float* iI  = (p_isFirst ? m_im1 : m_im2)->getPtr(c, i - hp2 + p) + j - wp2;
            const float* iT1 = (p_isFirst ? m_imStrobe11 : m_imStrobe21)->getPtr(c, i - hp2 + p) + j - wp2;
            const float* iT2 = (p_isFirst ? m_imStrobe12 : m_imStrobe22)->getPtr(c, i - hp2 + p) + j - wp2;
            const float* iK = m_params->window()->getPtr(0, p);
            int q = 0;

            //! SSE version
            for (; q < wp - 4; q += 4) {
              const __m128 xVal1 = _mm_loadu_ps(iI + q) - _mm_loadu_ps(iT1 + q);
              const __m128 xVal2 = _mm_loadu_ps(iI + q) - _mm_loadu_ps(iT2 + q);
              xD1 += _mm_loadu_ps(iK + q) * xVal1 * xVal1;
              xD2 += _mm_loadu_ps(iK + q) * xVal2 * xVal2;
            }

            //! Normal version
            for (; q < wp; q++) {
              const float dVal1 = iI[q] - iT1[q];
              const float dVal2 = iI[q] - iT2[q];
              d1 += iK[q] * dVal1 * dVal1;
              d2 += iK[q] * dVal2 * dVal2;
            }
          }
        }

        //! Get the normalized result
        float vDist1[4], vDist2[4];
        _mm_storeu_ps(vDist1, xD1);
        _mm_storeu_ps(vDist2, xD2);
        d1 = (d1 + vDist1[0] + vDist1[1] + vDist1[2] + vDist1[3]) / float(m_channels);
        d2 = (d2 + vDist2[0] + vDist2[1] + vDist2[2] + vDist2[3]) / float(m_channels);

        //! Update the mask
        iM[j] = iD2[j] - iD1[j] > std::max(d1, d2);
      }
    }
  }
}


//! Update Dmin and Dmax for both images.
void MSMW::updateWindowBoundary(
  const Image &p_window) {

  //! Parameters
  const float min1 = m_params->dmin1();
  const float max1 = m_params->dmax1();
  const float min2 = m_params->dmin2();
  const float max2 = m_params->dmax2();

  //! For convenience
  const int wp = (p_window.width () - 1) / 2;
  const int hp = (p_window.height() - 1) / 2;
  const int wi = m_width;
  const int hi = m_height;

  //! SSE convenience
  const __m128 x0    = _mm_setzero_ps();
  const __m128 xMin1 = _mm_set1_ps(min1);
  const __m128 xMax1 = _mm_set1_ps(max1);
  const __m128 xMin2 = _mm_set1_ps(min2);
  const __m128 xMax2 = _mm_set1_ps(max2);

#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
  for (int i = hp; i < hi - hp; i++) {
    const float* iM1 = m_imMask1->getPtr(0, i);
    const float* iM2 = m_imMask2->getPtr(0, i);

    float* oDmin1 = m_imDmin1->getPtr(0, i);
    float* oDmax1 = m_imDmax1->getPtr(0, i);
    float* oDmin2 = m_imDmin2->getPtr(0, i);
    float* oDmax2 = m_imDmax2->getPtr(0, i);

    int j = wp;

    //! SSE version
    for (; j < wi - wp - 4; j += 4) {

      //! Initialization
      const __m128 mask1 = _mm_cmpgt_ps(_mm_loadu_ps(iM1 + j), x0);
      const __m128 mask2 = _mm_cmpgt_ps(_mm_loadu_ps(iM2 + j), x0);

      const __m128 vMin1  = _mm_loadu_ps(oDmin1 + j);
      const __m128 vMax1  = _mm_loadu_ps(oDmax1 + j);
      const __m128 vMin2  = _mm_loadu_ps(oDmin2 + j);
      const __m128 vMax2  = _mm_loadu_ps(oDmax2 + j);

      __m128 mmin1 = vMin1;
      __m128 mmax1 = vMax1;
      __m128 mmin2 = vMin2;
      __m128 mmax2 = vMax2;

      //! Loop over the window
      for (int p = -hp; p <= hp; p++) {
        const float* iI1 = m_imDisp1->getPtr(0, i + p) + j;
        const float* iI2 = m_imDisp2->getPtr(0, i + p) + j;
        const float* pM1 = m_imMask1->getPtr(0, i + p) + j;
        const float* pM2 = m_imMask2->getPtr(0, i + p) + j;
        const float* iW  = p_window  .getPtr(0, p + hp);

        for (int q = -wp; q <= wp; q++) {
          if (iW[q + wp] > 0.f) {

            const __m128 cmp1 = _mm_cmpgt_ps(_mm_loadu_ps(pM1 + q), x0);
            const __m128 cmp2 = _mm_cmpgt_ps(_mm_loadu_ps(pM2 + q), x0);
            const __m128 xI1  = _mm_loadu_ps(iI1 + q);
            const __m128 xI2  = _mm_loadu_ps(iI2 + q);

            mmin1 = applyMask_ps(cmp1, _mm_min_ps(mmin1, xI1), mmin1);
            mmax1 = applyMask_ps(cmp1, _mm_max_ps(mmax1, xI1), mmax1);
            mmin2 = applyMask_ps(cmp2, _mm_min_ps(mmin2, xI2), mmin2);
            mmax2 = applyMask_ps(cmp2, _mm_max_ps(mmax2, xI2), mmax2);
          }
        }
      }

      //! Update the values
      _mm_storeu_ps(oDmin1 + j, applyMask_ps(mask1, mmin1, _mm_min_ps(xMin1, vMin1)));
      _mm_storeu_ps(oDmax1 + j, applyMask_ps(mask1, mmax1, _mm_max_ps(xMax1, vMax1)));
      _mm_storeu_ps(oDmin2 + j, applyMask_ps(mask2, mmin2, _mm_min_ps(xMin2, vMin2)));
      _mm_storeu_ps(oDmax2 + j, applyMask_ps(mask2, mmax2, _mm_max_ps(xMax2, vMax2)));
    }

    //! Normal version
    for (; j < wi -wp; j++) {

      //! First image
      if (iM1[j] > 0.0f) {

        //! Initialization
        float fmin = oDmin1[j];
        float fmax = oDmax1[j];

        //! Loop over the window
        for (int p = -hp; p <= hp; p++) {
          const float* iI = m_imDisp1->getPtr(0, i + p) + j;
          const float* pM = m_imMask1->getPtr(0, i + p) + j;
          const float* iW = p_window  .getPtr(0, p + hp);

          for (int q = -wp; q <= wp; q++) {
            if (iW[q + wp] > 0.f && pM[q] > 0.f) {
              fmin = std::min(iI[q], fmin);
              fmax = std::max(iI[q], fmax);
            }
          }
        }

        //! Update the values
        oDmin1[j] = fmin;
        oDmax1[j] = fmax;
      }
      else {
        oDmin1[j] = std::min(min1, oDmin1[j]);
        oDmax1[j] = std::max(max1, oDmax1[j]);
      }

      //! Second image
      if (iM2[j] > 0.0f) {

        //! Initialization
        float fmin = oDmin2[j];
        float fmax = oDmax2[j];

        //! Loop over the window
        for (int p = -hp; p <= hp; p++) {
          const float* iI = m_imDisp2->getPtr(0, i + p) + j;
          const float* pM = m_imMask2->getPtr(0, i + p) + j;
          const float* iW = p_window  .getPtr(0, p + hp);

          for (int q = -wp; q <= wp; q++) {
            if (iW[q + wp] > 0.f && pM[q] > 0.f) {
              fmin = std::min(iI[q], fmin);
              fmax = std::max(iI[q], fmax);
            }
          }
        }

        //! Update the values
        oDmin2[j] = fmin;
        oDmax2[j] = fmax;
      }
      else {
        oDmin2[j] = std::min(min2, oDmin2[j]);
        oDmax2[j] = std::max(max2, oDmax2[j]);
      }
    }
  }
}


//! Release the memory
void MSMW::releaseMemory() {

  //! Release the memory
  if (m_im1 != NULL) {
    delete m_im1;
    m_im1 = NULL;
  }
  if (m_im2 != NULL) {
    delete m_im2;
    m_im2 = NULL;
  }
  if (m_imDmin1 != NULL) {
    delete m_imDmin1;
    m_imDmin1 = NULL;
  }
  if (m_imDmax1 != NULL) {
    delete m_imDmax1;
    m_imDmax1 = NULL;
  }
  if (m_imDmin2 != NULL) {
    delete m_imDmin2;
    m_imDmin2 = NULL;
  }
  if (m_imDmax2 != NULL) {
    delete m_imDmax2;
    m_imDmax2 = NULL;
  }
  if (m_imDisp1 != NULL) {
    delete m_imDisp1;
    m_imDisp1 = NULL;
  }
  if (m_imDisp2 != NULL) {
    delete m_imDisp2;
    m_imDisp2 = NULL;
  }
  if (m_imDist1 != NULL) {
    delete m_imDist1;
    m_imDist1 = NULL;
  }
  if (m_imDist2 != NULL) {
    delete m_imDist2;
    m_imDist2 = NULL;
  }
  if (m_imMask1 != NULL) {
    delete m_imMask1;
    m_imMask1 = NULL;
  }
  if (m_imMask2 != NULL) {
    delete m_imMask2;
    m_imMask2 = NULL;
  }

  //! Release parameters
  if (m_params != NULL) {
    delete m_params;
    m_params = NULL;
  }

  //! Release windows
  if (m_windows != NULL) {
    delete[] m_windows;
    m_windows = NULL;
  }

  //! Release translated images
  if (m_imTrans1 != NULL) {
    delete m_imTrans1;
    m_imTrans1 = NULL;
  }
  if (m_imTrans2 != NULL) {
    delete m_imTrans2;
    m_imTrans2 = NULL;
  }
  if (m_imStrobe11 != NULL) {
    delete m_imStrobe11;
    m_imStrobe11 = NULL;
  }
  if (m_imStrobe12 != NULL) {
    delete m_imStrobe12;
    m_imStrobe12 = NULL;
  }
  if (m_imStrobe21 != NULL) {
    delete m_imStrobe21;
    m_imStrobe21 = NULL;
  }
  if (m_imStrobe22 != NULL) {
    delete m_imStrobe22;
    m_imStrobe22 = NULL;
  }
}
