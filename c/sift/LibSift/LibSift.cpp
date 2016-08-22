/**
 * @file LibSift.cpp
 *
 * @brief This class implements the SIFT method, based on its implementation
 *        initially done by I. Rey Otero, published on IPOL:
 *        "Anatomy of the SIFT Method."
 *        I. Rey Otero  and  M. Delbracio
 *        Image Processing Online, 2013.
 *        http://www.ipol.im/pub/algo/rd_anatomy_sift/
 *
 *        An IPOL demo is available at
 *        http://www.ipol.im/pub/demo/rd_anatomy_sift/
 *
 *        This implementation is identical to it, but rethought to increase its
 *        speed by an average factor of 4.
 *
 * @author (original version) I. Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 * @author (this version) Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <math.h>
#include <fstream>
#include <iostream>


//! Local includes
#include "LibSift.h"
#include "../Utilities/Memory.h"
#include "../Utilities/Utilities.h"


using namespace std;


//! Default constructor
Sift::Sift() :

  //! List of keypoints
  m_keyPoints((list<KeyPoint*>*) new list<KeyPoint*>(0)),
  
  //! Set of parameters
  m_params((Parameters*) new Parameters()),

  //! Image (gray version)
  m_im(NULL),
  m_width(0),
  m_height(0),

  //! Scalespaces
  m_s(NULL),
  m_d(NULL),
  m_sx(NULL),
  m_sy(NULL),

  //! Miscallaneous
  m_time((Time*) new Time()) {
}

//! Surcharged constructor.
Sift::Sift(
  const Parameters& p_params) :

  //! List of keypoints
  m_keyPoints((list<KeyPoint*>*) new list<KeyPoint*>(0)),

  //! Set of parameters
  m_params((Parameters*) new Parameters(p_params)),

  //! Image (gray version)
  m_im(NULL),
  m_width(0),
  m_height(0),

  //! Scalespaces
  m_s(NULL),
  m_d(NULL),
  m_sx(NULL),
  m_sy(NULL),

  //! Miscellaneous
  m_time((Time*) new Time()) {

}


//! Copy constructor.
Sift::Sift(
  const Sift& i_sift) :

  //! List of keypoints
  m_keyPoints((list<KeyPoint*>*) new list<KeyPoint*>(0)),

  //! Set of parameters
  m_params((Parameters*) new Parameters(*i_sift.m_params)),

  //! Image (gray version)
  m_im(i_sift.m_im == NULL ? NULL : (Image*) new Image(*i_sift.m_im)),
  m_width(i_sift.m_width),
  m_height(i_sift.m_height),

  //! Scalespaces
  m_s (i_sift.m_s  == NULL ? NULL : (ScaleSpace*) new ScaleSpace(*i_sift.m_s )),
  m_d (i_sift.m_d  == NULL ? NULL : (ScaleSpace*) new ScaleSpace(*i_sift.m_d )),
  m_sx(i_sift.m_sx == NULL ? NULL : (ScaleSpace*) new ScaleSpace(*i_sift.m_sx)),
  m_sy(i_sift.m_sy == NULL ? NULL : (ScaleSpace*) new ScaleSpace(*i_sift.m_sy)),

  //! Miscellaneous
  m_time((Time*) new Time(*i_sift.m_time)) {

  //! Copy the keypoints
  list<KeyPoint*>::iterator key = i_sift.m_keyPoints->begin();
  while (key != i_sift.m_keyPoints->end()) {
    m_keyPoints->push_back(new KeyPoint(**key));
    key++;
  }
}


//! Operator overload.
Sift& Sift::operator=(
  const Sift& i_sift) {

  if (&i_sift == this) {
    return *this;
  }

  releaseMemory();

  new (this) Sift(i_sift);
  return *this;
}


//! Default destructor
Sift::~Sift() {
  releaseMemory();
}


//! Compute the list of keypoints from the input image.
 void Sift::computeKeyPoints(
  const Image& i_im) {

  //! Initialize the scale space and the image
  this->init(i_im);

  //! Adapt the threshold to the scalespace discretization
  const float thresh = convertThreshold();

  //! Compute the scale space
  this->computeScaleSpace();
  this->computeDoG();

  //! Keypoint detection
  this->find3dDiscreteExtrema();

  //! Update the keypoints list
  this->discardKeyPointsWithLowResponse(0.8f * thresh);
  this->interpolateKeyPointsPosition();
  this->discardKeyPointsWithLowResponse(thresh);
  this->computeEdgeResponse();
  this->discardKeyPointsOnEdge();
  this->discardKeyPointsNearBorder();

  //! Precompute gradient scalespace
  this->computeScaleSpaceGradient();

  //! Keypoint description
  this->attributeKeyPointsOrientations();
  this->attributeKeyPointsDescriptors();
}


//! Initialization
void Sift::init(
  const Image& i_im) {

  //! Size of the input image
  m_width  = i_im.width ();
  m_height = i_im.height();

  //! Convert the image into grayscale
  m_im = new Image(m_width, m_height);
  i_im.rgb2gray(*m_im);

  //! Normalize the image into [0, 1[
  (*m_im) *= (1. / 256.);

  //! Determine the number of octaves.
  const size_t nbOct = getNbOctaves();

  //! Allocation of the scale-space structure
  m_s = new ScaleSpace();
  m_s->init(ScaleSpace::LOWE, nbOct, m_params->nbSpo(), m_width, m_height,
    m_params->deltaMin(), m_params->sigmaMin());
  m_d = new ScaleSpace();
  m_d->init(ScaleSpace::DOG , nbOct, m_params->nbSpo(), m_width, m_height,
    m_params->deltaMin(), m_params->sigmaMin());
  m_sx = new ScaleSpace(*m_s);
  m_sy = new ScaleSpace(*m_s);

  //! Print the time if desired
  if (m_params->verbose()) m_time->get_time(" - ScaleSpace allocation", 50);
}


//! Determine the number of octaves
size_t Sift::getNbOctaves() const {

  //! Minimal size (width or height) of images in the last octave
  const size_t hmin = 12;

  //! Size (min of width and height) of images in the first octave
  const size_t h0 = std::min(m_width, m_height) / m_params->deltaMin();

  //! Number of octaves
  return std::min(m_params->nbOct(), (size_t) (log(h0 / hmin) / M_LN2) + 1);
}


//! Adapt the threshold to the scalespace discretization.
float Sift::convertThreshold() const {

  //! Converting the threshold to make it consistent
  const float k_nspo = exp(M_LN2 / (float) m_params->nbSpo());
  const float k_3 =  exp(M_LN2 / (float) 3);
  return (k_nspo - 1) / (k_3 - 1) * m_params->dog();
}


//! Compute the Gaussian scalespace for SIFT.
void Sift::computeScaleSpace() {

  //! For convenience
  const float sigmaIn = m_params->sigmaIn();
  const size_t nbOctaves = m_s->getNbOctaves();

  //! Seed image
  const float sigmaMin = m_s->getOctave(0)->getSigma(0);
  const float deltaMin = m_s->getOctave(0)->getDelta();

  //! Loop over the octaves
  for (size_t n = 0; n < nbOctaves; n++) {

    //! For convenience
    Octave* oct = m_s->getOctave(n);

    //! Current parameters of the working octave
    const size_t nbScales = oct->getNbImages();
    const float delta = oct->getDelta();

    //! First image in the stack
    if (n == 0) {
      const float sigmaExtra = sqrt(sigmaMin * sigmaMin - sigmaIn * sigmaIn) / deltaMin;

      if (deltaMin < 1) {
        m_im->overSample(*oct->getImage(), deltaMin, 0, 0);
        oct->getImage()->applyGaussianBlur(*oct->getImage(), sigmaExtra, 0, 0);
      }
      else {
        m_im->applyGaussianBlur(*oct->getImage(), sigmaExtra, 0, 0);
      }
    }

    //! From previous octave
    else {
      m_s->getOctave(n - 1)->getImage()->subSample(*oct->getImage(), nbScales - 3, 0);
    }

    //! The rest of the image stack: add blur to previous image in the stack
    for (size_t k = 1; k < nbScales; k++) {
      const float sigmaPrev = oct->getSigma(k - 1);
      const float sigmaNext = oct->getSigma(k);
      const float sigmaExtra = sqrt(sigmaNext * sigmaNext - sigmaPrev * sigmaPrev) / delta;
      oct->getImage()->applyGaussianBlur(*oct->getImage(), sigmaExtra, k - 1, k);
    }
  }

  //! Print the elapsed time if wanted
  if (m_params->verbose()) m_time->get_time(" - Compute scale space", 50);
}


//! Compute the difference of Gaussians.
void Sift::computeDoG() {

  for (size_t o = 0; o < m_d->getNbOctaves(); o++) {
    const size_t ns = m_d->getOctave(o)->getNbImages();
    const size_t w  = m_d->getOctave(o)->getWidth();
    const size_t h  = m_d->getOctave(o)->getHeight();

    for (size_t s = 0; s < ns; s++){
      for (size_t i = 0; i < h; i++) {
        const float* iP = m_s->getOctave(o)->getPtr(s + 1, i);
        const float* iM = m_s->getOctave(o)->getPtr(s    , i);
        float* oI = m_d->getOctave(o)->getPtr(s, i);
        size_t j = 0;

#ifdef __AVX__
        //! AVX version
        for (; j < w - 8; j += 8) {
          _mm256_storeu_ps(oI + j, _mm256_loadu_ps(iP + j) - _mm256_loadu_ps(iM + j));
        }
#endif // __AVX__
#ifdef __SSE__
        //! SSE version
        for (; j < w - 4; j += 4) {
          _mm_storeu_ps(oI + j, _mm_loadu_ps(iP + j) - _mm_loadu_ps(iM + j));
        }
#endif // __SSE__

        //! Normal version
        for (;j < w; j++) {
          oI[j] = iP[j] - iM[j];
        }
      }
    }
  }

  //! Print the elapsed time if desired
  if (m_params->verbose()) m_time->get_time(" - Compute DoG", 50);
}


//! Find 3D discrete extrema.
void Sift::find3dDiscreteExtrema() {

  //! Loop over the octaves
  for (size_t o = 0; o < m_d->getNbOctaves(); o++) {
    const size_t ns = m_d->getOctave(o)->getNbImages();
    const size_t w = m_d->getOctave(o)->getWidth();
    const size_t h = m_d->getOctave(o)->getHeight();
    const float delta = m_d->getOctave(o)->getDelta();

    //! Loop through the samples of the image stack (one octave)
    for (size_t s = 1; s < ns - 1; s++) {
      for (size_t i = 1; i < h - 1; i++) {
        const float* iPT = m_d->getOctave(o)->getPtr(s - 1, i - 1);
        const float* iPC = m_d->getOctave(o)->getPtr(s - 1, i    );
        const float* iPB = m_d->getOctave(o)->getPtr(s - 1, i + 1);
        const float* iCT = m_d->getOctave(o)->getPtr(s    , i - 1);
        const float* iCC = m_d->getOctave(o)->getPtr(s    , i    );
        const float* iCB = m_d->getOctave(o)->getPtr(s    , i + 1);
        const float* iNT = m_d->getOctave(o)->getPtr(s + 1, i - 1);
        const float* iNC = m_d->getOctave(o)->getPtr(s + 1, i    );
        const float* iNB = m_d->getOctave(o)->getPtr(s + 1, i + 1);

        //! Allocation for the local extrema
        float* hasLocal = (float*) memalloc(16, w * sizeof(float));
        size_t j = 1;

#ifdef __AVX__
        //! AVX version
        const __m256 z0 = _mm256_set1_ps(0.f);
        const __m256 z1 = _mm256_set1_ps(1.f);
        for (; j < w - 1 - 8; j += 8) {
          const __m256 zVal = _mm256_loadu_ps(iCC + j);

          //! Look if the current pixel is a local minima
          __m256 zMin = _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPT + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPT + j    ), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPT + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPC + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPC + j    ), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPC + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPB + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPB + j    ), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPB + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCT + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCT + j    ), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCT + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCC + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCC + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCB + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCB + j    ), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCB + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNT + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNT + j    ), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNT + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNC + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNC + j    ), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNC + j + 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNB + j - 1), zVal, _CMP_GT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNB + j    ), zVal, _CMP_GT_OS),
                                      _mm256_cmp_ps(_mm256_loadu_ps(iNB + j + 1), zVal, _CMP_GT_OS))))))))))))))))))))))))));

          //! Look if the current pixel is a local maxima
          __m256 zMax = _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPT + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPT + j    ), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPT + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPC + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPC + j    ), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPC + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPB + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPB + j    ), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iPB + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCT + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCT + j    ), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCT + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCC + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCC + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCB + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCB + j    ), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iCB + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNT + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNT + j    ), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNT + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNC + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNC + j    ), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNC + j + 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNB + j - 1), zVal, _CMP_LT_OS),
                        _mm256_and_ps(_mm256_cmp_ps(_mm256_loadu_ps(iNB + j    ), zVal, _CMP_LT_OS),
                                      _mm256_cmp_ps(_mm256_loadu_ps(iNB + j + 1), zVal, _CMP_LT_OS))))))))))))))))))))))))));

          //! Keep track of it
          _mm256_storeu_ps(hasLocal + j, applyMask256_ps(_mm256_or_ps(zMin, zMax), z1, z0));
        }
#endif // __AVX__
#ifdef __SSE__
        //! SSE version
        const __m128 x0 = _mm_set1_ps(0.f);
        const __m128 x1 = _mm_set1_ps(1.f);
        for (; j < w - 1 - 4; j += 4) {
          const __m128 xVal = _mm_loadu_ps(iCC + j);

          //! Look if the current pixel is a local minima
          __m128 xMin = _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPT + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPT + j    ), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPT + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPC + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPC + j    ), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPC + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPB + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPB + j    ), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iPB + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCT + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCT + j    ), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCT + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCC + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCC + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCB + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCB + j    ), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iCB + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNT + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNT + j    ), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNT + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNC + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNC + j    ), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNC + j + 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNB + j - 1), xVal),
                        _mm_and_ps(_mm_cmpgt_ps(_mm_loadu_ps(iNB + j    ), xVal),
                                   _mm_cmpgt_ps(_mm_loadu_ps(iNB + j + 1), xVal))))))))))))))))))))))))));

          //! Look if the current pixel is a local maxima
          __m128 xMax = _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPT + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPT + j    ), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPT + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPC + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPC + j    ), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPC + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPB + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPB + j    ), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iPB + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCT + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCT + j    ), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCT + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCC + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCC + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCB + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCB + j    ), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iCB + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNT + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNT + j    ), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNT + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNC + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNC + j    ), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNC + j + 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNB + j - 1), xVal),
                        _mm_and_ps(_mm_cmplt_ps(_mm_loadu_ps(iNB + j    ), xVal),
                                   _mm_cmplt_ps(_mm_loadu_ps(iNB + j + 1), xVal))))))))))))))))))))))))));

          //! Keep track of it
          _mm_storeu_ps(hasLocal + j, applyMask_ps(_mm_or_ps(xMin, xMax), x1, x0));
        }
#endif // __SSE__

        //! Normal version
        for (; j < w - 1; j++) {
          const float cVal = iCC[j];

          //! Look if the current pixel is a local minima
          bool isLocalMin = true;
          isLocalMin &= (iPT[j - 1] > cVal) && (iPT[j] > cVal) && (iPT[j + 1] > cVal)
                     && (iPC[j - 1] > cVal) && (iPC[j] > cVal) && (iPC[j + 1] > cVal)
                     && (iPB[j - 1] > cVal) && (iPB[j] > cVal) && (iPB[j + 1] > cVal)
                     && (iCT[j - 1] > cVal) && (iCT[j] > cVal) && (iCT[j + 1] > cVal)
                     && (iCC[j - 1] > cVal) &&                    (iCC[j + 1] > cVal)
                     && (iCB[j - 1] > cVal) && (iCB[j] > cVal) && (iCB[j + 1] > cVal)
                     && (iNT[j - 1] > cVal) && (iNT[j] > cVal) && (iNT[j + 1] > cVal)
                     && (iNC[j - 1] > cVal) && (iNC[j] > cVal) && (iNC[j + 1] > cVal)
                     && (iNB[j - 1] > cVal) && (iNB[j] > cVal) && (iNB[j + 1] > cVal);

          //! Can skip max check if center point was determined to be a local min.
          bool isLocalMax = true;
          if (isLocalMin) {
            isLocalMax = false;
          }
          else {
            //! Check if the current pixel is a local maxima
            isLocalMax &= (iPT[j - 1] < cVal) && (iPT[j] < cVal) && (iPT[j + 1] < cVal)
                       && (iPC[j - 1] < cVal) && (iPC[j] < cVal) && (iPC[j + 1] < cVal)
                       && (iPB[j - 1] < cVal) && (iPB[j] < cVal) && (iPB[j + 1] < cVal)
                       && (iCT[j - 1] < cVal) && (iCT[j] < cVal) && (iCT[j + 1] < cVal)
                       && (iCC[j - 1] < cVal) &&                    (iCC[j + 1] < cVal)
                       && (iCB[j - 1] < cVal) && (iCB[j] < cVal) && (iCB[j + 1] < cVal)
                       && (iNT[j - 1] < cVal) && (iNT[j] < cVal) && (iNT[j + 1] < cVal)
                       && (iNC[j - 1] < cVal) && (iNC[j] < cVal) && (iNC[j + 1] < cVal)
                       && (iNB[j - 1] < cVal) && (iNB[j] < cVal) && (iNB[j + 1] < cVal);
          }

          //! Keep track of it
          hasLocal[j] = isLocalMin || isLocalMax;
        }

          //! If 3d discrete extrema, save a candidate keypoint
        for (size_t j = 1; j < w - 1; j++) {
          if (hasLocal[j]) {

            //! Create a new KeyPoint
            KeyPoint* key = new KeyPoint(m_params->nbOri(), m_params->nbHist(),
                                         m_params->nbBins());

            //! Update its values
            key->setI(i);
            key->setJ(j);
            key->setS(s);
            key->setO(o);
            key->setX(delta * float(i));
            key->setY(delta * float(j));
            key->setSigma(m_d->getOctave(o)->getSigma(s));
            key->setVal(iCC[j]);

            //! Add it to the list
            m_keyPoints->push_back(key);
          }
        }

        memfree(hasLocal);
      }
    }
  }

  //! Print the elapsed time if desired
  if (m_params->verbose()) m_time->get_time(" - Find extrema", 50);
}


//! Discard keypoints with low response.
void Sift::discardKeyPointsWithLowResponse(
  const float p_threshold) {

  //! For convenience
  std::list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Loop over all the keys
  while (key != m_keyPoints->end()) {
    if (fabs((*key)->getVal()) <= p_threshold) {

      //! Release memory
      delete *key;

      //! Remove it and go to the next key
      key = m_keyPoints->erase(key);
    }
    else {
      key++;
    }
  }

  //! Print the elapsed time if desired
  if (m_params->verbose()) m_time->get_time(" - Discard low response", 50);
}


//! Interpolate keypoints positions.
void Sift::interpolateKeyPointsPosition() {

  //! Parameters
  const float offMax = 0.6f;
  //! Ratio between two consecutive scales in the scalespace
  //! assuming the ratio is constant over all scales and over all octaves
  const float sigmaRatio = m_d->getOctave(0)->getSigma(1) /
                           m_d->getOctave(0)->getSigma(0);

  //! For convenience
  std::list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Loop over all the keypoints
  while (key != m_keyPoints->end()) {

    //! Loading keypoint and associated octave
    const size_t o    = (*key)->getO();
    const size_t s    = (*key)->getS();
    const int i       = (*key)->getI();
    const int j       = (*key)->getJ();
    const int w       = m_d->getOctave(o)->getWidth();
    const int h       = m_d->getOctave(o)->getHeight();
    //! WARNING this includes the auxiliary scales.
    const int ns      = m_d->getOctave(o)->getNbImages();
    const float delta = m_d->getOctave(o)->getDelta();

    //! Current value of i coordinate - at each interpolation
    int ic = i;
    int jc = j;
    int sc = s;
    size_t nIntrp = 0;
    bool isConv = false;
    float offX = 0.f;
    float offY = 0.f;
    float offS = 0.f;
    float val = (*key)->getVal();

    while (nIntrp < m_params->maxIter()) {

      //! Extrema interpolation via a quadratic function only if the detection
      //! is not too close to the border
      //! (so the discrete 3D Hessian is well defined)
      if ((0 < ic) && (ic < (h - 1)) && (0 < jc) && (jc < (w - 1))) {
        m_d->getOctave(o)->inverse3dTaylorSecondOrderExpansion(ic, jc, sc, offX, offY, offS, val);
      }
      else {
          isConv = false;
          offX = 5.0f;
          offY = 5.0f;
          offS = 5.0f;
      }

      //! Test if the quadratic model is consistent
      if (fabs(offX) < offMax && fabs(offY) < offMax && fabs(offS) < offMax) {
        isConv = true;
        break;
      }
      //! Move to another point
      else {
          //! Space...
          if ((offX > +offMax) && ((ic + 1) < (h - 1))) {ic++;}
          if ((offY > +offMax) && ((jc + 1) < (w - 1))) {jc++;}
          if ((offX < -offMax) && ((ic - 1) >  0     )) {ic--;}
          if ((offY < -offMax) && ((jc - 1) >  0     )) {jc--;}

          //! ... and scale.
          if ((offS > +offMax) && ((sc + 1) < (ns - 1))) {sc++;}
          if ((offS < -offMax) && ((sc - 1) >       0 )) {sc--;}
      }
      nIntrp++;
    }

    //! Keep the key
    if (isConv) {

      //! Update the key
      (*key)->setX((ic + offX) * delta);
      (*key)->setY((jc + offY) * delta);
      (*key)->setI(ic);
      (*key)->setJ(jc);
      (*key)->setS(sc);
      (*key)->setSigma(m_d->getOctave(o)->getSigma(sc) * pow(sigmaRatio, offS));
      (*key)->setVal(val);
      key++;
    }
    //! Remove the key
    else {

      //! Release memory
      delete *key;

      //! Remove it from the list and go to the next key
      key = m_keyPoints->erase(key);
    }
  }

  //! Print elapsed time if desired
  if (m_params->verbose()) m_time->get_time(" - Interpolate position", 50);
}


//! Compute edge response.
void Sift::computeEdgeResponse() {

  //! For convenience
  std::list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Loop over all the keypoints
  while (key != m_keyPoints->end()) {

    //! Loading keypoint and associated octave
    const int o = (*key)->getO();
    const int s = (*key)->getS();
    const int i = (*key)->getI();
    const int j = (*key)->getJ();
    const float* iT = m_d->getOctave(o)->getPtr(s, i - 1);
    const float* iC = m_d->getOctave(o)->getPtr(s, i    );
    const float* iB = m_d->getOctave(o)->getPtr(s, i + 1);

    //! Compute the 2D Hessian at pixel (i, j)
    const float hXX = iT[j    ] + iB[j    ] - 2 * iC[j];
    const float hYY = iC[j + 1] + iC[j - 1] - 2 * iC[j];
    const float hXY = ((iB[j + 1] - iB[j - 1]) - (iT[j + 1] - iT[j - 1])) / 4.f;

    //! Harris and Stephen Edge response computed on the DoG operator
    (*key++)->setEdgeResp((hXX + hYY) * (hXX + hYY) / (hXX * hYY - hXY * hXY));
  }

  //! Print elapsed time if desired
  if (m_params->verbose()) m_time->get_time(" - Compute edge response", 50);
}


//! Discard keypoints on edge.
void Sift::discardKeyPointsOnEdge() {

  //! For convenience
  const float threshold = (m_params->edge() + 1) * (m_params->edge() + 1) / m_params->edge();
  std::list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Loop over the keypoints
  while (key != m_keyPoints->end()) {

    //! Check if it isn't accepted
    if (fabs((*key)->getEdgeResp()) > threshold) {

      //! Release memory
      delete *key;

      //! Remove it from the list and go to the next key
      key = m_keyPoints->erase(key);
    }
    //! Keep it and move on
    else{
      key++;
    }
  }

  //! Print elapsed time if desired.
  if (m_params->verbose()) m_time->get_time(" - Discard on edge", 50);
}


//! Discard keypoints near the border.
void Sift::discardKeyPointsNearBorder() {

  //! For convenience
  const float h = m_height;
  const float w = m_width ;
  const float lambda = 1.f;
  std::list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Loop over the keypoints
  while (key != m_keyPoints->end()) {

    //! Get the current key
    const float x = (*key)->getX();
    const float y = (*key)->getY();
    const float ls = lambda * (*key)->getSigma();

    //! Check if it isn't accepted then remove it from the list
    if ((x <= ls) || (x + ls >= h) || (y <= ls) || (y + ls >= w)) {

      //! Release memory
      delete *key;

      //! Remove it from the list and go to the next key
      key = m_keyPoints->erase(key);
    }
    //! Otherwise keep it and move on.
    else {
      key++;
    }
  }

  //! Print elapsed time if desired.
  if (m_params->verbose()) m_time->get_time(" - Discard near border", 50);
}


//! Compute the scalespace gradient.
void Sift::computeScaleSpaceGradient() {

  //! For each octave
  for (size_t o = 0; o < m_s->getNbOctaves(); o++) {
    const size_t nbScales = m_s->getOctave(o)->getNbImages();

    //! For each scale
    for (size_t s = 0; s < nbScales; s++) {
      m_s->getOctave(o)->getImage()->computeGradient(
        *m_sx->getOctave(o)->getImage(),
        *m_sy->getOctave(o)->getImage(), s, s);
    }
  }

  //! Print elapsed time if desired.
  if (m_params->verbose()) m_time->get_time(" - Compute gradient", 50);
}


//! Attribute the orientation for all the keypoints.
void Sift::attributeKeyPointsOrientations() {

  //! For convenience
  const size_t nbBins = m_params->nbBins();
  std::list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Loop over the keypoints
  while (key != m_keyPoints->end()) {

    //! Accumulate gradient orientation histogram
    this->accumulateOrientationHistogram(*key);

    //! Extract principal orientation
    float* principalOrientations = (float*) memalloc(16, nbBins * sizeof(float));
    const size_t nbPrOri = (*key)->extractPrincipalOrientations(nbBins,
      m_params->t(), principalOrientations);

    //! Updating keypoints and save in new list
    for (size_t n = 0; n < nbPrOri; n++) {
      KeyPoint* newKey = new KeyPoint(**key);
      newKey->setTheta(principalOrientations[n]);
      m_keyPoints->insert(key, newKey);
    }

    //! Release memory
    memfree(principalOrientations);
    delete *key;

    //! Remove this useless keypoint from the list and go to the next key
    key = m_keyPoints->erase(key);
  }

  //! Print elapsed time if desired.
  if (m_params->verbose()) m_time->get_time(" - Attribute orientations", 50);
}


//! Attribute the descriptors for all the keypoints.
void Sift::attributeKeyPointsDescriptors() {

  //! Parameters
  const size_t nbDescr = m_params->nbHist() * m_params->nbHist() * m_params->nbOri();

  //! For convenience
  std::list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Loop over the keypoints
  for (; key != m_keyPoints->end(); key++) {

    //! Compute descriptor representation
    (*key)->extractFeatureVector(m_sx, m_sy, *m_params);

    //! Threshold and quantization of the descriptor
    (*key)->thresholdAndQuantizeFeatureVector(nbDescr);
  }

  //! Print elapsed time if desired
  if (m_params->verbose()) m_time->get_time(" - Attribute descriptors", 50);
}


//! Accumulate orientation histogram
void Sift::accumulateOrientationHistogram(
  KeyPoint* io_key) {

  //! For convenience
  const size_t nbins = m_params->nbBins();
  const size_t o = io_key->getO();
  const size_t s = io_key->getS();
  const float delta = m_sx->getOctave(o)->getDelta();
  const float lambdaOri = m_params->lambdaOri();
  const float xk = io_key->getX() / delta;
  const float yk = io_key->getY() / delta;
  const float sk = io_key->getSigma() / delta;
  float* oH = io_key->getPtrHist();
  const int w = m_sx->getOctave(o)->getWidth();
  const int h = m_sx->getOctave(o)->getHeight();

  //! Initialize output vector
  for (size_t n = 0; n < nbins; n++) {
    oH[n] = 0.f;
  }

  //! Contributing pixels are inside a patch [siMin;siMax] X [sjMin;sjMax]
  //! of width w 6*lambda_ori*sigma_key (=9*sigma_key)
  const float R = 3 * lambdaOri * sk;
  const size_t siMin = (size_t) std::max(0, (int)(xk - R + 0.5));
  const size_t sjMin = (size_t) std::max(0, (int)(yk - R + 0.5));
  const size_t siMax = (size_t) std::min((int)(xk + R + 0.5), h - 1);
  const size_t sjMax = (size_t) std::min((int)(yk + R + 0.5), w - 1);
  const size_t N = (sjMax - sjMin + 1) * (siMax - siMin + 1);

  //! For speed up
  float* vg = (float*) memalloc(16, N * sizeof(float));
  float* vm = (float*) memalloc(16, N * sizeof(float));

  //! For convenience
#ifdef __SSE__
  const __m128 xNorm = _mm_set1_ps(1.f / (2.f * lambdaOri * lambdaOri));
  const __m128 yk128 = _mm_set1_ps(yk);
  const __m128 sk128 = _mm_set1_ps(sk);
  const __m128 x2PI  = _mm_set1_ps(2 * M_PI);
#endif // __SSE__
#ifdef __AVX__
  const __m256 zNorm = _mm256_set1_ps(1.f / (2.f * lambdaOri * lambdaOri));
  const __m256 yk256 = _mm256_set1_ps(yk);
  const __m256 sk256 = _mm256_set1_ps(sk);
  const __m256 z2PI  = _mm256_set1_ps(2 * M_PI);
#endif // __AVX__

  //! For each pixel inside the patch.
  for (size_t si = siMin; si <= siMax; si++) {
    const float* iX = m_sx->getOctave(o)->getPtr(s, si);
    const float* iY = m_sy->getOctave(o)->getPtr(s, si);
    size_t sj = sjMin;
    size_t n = (si - siMin) * (sjMax - sjMin + 1);

#ifdef __AVX__
    //! AVX version
    for (; sj <= sjMax - 8; sj += 8, n += 8) {

      //! Compute pixel coordinates (sX, sY) on keypoint's invariant referential
      const __m256 sX = _mm256_set1_ps((si - xk) / sk);
      const __m256 sY = (_mm256_set_ps(sj + 7, sj + 6, sj + 5, sj + 4, sj + 3,
                         sj + 2, sj + 1, sj) - yk256) / sk256;

      //! Gradient orientation (theta)
      const __m256 dx  = _mm256_loadu_ps(iX + sj);
      const __m256 dy  = _mm256_loadu_ps(iY + sj);
      const __m256 ori = modulus_256(atan2_256(dy, dx), z2PI);

      //! Gradient magnitude with Gaussian weighting
      _mm256_storeu_ps(vm + n, hypot_256(dx, dy) * exp_256(-(sX * sX + sY * sY) * zNorm));

      //! Determine the bin index in the circular histogram
      _mm256_storeu_ps(vg + n, ori_to_bin_256(ori, nbins));
    }
#endif // __SSE__
#ifdef __SSE__
    //! SSE version
    for (; sj <= sjMax - 4; sj += 4, n += 4) {

      //! Compute pixel coordinates (sX, sY) on keypoint's invariant referential
      const __m128 sX = _mm_set1_ps((si - xk) / sk);
      const __m128 sY = (_mm_set_ps(sj + 3, sj + 2, sj + 1, sj) - yk128) / sk128;

      //! Gradient orientation (theta)
      const __m128 dx  = _mm_loadu_ps(iX + sj);
      const __m128 dy  = _mm_loadu_ps(iY + sj);
      const __m128 ori = modulus_128(atan2_128(dy, dx), x2PI);

      //! Gradient magnitude with Gaussian weighting
      _mm_storeu_ps(vm + n, hypot_128(dx, dy) * exp_128(-(sX * sX + sY * sY) * xNorm));

      //! Determine the bin index in the circular histogram
      _mm_storeu_ps(vg + n, ori_to_bin_128(ori, nbins));
    }
#endif // __SSE__

    //! Normal version
    for (; sj <= sjMax; sj++, n++) {

      //! Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
      const float sX = (si - xk) / sk;
      const float sY = (sj - yk) / sk;

      //! Gradient orientation (theta)
      const float dx = iX[sj];
      const float dy = iY[sj];
      const float ori = mod(atan2(dy, dx), 2 * M_PI);

      //! Gradient magnitude with Gaussian weighting
      const float r2 = sX * sX + sY * sY;
      vm[n] = hypot(dx, dy) * exp(-r2 / (2 * lambdaOri * lambdaOri));

      /// Determine the bin index in the circular histogram
      vg[n] = ori2bin(ori, nbins);
    }
  }

  /// Add the contribution to the orientation histogram
  for (size_t n = 0; n < N; n++) {
    oH[int(vg[n])] += vm[n];
  }

  //! Release memory
  memfree(vg);
  memfree(vm);
}


//! Get the keypoints.
void Sift::getKeyPoints(
  std::vector<KeyPoint>& o_keyPoints) const {

  //! Initialization
  o_keyPoints.clear();

  //! For convenience
  list<KeyPoint*>::iterator key = m_keyPoints->begin();

  //! Copy the list of keypoints
  for (; key != m_keyPoints->end(); key++) {
    KeyPoint newKey(**key);
    o_keyPoints.push_back(newKey);
  }
}


//! Write the keypoints into a file.
void Sift::writeKeyPoints(
  std::string const& i_fileName) const {

  //! Open the file
  ofstream myFile;
  myFile.open(i_fileName.c_str(), ios::out | ios::trunc);

  //! Check that the file has been correctly opened
  if (!myFile.is_open()) {
    cout << "Can't open the file " << i_fileName << ". Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Loop over the keys
  list<KeyPoint*>::iterator key = m_keyPoints->begin();
  for (; key != m_keyPoints->end(); key++) {

    //! For convenience
    const size_t nbDesc = (*key)->getNbHist() * (*key)->getNbHist() *
                          (*key)->getNbOri();

    //! Write the coordinates
    myFile << (*key)->getY() << " " << (*key)->getX() << " ";

    //! Write the scale and orientations
    myFile << (*key)->getSigma() << " " << (*key)->getTheta() << " ";

    //! Write the descriptors
    for (size_t n = 0; n < nbDesc; n++) {
      myFile << (*key)->getDescr(n) << " ";
    }
    myFile << endl;
  }

  //! Close the file
  myFile.close();
}


//! Release memory.
void Sift::releaseMemory() {

  //! Parameters
  if (m_params != NULL) {
    delete m_params;
    m_params = NULL;
  }

  //! List of keypoints
  if (m_keyPoints != NULL) {
    list<KeyPoint*>::iterator key = m_keyPoints->begin();
    while (key != m_keyPoints->end()) {
      delete *key;
      key++;
    }
    delete m_keyPoints;
    m_keyPoints = NULL;
  }

  //! Image (gray version)
  if (m_im != NULL) {
    delete m_im;
    m_im = NULL;
  }

  //! Scalespaces
  if (m_s != NULL) {
    delete m_s;
    m_s = NULL;
  }
  if (m_d != NULL) {
    delete m_d;
    m_d = NULL;
  }
  if (m_sx != NULL) {
    delete m_sx;
    m_sx = NULL;
  }
  if (m_sy != NULL) {
    delete m_sy;
    m_sy = NULL;
  }

  //! Miscellaneous
  if (m_time != NULL) {
    delete m_time;
    m_time = NULL;
  }
}
