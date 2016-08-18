/**
 * @file Octave.cpp
 *
 * @brief Class to handle the octaves.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <iostream>


//! Local includes
#include "Octave.h"
#include "../Utilities/Memory.h"


using namespace std;


//! Default constructor
Octave::Octave() :
  m_delta   (0.f),
  m_width   (0),
  m_height  (0),
  m_nbImages(0),
  m_sigmas  (NULL),
  m_image   (NULL) {
}


//! Copy constructor.
Octave::Octave(
  const Octave& i_octave) :
  m_delta   (i_octave.m_delta),
  m_width   (i_octave.m_width),
  m_height  (i_octave.m_height),
  m_nbImages(i_octave.m_nbImages),
  m_sigmas  ((float*) memalloc(16, m_nbImages * sizeof(float))),
  m_image   (new Image(m_width, m_height, m_nbImages)) {

  //! Copy the data
  for (size_t n = 0; n < m_nbImages; n++) {
    m_sigmas[n] = i_octave.m_sigmas[n];
  }
}


//! Surcharged constructor.
Octave::Octave(
  const float i_delta,
  const size_t i_width,
  const size_t i_height,
  const size_t i_nbImages,
  const float* i_sigmas) :

  m_delta   (i_delta),
  m_width   (i_width),
  m_height  (i_height),
  m_nbImages(i_nbImages),
  m_sigmas  ((float*) memalloc(16, m_nbImages * sizeof(float))),
  m_image   (new Image(m_width, m_height, m_nbImages)) {

  //! Copy the data
  for (size_t n = 0; n < m_nbImages; n++) {
    m_sigmas[n] = i_sigmas[n];
  }
}


//! Operator overload.
Octave& Octave::operator=(
  const Octave& i_octave) {

  if (&i_octave == this) {
    return *this;
  }

  releaseMemory();

  new (this) Octave(i_octave);
  return *this;
}


//! Default destructor
Octave::~Octave() {
  releaseMemory();
}


//! Extrema interpolation via a quadratic function.
void Octave::inverse3dTaylorSecondOrderExpansion(
  const int i,
  const int j,
  const int s,
  float& o_di,
  float& o_dj,
  float& o_ds,
  float& o_val) {

  //! For convenience
  const float* iPT = m_image->getPtr(s - 1, i - 1);
  const float* iPC = m_image->getPtr(s - 1, i    );
  const float* iPB = m_image->getPtr(s - 1, i + 1);
  const float* iCT = m_image->getPtr(s    , i - 1);
  const float* iCC = m_image->getPtr(s    , i    );
  const float* iCB = m_image->getPtr(s    , i + 1);
  const float* iNT = m_image->getPtr(s + 1, i - 1);
  const float* iNC = m_image->getPtr(s + 1, i    );
  const float* iNB = m_image->getPtr(s + 1, i + 1);

  //! Compute the 3D Hessian at pixel (i, j) of the image s. Finite difference scheme.
  const float hXX = iCT[j    ] + iCB[j    ] - 2 * iCC[j];
  const float hYY = iCC[j + 1] + iCC[j - 1] - 2 * iCC[j];
  const float hSS = iPC[j    ] + iNC[j    ] - 2 * iCC[j];
  const float hXY = 0.25 * ((iCB[j + 1] - iCB[j - 1]) - (iCT[j + 1] - iCT[j - 1]));
  const float hXS = 0.25 * ((iNB[j    ] - iNT[j    ]) - (iPB[j    ] - iPT[j    ]));
  const float hYS = 0.25 * ((iNC[j + 1] - iNC[j - 1]) - (iPC[j + 1] - iPC[j - 1]));

  //! Compute the 3D gradient at pixel (i, j) for the image s.
  const float gX = 0.5 * (iCB[j    ] - iCT[j    ]);
  const float gY = 0.5 * (iCC[j + 1] - iCC[j - 1]);
  const float gS = 0.5 * (iNC[j    ] - iPC[j    ]);

  //! Inverse the Hessian - Fitting a quadratic function
  const float det = hXX * hYY * hSS - hXX * hYS * hYS - hXY * hXY * hSS
                  + 2 * hXY * hXS * hYS - hXS * hXS * hYY;
  const float aa = (hYY * hSS - hYS * hYS) / det;
  const float ab = (hXS * hYS - hXY * hSS) / det;
  const float ac = (hXY * hYS - hXS * hYY) / det;
  const float bb = (hXX * hSS - hXS * hXS) / det;
  const float bc = (hXY * hXS - hXX * hYS) / det;
  const float cc = (hXX * hYY - hXY * hXY) / det;

  //! Offset in position (output value)
  o_di = -aa * gX - ab * gY - ac * gS;
  o_dj = -ab * gX - bb * gY - bc * gS;
  o_ds = -ac * gX - bc * gY - cc * gS;

  //! Compute the DoG value offset (output value)
  o_val = iCC[j] + 0.5 * (gX * o_di + gY * o_dj + gS * o_ds);
}


//! Release memory.
void Octave::releaseMemory() {

  memfree(m_sigmas);
  if (m_image != NULL) {
    delete m_image;
    m_image = NULL;
  }
}









