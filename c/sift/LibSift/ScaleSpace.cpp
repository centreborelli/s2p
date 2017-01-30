/**
 * @file ScaleSpace.cpp
 *
 * @brief Class to handle the ScaleSpace.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <iostream>
#include <cmath>


//! Local includes
#include "ScaleSpace.h"
#include "../Utilities/Memory.h"


using namespace std;


//! Default constructor
ScaleSpace::ScaleSpace() :
  m_nbOctaves(0),
  m_octaves  (NULL) {
}


//! Copy constructor.
ScaleSpace::ScaleSpace(
  const ScaleSpace& i_scaleSpace) :
  m_nbOctaves(i_scaleSpace.m_nbOctaves),
  m_octaves  ((Octave**) memalloc(16, m_nbOctaves * sizeof(Octave*))) {

  //! Copy the data
  for (size_t n = 0; n < m_nbOctaves; n++) {
    m_octaves[n] = new Octave(*i_scaleSpace.m_octaves[n]);
  }
}


//! Operator overload.
ScaleSpace& ScaleSpace::operator=(
  const ScaleSpace& i_scaleSpace) {

  if (&i_scaleSpace == this) {
    return *this;
  }

  releaseMemory();

  new (this) ScaleSpace(i_scaleSpace);
  return *this;
}


//! Default destructor
ScaleSpace::~ScaleSpace() {
  releaseMemory();
}


//! Initialization. Depends on the wanted space.
void ScaleSpace::init(
  ScaleSpace::Space const& p_space,
  const size_t p_nbOctaves,
  const size_t p_nbScales,
  const size_t p_width,
  const size_t p_height,
  const float p_delta,
  const float p_sigma) {

  //! Memory allocation
  m_nbOctaves = p_nbOctaves;
  m_octaves = (Octave**) memalloc(16, m_nbOctaves * sizeof(Octave*));

  //! Initialization for DoG
  if (p_space == ScaleSpace::DOG) {

    //! Sizes initialization
    const size_t h0 = (size_t) (float(p_height) / p_delta);
    const size_t w0 = (size_t) (float(p_width ) / p_delta);

    //! Initialize each octave
    for (size_t o = 0; o < m_nbOctaves; o++) {
      const float dn  = p_delta * float(1 << o);
      const size_t hn = h0 / (1 << o);
      const size_t wn = w0 / (1 << o);

      //! Sigmas initialization.
      float* sn = (float*) memalloc(16, (p_nbScales + 2) * sizeof(float));
      for (size_t s = 0; s < p_nbScales + 2; s++) {
        sn[s] = dn / p_delta * p_sigma * pow(2, float(s) / float(p_nbScales));
      }

      //! Create a new octave
      m_octaves[o] = new Octave(dn, wn, hn, p_nbScales + 2, sn);

      //! Release memory
      memfree(sn);
    }
  }

  //! Initialization for LOWE
  else if (p_space == ScaleSpace::LOWE) {

    //! Sizes initialization
    const size_t h0 = (size_t) (float(p_height) / p_delta);
    const size_t w0 = (size_t) (float(p_width ) / p_delta);

    //! Initialize each octave
    for (size_t o = 0; o < m_nbOctaves; o++) {
      const float dn  = p_delta * float(1 << o);
      const size_t hn = h0 / (1 << o);
      const size_t wn = w0 / (1 << o);

      //! Sigmas initialization. Warning, 3 extra images int the stack
      //! 1 for DoG computation and 2 for 3D discrete extrema definition.
      float* sn = (float*) memalloc(16, (p_nbScales + 3) * sizeof(float));
      for (size_t s = 0; s < p_nbScales + 3; s++) {
        sn[s] = dn / p_delta * p_sigma * pow(2, float(s) / float(p_nbScales));
      }

      //! Create a new octave
      m_octaves[o] = new Octave(dn, wn, hn, p_nbScales + 3, sn);

      //! Release memory
      memfree(sn);
    }
  }

  //! Unknown type of space
  else {
    cout << "ScaleSpace::init: unknown provided type of space. Abort." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Release memory.
void ScaleSpace::releaseMemory() {

  if (m_octaves != NULL) {
    for (size_t n = 0; n < m_nbOctaves; n++) {
      delete m_octaves[n];
    }
  }
  memfree(m_octaves);
}









