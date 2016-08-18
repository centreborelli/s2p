#ifndef SCALESPACE_H_INCLUDED
#define SCALESPACE_H_INCLUDED


//! Global includes


//! Local includes
#include "Octave.h"


/**
 * @brief A ScaleSpace contains multiple octaves.
 **/
class ScaleSpace {

  //! Methods
  public:
    /**
     * @brief Differents type of initialization.
     **/
    enum Space {
      DOG,
      LOWE
    };


    /**
     * @brief Default constructor.
     **/
    ScaleSpace();


    /**
     * @brief Copy constructor.
     **/
    ScaleSpace(
      const ScaleSpace& i_scaleSpace);


    /**
     * @brief Operator overload.
     **/
    ScaleSpace& operator=(
      const ScaleSpace& i_scaleSpace);


    /**
     *@brief Default destructor.
     **/
    ~ScaleSpace();


    /**
     * @brief Initialization. Depends on the wanted space.
     **/
    void init(
      ScaleSpace::Space const& p_space,
      const size_t p_nbOctaves,
      const size_t p_nbScales = 0,
      const size_t p_width = 0,
      const size_t p_height = 0,
      const float p_delta = 0,
      const float p_sigma = 0);


    /**
     * @brief Getters.
     **/
    size_t getNbOctaves() const {return m_nbOctaves;}
    Octave* getOctave(const size_t p_n) const {return m_octaves[p_n];}


  private:

    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


  //! Data members
  private:
    size_t m_nbOctaves;  // number of octaves
    Octave** m_octaves; // octaves
};
#else
class Octave;

#endif // SCALESPACE_H_INCLUDED
