#ifndef OCTAVE_H_INCLUDED
#define OCTAVE_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief An octave is a stack of images sharing same w, h and intersample
          distance. Each image of the stack has a level of blur.
 **/
class Octave {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    Octave();


    /**
     * @brief Copy constructor.
     **/
    Octave(
      const Octave& i_octave);


    /**
     * @brief Surcharged constructor.
     **/
    Octave(
      const float i_delta,
      const size_t i_width,
      const size_t i_height,
      const size_t i_nbImages,
      const float* i_sigmas);


    /**
     * @brief Operator overload.
     **/
    Octave& operator=(
      const Octave& i_octave);


    /**
     * @brief Default destructor.
     **/
    ~Octave();


    /**
     * @brief Getters.
     **/
    float getDelta    () const {return m_delta   ;}
    size_t getWidth   () const {return m_width   ;}
    size_t getHeight  () const {return m_height  ;}
    size_t getNbImages() const {return m_nbImages;}
    float getSigma(const size_t p_n) const {return m_sigmas[p_n];}
    float* getPtr(
      const int p_c = 0,
      const int p_i = 0) const {return m_image->getPtr(p_c, p_i);}
    Image* getImage() const {return m_image;}


    /**
     * @brief Extrema interpolation via a quadratic function.
     **/
    void inverse3dTaylorSecondOrderExpansion(
      const int i,
      const int j,
      const int s,
      float& o_di,
      float& o_dj,
      float& o_ds,
      float& o_val);


  private:


    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


  //! Data members
  private:
    float m_delta;    // sampling rate in this octave
    size_t m_width;   // the same for all images in the stack
    size_t m_height;  // the same for all images in the stack
    size_t m_nbImages;// number of images in the stack
    float* m_sigmas;  // intrinsic levels of blur
    Image* m_image; // Stack of images indexed from fine to coarse.
};
#else
class Octave;

#endif // OCTAVE_H_INCLUDED
