#ifndef LIBIMAGES_H_INCLUDED
#define LIBIMAGES_H_INCLUDED

//! Global includes
#include <stdlib.h>
#include <string>
#include <vector>
#include <xmmintrin.h>
#include <x86intrin.h>

//! Local includes


/**
 * @brief Small Image class, which contains basic manipulations.
 *
 **/
class Image {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    Image();


    /**
     * @brief Surcharged constructor with a specified size.
     *
     * @param i_width, i_height, i_channels : size of the wanted image;
     * @param i_border: size of the border to add around the image.
     **/
    Image(
      const size_t i_width,
      const size_t i_height,
      const size_t i_channels = 1,
      const size_t i_border = 0);


    /**
     * @brief Copy constructor.
     **/
    Image(
      const Image& i_im);


    Image(const float *i_ptr, const size_t i_width, const size_t i_height, const size_t i_channels);

    /**
     * @brief Operator overload.
     **/
    Image& operator=(
      const Image& i_im);


    /**
     * @brief Operator overload.
     **/
    Image& operator=(
      const float p_value);
    void operator*=(
      const float p_value);

    /**
     * @brief Getter
     **/
    size_t width    () const {return m_width    ;}
    size_t height   () const {return m_height   ;}
    size_t size     () const {return m_size     ;}
    size_t channels () const {return m_channels ;}
    size_t border   () const {return m_border   ;}
    size_t nHeight  (const size_t p_n) const {return m_heights[p_n];}


    /**
     * @brief Initialize an image full of zeros.
     **/
    void init(
      const size_t i_width,
      const size_t i_height,
      const size_t i_channels = 1,
      const size_t i_border = 0) {
      this->releaseMemory();
      new (this) Image(i_width, i_height, i_channels, i_border);}

    void init(
      const Image& i_im) {
      this->init(i_im.m_width, i_im.m_height, i_im.m_channels, i_im.m_border);}


    //! Default destructor
    ~Image();


    /**
     * @brief Get the pointer to the channel c, row i.
     *
     * @param p_c : current channel to adress,
     * @param p_i : current row to adress.
     *
     * @return the pointer of I(c, i, 0);
     **/
    float* getPtr(
      const int p_c = 0,
      const int p_i = 0) const {
      return m_ptr + (p_c * int(m_height + 2 * m_border) + p_i + (int) m_border)
                          * int(m_width  + 2 * m_border) + (int) m_border;
    }


    /**
     * @brief Convert an image from rgb to grayscale.
     **/
    void rgb2gray(
      Image& o_im) const;


    /**
     * @brief Interpolate the image with a bilinear model.
     *
     * @param o_im: output image of size (width / delta, height / delta).
     * @param p_delta: inter-pixel distance in the output image;
     * @param p_ic: working channel for the current image;
     * @param p_oc: working channel for the output image.
     **/
    void overSample(
      Image& o_im,
      const float p_delta,
      const size_t p_ic,
      const size_t p_oc) const;


    /**
     * @brief Image subsampling by a factor 2.
     *
     * @param o_im: output image (size must be allocated).
     * @param p_ic: working channel for the current image;
     * @param p_oc: working channel for the output image.
     **/
    void subSample(
      Image& o_im,
      const size_t p_ic,
      const size_t p_oc) const;


    /**
     * @brief Apply a Gaussian blur over an image.
     *
     * @param o_im: output image (size must be allocated);
     * @param p_ic: working channel for the current image;
     * @param p_oc: working channel for the output image.
     **/
    void applyGaussianBlur(
      Image& o_im,
      const float p_sigma,
      const size_t p_ic,
      const size_t p_oc) const;


    /**
     * @brief Compute image gradient via symmetric finite difference schemes.
     *        Image extension :  symmetrization at x=-1/2.
     *
     * @param o_gradX: output horizontal gradient (size must be allocated);
     * @param o_gradY: output vertical gradient (size must be allocated);
     * @param p_ic: working channel for the current image;
     * @param p_oc: working channel for the output images.
     **/
    void computeGradient(
      Image& o_gradX,
      Image& o_gradY,
      const size_t p_ic,
      const size_t p_oc) const;

  private:
    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


  //! Data members
  private:

    //! Size
    size_t m_width;
    size_t m_height;
    size_t m_channels;
    size_t m_border;
    size_t m_size;

    //! Miscalleneous
    size_t m_nbThreads;

    //! Pointers
    float* m_ptr;
    size_t* m_heights;
};
#else
class Image;

#endif // LIBIMAGES_H_INCLUDED
