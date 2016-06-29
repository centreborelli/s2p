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
      const size_t i_channels,
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
     * @brief Getter
     **/
    size_t width    () const {return m_width    ;}
    size_t height   () const {return m_height   ;}
    size_t channels () const {return m_channels ;}
    size_t border   () const {return m_border   ;}
    size_t nHeight  (const size_t p_n) const {return m_heights[p_n];}


    /**
     * @brief Initialize an image full of zeros.
     **/
    void init(
      const size_t i_width,
      const size_t i_height,
      const size_t i_channels,
      const size_t i_border = 0) {
      this->releaseMemory();
      new (this) Image(i_width, i_height, i_channels, i_border);}

    void init(
      const Image& i_im) {
      this->init(i_im.m_width, i_im.m_height, i_im.m_channels, i_im.m_border);}


    //! Default destructor
    ~Image();


    /**
     * @brief Generic read of an image. Call readPng or readTiff according to
     *        the extension of the input name of the image to read.
     **/
    void read(
      const char* p_name,
      const size_t i_border = 0);


    /**
     * @brief Generic write of an image. Call writePng or writeTiff according to
     *        the extension of the input name of the image to write.
     *
     * @param p_name : path to the image to save;
     * @param p_quad: if true, apply a zoom-in by duplication of 1 pixel into 4.
     *
     * @return none.
     **/
     void write(
      const char* p_name,
      const bool p_quad = false) const;


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
     * @brief Add a border to the current image.
     *
     * @param i_border: size of the border to add.
     *
     * @return none.
     **/
    void addBorder(
      const size_t i_border);


    /**
     * @brief Set a value for the border of the current image.
     *
     * @param p_value : value for the initialization.
     *
     * @return none.
     **/
    void setBorder(
      const float p_value);


    /**
     * @brief Convolve the image with a Gaussian of width sigma and store result
     *        back in the image. This routine creates the Gaussian kernel, and
     *        then applies it sequentially in both horizontal and vertical
     *        directions.
     **/
    void convolveGaussian(
      const float p_sigma);


  private:
    /**
     * @brief Copy the inner image of the current image to another image.
     *
     * @param o_im : will recieve the inner image of this, without modifying its
     *               borders.
     *
     * @return none.
     **/
    void copyInner(
      Image& o_im,
      const int p_tid =-1) const;


    /**
     * @brief Read an image via the LibGDAL library. Will exit the main program
     *        in case of problem.
     *
     * @param i_name : path to the image which will be filled into p_ptr;
     * @param i_border : size of the border to add around the image (will be
     *                   initialized full of zeros.
     *
     * @return none.
     **/
    void readGDAL(
      const char* p_name,
      const size_t i_border = 0);


    /**
     * @brief Write an image via the LibGdal library. Will exit the main problem
     *        in case of problem.
     *
     * @param i_name: path to the image to save;
     * @param p_quad: if true, apply a zoom-in by duplication of 1 pixel into 4.
     **/
    void writeGDAL(
      const char* p_name,
      const bool p_quad = false,
      const char* pszFormat = "GTiff") const;


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
