/**
 * @brief Class to deal with images, and do conversion from/to cv::Mat.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

//! Global includes
#include <fftw3.h>
#include <tiffio.h>
#include <iostream>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <xmmintrin.h>
#include <x86intrin.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>

//! Local includes
#include "LibImages.h"
#include "../Utilities/Utilities.h"
#include "../Utilities/Memory.h"


using namespace std;


//! Default constructor
Image::Image() :

  //! Size
  m_width   (0),
  m_height  (0),
  m_channels(0),
  m_border  (0),
  m_size    (0),

  //! Miscalleneous
#ifdef _OPENMP
  m_nbThreads(omp_get_max_threads()),
#else
  m_nbThreads(1),
#endif // _OPENMP

  //! Pointers
  m_ptr     (NULL),
  m_heights (NULL) {

}


//! Surcharged constructor with a specified size.
Image::Image(
  const size_t i_width,
  const size_t i_height,
  const size_t i_channels,
  const size_t i_border) :

  //! Size
  m_width   (i_width),
  m_height  (i_height),
  m_channels(i_channels),
  m_border  (i_border),
  m_size    (i_channels * (i_width + 2 * i_border) * (i_height + 2 * i_border)),

  //! Miscalleneous
#ifdef _OPENMP
  m_nbThreads(omp_get_max_threads()),
#else
  m_nbThreads(1),
#endif // _OPENM

  //! Pointers
  m_ptr     ((float*) memalloc(16, m_size * sizeof(float))),
  m_heights ((size_t*) malloc((m_nbThreads + 1) * sizeof(size_t))) {

  //! Fill the image with zeros
  size_t k = 0;

#ifdef USE_AVX
  //! AVX version
  for (; k < m_size - 8; k += 8) {
    _mm256_storeu_ps(m_ptr + k, _mm256_setzero_ps());
  }
#else
#ifdef USE_SSE
  //! SSE version
  for (; k < m_size - 4; k += 4) {
    _mm_storeu_ps(m_ptr + k, _mm_setzero_ps());
  }
#endif // USE_SSE
#endif // USE_AVX

  //! Normal version
  for (; k < m_size; k++) {
    m_ptr[k] = float(0);
  }

  //! Initialize the beginning and end for each slices of the image
  initializeHeights(m_height, m_heights, m_nbThreads);
}


//! Copy constructor.
Image::Image(
  const Image& i_im) :

  //! Size
  m_width   (i_im.m_width),
  m_height  (i_im.m_height),
  m_channels(i_im.m_channels),
  m_border  (i_im.m_border),
  m_size    (i_im.m_size),

  //! Miscalleneous
  m_nbThreads(i_im.m_nbThreads),

  //! Pointers
  m_ptr     ((float*) memalloc(16, m_size * sizeof(float))),
  m_heights ((size_t*) malloc((m_nbThreads + 1) * sizeof(size_t))) {

  //! Copy the data
  size_t k = 0;
#ifdef USE_AVX
  //! AVX version
  for (; k < m_size - 8; k += 8) {
    _mm256_storeu_ps(m_ptr + k, _mm256_loadu_ps(i_im.m_ptr + k));
  }
#else
#ifdef USE_SSE
  //! SSE version
  for (; k < m_size - 4; k += 4) {
    _mm_storeu_ps(m_ptr + k, _mm_loadu_ps(i_im.m_ptr + k));
  }
#endif // USE_SSE
#endif // USE_AVX

  //! Normal version
  for (; k < m_size; k++) {
    m_ptr[k] = i_im.m_ptr[k];
  }

  //! Copy the slices
  for (size_t k = 0; k < m_nbThreads + 1; k++) {
    m_heights[k] = i_im.m_heights[k];
  }
}


//! Default destructor
Image::~Image() {

  //! Release the memory
  releaseMemory();
}


//! Operator overload.
Image& Image::operator=(
  const Image& i_im) {

  if (&i_im == this) {
    return *this;
  }

  releaseMemory();

  new (this) Image(i_im);
  return *this;
}


//! Operator overload.
Image& Image::operator=(
  const float p_value) {

  for (size_t k = 0; k < m_size; k++) {
    m_ptr[k] = p_value;
  }

  return *this;
}


void Image::operator*=(
  const float p_value) {

  for (size_t k = 0; k < m_size; k++) {
    m_ptr[k] *= p_value;
  }
}


//! Generic read of an image.
void Image::read(
  const char* p_name,
  const size_t i_border) {

  //! Parameters check
  if (NULL == p_name) {
    cout << "Null name" << endl;
    exit(EXIT_FAILURE);
  }
  //! Check the extension
  string ext = p_name;
  const int pos = ext.find_last_of(".");
  ext = ext.substr(pos + 1, ext.size());

  //! Call the right function
  if (ext == "png") {
    this->readPng(p_name, i_border);
  }
  else if (ext == "tif" || ext == "tiff") {
    this->readTiff(p_name, i_border);
  }
  else {
    cout << "Extension " << ext << " not known. Abort." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Read png-only image.
void Image::readPng(
  const char* p_name,
  const size_t i_border) {

  //! Initialization
  const size_t png_sig_len = 4;
  int transform = PNG_TRANSFORM_IDENTITY;

  //! Open the PNG input file
  FILE *volatile fp = NULL;
  if (0 == strcmp(p_name, "-")) {
    fp = stdin;
  }
  else if (NULL == (fp = fopen(p_name, "rb"))) {
    cout << "Can't open the file " << p_name << "." << endl;
    exit(EXIT_FAILURE);
  }

  //! Read in some of the signature bytes and check this signature
  png_byte png_sig[4];
  if ((png_sig_len != fread(png_sig, 1, png_sig_len, fp))
    || 0 != png_sig_cmp(png_sig, (png_size_t) 0, png_sig_len)) {
    cout << "Bad signature." << endl;
    (void) fclose(fp);
    exit(EXIT_FAILURE);
  }

  //! Create and initialize the png_struct with the default stderr and error
  //! handling
  png_structp png_ptr;
  if (NULL == (png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL))) {
    cout << "Can't initialize the png_struct." << endl;
    (void) fclose(fp);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    exit(EXIT_FAILURE);
  }

  //! Allocate/initialize the memory for image information
  png_infop info_ptr;
  if (NULL == (info_ptr = png_create_info_struct(png_ptr))) {
    cout << "Can't allocate the memory for image information." << endl;
    (void) fclose(fp);
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    exit(EXIT_FAILURE);
  }

  //! Set error handling
  if (0 != setjmp(png_jmpbuf(png_ptr))) {
    cout << "Can't set the error handling." << endl;
    (void) fclose(fp);
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    exit(EXIT_FAILURE);
  }

  //! Set up the input control using standard C streams
  png_init_io(png_ptr, fp);

  //! let libpng know that some bytes have been read
  png_set_sig_bytes(png_ptr, png_sig_len);

  //! Set the read filter transforms, to get 8bit RGB whatever the original
  //! file may contain:
  //! PNG_TRANSFORM_STRIP_16      strip 16-bit samples to 8 bits
  //! PNG_TRANSFORM_PACKING       expand 1, 2 and 4-bit
  //!                             samples to bytes
  transform |= (PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING);

  //! Read in the entire image at once
  png_read_png(png_ptr, info_ptr, transform, NULL);

  //! Get image informations
  m_width    = (size_t) png_get_image_width (png_ptr, info_ptr);
  m_height   = (size_t) png_get_image_height(png_ptr, info_ptr);
  m_channels = (size_t) png_get_channels    (png_ptr, info_ptr);
  png_bytepp row_pointers;
  row_pointers = png_get_rows(png_ptr, info_ptr);

  //! Allocate the image
  new (this) Image(m_width, m_height, m_channels, i_border);

  //! Deinterlace and convert png RGB RGB RGB 8bit to RRR GGG BBB
  //! the image is deinterlaced layer after layer
  //! this generic loop also works for one single channel
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < m_height; i++) {
      float* oI = this->getPtr(c, i);
      png_bytep row_ptr = row_pointers[i] + c;

      for (size_t j = 0; j < m_width; j++) {
        oI[j] = (float) *row_ptr;
        row_ptr += m_channels;
      }
    }
  }

  //! Test if image is really a color image and exclude the alpha channel
  if (m_channels == 2) {
    cout << "Warning: this image has only 2 channels. Only one will be kept\n";
    m_channels = 1;
  }
  else if (m_channels >= 3) {
    const size_t wh = m_width * m_height;
    size_t k = 0;
    const float* iR = this->getPtr(0, 0);
    const float* iG = this->getPtr(1, 0);
    const float* iB = this->getPtr(2, 0);
    while (k < wh && iR[k] == iG[k] && iR[k] == iB[k]) {
      k++;
    }
    m_channels = (k == wh ? 1 : 3);
  }

  //! clean up and memfree any memory allocated, close the file
  (void) fclose(fp);
  png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
}


//! Read a tiff-only image.
void Image::readTiff(
  const char* p_name,
  const size_t i_border) {

  //! Open the tiff file
  TIFF *tif = TIFFOpen(p_name, "r");

  //! Check that the file has been open
  if (!tif) {
    cout << "Unable to read TIFF file " << p_name << ". Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Initialization
  uint32 w = 0, h = 0;
  uint16 spp = 0, bps = 0, fmt = 0;

  //! Get the metadata
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
  TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &fmt);

  //! Check the metadata
  if (spp != 1 || bps != (uint16) sizeof(float) * 8 ||
      fmt != SAMPLEFORMAT_IEEEFP) {
    cout << "readTiff: metadata non-conform. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Allocate the image
  new (this) Image(w, h, 1, i_border);

  //! Read the values
  for (size_t i = 0; i < m_height; i++) {
    float* oI = this->getPtr(0, i);
    if (TIFFReadScanline(tif, oI, i, 0) < 0) {
      cout << "readTiff: error reading row " << i << endl;
      exit(EXIT_FAILURE);
    }
  }

  //! Close the file
  TIFFClose(tif);
}


//! Generic write image.
void Image::write(
  const char* p_name,
  const bool p_quad) const {

  //! parameters check
  if (0 >= m_width || 0 >= m_height || 0 >= m_channels) {
    cout << "Inconsistent size." << endl;
    exit(EXIT_FAILURE);
  }

  //! Parameters check
  if (NULL == p_name) {
    cout << "Null name" << endl;
    exit(EXIT_FAILURE);
  }
  //! Check the extension
  string ext = p_name;
  const int pos = ext.find_last_of(".");
  ext = ext.substr(pos + 1, ext.size());

  //! Call the right function
  if (ext == "png") {
    this->writePng(p_name, p_quad);
  }
  else if (ext == "tif" || ext == "tiff") {
    this->writeTiff(p_name, p_quad);
  }
  else {
    cout << "Extension " << ext << " not known. Abort." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Write an image via the libpng write function
void Image::writePng(
  const std::string &p_name,
  const bool p_quad) const {

  //! open the PNG output file
  FILE *volatile fp;
  if (0 == strcmp(p_name.c_str(), "-")) {
    fp = stdout;
  }
  else if (NULL == (fp = fopen(p_name.c_str(), "wb"))) {
    exit(EXIT_FAILURE);
  }

  //! For convenience
  const size_t w = p_quad ? 2 * m_width  : m_width ;
  const size_t h = p_quad ? 2 * m_height : m_height;

  //! allocate the interlaced array and row pointers
  const size_t whc = w* h * m_channels;
  png_byte *idata = NULL;
  if (NULL == (idata = (png_byte *) malloc(whc * sizeof(png_byte)))) {
    (void) fclose(fp);
    exit(EXIT_FAILURE);
  }

  png_bytep *row_pointers = NULL;
  if (NULL == (row_pointers = (png_bytep *) malloc(h * sizeof(png_bytep)))) {
    (void) fclose(fp);
    memfree(idata);
    exit(EXIT_FAILURE);
  }

  //! create and initialize the png_struct with the default stderr and error
  //! handling
  png_structp png_ptr;
  if (NULL == (png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL))) {
    (void) fclose(fp);
    memfree(idata);
    memfree(row_pointers);
    exit(EXIT_FAILURE);
  }

  //! allocate/initialize the memory for image information
  png_infop info_ptr;
  if (NULL == (info_ptr = png_create_info_struct(png_ptr))) {
    (void) fclose(fp);
    memfree(idata);
    memfree(row_pointers);
    png_destroy_write_struct(&png_ptr, NULL);
    exit(EXIT_FAILURE);
  }

  //! set error handling
  if (0 != setjmp(png_jmpbuf(png_ptr))) {
    (void) fclose(fp);
    memfree(idata);
    memfree(row_pointers);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    exit(EXIT_FAILURE);
  }

  //! set up the input control using standard C streams
  png_init_io(png_ptr, fp);

  //! set image informations
  const png_byte bit_depth = 8;
  int color_type;
  switch (m_channels) {
    case 1:
      color_type = PNG_COLOR_TYPE_GRAY;
      break;
    case 2:
      color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
      break;
    case 3:
      color_type = PNG_COLOR_TYPE_RGB;
      break;
    case 4:
      color_type = PNG_COLOR_TYPE_RGB_ALPHA;
      break;
    default:
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      memfree(row_pointers);
      memfree(idata);
      (void) fclose(fp);
      exit(EXIT_FAILURE);
  }
  const int interlace = PNG_INTERLACE_ADAM7;
  const int compression = PNG_COMPRESSION_TYPE_BASE;
  const int filter = PNG_FILTER_TYPE_BASE;

  //! set image header
  png_set_IHDR(png_ptr, info_ptr, (png_uint_32) w, (png_uint_32) h,
    bit_depth, color_type, interlace, compression, filter);
  png_write_info(png_ptr, info_ptr);

  //!interlace and convert RRR GGG BBB to RGB RGB RGB
  //! the image is interlaced layer after layer
  //! this involves more memory exchange, but allows a generic loop
  for (size_t c = 0; c < m_channels; c++) {
    png_byte *idata_ptr = idata + c;

    for (size_t i = 0; i < h; i++) {
      const float* iI = this->getPtr(c, p_quad ? i / 2 : i);

      for (size_t j = 0; j < w; j++) {
        const int tmp = (int) floor((float) iI[p_quad ? j / 2 : j] + 0.5f);
        *idata_ptr = (png_byte) std::min(255, std::max(0, tmp));
        idata_ptr += m_channels;
      }
    }
  }

  //! set row pointers
  for (size_t i = 0; i < h; i++) {
    row_pointers[i] = idata + (size_t) (m_channels * w * i);
  }

  //! write out the entire image and end it
  png_write_image(png_ptr, row_pointers);
  png_write_end(png_ptr, info_ptr);

  //! clean up and memfree any memory allocated, close the file
  (void) fclose(fp);
  memfree(idata);
  memfree(row_pointers);
  png_destroy_write_struct(&png_ptr, &info_ptr);
}


//! Write an image via the libtiff write function
void Image::writeTiff(
  const std::string &p_name,
  const bool p_quad) const {

  //! Open the file
  TIFF *tif = TIFFOpen(p_name.c_str(), "w");

  //! Check that the file has been created
  if (!tif) {
    cout << "Unable to write TIFF file " << p_name << endl;
    exit(EXIT_FAILURE);
  }

  //! For convenience
  const size_t h = m_height * (p_quad ? 2 : 1);
  const size_t w = m_width  * (p_quad ? 2 : 1);
  uint32 rowsperstrip;
  float* line = (float*) memalloc(16, (p_quad ? w : m_width) * sizeof(float));

  //! Header of the tiff file
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32) w);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32) h);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, (uint16) m_channels);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, (uint16) sizeof(float) * 8);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  rowsperstrip = TIFFDefaultStripSize(tif, (uint32) h);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

  //! Write the file data
  for (size_t c = 0, ok = 1; ok && c < m_channels; c++) {
    for (size_t i = 0; ok && i < h; i++) {
      const float* iI = this->getPtr(c, p_quad ? i / 2 : i);

      //! Copy the line
      if (p_quad) {
        for (size_t j = 0; j < m_width; j++) {
          line[2 * j + 0] = line[2 * j + 1] = iI[j];
        }
      }
      else {
        for (size_t j = 0; j < m_width; j++) {
          line[j] = iI[j];
        }
      }

      //! Write the line
      if (TIFFWriteScanline(tif, line, (uint32) i, (tsample_t) c) < 0) {
        cout << "WriteTIFF: error writing row " << i << endl;
        ok = 0;
      }
    }
  }

  //! Release memory
  memfree(line);

  //! Close the file
  TIFFClose(tif);
}

//! Add a border to the current image.
void Image::addBorder(
  const size_t i_border) {

  //! Temporary image with border
  Image tmp(m_width, m_height, m_channels, i_border);

  //! Copy the inner image into it
  this->copyInner(tmp);

  //! Retrieve the bigger image
  *this = tmp;
}


//! Set a value for the border of the current image.
void Image::setBorder(
  const float p_value) {

#ifdef USE_AVX
  const __m256 xValue = _mm256_set1_ps(p_value);
#else
#ifdef USE_SSE
  const __m128 xValue = _mm_set1_ps(p_value);
#endif // USE_SSE
#endif // USE_AVX

  //! For every channels
  for (size_t c = 0; c < m_channels; c++) {

    //! Top and bottom borders
    for (int n = 0; n < int(m_border); n++) {
      float* iT = this->getPtr(c, - n - 1);
      float* iB = this->getPtr(c, m_height + n);
      int j = -int(m_border);

#ifdef USE_AVX
      //! AVX version
      for (; j < int(m_width + m_border) - 8; j += 8) {
        _mm256_storeu_ps(iT + j, xValue);
        _mm256_storeu_ps(iB + j, xValue);
      }
#else
#ifdef USE_SSE
      //! SSE version
      for (; j < int(m_width + m_border) - 4; j += 4) {
        _mm_storeu_ps(iT + j, xValue);
        _mm_storeu_ps(iB + j, xValue);
      }
#endif // USE_SSE
#endif // USE_AVX

      //! Normal version
      for (; j < int(m_width + m_border); j++) {
        iT[j] = p_value;
        iB[j] = p_value;
      }
    }

    //! Left and right borders
    for (int i = -int(m_border); i < int(m_height + m_border); i++) {
      float* iT = this->getPtr(c, i);

      for (int n = 0; n < int(m_border); n++) {
        iT[-n - 1] = p_value;
        iT[m_width + n] = p_value;
      }
    }
  }
}


//! Copy the inner image of the current image to another image.
void Image::copyInner(
  Image& o_im,
  const int p_tid) const {

  //! Check the size of the output image.
  if (o_im.width   () != m_width   ||
      o_im.height  () != m_height  ||
      o_im.channels() != m_channels) {

    o_im.init(m_width, m_height, m_channels, o_im.border());
  }

  //! For convenience
  const size_t hBegin = p_tid >= 0 ? m_heights[p_tid] : 0;
  const size_t hEnd   = p_tid >= 0 ? m_heights[p_tid + 1] : m_height;

  //! Copy the values
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = hBegin; i < hEnd; i++) {
      const float* iI = this->getPtr(c, i);
      float* oI = o_im.getPtr(c, i);
      size_t j = 0;

#ifdef USE_AVX
      //! AVX version
      for (; j < m_width - 8; j += 8) {
        _mm256_storeu_ps(oI + j, _mm256_loadu_ps(iI + j));
      }
#else
#ifdef USE_SSE
      //! AVX version
      for (; j < m_width - 4; j += 4) {
        _mm_storeu_ps(oI + j, _mm_loadu_ps(iI + j));
      }
#endif // USE_SSE
#endif // USE_AVX

      //! Normal version
      for (; j < m_width; j++) {
        oI[j] = iI[j];
      }
    }
  }
}


//! In-place convolution with a Gaussian.
void Image::convolveGaussian(
  const float p_sigma) {

  //! For convenience
  const int h = m_height;
  const int w = m_width;

  //! Kernel initialization
  const int kSize = 4 * ceilf(p_sigma) + 1;
  const int kHalf = (kSize - 1) / 2;
  float *kernel = (float*) memalloc(16, kSize * sizeof(float));

  //! Special case
  if (kSize == 1) {
    kernel[0] = 1.f;
  }
  else {
    const int hSize = (kSize - 1) / 2;

    //! Build the kernel
    for (int k = hSize; k < kSize; k++) {
      const float val = float(k - hSize) / p_sigma;
      kernel[k] = kernel[kSize - 1 - k] = (float) exp(-0.5 * val * val);
    }

    float sum = 0.f;
    for (int k = 0; k < kSize; k++) {
      sum += kernel[k];
    }

    //! Normalization
    for (int k = 0; k < kSize; k++) {
      kernel[k] /= sum;
    }
  }

  //! Loop over the channels
  for (size_t c = 0; c < m_channels; c++) {

    //! Buffer allocation
    float* line = (float*) memalloc(16, (2 * kHalf + w) * sizeof(float));

    //! Horizontal convolution
    for (int i = 0; i < h; i++) {
      float* iI = this->getPtr(c, i);

      //! Copy the line into a buffer
      for (int j = 0; j < kHalf; j++) {
        //line[j] = iI[0];
        line[j] = iI[kHalf - 1 - j];
        //line[kHalf + w + j] = iI[w - 1];
        line[kHalf + w + j] = iI[w - 1 - j];
      }
      for (int j = 0; j < w; j++) {
        line[kHalf + j] = iI[j];
      }

      //! Apply the convolution - SSE version
      int j = 0;
      for (; j < w - 4; j += 4) {
        __m128 value = _mm_setzero_ps();

        //! Unroll the loops
        int k = 0;
        while (k + 4 < kSize) {
          value += _mm_set1_ps(kernel[k + 0]) * _mm_loadu_ps(line + j + k + 0) +
                   _mm_set1_ps(kernel[k + 1]) * _mm_loadu_ps(line + j + k + 1) +
                   _mm_set1_ps(kernel[k + 2]) * _mm_loadu_ps(line + j + k + 2) +
                   _mm_set1_ps(kernel[k + 3]) * _mm_loadu_ps(line + j + k + 3);
          k += 4;
        }
        for (; k < kSize; k++) {
          value += _mm_set1_ps(kernel[k]) * _mm_loadu_ps(line + j + k);
        }
        _mm_storeu_ps(iI + j, value);
      }

      //! Normal version
      for (; j < w; j++) {
        float value = 0.f;

        //! Apply the convolution
        for (int k = 0; k < kSize; k++) {
          value += kernel[k] * line[j + k];
        }

        iI[j] = value;
      }
    }

    //! Release memory
    memfree(line);

    //! Buffer allocation
    __m128* xCol = (__m128*) memalloc(16, (2 * kHalf + h) * sizeof(__m128));

    //! Vertical convolution - SSE version
    int j = 0;
    for (; j < w - 4; j += 4) {

      //! Copy the line into a buffer
      for (int i = 0; i < kHalf; i++) {
        //xCol[i] = _mm_loadu_ps(this->getPtr(c, 0) + j);
        xCol[i] = _mm_loadu_ps(this->getPtr(c, kHalf - 1 - i) + j);
        //xCol[kHalf + h + i] = _mm_loadu_ps(this->getPtr(c, h - 1) + j);
        xCol[kHalf + h + i] = _mm_loadu_ps(this->getPtr(c, h - 1 - i) + j);
      }
      for (int i = 0; i < h; i++) {
        xCol[kHalf + i] = _mm_loadu_ps(this->getPtr(c, i) + j);
      }

      //! Apply the convolution
      for (int i = 0; i < h; i++) {
        __m128 value = _mm_setzero_ps();

        //! Unroll the loop
        int k = 0;
        while (k + 4 < kSize) {
          value += _mm_set1_ps(kernel[k + 0]) * xCol[i + k + 0] +
                   _mm_set1_ps(kernel[k + 1]) * xCol[i + k + 1] +
                   _mm_set1_ps(kernel[k + 2]) * xCol[i + k + 2] +
                   _mm_set1_ps(kernel[k + 3]) * xCol[i + k + 3];
          k += 4;
        }
        for (; k < kSize; k++) {
          value += _mm_set1_ps(kernel[k]) * xCol[i + k];
        }

        xCol[i] = value;
      }

      //! Get the value
      for (int i = 0; i < h; i++) {
        _mm_storeu_ps(this->getPtr(c, i) + j, xCol[i]);
      }
    }

    //! Release memory
    memfree(xCol);

    //! Buffer allocation
    float* col = (float*) memalloc(16, (2 * kHalf + h) * sizeof(float));

    //! Vertical convolution - normal version
    for (; j < w; j++) {

      //! Copy the line into a buffer
      for (int i = 0; i < kHalf; i++) {
        //col[i] = this->getPtr(c, 0)[j];
        col[i] = this->getPtr(c, kHalf - 1 - i)[j];
        //col[kHalf + h + i] = this->getPtr(c, h - 1)[j];
        col[kHalf + h + i] = this->getPtr(c, h - 1 - i)[j];
      }
      for (int i = 0; i < h; i++) {
        col[kHalf + i] = this->getPtr(c, i)[j];
      }

      //! Apply the convolution
      for (int i = 0; i < h; i++) {
        float value = 0.f;

        //! Unroll the loop
        int k = 0;
        while (k + 4 < kSize) {
          value += kernel[k] * col[i + k] + kernel[k + 1] * col[i + k + 1] +
            kernel[k + 2] * col[i + k + 2] + kernel[k + 3] * col[i + k + 3];
          k += 4;
        }
        for (; k < kSize; k++) {
          value += kernel[k] * col[i + k];
        }

        //! Get the value
        this->getPtr(c, i)[j] = value;
      }
    }

    //! Release memory
    memfree(col);
  }

  //! Release memory
  memfree(kernel);
}


//! Set to 0 all NaN and Inf values.
void Image::removeNaN() {

  //! Loop over the pixels
  for (size_t k = 0; k < m_size; k++) {
    m_ptr[k] = isfinite(m_ptr[k]) ? m_ptr[k] : 0.f;
  }
}


//! Apply an in-place sub-pixelian translation horizontaly.
void Image::translate(
  const float p_delta) {

  //! Add symmetrical boundaries
  Image im(2 * m_width, m_height, m_channels);
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < m_height; i++) {
      const float* iI = this->getPtr(c, i);
      float* oI = im.getPtr(c, i);

      for (size_t j = 0; j < m_width; j++) {
        oI[j] = oI[2 * m_width - 1 - j] = iI[j];
      }
    }
  }

  //! Apply the translation
  for (size_t c = 0; c < m_channels; c++) {
    im.horizontalShear(p_delta, c);
  }

  //! Get rid off the boundaries
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < m_height; i++) {
      const float* iI = im.getPtr(c, i);
      float* oI = this->getPtr(c, i);

      for (size_t j = 0; j < m_width; j++) {
        oI[j] = iI[j];
      }
    }
  }
}


//! Apply an in-place horizontal shear for the current channel.
void Image::horizontalShear(
  const float p_delta,
  const size_t p_currentChnl) {

  //! For convenience
  const size_t h   = m_height;
  const size_t w   = m_width;
  const size_t w2  = w >> 1;
  const size_t odd = (w & 1) != 0;

  //! Compute the real amount of translation
  const double translation = -(p_delta) * 2.0 * M_PI;

  //! Allocation for the threads
  fftwf_complex** avec = (fftwf_complex**) memalloc(16, m_nbThreads * sizeof(fftwf_complex*));
  fftwf_complex** afft = (fftwf_complex**) memalloc(16, m_nbThreads * sizeof(fftwf_complex*));
  fftwf_plan* pfwd = (fftwf_plan*) memalloc(16, m_nbThreads * sizeof(fftwf_plan));
  fftwf_plan* pbwd = (fftwf_plan*) memalloc(16, m_nbThreads * sizeof(fftwf_plan));

  //! Create temporary array to compute the fourier transform of a line
  //! Create FFTW forward and backward plans
  for (size_t n = 0; n < m_nbThreads; n++) {
    avec[n] = (fftwf_complex*) fftwf_alloc_complex(w * sizeof(fftwf_complex));
    afft[n] = (fftwf_complex*) fftwf_alloc_complex(w * sizeof(fftwf_complex));
    pfwd[n] = fftwf_plan_dft_1d(w, avec[n], afft[n], FFTW_BACKWARD, FFTW_ESTIMATE);
    pbwd[n] = fftwf_plan_dft_1d(w, afft[n], avec[n], FFTW_FORWARD , FFTW_ESTIMATE);
  }


#ifdef _OPENMP
#pragma omp parallel
{
#endif // _OPENMP

  //! Loop over the lines
#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif // _OPENMP
  for (size_t i = 0; i < h; i++) {

#ifdef _OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // _OPENMP
    fftwf_complex* vec = avec[tid];
    fftwf_complex* fft = afft[tid];

    float* iI = this->getPtr(p_currentChnl, i);

    //! Copy the input real data into working complex data
    for (size_t j = 0; j < w; j++) {
      vec[j][0] = iI[j];
      vec[j][1] = 0.f;
    }

    //! Compute the FFT of the current line of the image
    fftwf_execute(pfwd[tid]);

    //! Apply the translation
    for (size_t j = 1; j < w2 + odd; j++) {
      const float fcos = (float) cos(j * translation / w);
      const float fsin = (float) sin(j * translation / w);
      const float re = fft[j][0];
      const float im = fft[j][1];

      fft[    j][0] = fcos * re - fsin * im;
      fft[    j][1] = fsin * re + fcos * im;
      fft[w - j][0] = fcos * re - fsin * im;
      fft[w - j][1] =-fsin * re - fcos * im;;
    }

    if (odd == 0) {
      fft[w2][0] = cos(translation * 0.5) * fft[w2][0] - sin(translation * 0.5) * fft[w2][1];
      fft[w2][1] = 0.f;
    }

    //! Compute the inverse FFT of the current line
    fftwf_execute(pbwd[tid]);

    //! Get the normalized result
    for (size_t j = 0; j < w; j++) {
      iI[j] = vec[j][0] / float(w);
    }
  }
#ifdef _OPENMP
}
#endif // _OPENMP

  //! Release memory
  for (size_t n = 0; n < m_nbThreads; n++) {
    fftwf_free(avec[n]);
    fftwf_free(afft[n]);
    fftwf_destroy_plan(pfwd[n]);
    fftwf_destroy_plan(pbwd[n]);
  }
  memfree(pfwd);
  memfree(pbwd);
  memfree(avec);
  memfree(afft);

  //! FFTW cleanup
  fftwf_cleanup();
}


//! Sub-sample an image by a factor 2 after applying a Gaussian convolution.
void Image::subSample(
  Image& o_im) const {

  //! First, apply a Gaussian convolution
  Image imG = *this;
  imG.convolveGaussian(0.8f);

  //! Size initialization
  const size_t ws = (size_t) floor(float(m_width ) / 2.f);
  const size_t hs = (size_t) floor(float(m_height) / 2.f);
  o_im.init(ws, hs, m_channels, m_border);

  //! Then, sub-sample
  for (size_t c = 0; c < m_channels; c++) {
    for (size_t i = 0; i < hs; i++) {
      const float* iI0 = imG.getPtr(c, 2 * i);
      const float* iI1 = imG.getPtr(c, 2 * i + 1);
      float* oI = o_im.getPtr(c, i);

      for (size_t j = 0; j < ws; j++) {
        if (2 * i + 2 < m_height && 2 * j + 2 < m_width) {
          oI[j] = (iI0[2 * j] + iI0[2 * j + 1] + iI1[2 * j] + iI1[2 * j + 1]) / 4.f;
        }
        else {
          oI[j] = iI0[2 * j];
        }
      }
    }
  }
}


//! Up-sample an image by a factor 2.
void Image::upSample(
  Image& o_im) const {

  //! Size initialization
  const size_t wo = (size_t) rintf(2.f * float(m_width ));
  const size_t ho = (size_t) rintf(2.f * float(m_height));
  o_im.init(wo, ho, m_channels, m_border);

  //! Up sample the image by closest neighbour
  for (size_t c = 0; c < m_channels; c++) {
    for(size_t i = 0; i < ho; i++) {
      const float* iI = this->getPtr(c, int(floor(i * 0.5)));
      float* oI = o_im.getPtr(c, i);

      for(size_t j = 0; j < wo; j++) {
        oI[j] = iI[int(floor(j * 0.5))];
      }
    }
  }
}


//! Normalize (L1 norm) an image.
void Image::normalizeL1() {

  //! Initialization
  float sum = 0.f;

  //! Compute the sum of the pixels for all channels
  for (size_t k = 0; k < m_size; k++) {
    sum += m_ptr[k];
  }

  //! Normalize it
  if (sum != 0.f) {
    sum = 1.f / sum;

    for (size_t k = 0; k < m_size; k++) {
      m_ptr[k] *= sum;
    }
  }
}


//! Release the memory
void Image::releaseMemory() {

  //! Release the memory
  if (m_ptr != NULL) {
    memfree(m_ptr);
    m_ptr = NULL;
  }

  if (m_heights != NULL) {
    memfree(m_heights);
    m_heights = NULL;
  }
}












