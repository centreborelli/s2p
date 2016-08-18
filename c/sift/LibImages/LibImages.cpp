/**
 * @brief Class to deal with images, and do conversion from/to cv::Mat.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

//! Global includes
#include <tiffio.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <xmmintrin.h>
#include <x86intrin.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>

//! Local includes
#include "LibImages.h"
#include "../LibSSE/LibSSE.h"
#include "../Utilities/Memory.h"
#include "../Utilities/Utilities.h"


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
  m_ptr     ((float*) memalloc(32, m_size * sizeof(float))),
  m_heights ((size_t*) malloc((m_nbThreads + 1) * sizeof(size_t))) {

  //! Fill the image with zeros
  size_t k = 0;

#ifdef __AVX__
  //! AVX version
  for (; k < m_size - 8; k += 8) {
    _mm256_storeu_ps(m_ptr + k, _mm256_setzero_ps());
  }
#endif // __AVX__
#ifdef __SSE__
  //! SSE version
  for (; k < m_size - 4; k += 4) {
    _mm_storeu_ps(m_ptr + k, _mm_setzero_ps());
  }
#endif // __SSE__

  //! Normal version
  for (; k < m_size; k++) {
    m_ptr[k] = float(0);
  }

  //! Initialize the beginning and end for each slices of the image
  initializeHeights(m_height, m_heights, m_nbThreads);
}


//! Surcharged constructor with a specific size and pointer to data.
Image::Image(
    const float* i_ptr,
    const size_t i_width,
    const size_t i_height,
    const size_t i_channels):

    //! Size
    m_width   (i_width),
    m_height  (i_height),
    m_channels(i_channels),
    m_border  (0),
    m_size    (i_channels * i_width * i_height),

    //! Miscalleneous
#ifdef _OPENMP
    m_nbThreads(omp_get_max_threads()),
#else
    m_nbThreads(1),
#endif // _OPENMP

    //! Pointers
    m_ptr     ((float*) memalloc(16, m_size * sizeof(float))),
    m_heights ((size_t*) malloc((m_nbThreads + 1) * sizeof(size_t))) {

    //! Copy the data
    size_t k = 0;
#ifdef USE_AVX
    //! AVX version
    for (; k < m_size - 8; k += 8) {
        _mm256_storeu_ps(m_ptr + k, _mm256_loadu_ps(i_ptr + k));
    }
#else
#ifdef USE_SSE
    //! SSE version
    for (; k < m_size - 4; k += 4) {
        _mm_storeu_ps(m_ptr + k, _mm_loadu_ps(i_ptr + k));
    }
#endif // USE_SSE
#endif // USE_AVX

    //! Normal version
    for (; k < m_size; k++)
        m_ptr[k] = i_ptr[k];

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
#ifdef __AVX__
  //! AVX version
  for (; k < m_size - 8; k += 8) {
    _mm256_storeu_ps(m_ptr + k, _mm256_loadu_ps(i_im.m_ptr + k));
  }
#endif // __AVX__
#ifdef __SSE__
  //! SSE version
  for (; k < m_size - 4; k += 4) {
    _mm_storeu_ps(m_ptr + k, _mm_loadu_ps(i_im.m_ptr + k));
  }
#endif // __SSE__

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


//! Convert an image from rgb to grayscale.
void Image::rgb2gray(
  Image& o_im) const {

  //! Check the number of channels
  if (m_channels == 1) {
    o_im = *this;
    return;
  }
  else if (m_channels != 3) {
    cout << "rgb2gray: the input image must have 1 or 3 channels. ";
    cout << m_channels << " are provided. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Size initialization
  o_im.init(m_width, m_height, 1, m_border);

  //! Convert to gray
  for (size_t i = 0; i < m_height; i++) {

    //! Get the pointer to the three color channels
    float* iR = this->getPtr(0, i);
    float* iG = this->getPtr(1, i);
    float* iB = this->getPtr(2, i);
    float* oI = o_im.getPtr(0, i);

    for (size_t j = 0; j < m_width; j++) {
      oI[j] = 0.212639005871510 * iR[j] +
              0.715168678767756 * iG[j] +
              0.072192315360734 * iB[j];
    }
  }
}


//! Interpolate the image with a bilinear model.
void Image::overSample(
  Image& o_im,
  const float p_delta,
  const size_t p_ic,
  const size_t p_oc) const {

  //! For convenience
  const size_t hi = m_height;
  const size_t wi = m_width ;
  const size_t ho = (size_t) floor(float(m_height) / p_delta);
  const size_t wo = (size_t) floor(float(m_width ) / p_delta);

  //! Check size
  if (o_im.width() != wo || o_im.height() != ho) {
    cout << "Images::overSample: output size not consistant." << endl;
    exit(EXIT_FAILURE);
  }

  //! Loop over the channels
  for(size_t i = 0; i < ho; i++) {

    //! New coordinates
    const float x = i * p_delta;
    const size_t dx = (size_t)x;

    //! Image extension by symmetrization
    const size_t im = dx >= hi ? 2 * hi - 1 - dx : dx;
    const size_t ip = dx + 1 >= hi ? 2 * hi - 1 - dx - 1 : dx + 1;

    //! Pointers to the corresponding rows
    const float* pI = this->getPtr(p_ic, ip);
    const float* mI = this->getPtr(p_ic, im);
    float* oI = o_im.getPtr(p_oc, i);

    //! Loop over the columns
    for(size_t j = 0; j < wo; j++) {

      //! New coordinates
      const float y = j * p_delta;
      const size_t dy = (size_t)y;

      //! Image extension by symmetrization
      const size_t jm = dy >= wi ? 2 * wi - 1 - dy : dy;
      const size_t jp = dy + 1 >= wi ? 2 * wi - 1 - dy - 1 : dy + 1;

      //! Sampling coefficients
      const float fx = x - floor(x);
      const float fy = y - floor(y);

      //! Get the result
      oI[j] =           fx  * (fy * pI[jp] + (1.f - fy) * pI[jm])
               + (1.f - fx) * (fy * mI[jp] + (1.f - fy) * mI[jm]);
    }
  }
}


//! Image subsampling by a factor 2.
void Image::subSample(
  Image& o_im,
  const size_t p_ic,
  const size_t p_oc) const {

  //! For convenience
  const size_t ho = m_height / 2;
  const size_t wo = m_width  / 2;

  //! Check the size
  if (o_im.width() != wo || o_im.height() != ho) {
    cout << "subSample: sizes aren't consistant. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Sub sample
  for (size_t i = 0; i < ho; i++) {
    const float* iI = this->getPtr(p_ic, 2 * i);
    float* oI = o_im.getPtr(p_oc, i);

    for (size_t j = 0; j < wo; j++) {
      oI[j] = iI[2 * j];
    }
  }
}


//! Apply a Gaussian blur over an image.
void Image::applyGaussianBlur(
  Image& o_im,
  const float p_sigma,
  const size_t p_ic,
  const size_t p_oc) const {

  //! For convenience
  const int h = m_height;
  const int w = m_width;

  //! Check size
  if (o_im.width() != m_width || o_im.height() != m_height) {
    cout << "applyGaussianBlur: sizes aren't consistant. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Kernel initialization
  const int kHalf = (int) ceil(4 * p_sigma);
  float* ker = (float*) memalloc(16, (kHalf + 1) * sizeof(float));
  ker[0] = 1.f;
  if (p_sigma > 0) {
    float sum = ker[0];
    for (int k = 1; k <= kHalf; k++) {
      ker[k] = exp(-0.5*(float)k*(float)k/p_sigma/p_sigma);
      sum += 2 * ker[k];
    }
    for (int k = 0; k <= kHalf; k++) {
      ker[k] /= sum;
    }
  }
  else {
    for(int k = 1; k <= kHalf; k++){
      ker[k] = 0.f;
    }
  }

  //! Vertical convolution
  int j = 0;

#ifdef __AVX__
  //! AVX version
  //! Column buffer
  __m256* zcol = (__m256*) memalloc(32, (h + 2 * kHalf) * sizeof(__m256));
  for (; j < w - 8; j += 8) {

    //! Fill the buffer with the current column - first out-of-domain pixels
    for (int k = 0; k < std::min(kHalf, h); k++) {
      zcol[kHalf - k - 1] = _mm256_loadu_ps(this->getPtr(p_ic, k) + j);
      zcol[h + kHalf + k] = _mm256_loadu_ps(this->getPtr(p_ic, h - k - 1) + j);
    }
    //! Then in-domain pixels
    for (int i = 0; i < h; i++) {
      zcol[kHalf + i] = _mm256_loadu_ps(this->getPtr(p_ic, i) + j);
    }

    //! Compute the convolution
    for (int i = 0; i < h; i++) {
      __m256 zsum = _mm256_set1_ps(ker[0]) * zcol[i + kHalf];

      //! Unroll the loop
      int k = 1;
      for (; k <= kHalf - 4; k += 4) {
        zsum += _mm256_set1_ps(ker[k + 0]) * (zcol[i - k + kHalf - 0] +
                                              zcol[i + k + kHalf + 0]);
        zsum += _mm256_set1_ps(ker[k + 1]) * (zcol[i - k + kHalf - 1] +
                                              zcol[i + k + kHalf + 1]);
        zsum += _mm256_set1_ps(ker[k + 2]) * (zcol[i - k + kHalf - 2] +
                                              zcol[i + k + kHalf + 2]);
        zsum += _mm256_set1_ps(ker[k + 3]) * (zcol[i - k + kHalf - 3] +
                                              zcol[i + k + kHalf + 3]);
      }

      //! Finish the loop
      for (; k <= kHalf; k++) {
        zsum += _mm256_set1_ps(ker[k]) * (zcol[i - k + kHalf] +
                                          zcol[i + k + kHalf]);
      }

      //! Get the result
      _mm256_storeu_ps(o_im.getPtr(p_oc, i) + j, zsum);
    }
  }

  //! Release the column buffer
  memfree(zcol);
#endif // __AVX__
#ifdef __SSE__
  //! SSE version
  //! Column buffer
  __m128* xcol = (__m128*) memalloc(16, (h + 2 * kHalf) * sizeof(__m128));
  for (; j < w - 4; j += 4) {

    //! Fill the buffer with the current column - first out-of-domain pixels
    for (int k = 0; k < std::min(kHalf, h); k++) {
      xcol[kHalf - k - 1] = _mm_loadu_ps(this->getPtr(p_ic, k) + j);
      xcol[h + kHalf + k] = _mm_loadu_ps(this->getPtr(p_ic, h - k - 1) + j);
    }
    //! Then in-domain pixels
    for (int i = 0; i < h; i++) {
      xcol[kHalf + i] = _mm_loadu_ps(this->getPtr(p_ic, i) + j);
    }

    //! Compute the convolution
    for (int i = 0; i < h; i++) {
      __m128 xsum = _mm_set1_ps(ker[0]) * xcol[i + kHalf];

      //! Unroll the loop
      int k = 1;
      for (; k <= kHalf - 4; k += 4) {
        xsum += _mm_set1_ps(ker[k + 0]) * (xcol[i - k + kHalf - 0] +
                                           xcol[i + k + kHalf + 0]);
        xsum += _mm_set1_ps(ker[k + 1]) * (xcol[i - k + kHalf - 1] +
                                           xcol[i + k + kHalf + 1]);
        xsum += _mm_set1_ps(ker[k + 2]) * (xcol[i - k + kHalf - 2] +
                                           xcol[i + k + kHalf + 2]);
        xsum += _mm_set1_ps(ker[k + 3]) * (xcol[i - k + kHalf - 3] +
                                           xcol[i + k + kHalf + 3]);
      }

      //! Finish the loop
      for (; k <= kHalf; k++) {
        xsum += _mm_set1_ps(ker[k]) * (xcol[i - k + kHalf] +
                                       xcol[i + k + kHalf]);
      }

      //! Get the result
      _mm_storeu_ps(o_im.getPtr(p_oc, i) + j, xsum);
    }
  }

  //! Release the column buffer
  memfree(xcol);
#endif // __SSE__

  //! Normal version
  //! Column buffer
  float* col = (float*) memalloc(16, (h + 2 * kHalf) * sizeof(float));
  for (; j < w; j++) {

    //! Fill the buffer with the current column - first out-of-domain pixels
    for (int k = 0; k < std::min(kHalf, h); k++) {
      col[kHalf - k - 1] = this->getPtr(p_ic, k)[j];
      col[h + kHalf + k] = this->getPtr(p_ic, h - k - 1)[j];
    }
    //! Then in-domain pixels
    for (int i = 0; i < h; i++) {
      col[kHalf + i] = this->getPtr(p_ic, i)[j];
    }

    //! Compute the convolution
    for (int i = 0; i < h; i++) {
      float sum = ker[0] * col[i + kHalf];

      //! Unroll the loop
      int k = 1;
      for (; k <= kHalf - 4; k += 4) {
        sum += ker[k + 0] * (col[i - k + kHalf - 0] + col[i + k + kHalf + 0]);
        sum += ker[k + 1] * (col[i - k + kHalf - 1] + col[i + k + kHalf + 1]);
        sum += ker[k + 2] * (col[i - k + kHalf - 2] + col[i + k + kHalf + 2]);
        sum += ker[k + 3] * (col[i - k + kHalf - 3] + col[i + k + kHalf + 3]);
      }

      //! Finish the loop
      for (; k <= kHalf; k++) {
        sum += ker[k] * (col[i - k + kHalf] + col[i + k + kHalf]);
      }

      //! Get the result
      o_im.getPtr(p_oc, i)[j] = sum;
    }
  }

  //! Release the column buffer
  memfree(col);

  //! Line buffer
  float* line = (float*) memalloc(16, (w + 2 * kHalf) * sizeof(float));
  for (int i = 0; i < h; i++) {
    float* oI = o_im.getPtr(p_oc, i);

    //! Fill the buffer for the current line - first out-of-domain pixels
    for (int k = 0; k < std::min(kHalf, w); k++) {
      line[kHalf - k - 1] = oI[k];
      line[w + kHalf + k] = oI[w - k - 1];
    }
    //! Then in-domain pixels
    for (int j = 0; j < w; j++) {
      line[kHalf + j] = oI[j];
    }

    //! Compute the convolution
    int j = 0;
    const float* iL = line + kHalf;

#ifdef __AVX__
    //! AVX version
    for (; j < w - 8; j += 8) {
      __m256 xsum = _mm256_set1_ps(ker[0]) * _mm256_loadu_ps(iL + j);

      //! Unroll the loop
      int k = 1;
      for (; k <= kHalf - 4; k += 4) {
        xsum += _mm256_set1_ps(ker[k + 0]) * (_mm256_loadu_ps(iL + j - k - 0)
                                           + _mm256_loadu_ps(iL + j + k + 0));
        xsum += _mm256_set1_ps(ker[k + 1]) * (_mm256_loadu_ps(iL + j - k - 1)
                                           + _mm256_loadu_ps(iL + j + k + 1));
        xsum += _mm256_set1_ps(ker[k + 2]) * (_mm256_loadu_ps(iL + j - k - 2)
                                           + _mm256_loadu_ps(iL + j + k + 2));
        xsum += _mm256_set1_ps(ker[k + 3]) * (_mm256_loadu_ps(iL + j - k - 3)
                                           + _mm256_loadu_ps(iL + j + k + 3));
      }

      //! Finish the loop
      for (; k <= kHalf; k++) {
        xsum += _mm256_set1_ps(ker[k]) * (_mm256_loadu_ps(iL + j - k)
                                        + _mm256_loadu_ps(iL + j + k));
      }

      //! Get the result
      _mm256_storeu_ps(line + j, xsum);
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    for (; j < w - 4; j += 4) {
      __m128 xsum = _mm_set1_ps(ker[0]) * _mm_loadu_ps(iL + j);

      //! Unroll the loop
      int k = 1;
      for (; k <= kHalf - 4; k += 4) {
        xsum += _mm_set1_ps(ker[k + 0]) * (_mm_loadu_ps(iL + j - k - 0)
                                         + _mm_loadu_ps(iL + j + k + 0));
        xsum += _mm_set1_ps(ker[k + 1]) * (_mm_loadu_ps(iL + j - k - 1)
                                         + _mm_loadu_ps(iL + j + k + 1));
        xsum += _mm_set1_ps(ker[k + 2]) * (_mm_loadu_ps(iL + j - k - 2)
                                         + _mm_loadu_ps(iL + j + k + 2));
        xsum += _mm_set1_ps(ker[k + 3]) * (_mm_loadu_ps(iL + j - k - 3)
                                         + _mm_loadu_ps(iL + j + k + 3));
      }

      //! Finish the loop
      for (; k <= kHalf; k++) {
        xsum += _mm_set1_ps(ker[k]) * (_mm_loadu_ps(iL + j - k)
                                     + _mm_loadu_ps(iL + j + k));
      }

      //! Get the result
      _mm_storeu_ps(line + j, xsum);
    }
#endif // __SSE__

    //! Normal version
    for (; j < w; j++) {
      float sum = ker[0] * iL[j];

      //! Unroll the loop
      int k = 1;
      for (; k <= kHalf - 4; k += 4) {
        /*sum += ker[k + 0] * (line[j - k + kHalf - 0] + line[j + k + kHalf + 0]);
        sum += ker[k + 1] * (line[j - k + kHalf - 1] + line[j + k + kHalf + 1]);
        sum += ker[k + 2] * (line[j - k + kHalf - 2] + line[j + k + kHalf + 2]);
        sum += ker[k + 3] * (line[j - k + kHalf - 3] + line[j + k + kHalf + 3]);*/
        sum += ker[k + 0] * (iL[j - k - 0] + iL[j + k + 0]);
        sum += ker[k + 1] * (iL[j - k - 1] + iL[j + k + 1]);
        sum += ker[k + 2] * (iL[j - k - 2] + iL[j + k + 2]);
        sum += ker[k + 3] * (iL[j - k - 3] + iL[j + k + 3]);
      }

      //! Finish the loop
      for (; k <= kHalf; k++) {
        //sum += ker[k] * (line[j - k + kHalf] + line[j + k + kHalf]);
        sum += ker[k] * (iL[j - k] + iL[j + k]);
      }

      //! Get the result
      line[j] = sum;
    }

    //! Put the line into the output image
    for (size_t j = 0; j < m_width; j++) {
      oI[j] = line[j];
    }
  }

  //! Release memory for the buffer
  memfree(line);

  //! Release memory of the kernel
  memfree(ker);
}


//! Compute image gradient via symmetric finite difference schemes.
//! Image extension :  symmetrization at x=-1/2.
void Image::computeGradient(
  Image& o_gradX,
  Image& o_gradY,
  const size_t p_ic,
  const size_t p_oc) const {

  //! Check sizes
  if (o_gradX.width() != m_width || o_gradX.height() != m_height ||
      o_gradY.width() != m_width || o_gradY.height() != m_height) {
    cout << "computeGradient: sizes aren't consistant. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  //! Horizontal gradient
  for (size_t i = 0; i < m_height; i++) {
    const float* iI = this->getPtr(p_ic, i);
    float* oY = o_gradY.getPtr(p_oc, i);

    //! Inner image
    size_t j = 1;
#ifdef __AVX__
    //! AVX version
    const __m256 z2 = _mm256_set1_ps(0.5f);
    for (; j < m_width - 1 - 8; j += 8) {
      _mm256_storeu_ps(oY + j, z2 * (_mm256_loadu_ps(iI + j + 1) - _mm256_loadu_ps(iI + j - 1)));
    }
#endif // __AVX__
#ifdef __SSE__
    //! SSE version
    const __m128 x2 = _mm_set1_ps(0.5f);
    for (; j < m_width - 1 - 4; j += 4) {
      _mm_storeu_ps(oY + j, x2 * (_mm_loadu_ps(iI + j + 1) - _mm_loadu_ps(iI + j - 1)));
    }
#endif // __SSE__

    //! Normal version
    for (; j < m_width - 1; j++) {
      oY[j] = (iI[j + 1] - iI[j - 1]) * 0.5f;
    }

    //! Pixels on border
    oY[0] = iI[1] - iI[0];
    oY[m_width - 1] = iI[m_width - 1] - iI[m_width - 2];
  }

  //! Vertical gradient - inner image
  for (size_t i = 1; i < m_height - 1; i++) {
    const float* iT = this->getPtr(p_ic, i - 1);
    const float* iB = this->getPtr(p_ic, i + 1);
    float* oX = o_gradX.getPtr(p_oc, i);
    size_t j = 0;

#ifdef __AVX__
    const __m256 z2 = _mm256_set1_ps(0.5f);
    for (; j < m_width - 8; j += 8) {
      _mm256_storeu_ps(oX + j, z2 * (_mm256_loadu_ps(iB + j) - _mm256_loadu_ps(iT + j)));
    }
#endif // __AVX__
#ifdef __SSE__
    const __m128 x2 = _mm_set1_ps(0.5f);
    for (; j < m_width - 8; j += 8) {
      _mm_storeu_ps(oX + j, x2 * (_mm_loadu_ps(iB + j) - _mm_loadu_ps(iT + j)));
    }
#endif // __SSE__
    for (; j < m_width; j++) {
      oX[j] = (iB[j] - iT[j]) * 0.5f;
    }
  }

  //! Pixels on border
  const float* i0 = this->getPtr(p_ic, 0);
  const float* i1 = this->getPtr(p_ic, 1);
  const float* i2 = this->getPtr(p_ic, m_height - 2);
  const float* i3 = this->getPtr(p_ic, m_height - 1);
  float* oX0 = o_gradX.getPtr(p_oc, 0);
  float* oXh = o_gradX.getPtr(p_oc, m_height - 1);
  for (size_t j = 0; j < m_width; j++) {
    oX0[j] = i1[j] - i0[j];
    oXh[j] = i3[j] - i2[j];
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












