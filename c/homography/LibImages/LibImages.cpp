/**
 * @brief Class to deal with images.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com> (original version)
 * @author Carlo de Franchis <carlodef@gmail.com> (modified version)
 **/

//! Global includes
#include <limits.h>
#include <inttypes.h>
#include <gdal_priv.h>
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
#endif // _OPENMP

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
  if (ext == "tif" || ext == "tiff" || ext == "png") {
    this->readGDAL(p_name, i_border);
  }
  else {
    cout << "Extension " << ext << " not known. Abort." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Read a mono-channel uint16 tiff (or png) image.
void Image::readGDAL(
  const char* p_name,
  const size_t i_border) {

  GDALAllRegister();
  GDALDataset *poDataset = (GDALDataset *) GDALOpen(p_name, GA_ReadOnly);

  if (poDataset == NULL) {
    fprintf(stderr, "ERROR: can't open %s\n", p_name);
    exit(EXIT_FAILURE);
  }

  int w = poDataset->GetRasterXSize();
  int h = poDataset->GetRasterYSize();

  // read the needed ROI in the input image
  GDALRasterBand *poBand = poDataset->GetRasterBand(1);
  float *roi = (float *) CPLMalloc(sizeof(float)*w*h);
  int errorRasterIO = poBand->RasterIO(GF_Read, 0, 0, w, h, roi, w, h, GDT_Float32, 0, 0);
  if (errorRasterIO != CPLE_None)
    fprintf(stderr, "errorRasterIO = %d\n", errorRasterIO);
  GDALClose((GDALDatasetH) poDataset);

  //! Allocate the image
  new (this) Image(roi, (const size_t) w, (const size_t) h, 1);
  this->m_border = i_border;

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
  if (ext == "tif" || ext == "tiff" || ext == "png") {
    if (ext == "png") {
      this->writeGDAL(p_name, p_quad, "PNG");
    }
    else {
      this->writeGDAL(p_name, p_quad, "GTiff");
    }
  }
  else {
    cout << "Extension " << ext << " not known. Abort." << endl;
    exit(EXIT_FAILURE);
  }
}


//! Write an image via the libgdal write function
void Image::writeGDAL(
  const char *p_name,
  const bool p_quad,
  const char *pszFormat) const {

  //! For convenience
  const size_t h = m_height * (p_quad ? 2 : 1);
  const size_t w = m_width  * (p_quad ? 2 : 1);

  GDALDriver *poDriver;
  GDALDataset *poDstDS;
  char **papszOptions = NULL;

  if (strcmp(pszFormat,"GTiff") == 0){
      poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
      poDstDS = poDriver->Create( p_name, w, h, m_channels, GDT_Float32, papszOptions );
    }
  else if (strcmp(pszFormat,"PNG") == 0) {
      poDriver = GetGDALDriverManager()->GetDriverByName("MEM");
      poDstDS = poDriver->Create("", w, h, m_channels, GDT_Float32, papszOptions);
    }
  else {
    cout << "Driver not know. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  GDALRasterBand *poBand;

  float* line = (float*) memalloc(16, (p_quad ? w : m_width) * sizeof(float));

  //! Write the file data
  for (size_t c = 0; c < m_channels; c++) {
    poBand = poDstDS->GetRasterBand(c+1);
    for (size_t i = 0; i < h; i++) {
      const float* iI = this->getPtr(c, p_quad ? i / 2 : i);

      //! Copy the line
      if (p_quad)
        for (size_t j = 0; j < m_width; j++)
          line[2 * j + 0] = line[2 * j + 1] = iI[j];
      else
        for (size_t j = 0; j < m_width; j++)
          line[j] = iI[j];

      CPLErr ee = poBand->RasterIO(GF_Write, 0, i, w, 1, line, w, 1, GDT_Float32, 0, 0);
    }
  }

  //! Release memory
  memfree(line);

  // CreateCopy for PNG format
  if (strcmp(pszFormat,"PNG") == 0) {
    GDALDriver      *pDriver;
    GDALDataset     *pDs;
    pDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    pDs = pDriver->CreateCopy( p_name, poDstDS, FALSE, NULL, NULL, NULL);
    GDALClose( (GDALDatasetH) pDs );
  }

  /* Once we're done, close properly the dataset */
  GDALClose( (GDALDatasetH) poDstDS );
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

  //! The kernel is truncated at 4 sigmas from center.
  int kSize = (int) (8.f * p_sigma + 1.f);

  //! The kernel size should be odd
  kSize = kSize > 3 ? kSize + 1 - kSize % 2 : 3;
  const int kHalf = kSize / 2;

  //! Allocate the kernel
  float* kernel = (float*) memalloc(16, kSize * sizeof(float));

  //! Compute the kernel
  float sum = 0.f;
  for (int i = 0; i < kSize; i++) {
    const float x = (float) (i - kHalf);
    kernel[i] = exp(-x * x / (2.f * p_sigma * p_sigma));
    sum += kernel[i];
  }

  //! Normalize the kernel
  for (int i = 0; i < kSize; i++) {
    kernel[i] /= sum;
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
        line[j] = iI[0];
        line[kHalf + w + j] = iI[w - 1];
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
        xCol[i] = _mm_loadu_ps(this->getPtr(c, 0) + j);
        xCol[kHalf + h + i] = _mm_loadu_ps(this->getPtr(c, h - 1) + j);
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
        col[i] = this->getPtr(c, 0)[j];
        col[kHalf + h + i] = this->getPtr(c, h - 1)[j];
      }
      for (int i = 0; i < h; i++) {
        col[kHalf + i] = this->getPtr(c, i)[j];
      }

      //! Apply the convolution
      //ConvBufferFast(col, kernel, h, kSize);
      for (int i = 0; i < h; i++) {
        float value = 0.f;

        //! Unroll the loop
        int k = 0;
        while (k + 4 < kSize) {
          value += kernel[k + 0] * col[i + k + 0] +
                   kernel[k + 1] * col[i + k + 1] +
                   kernel[k + 2] * col[i + k + 2] +
                   kernel[k + 3] * col[i + k + 3];
          k += 4;
        }
        for (; k < kSize; k++) {
          value += kernel[k] * col[i + k];
        }

        //col[i] = value;
        //! Get the value
        this->getPtr(c, i)[j] = col[i];
      }

      //! Get the value
      /*for (int i = 0; i < h; i++) {
        this->getPtr(c, i)[j] = col[i];
      }*/
    }

    //! Release memory
    memfree(col);
  }

  //! Release memory
  memfree(kernel);
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
