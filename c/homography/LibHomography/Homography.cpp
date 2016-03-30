/**
 * @file Homography.cpp
 *
 * @brief Class to apply an homography to an image with anti-aliasing filter.
 *
 * @author Pascal Monasse (original version),
 * @author Gabriele Facciolo (modified version),
 * @author Marc Lebrun (final version) <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <cmath>
#include <cfloat>


//! Local includes
#include "Homography.h"
#include "Splines.h"
#include "../Utilities/Time.h"
#include "../Utilities/Memory.h"


using namespace std;


//! Apply the homography to the input image.
void runHomography(
  const Image& i_im,
  const double i_mat[9],
  Image& o_im,
  const Parameters& p_params) {

  //! Check the new size
  bool adjustSize = p_params.adjustSize();
  if (p_params.oWidth() > 0 && p_params.oHeight() > 0) {
    o_im.init(p_params.oWidth(), p_params.oHeight(), i_im.channels());
    adjustSize = false;
  }
  else {
    o_im.init(i_im.width(), i_im.height(), i_im.channels());
  }

  //! Apply the homography
  mapImage(i_im, i_mat, o_im, p_params.verbose(), adjustSize);
}


//! Map the image. May be recursive if zoom != 1.
void mapImage(
  const Image& i_im,
  const double i_mat[9],
  Image& o_im,
  const bool p_verbose,
  const bool p_adjustSize,
  const bool p_useAntiAliasing) {

  //! Time test
  Time time;

  //! Check the zoomOut
  const float zoomOut = p_useAntiAliasing ? getMinZoomOut(i_mat, i_im.width(), i_im.height()) : 1.f;

  //! Adjust the size of the image if needed
  int offsetI = 0, offsetJ = 0;
  if (p_adjustSize) {
    adjustImage(i_im, i_mat, o_im, offsetI, offsetJ);
    if (p_verbose) time.get_time(" - adjustImage");
  }

  //! Temporary image if the zoom has to be applied
  Image imTmp;

  //! Temporary homography if the zoom has to be applied
  double matZ[9] = {0.0};

  //! Check for the zoom
  const bool useZ = zoomOut < 1.f;
  if(zoomOut < 1.0f) {

    //! Get the zoom factor
    const float zoomIn = 1.0f / zoomOut;

    //! Compute the new size of the zoomed image
    const int wZoom = (int) std::ceil(o_im.width () * zoomIn * 1.5);
    const int hZoom = (int) std::ceil(o_im.height() * zoomIn * 1.5);

    //! Allocate the new image
    imTmp.init(wZoom, hZoom, i_im.channels());

    //! Compute the new homography associated to the zoom value
    matZ[0] = zoomIn * i_mat[0]; matZ[1] = zoomIn * i_mat[1]; matZ[2] = zoomIn * i_mat[2];
    matZ[3] = zoomIn * i_mat[3]; matZ[4] = zoomIn * i_mat[4]; matZ[5] = zoomIn * i_mat[5];
    matZ[6] =          i_mat[6]; matZ[7] =          i_mat[7]; matZ[8] =          i_mat[8];

    //! Call recursively the function
    mapImage(i_im, matZ, imTmp, false, false);
    if (p_verbose) time.get_time(" - recursive call");

    //! Apply a Gaussian convolution
    const float sigma = 0.8f * sqrt(zoomIn * zoomIn - 1.f);
    imTmp.convolveGaussian(sigma);
    if (p_verbose) time.get_time(" - gaussianConvolution");

    //! Update the homography for the next part.
    matZ[0] = zoomOut; matZ[1] = 0.0    ; matZ[2] = 0.0;
    matZ[3] = 0.0    ; matZ[4] = zoomOut; matZ[5] = 0.0;
    matZ[6] = 0.0    ; matZ[7] = 0.0    ; matZ[8] = 1.0;
  }

  //! Temporary image
  if (!useZ) {
    imTmp = i_im;
  }

  //! Prepare the image for cardinal spline interpolation
  prepareSpline(imTmp);
  if (p_verbose) time.get_time(" - prepareSpline");

  //! Invert the input homography
  double matInv[9];
  invert(useZ ? matZ : i_mat, matInv);

  //! Apply the homography
  for (int i = 0; i < (int) o_im.height(); i++) {
    for (int j = 0; j < (int) o_im.width(); j++) {

      //! Get the inverse position of the current pixel
      double x = j + offsetJ, y = i + offsetI;
      apply(matInv, x, y);

      //! Interpolate the pixels
      interpolateSpline(imTmp, o_im, x + 0.5, y + 0.5, i, j);
    }
  }
  if (p_verbose) time.get_time(" - interpolateSpline");

  //! If a zoom has been used, put the NaN mask back
  if (useZ) {

    //! NaN
    const float NaN = sqrtf(-1.f);

    //! Invert the original homography
    double matInv[9] = {0.0};
    invert(i_mat, matInv);

    //! Put the mask back
    for (size_t c = 0; c < o_im.channels(); c++) {
      for (size_t i = 0; i < o_im.height(); i++) {
        float* oI = o_im.getPtr(c, i);

        for (size_t j = 0; j < o_im.width(); j++) {

          //! Get the inverse position of the current pixel
          double x = j + offsetJ, y = i + offsetI;
          apply(matInv, x, y);

          //! Check its position
          if (x < 0 || x >= i_im.width() || y < 0 || y >= i_im.height()) {
            oI[j] = NaN;
          }
        }
      }
    }
    if (p_verbose) time.get_time(" - put the mask back");
  }
}


/**
 * @brief Min singular value of the Jacobian of H at point (x, y).
 **/
float getMinSVJacob(
  const double i_mat[9],
  const double i_x,
  const double i_y) {

  //! Get the point after applying the homography
  const double X[3] = {i_mat[0] * i_x + i_mat[1] * i_y + i_mat[2],
                       i_mat[3] * i_x + i_mat[4] * i_y + i_mat[5],
                       i_mat[6] * i_x + i_mat[7] * i_y + i_mat[8]};

  //! For convenience
  const double z = 1.0 / X[2];
  const double x = z   * X[0];
  const double y = z   * X[1];

  //! Get the Jacobian: J = (a b, c d)
  const double a = z * (i_mat[0]- i_mat[6] * x);
  const double b = z * (i_mat[1]- i_mat[7] * x);
  const double c = z * (i_mat[3]- i_mat[6] * y);
  const double d = z * (i_mat[4]- i_mat[7] * y);

  //! Return the smallest singular value
  return sqrt(0.5 * (a * a + b * b + c * c + d * d -
         sqrt((a * a + b * b - c * c - d * d) * (a * a + b * b - c * c - d * d)
              + 4.0 * (a * c + b * d) * (a * c + b * d))));
}


//! Compute min linear compression of H in rectangle [0,w]x[0,h].
float getMinZoomOut(
  const double i_mat[9],
  const size_t p_w,
  const size_t p_h) {

  return min(1.f,
         min(getMinSVJacob(i_mat, 0, 0),
         min(getMinSVJacob(i_mat, p_w, 0),
         min(getMinSVJacob(i_mat, 0, p_h), getMinSVJacob(i_mat, p_w, p_h)))));
}


//! Adjust the size of the output image according to the input image size and
//! the homography that will be applied on it.
void adjustImage(
  const Image& i_im,
  const double i_mat[9],
  Image& o_im,
  int &o_offsetI,
  int &o_offsetJ) {

  //! Initialization of the rectangle
  float x1 = FLT_MAX, y1 = FLT_MAX, x2 = -FLT_MAX, y2 = -FLT_MAX;

  //! Add image of rectangle (0,0,\a w,\a h) by \a i_mat to the size.
  //! Top left
  double x = 0.0, y = 0.0;
  apply(i_mat, x, y);
  if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
  if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;

  //! Top right
  x = double(i_im.width()) - 1.0; y = 0.0;
  apply(i_mat, x, y);
  if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
  if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;

  //! Bottom right
  x = double(i_im.width()) - 1.0; y = double(i_im.height()) - 1.0;
  apply(i_mat, x, y);
  if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
  if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;

  //! Bottom left
  x = 0.0; y = double(i_im.height()) - 1.0;
  apply(i_mat, x, y);
  if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
  if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;

  //! Get the new border of the image
  const int w = (x2 < x1)? 0: int(x2 - x1 + 1);
  const int h = (y2 < y1)? 0: int(y2 - y1 + 1);
  o_offsetI = int(y1);
  o_offsetJ = int(x1);

  //! Resize the output image
  o_im.init(w, h, i_im.channels());
}


//! Apply the homography to (x, y, 1).
void apply(
  const double i_mat[9],
  double& io_x,
  double& io_y) {

  //! Product between the homography and the vector (x, y, 1)
  const double x = i_mat[0] * io_x + i_mat[1] * io_y + i_mat[2];
  const double y = i_mat[3] * io_x + i_mat[4] * io_y + i_mat[5];
  const double z = i_mat[6] * io_x + i_mat[7] * io_y + i_mat[8];

  //! Get the normalized result
  io_x = x / z;
  io_y = y / z;
}


//! Invert an homography.
void invert(
  const double i_mat[9],
  double o_mat[9]) {

  //! Get the determinant
  const double det = 1.0 / (i_mat[0] * i_mat[4] * i_mat[8] -
                            i_mat[0] * i_mat[5] * i_mat[7] -
                            i_mat[1] * i_mat[3] * i_mat[8] +
                            i_mat[1] * i_mat[5] * i_mat[6] +
                            i_mat[2] * i_mat[3] * i_mat[7] -
                            i_mat[2] * i_mat[4] * i_mat[6]);

  //! Invert the matrix
  o_mat[0] = det * (i_mat[4] * i_mat[8] - i_mat[5] * i_mat[7]);
  o_mat[1] = det * (i_mat[2] * i_mat[7] - i_mat[1] * i_mat[8]);
  o_mat[2] = det * (i_mat[1] * i_mat[5] - i_mat[2] * i_mat[4]);
  o_mat[3] = det * (i_mat[5] * i_mat[6] - i_mat[3] * i_mat[8]);
  o_mat[4] = det * (i_mat[0] * i_mat[8] - i_mat[2] * i_mat[6]);
  o_mat[5] = det * (i_mat[2] * i_mat[3] - i_mat[0] * i_mat[5]);
  o_mat[6] = det * (i_mat[3] * i_mat[7] - i_mat[4] * i_mat[6]);
  o_mat[7] = det * (i_mat[1] * i_mat[6] - i_mat[0] * i_mat[7]);
  o_mat[8] = det * (i_mat[0] * i_mat[4] - i_mat[1] * i_mat[3]);
}
