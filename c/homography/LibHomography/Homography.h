#ifndef HOMOGRAPHY_H_INCLUDED
#define HOMOGRAPHY_H_INCLUDED


//! Global includes


//! Local includes
#include "../Utilities/Parameters.h"
#include "../LibImages/LibImages.h"

/**
 * @brief new function.
 **/
void runHomography(
  const Image& i_im,
  const double i_mat[9],
  Image& o_im,
  const Parameters& p_params);


/**
 * @brief Map the image. May be recursive if zoom != 1.
 **/
void mapImage(
  const Image& i_im,
  const double i_mat[9],
  Image& o_im,
  const bool p_verbose,
  const bool p_adjustSize,
  const bool p_useAntiAliasing = true);


/**
 * @brief Adjust the size of the output image according to the input image size
 *        and the homography that will be applied on it.
 **/
void adjustImage(
  const Image& i_im,
  const double i_mat[9],
  Image& o_im,
  int &o_offsetI,
  int &o_offsetJ);


/**
 * @brief Apply the homography to (x, y, 1).
 **/
void apply(
  const double i_mat[9],
  double& io_x,
  double& io_y);


/**
 * @brief Invert an homography.
 **/
void invert(
  const double i_mat[9],
  double o_mat[9]);


/**
 * @brief Compute min linear compression of H in rectangle [0,w]x[0,h].
 **/
float getMinZoomOut(
  const double i_mat[9],
  const size_t p_w,
  const size_t p_h);


/**
 * @brief Min singular value of the Jacobian of H at point (x, y).
 **/
float getMinSVJacob(
  const double i_mat[9],
  const double i_x,
  const double i_y);


#endif // HOMOGRAPHY_H_INCLUDED
