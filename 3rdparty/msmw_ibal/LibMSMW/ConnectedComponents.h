#ifndef CONNECTEDCOMPONENTS_H_INCLUDED
#define CONNECTEDCOMPONENTS_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief Structure that contains pixel value and position.
 **/
struct pix_new{
  int i;
  int j;
  float v;
};


/**
 * @brief Connected components of positive pixels of the image rep (8 connected)
 *
 * @param io_rep: repository of index position;
 * @param p_width, p_height: size of io_rep 2D array stored in 1D.
 **/
void updatePositiveConnectedComponent(
  int *io_rep,
  const int p_width,
  const int p_height);


/**
 * @brief Find index such as i_array[index] = index and update i_array[p_value].
 *
 * @param i_array: array of indexes that will be updated;
 * @param p_value: index.
 **/
int findCC(
  int *io_array,
  const int p_value);


/**
 * @brief Join connected components.
 *
 * @param io_array: array of index to be updated;
 * @param p_a: first index,
 * @param p_b: second index.
 **/
void joinCC(
  int *io_array,
  const int p_a,
  const int p_b);


/**
 * @brief Initialize index of an image.
 *
 * @param i_im: input image;
 * @param o_rep: output index repository. Must be allocated.
 **/
void initRep(
  const Image& i_im,
  int* o_rep);


/**
 * @brief Remove small regions of an image based on connected components.
 **/
void removeRegions(
  Image& io_im,
  int* io_rep,
  const size_t p_minArea);


#endif // CONNECTEDCOMPONENTS_H_INCLUDED
