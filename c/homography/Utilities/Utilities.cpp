/**
 * @file Utilites.cpp
 *
 * @brief Miscellaneous functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <png.h>
#include <tiffio.h>
#include <omp.h>


//! Local includes
#include "Utilities.h"


using namespace std;


//! Compute an amount of time elapsed.
void getTime(
  struct timespec &io_start,
  const char* p_name) {

  //! Check the current time
  struct timespec finish;
  clock_gettime(CLOCK_MONOTONIC, &finish);

  //! Compute the elapsed time
  double elapsed = (finish.tv_sec - io_start.tv_sec) * 1000;
  elapsed += (finish.tv_nsec - io_start.tv_nsec) / 1000000.0;

  //! Print the result
  cout << p_name << ": (ms) = " << elapsed << endl;

  //! Start a new timer
  clock_gettime(CLOCK_MONOTONIC, &io_start);
}


//! Return the optimal cut of the lines of an image according to the number of
//! threads available.
void initializeHeights(
  const size_t i_height,
  size_t* o_heights,
  const size_t p_nbThreads) {

  //! Get the initial number of lines for each thread to deal with.
  const size_t nb = p_nbThreads;
  const size_t h = ceil(float(i_height) / float(nb));

  //! Initialize the number of lines
  size_t* lines = (size_t*) malloc(nb * sizeof(size_t));
  for (size_t n = 0; n < nb; n++) {
    lines[n] = h;
  }

  //! Remove one by one the lines that are not needed
  int nCur = nb - 1;
  int k = int(h * nb) - int(i_height);
  while (k > 0) {
    lines[nCur]--;
    if (nCur > 0) {
      nCur--;
    }
    else {
      nCur = nb - 1;
    }
    k--;
  }

  //! Get the lines position
  o_heights[0] = 0;
  for (size_t n = 0; n < nb; n++) {
    o_heights[n + 1] = o_heights[n] + lines[n];
  }

  //! Release memory
  free(lines);
}


//! Read an homography.
void readHomography(
  const char* i_fileName,
  double o_mat[9]) {

  //! Open the file
  ifstream file(i_fileName);

  //! Check that the file has been opened
  if (!file) {
    cout << "Can't read the file: " << i_fileName << endl;
    exit(EXIT_FAILURE);
  }

  //! Read the 9 values in the 3 lines
  for (size_t l = 0; l < 3; l++) {

    //! Get the line
    string line;
    getline(file, line);

    //! Get the values
    istringstream iss(line);

    //! First line begins with '[' character.
    if (l == 0) {
      char c;
      iss >> c >> o_mat[3 * l] >> o_mat[3 * l + 1] >> o_mat[3 * l + 2];
    }
    else {
      iss >> o_mat[3 * l] >> o_mat[3 * l + 1] >> o_mat[3 * l + 2];
    }
  }

  //! Close the file
  file.close();
}













