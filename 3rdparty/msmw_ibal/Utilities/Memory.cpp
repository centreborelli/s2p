/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file Memory.cpp
 * @brief Generic function for allocation and freeing of aligned memory.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <cstdlib>
#include <iostream>
#ifdef _WIN32
#include <malloc.h>
#endif


//! Local includes
#include "Memory.h"


//! Memory aligned allocation.
void* memalloc(
  const size_t p_pad,
  const size_t p_size) {

#ifdef _WIN32
  return _aligned_malloc(p_size, p_pad);
#else
  void* ptr;
  if (posix_memalign(&ptr, p_pad, p_size) != 0) {
    std::cout << "Error during the allocation." << std::endl;
    exit(EXIT_FAILURE);
  }
  return ptr;
#endif
}


//! Release the memory allocated with memalloc.
void memfree(
  void* i_ptr) {

  if (i_ptr != NULL) {
#ifdef _WIN32
    _aligned_free(i_ptr)
#else
    free(i_ptr);
#endif
    i_ptr = NULL;
  }
}
