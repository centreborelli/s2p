/*
 * Copyright (c) 2014, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED


/**
 * @brief Memory aligned allocation.
 **/
void* memalloc(
  const size_t p_pad,
  const size_t p_size);


/**
 * @brief Release the memory allocated with memalloc.
 **/
void memfree(
  void* i_ptr);



#endif // MEMORY_H_INCLUDED
