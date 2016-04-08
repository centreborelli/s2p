/*
IPOL SIFT
Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan 
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20140911 (September 11th, 2014)

This C ANSI source code is related to the IPOL publication

    [1] "Anatomy of the SIFT Method." 
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2013.
        http://www.ipol.im/pub/algo/rd_anatomy_sift/

An IPOL demo is available at
        http://www.ipol.im/pub/demo/rd_anatomy_sift/





== Patent Warning and License =================================================

The SIFT method is patented 

    [2] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89
  
 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.


This file also implements the colormap published in the paper
   [3] "Diverging Color Maps for Scientific Visualization."
        Kenneth Moreland
        International Symposium on Visual Computing 2009



*/


#ifndef _SIFT_IO_SCALESPACE_H_
#define _SIFT_IO_SCALESPACE_H_

#include "lib_scalespace.h"

void print_sift_scalespace_gray(const struct sift_scalespace* scalespace, const char* basename);
void print_sift_scalespace_gray_nearestneighbor(const struct sift_scalespace* scalespace, const char* basename);
void print_sift_scalespace_rgb(const struct sift_scalespace* scalespace, const char* basename);
void print_sift_scalespace_rgb_nointerp(const struct sift_scalespace* scalespace, const char* basename);

#endif
