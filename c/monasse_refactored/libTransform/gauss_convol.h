#ifndef GAUSSCONVOL_H
#define GAUSSCONVOL_H

#include "libLWImage/LWImage.h"

typedef float flnum;

void gauss_convol(LWImage<flnum>& image, flnum sigma);

#endif
