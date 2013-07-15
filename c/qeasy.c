#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

#define xmalloc malloc

int main(int c, char *v[])
{
	if (c != 5 && c != 4 && c != 3) {
		fprintf(stderr,"usage:\n\t%s black white  [in [out]]\n", *v);
		//                         0 1     2       3   4
		return EXIT_FAILURE;
	}
	float black = atof(v[1]);
	float white = atof(v[2]);
	char *in = c > 3 ? v[3] : "-";
	char *out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	uint8_t *y = xmalloc(w*h*pd);
	for (int i = 0; i < w*h*pd; i++) {
		float g = x[i];
		g = floor(255 * (g - black)/(white - black));
		if (g < 0) g = 0;
		if (g > 255) g = 255;
		y[i] = g;
	}
	iio_save_image_uint8_vec(out, y, w, h, pd);
	return EXIT_SUCCESS;
}
