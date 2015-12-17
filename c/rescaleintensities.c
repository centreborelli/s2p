#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

#define xmalloc malloc


int main(int c, char *v[])
{
	if (c != 5 ) {
		fprintf(stderr, "usage:\n\t%s in out min max\n", *v);
		//                          0  1  2   3   4
		return EXIT_FAILURE;
	}
	char *in = v[1];
	char *out = v[2];
    unsigned long rmin = strtoul(v[3],NULL,0);
    unsigned long rmax = strtoul(v[4],NULL,0);

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);

	uint8_t *y = xmalloc(w*h*pd);
	for (int i = 0; i < w*h*pd; i++) {
		float g = x[i];
		g = floor(255 * (g - rmin)/(rmax - rmin));
		if (g < 0) g = 0;
		if (g > 255) g = 255;
		y[i] = g;
	}
    
	iio_save_image_uint8_vec(out, y, w, h, pd);
    
	return EXIT_SUCCESS;
}
