#include <stdio.h> // only for "fprintf"
#include <stdlib.h> // only for "free"
#include <stdint.h>
#include "iio.h"

// read an image in any format from STDIN and write a ppm to STDOUT
int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s infile outfile\n", *v);
		return 1;
	}
	int w, h, pixeldim;
	//uint16_t *x = iio_read_image_uint16_vec(v[1], &w, &h, &pixeldim);
	float *x = iio_read_image_float_vec(v[1], &w, &h, &pixeldim);
	if (!x) {
		fprintf(stderr, "failed to read an image from file "
				"\"%s\"\n", v[1]);
		return 1;
	}
	fprintf(stderr, "got a %dx%d image with %d channels\n", w, h, pixeldim);
	for (int i = 0; i < pixeldim; i++)
		fprintf(stderr, "p0: %g\n", (float)x[i]);
	//iio_save_image_uint16_vec(v[2], x, w, h, pixeldim);
	iio_save_image_float_vec(v[2], x, w, h, pixeldim);
	free(x);
	return 0;
}
