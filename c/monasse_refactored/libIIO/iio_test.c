#include "iio.h"
#include <stdio.h> // only for "fprintf"
#include <stdlib.h> // only for "free"

// read an image in any format from STDIN and write a ppm to STDOUT
int main(int c, char *v[])
{
	int w, h, pixeldim;
	float *x = iio_read_image_float_vec("-", &w, &h, &pixeldim);
	fprintf(stderr, "got a %dx%d image with %d channels\n", w, h, pixeldim);
	iio_save_image_float_vec("-", x, w, h, pixeldim);
	free(x);
	return 0;
}
