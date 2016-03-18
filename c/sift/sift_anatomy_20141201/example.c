#include <stdlib.h>
#include "lib_sift.h"

int main(void)
{
	// create input image
	int w = 512;
	int h = 512;
	float *x = (float*)malloc(w*h*sizeof(float));
	for (int i = 0; i < w*h; i++)
		x[i] = rand();

	// compute sift keypoints
	int n;
	struct sift_keypoint_std *k = sift_compute_points(x, w, h, &n);

	// write to standard output
	sift_write_to_file("/dev/stdout", k, n);

	// cleanup
	free(k);
	free(x);
	return 0;
}
