#include <stdlib.h>
#include <stdio.h>
#include "lib_sift.h"
#include "io_png.h"


int main(int argc, char **argv)
{
    if(argc != 2){
        fprintf(stderr, "usage:\n./sift_basic image.png\n");
        return -1;
    }

	// Loading image
	size_t w, h;
    float* x = io_png_read_f32_gray(argv[1], &w, &h);
    if(!x)
        fatal_error("File \"%s\" not found.", argv[1]);
    for(int i=0; i < w*h; i++)
        x[i] /=256.;

	// compute sift keypoints
	int n;
	struct sift_keypoint_std *k = sift_compute_features(x, w, h, &n);

	// write to standard output
	sift_write_to_file("/dev/stdout", k, n);

	// cleanup
	free(k);
	free(x);
	return 0;
}
