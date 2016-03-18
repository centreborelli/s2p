#include <stdlib.h>

#include "fancy_image.h"
#include "sift_anatomy_20141201/lib_sift.h"
#include "iio.h"

int main(int c, char *v[])
{
    // process input arguments
    if (c != 2 && c != 6) {
    	fprintf(stderr, "usage:\n\t%s file.tif [x y w h]\n", *v);
        //                          0 1         2 3 4 5
    	return 1;
    }

    // open the image
    struct fancy_image *fimg = fancy_image_open(v[1], "");

    // define the rectangular region of interest (roi)
    int x, y, w, h;
    if (c == 2) {
        x = 0;
        y = 0;
        w = fimg->w;
        h = fimg->h;
    } else {
        x = atoi(v[2]);
        y = atoi(v[3]);
        w = atoi(v[4]);
        h = atoi(v[5]);
    }

    // read the roi in the input image
    float *roi = (float*) malloc(w * h * sizeof(float));
    fancy_image_fill_rectangle_float_split(roi, w, h, fimg, 0, x, y);

    // write roi (debug)
    //iio_save_image_float("/tmp/roi.tif", roi, w, h);

    // compute sift keypoints
    int n;
    struct sift_keypoint_std *k = sift_compute_features(roi, w, h, &n);

    // add (x, y) offset to keypoints coordinates
    if (x != 0 || y != 0)
        for (int i = 0; i < n; i++) {
            k[i].x += x; // TODO: check that there is no x/y swap with Ives convention...
            k[i].y += y;
        }

    // write to standard output
    sift_write_to_file("/dev/stdout", k, n);

    // cleanup
    free(k);
    free(roi);
    fancy_image_close(fimg);
    return 0;
}
