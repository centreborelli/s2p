#include <stdlib.h>
#include <limits.h>

#include "iio.h"
#include "pickopt.h"
#include "fancy_image.h"
#include "sift_anatomy_20141201/lib_sift_anatomy.h"


void print_help(char *v[])
{
    fprintf(stderr, "usage:\n\t%s file.tif [x y w h] [--max-nb-pts n] [-b]"
    //                          0 1         2 3 4 5
            " [--thresh-dog t]\n", *v);
}


int main(int c, char *v[])
{
    // process input arguments
    if (c < 2 || c > 11) {
        print_help(v);
    	return 1;
    }

    // optional arguments
    char s[100];
    sprintf(s, "%d", INT_MAX);
    char *tmp1 = pick_option(&c, &v, "-max-nb-pts", s);
    int max_nb_pts = atoi(tmp1);

    char *tmp2 = pick_option(&c, &v, "-thresh-dog", "0.0133");
    float thresh_dog = atoi(tmp2);

    bool binary = (bool) pick_option(&c, &v, "b", NULL);

    // open the image
    struct fancy_image *fimg = fancy_image_open(v[1], "");

    // define the rectangular region of interest (roi)
    int x, y, w, h;
    if (c == 2) {
        x = 0;
        y = 0;
        w = fimg->w;
        h = fimg->h;
    } else if (c == 6) {
        x = atoi(v[2]);
        y = atoi(v[3]);
        w = atoi(v[4]);
        h = atoi(v[5]);
    } else {
        print_help(v);
        return 1;
    }

    // read the roi in the input image
    float *roi = (float*) malloc(w * h * sizeof(float));
    fancy_image_fill_rectangle_float_split(roi, w, h, fimg, 0, x, y);

    // prepare sift parameters
    struct sift_parameters* p = sift_assign_default_parameters();
    p->C_DoG = thresh_dog;

    // compute sift keypoints
    struct sift_scalespace **ss = malloc(4 * sizeof(struct sift_scalespace*));
    struct sift_keypoints **kk = malloc(6 * sizeof(struct sift_keypoints*));
    for (int i = 0; i < 6; i++)
        kk[i] = sift_malloc_keypoints();
    struct sift_keypoints* kpts = sift_anatomy(roi, w, h, p, ss, kk);

    // add (x, y) offset to keypoints coordinates
    if (x != 0 || y != 0)
        for (int i = 0; i < kpts->size; i++) {
            kpts->list[i]->x += y;  // in Ives' conventions x is the row index
            kpts->list[i]->y += x;
        }

    // write to standard output
    fprintf_keypoints(stdout, kpts, max_nb_pts, binary, 1);

    // debug
    iio_save_image_float("/tmp/roi.tif", roi, w, h);

    // cleanup
    free(kpts);
    free(roi);
    fancy_image_close(fimg);
    return 0;
}
