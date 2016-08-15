/* Copyright (C) 2016, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "abstract_dsf.c"

// connected components of positive pixels of the image rep
static void connected_component_filter(int *rep, int w, int h, float *in, float intensity_threshold)
{
   adsf_begin(rep, w*h);

   // remove from dsf pixels with NANs in input
   for (int i = 0; i < w*h; i++)
      if (isnan(in[i]))
	      rep[i] = -1;

   for (int j = 0; j < h - 1; j++)
      for (int i = 0; i < w - 1; i++)
      {
         int p0 = j*w + i;
         int p1 = j*w + i+1;
         int p2  = (j+1)*w + i;
         if (rep[p0] >= 0 && rep[p1] >= 0 && 
             fabs(in[p0] - in[p1]) < intensity_threshold)
            adsf_union(rep, w*h, p0, p1); 
         if (rep[p0] >= 0 && rep[p2] >= 0 && 
             fabs(in[p0] - in[p2]) < intensity_threshold)
            adsf_union(rep, w*h, p0, p2);

      }
   for (int i = 0; i < w*h; i++)
      if (rep[i] >= 0)
         rep[i] = adsf_find(rep, w*h, i);
}



// removes small connected components of non-nan values in the input
int remove_small_cc(int w, int h, float *in, float *out, int minarea, float intensity_threshold)
{
	// initialize the output
	for (int i = 0; i < w*h; i++)
		out[i] = in[i];

	// identify the connected components of non nan values
	int *rep = (int*) malloc(w * h * sizeof*rep);
	connected_component_filter(rep, w, h, in, intensity_threshold);

	// measure the area of connected components
	int *area = (int*) malloc(w * h * sizeof*area);
	for (int i = 0; i < w*h; i++)
		area[i] = 0;
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0)
			area[rep[i]] += 1;

	// set to NAN the pixels whose connected component is too small
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0 && area[rep[i]] <= minarea)
			out[i] = NAN;

	// count how many connected components have been kept
	int remaining_cc = 0;
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0 && area[rep[i]] > minarea)
			remaining_cc += 1;

	// cleanup and exit
	free(rep);
	free(area);
	return remaining_cc;
}


#ifndef SKIP_MAIN


#include "iio.h"
#include "pickopt.c"


int main(int argc, char **argv)
{

   if (argc < 3) {
      fprintf(stderr, " %s img out area [intensity_threshold (default INF)]\n", argv[0]);
      fprintf(stderr, "   remove connected compotents of img smaller than area\n"
                      "   if intensity difference between two neighboring pixels of img\n"
                      "   is larger than the threshold then the pixels are considered disjoint\n");
      return 1;
   }

   char *img_file = argv[1];
   char *out_file = argv[2];
   int area = atoi(argv[3]);
   float intensity_threshold = argc > 4 ? atof(argv[4]): INFINITY;

   int w,h,nch;

   float *img = iio_read_image_float_split(img_file, &w, &h, &nch);
   float *out = calloc(w*h*nch, sizeof(*out));
   int n_cc = remove_small_cc(w, h, img, out, area, intensity_threshold);
   fprintf(stderr, "remaining components %d\n", n_cc);


	/* write the results */
   iio_save_image_float_split(out_file, out,w,h,nch);
   free(img);
   free(out);
}


#endif
