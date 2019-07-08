/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
   name = {fconvol};
   version = {"1.2"};
   author = {"Jacques Froment"};
   function = {"2D-direct convolution of a fimage"};
   usage = {
    in->in           "Input fimage",
    filter->filtre   "convolution filter (fimage)",
    out<-out         "Output fimage"
   };
*/
/*----------------------------------------------------------------------
 v1.1: fixed kmin/lmin/dyS bugs (L.Moisan)
 v1.2: fixed mwcommand syntax (JF)
----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "iio.h"


// dct type II symmetry at the boundary
   int p_sym(int nx, int ny, int x, int y) { 
      if(x < 0)
         x = -x-1;
      if(y < 0)
         y = -y-1;
      if(x >= nx)
         x = -x+2*nx-1;
      if(y >= ny)
         y = -y+2*ny-1;
      return x+nx*y;
   }

// OPERATIONS
float op_min(float* data, int N){
   float vmin=INFINITY; 
   int i;
   for(i=0;i<N;i++) 
      vmin = fmin(vmin,data[i]);
   return vmin;
}
float op_max(float* data, int N){
   float vmax=-INFINITY; 
   int i;
   for(i=0;i<N;i++) 
      vmax = fmax(vmax,data[i]);
   return vmax;
}
float op_average(float* data, int N){
   float vnum=0, vden=0;
   int i;
   for(i=0;i<N;i++) {
      vnum += data[i];
      vden ++;
   }
   return vnum/vden;
}
int compare_floats(const void *a, const void *b)
{
   const float *da = (const float *) a;
   const float *db = (const float *) b;
   return (*da > *db) - (*da < *db);
}
#ifndef ODDP
#define ODDP(x) ((x)&1)
#endif
#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif
#define STATISTIC_MEDIAN_BIAS 0
float op_median(float* data, int N){
   int i;
   float median;
   if (N == 1)
   {
      return data[0];
   }
   // this will modify data
   qsort(data, N, sizeof*data, compare_floats);
   median = data[N/2-1];
   if (EVENP(N))
   {
      int mtype = STATISTIC_MEDIAN_BIAS;
      switch(mtype)
      {
         case -1: break;
         case 0: median += data[N/2]; median /=2; break;
         case 1: median  = data[N/2]; break;
         default: fprintf(stderr,"bad STATISTIC_MEDIAN_BIAS %d", mtype);
      }
   }
   return median;
}
int randombounds(int a, int b)
{
   if (b < a)
      return randombounds(b, a);
   if (b == a)
      return b;
   return a + rand()%(b - a + 1);
}
float op_random(float* data, int N){
   return data[randombounds(0, N-1)];
}
float op_rank(float* data, int N, float center_value){
   int i;
   // this will modify data
   qsort(data, N, sizeof*data, compare_floats);

   // scan until finding the right position
   for(i=0;i<N;i++) 
      if( data[i] >= center_value ) return i;

   // this should never happen
   return -1;
}

// apply a morphologic operation to I (of size nc x nr), using the structuring element S (size se_nc x se_nr)
// The shape of the structuring element is given by the non-zero pixels of S, the center of the structuring element is (ctr_c,ctr_r).
// The operation is indicated with the string opetarion : min, max, median, average, random
// The boundaries of the image are symmetrized, and all the arrays are stored in row major order.
void morphoop(float *ptrI, int nc, int nr,  float *ptrS, int se_nc, int se_nr, int ctr_c, int ctr_r, char* operation, float *ptrO)
{
   int n, m, k, l, kmax, kmin, lmax, lmin;
   float pdata[se_nr*se_nc], S;

   kmin = -ctr_r;         // rows
   kmax = se_nr-1-ctr_r;

   lmin = -ctr_c;         // columns
   lmax = se_nc-1-ctr_c;

   // symmetric boundaries
   for (n = 0; n < nr; n++) // rows
      for (m = 0; m < nc; m++)     // columns
      {
         // scan the structuring element and store the values in pdata
         int i=0;
         for (k = kmin; k <= kmax; k++) // rows
            for (l = lmin; l <= lmax; l++)    // columns
               if(ptrS[se_nc*(k-kmin) + (l-lmin)]) {
                  float v = ptrI[ p_sym(nc, nr, (m + l), (n + k))];
                  if (!isnan(v)) pdata[i++] = v;
               }

         if(strcmp(operation,"min")==0) {
            S = op_min(pdata, i);
         } else if(strcmp(operation,"max")==0) {
            S = op_max(pdata, i);
         } else if(strcmp(operation,"median")==0) {
            S = op_median(pdata, i);
         } else if(strcmp(operation,"average")==0) {
            S = op_average(pdata, i);
         } else if(strcmp(operation,"random")==0) {
            S = op_random(pdata, i);
         } else if(strcmp(operation,"rank")==0) {
            S = op_rank(pdata, i, ptrI[ p_sym(nc, nr, m, n) ]);
         } else {
            fprintf(stderr,"unknown operation: %s", operation);
            exit(1);
         }

	 if (i!=0)
	   *ptrO++ = (float) S;
	 else
	   *ptrO++ = NAN;
      }
}


int main (int argc, char **argv)
{
   /* ppatameter parsing - parameters*/
   if(argc<5) 
   {
      fprintf (stderr, "too few parameters\n");
      fprintf (stderr, "apply a 'morphologic' operation over the square structuring element of size sz\n");
      fprintf (stderr, "   usage: %s input {min,max,median,average,random,rank} sz out [ struc_elem ]\n",argv[0]);
      return 1;
   }


   // input
   int nc,nr,nch;
   float *in = iio_read_image_float_split(argv[1], &nc, &nr, &nch);

   // operation
   char *operation = argv[2];

   // construct structuring element
   int se_nr, se_nc;
   se_nr = se_nc = atoi(argv[3]);
   float sse[se_nr*se_nc];  // static memory
   float *se = sse;
   int i;
   for (i=0;i<se_nr*se_nc;i++) se[i] = 1;
   int ctr_c = se_nc /2; // the center
   int ctr_r = se_nr /2;

   // if a structuring element is passed use it insted of the one constructed
   if (argc >5) {
      se = iio_read_image_float(argv[5], &se_nc, &se_nr);
      ctr_c = se_nc /2; // the center
      ctr_r = se_nr /2;
   }

   // allocate output
   float *out = malloc(nc*nr*nch*sizeof*out);

   // call the algorithm
   for (int i=0;i<nch;i++)
      morphoop(in + nc*nr*i,nc,nr, se, se_nc, se_nr, ctr_c, ctr_r, operation, out + nc*nr*i);


   iio_write_image_float_split(argv[4], out, nc, nr, nch);

   free(in);
   free(out);

   return 0;
}



