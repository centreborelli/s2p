// Copyright 2014, Gabriele Facciolo <gabriele.facciolo@upf.edu>
// All rights reserved.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <assert.h>
//extern "C"{
//#include "iio.h"
//}
#include "vvector.h"



	/*   dim=3
    *  a x + b y + c  
    * */
float eval_model(int dim, float *data, double *model)  {
	return (model[0]*data[0] + model[1]*data[1] + model[2]);
}

float generate_model ( int dim, float *data, int n, double *result) {

   double A[3][3]={{0,0,0},{0,0,0},{0,0,0}};
   double b[3]={0,0,0};
      
   for(int ii=0;ii<n;ii++) {
      for (int i=0;i<3;i++) {
         double d[4] = {data[dim*ii+0], data[dim*ii+1], 1, data[dim*ii+2]};
         for (int j=0;j<3;j++) {
            A[i][j] += d[i]*d[j];
         }
         b[i] += d[i]*d[3];
      }
   }

   double iA[3][3];
   double detA;
   INVERT_3X3(iA, detA, A);

//     for (int i=0;i<3;i++) 
//         printf("%f %f %f    %f\n", A[i][0],A[i][1],A[i][2],b[i]);
//         printf("\n");
//     for (int i=0;i<3;i++) 
//         printf("%f %f %f    %f\n", iA[i][0],iA[i][1],iA[i][2],b[i]);
//         printf("\n");

   if(fabs(detA)>0.00001) {
      MAT_DOT_VEC_3X3(result, iA,b);
   } else {
      // return just the average
		result[0]=0;
		result[1]=0;
		result[2]=b[2]/n;
	}

}

void TEST_LSE()
{
   float data[3*6]={
      510,700,1,
      510,710,1,
      500,710,1,
      501,703,1,
      502,702,1,
      503,701,1
   };
   double model[3];

      generate_model ( 3, data, 6, model);

      printf("%f %f %f\n", model[0],model[1], model[2]);
      for(int i = 0;i<6;i++){
         printf("%f\n", fabs(eval_model(3, &(data[i*3]), model) - data[i*3+2]));
      }
}

static int dsf_find(int *t, int a)
{
   if (a != t[a])
      t[a] = dsf_find(t, t[a]);
   return t[a];
}

static int dsf_make_link(int *t, int a, int b)
{
   if (a < b) { // arbitrary choice
      t[b] = a;
      return a;
   } else {
      t[a] = b;
      return b;
   }
}

static int dsf_join(int *t, int a, int b)
{
   a = dsf_find(t, a);
   b = dsf_find(t, b);
   if (a != b)
      b = dsf_make_link(t, a, b);
   return b;
}



// connected components of positive pixels of the image rep
static void positive_connected_component_filter(int *rep, int w, int h)
{
   for (int i = 0; i < w*h; i++)
      if (rep[i] >= 0)
         rep[i] = i;
   for (int j = 0; j < h - 1; j++)
      for (int i = 0; i < w - 1; i++)
      {
         int p0 = j*w + i;
         int p1 = j*w + i+1;
         int p2  = (j+1)*w + i;
         if (rep[p0] >= 0 && rep[p1] >= 0)
            dsf_join(rep, p0, p1);
         if (rep[p0] >= 0 && rep[p2] >= 0)
            dsf_join(rep, p0, p2);

      }
   for (int i = 0; i < w*h; i++)
      if (rep[i] >= 0)
         rep[i] = dsf_find(rep, i);
}

// connected components of positive pixels of the image rep (8 connected)
static void positive_connected_component_filter_8(int *rep, int w, int h)
{
   for (int i = 0; i < w*h; i++)
      if (rep[i] >= 0)
         rep[i] = i;
   for (int j = 0; j < h - 1; j++)
      for (int i = 0; i < w; i++)
      {
         int p0 = j*w + i;
         int p1 = j*w + i+1;
         int p2 = (j+1)*w + i;
         int p3 = (j+1)*w + i+1;
         int p4 = (j+1)*w + i-1;
         if (i<w && rep[p0] >= 0 && rep[p1] >= 0)
            dsf_join(rep, p0, p1);
         if (rep[p0] >= 0 && rep[p2] >= 0)
            dsf_join(rep, p0, p2);
         if (i<w && rep[p0] >= 0 && rep[p3] >= 0)
            dsf_join(rep, p0, p3);
         if (i>0 && rep[p0] >= 0 && rep[p4] >= 0)
            dsf_join(rep, p0, p4);

      }
   for (int i = 0; i < w*h; i++)
      if (rep[i] >= 0)
         rep[i] = dsf_find(rep, i);
}


struct pix{
   int x;
   int y;
   float v;
};

void extract_cc_boundary_at_pix(int w,int h,int i,int j,int *rep,float *img, std::vector<struct pix > *ccboundary)
{
   int dx[] = {-1,0,1,0};
   int dy[] = {0,-1,0,1};

   int idx = j*w + i;
   int ccidx = rep[idx];
   // store the boundary values
   for(int u=0;u<4;u++) 
   {
      int ii= i+dx[u];
      int jj= j+dy[u];
      if(ii>=0 && jj>=0 && ii<w && jj<h) {
         int iidx = ii+jj*w;
         if( rep[iidx] != ccidx) {
            struct pix pp;
            pp.x=ii;
            pp.y=jj;
            pp.v=img[iidx];
            ccboundary->push_back(pp);
         }
      }
   }
}



int cc_postprocess(int w, int h, float *img, float *msk, float *out, float *outmsk, float MAXERR, float MAXDIFF, int LR_REVERSE) {

   // initialize the image to 
   std::vector< int > rep(w*h); 
   for (int i=0;i<w*h;i++) 
      rep[i] = msk[i] > 0 ? -1 : 1;

   // identify the connected components of positive values
   positive_connected_component_filter(&(rep.front()), w, h);

   // count connected components
   int n_cc = 0;
   std::vector< int > cc_ids(w*h); 
   for (int i=0;i<w*h;i++) {
      if (rep[i] == i) {
         cc_ids[n_cc] = i;
         n_cc += 1;
      }
   }

   // store the cc_pixels  and cc_boundary in vectors of vectors
   std::vector< std::vector< struct pix > > cc_pixels(w*h);
   std::vector< std::vector< struct pix > > cc_boundary(w*h);

   for (int j = 0; j < h; j++)
      for (int i = 0; i < w; i++)
      {
         int idx = j*w + i;
         int cc = rep[idx];
         if (cc>=0) { // if it is a region 
            struct pix pp;
            pp.x=i;
            pp.y=j;
            pp.v=img[idx];
            cc_pixels[cc].push_back(pp);
            extract_cc_boundary_at_pix(w, h, i, j, &(rep.front()), img, &(cc_boundary[cc]));
         }
      }
   //printf("%d\n", n_cc);

   for (int i=0;i<w*h;i++) {
      out[i] = img[i];
      outmsk[i] = msk[i];
   }

   // for each connected component
   std::vector< float > data(3*w*h);
   for(int t =0;t < n_cc; t++) 
   {
      int cc = cc_ids[t];
      std::vector< struct pix > pixels   =  cc_pixels[cc];
      std::vector< struct pix > boundary =  cc_boundary[cc];
      int l=boundary.size();

      // copy the data at the boundary and compute the affine model
      double model[3];
      for(int i = 0;i<l;i++) {
         data[i*3+0]=(float)boundary[i].x;
         data[i*3+1]=(float)boundary[i].y;
         data[i*3+2]=(float)boundary[i].v;
      }
      generate_model ( 3, &(data.front()), l, model);

      // compute stuff
      float maxerr = 0;
      float sum = 0;
      float vmax = -INFINITY;
      float vmin = INFINITY;
      for(int i = 0;i<l;i++){
         sum += boundary[i].v;
         vmax = fmax(boundary[i].v,vmax);
         vmin = fmin(boundary[i].v,vmin);

         float err =  eval_model(3, &(data[i*3]), model) - data[i*3+2];
         maxerr = fmax(maxerr, fabs(err));
      }
      sum/=l;

      // interpolate the region if necessary
      if(LR_REVERSE*model[0] < 1.0 )      //      if(model[0] > -1.0 )
      if(maxerr< MAXERR && vmax-vmin< MAXDIFF ) {
         for(int i = 0;i<pixels.size();i++) {
            int idx = pixels[i].x + pixels[i].y*w;
            float vv[3] = {pixels[i].x, pixels[i].y, 0};
            out[idx] = eval_model(3, vv, model);
            outmsk[idx] = 254;
         }
         //      printf("id:%d sz:%d per:%d vmax-vmin:%f avg:%f maxerr:%f\n", 
         //            cc, (int) pixels.size(), (int)boundary.size(), vmax-vmin, sum, maxerr);
      }
   }
   return n_cc;
}




// remove small connected components (positive valued) of the mask
int cc_grain(int w, int h, const float *msk, float *outmsk, int minarea) {

   // initialize the rep image to  -1 or 1
   std::vector< int > rep(w*h); 
   for (int i=0;i<w*h;i++) 
      rep[i] = msk[i] > 0 ? 1 : -1;

   // identify the connected components of positive values
   positive_connected_component_filter_8(&(rep.front()), w, h);

   // count connected components
   int n_cc = 0;
   std::vector< int > cc_ids(w*h); 
   for (int i=0;i<w*h;i++) {
      if (rep[i] == i) {
         cc_ids[n_cc] = i;
         n_cc += 1;
      }
   }

   // store the cc_pixels
   std::vector< std::vector< struct pix > > cc_pixels(w*h);

   for (int j = 0; j < h; j++)
      for (int i = 0; i < w; i++)
      {
         int idx = j*w + i;
         int cc = rep[idx];
         if (cc >= 0) { // if it is a region 
            struct pix pp;
            pp.x=i;
            pp.y=j;
            pp.v=msk[idx];
            cc_pixels[cc].push_back(pp);
         }
      }

   for (int i=0;i<w*h;i++) {
      outmsk[i] = msk[i];
   }

   // for each connected component
   int removed_cc = 0;
   for(int t = 0; t < n_cc; t++) 
   {
      int cc = cc_ids[t];
      std::vector< struct pix > pixels = cc_pixels[cc];

      // remove small regions 
      if (pixels.size() < minarea) {
         for(int i = 0;i<pixels.size();i++) {
            int idx = pixels[i].x + pixels[i].y*w;
            outmsk[idx] = 0;
         } 
         removed_cc++;
      }

   }
   return removed_cc;
}


#include "smartparameter.h"
SMART_PARAMETER(MAXERR,2.0)
SMART_PARAMETER(MAXDIFF,5.0)
SMART_PARAMETER(LR_REVERSE,1)





///***************************/
//
//// c: pointer to original argc
//// v: pointer to original argv
//// o: option name after hyphen
//// d: default value (if NULL, the option takes no argument)
//static char *pick_option(int *c, char ***v, char *o, char *d)
//{
//   int argc = *c;
//   char **argv = *v;
//   int id = d ? 1 : 0;
//   for (int i = 0; i < argc - id; i++)
//      if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, o)) {
//         char *r = argv[i + id] + 1 - id;
//         *c -= id + 1;
//         for (int j = i; j < argc - id; j++)
//            (*v)[j] = (*v)[j + id + 1];
//         return r;
//      }
//   return d;
//}
//
//
//int main(int argc, char **argv)
//{
//
//   if (argc < 3) {
//      fprintf(stderr, "[MAXERR=%.2f MAXDIFF=%.2f LR_REVERSE={1(LR),0(disabled),-1(RL)}] %s img mask out [outmsk]\n", MAXERR(), MAXDIFF(), argv[0]);
//      return 1;
//   }
//
//   char *img_file = argv[1];
//   char *msk_file = argv[2];
//   char *out_file = argv[3];
//   char *outmsk_file = argc>3 ? argv[4]: NULL;
//
//   int w,h,nch;
//
//   float *img = iio_read_image_float_split(img_file, &w, &h, &nch);
//   float *msk = iio_read_image_float_split(msk_file, &w, &h, &nch);
//   float *out = (float*) calloc(w*h*nch, sizeof(*out));
//   float *outmsk = (float*) calloc(w*h*nch, sizeof(*outmsk));
//   nch=1;
//
//   
//   int n_cc = cc_postprocess(w, h, img, msk, out, outmsk, MAXERR(), MAXDIFF(), LR_REVERSE());
//   printf("%d\n", n_cc);
//
////   // initialize the image to 
////   std::vector< int > rep(w*h); 
////   for (int i=0;i<w*h;i++) 
////      rep[i] = msk[i] > 0 ? -1 : 1;
////
////   // identify the connected components of positive values
////   positive_connected_component_filter(&(rep.front()), w, h);
////
////   // count connected components
////   int n_cc = 0;
////   std::vector< int > cc_ids(w*h); 
////   for (int i=0;i<w*h;i++) {
////      if (rep[i] == i) {
////         cc_ids[n_cc] = i;
////         n_cc += 1;
////      }
////   }
////
////   // store the cc_pixels  and cc_boundary in vectors of vectors
////   std::vector< std::vector< struct pix > > cc_pixels(w*h);
////   std::vector< std::vector< struct pix > > cc_boundary(w*h);
////
////   for (int j = 0; j < h; j++)
////      for (int i = 0; i < w; i++)
////      {
////         int idx = j*w + i;
////         int cc = rep[idx];
////         if (cc>=0) { // if it is a region 
////            struct pix pp;
////            pp.x=i;
////            pp.y=j;
////            pp.v=img[idx];
////            cc_pixels[cc].push_back(pp);
////            extract_cc_boundary_at_pix(w, h, i, j, &(rep.front()), img, &(cc_boundary[cc]));
////         }
////      }
////   printf("%d\n", n_cc);
////
////   for (int i=0;i<w*h;i++) {
////      out[i] = img[i];
////      outmsk[i] = msk[i];
////   }
////
////   // for each connected component
////   std::vector< float > data(3*w*h);
////   for(int t =0;t < n_cc; t++) 
////   {
////      int cc = cc_ids[t];
////      std::vector< struct pix > pixels   =  cc_pixels[cc];
////      std::vector< struct pix > boundary =  cc_boundary[cc];
////      int l=boundary.size();
////
////      // copy the data at the boundary and compute the affine model
////      double model[3];
////      for(int i = 0;i<l;i++) {
////         data[i*3+0]=(float)boundary[i].x;
////         data[i*3+1]=(float)boundary[i].y;
////         data[i*3+2]=(float)boundary[i].v;
////      }
////      generate_model ( 3, &(data.front()), l, model);
////
////      // compute stuff
////      float maxerr = 0;
////      float sum = 0;
////      float vmax = -INFINITY;
////      float vmin = INFINITY;
////      for(int i = 0;i<l;i++){
////         sum += boundary[i].v;
////         vmax = fmax(boundary[i].v,vmax);
////         vmin = fmin(boundary[i].v,vmin);
////
////         float err =  eval_model(3, &(data[i*3]), model) - data[i*3+2];
////         maxerr = fmax(maxerr, fabs(err));
////      }
////      sum/=l;
////
////      // interpolate the region if necessary
////      if(LR_REVERSE()*model[0] < 1.0 )      //      if(model[0] > -1.0 )
////      if(maxerr< MAXERR() && vmax-vmin< MAXDIFF() ) {
////         for(int i = 0;i<pixels.size();i++) {
////            int idx = pixels[i].x + pixels[i].y*w;
////            float vv[3] = {pixels[i].x, pixels[i].y, 0};
////            out[idx] = eval_model(3, vv, model);
////            outmsk[idx] = 254;
////         }
////         //      printf("id:%d sz:%d per:%d vmax-vmin:%f avg:%f maxerr:%f\n", 
////         //            cc, (int) pixels.size(), (int)boundary.size(), vmax-vmin, sum, maxerr);
////      }
////   }
//
//	/* write the results */
//   iio_save_image_float_split(out_file, out,w,h,nch);
//   if(outmsk_file)
//      iio_save_image_float_split(outmsk_file, outmsk,w,h,nch);
//   free(out);
//   free(outmsk);
//   free(msk);
//   free(img);
//}
//

// /***************************/
// 
// // c: pointer to original argc
// // v: pointer to original argv
// // o: option name after hyphen
// // d: default value (if NULL, the option takes no argument)
// static char *pick_option(int *c, char ***v, char *o, char *d)
// {
//    int argc = *c;
//    char **argv = *v;
//    int id = d ? 1 : 0;
//    for (int i = 0; i < argc - id; i++)
//       if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, o)) {
//          char *r = argv[i + id] + 1 - id;
//          *c -= id + 1;
//          for (int j = i; j < argc - id; j++)
//             (*v)[j] = (*v)[j + id + 1];
//          return r;
//       }
//    return d;
// }
// extern "C"{
// #include "iio.h"
// }
// 
// 
// int main(int argc, char **argv)
// {
// 
//    if (argc < 2) {
//       fprintf(stderr, "[MAXERR=%.2f MAXDIFF=%.2f LR_REVERSE={1(LR),0(disabled),-1(RL)}] %s mask out [minsz]\n", MAXERR(), MAXDIFF(), argv[0]);
//       return 1;
//    }
// 
//    char *img_file = argv[1];
//    char *out_file = argv[2];
//    int minsz = argc>2 ? atoi(argv[3]): 25;
// 
//    int w,h,nch;
// 
//    float *img = iio_read_image_float_split(img_file, &w, &h, &nch);
//    float *out = (float*) calloc(w*h*nch, sizeof(*out));
//    nch=1;
// 
//    
//    int n_cc = cc_grain(w, h, img, out, minsz); 
//    printf("%d\n", n_cc);
// 
// //   // initialize the image to 
// //   std::vector< int > rep(w*h); 
// //   for (int i=0;i<w*h;i++) 
// //      rep[i] = msk[i] > 0 ? -1 : 1;
// //
// //   // identify the connected components of positive values
// //   positive_connected_component_filter(&(rep.front()), w, h);
// //
// //   // count connected components
// //   int n_cc = 0;
// //   std::vector< int > cc_ids(w*h); 
// //   for (int i=0;i<w*h;i++) {
// //      if (rep[i] == i) {
// //         cc_ids[n_cc] = i;
// //         n_cc += 1;
// //      }
// //   }
// //
// //   // store the cc_pixels  and cc_boundary in vectors of vectors
// //   std::vector< std::vector< struct pix > > cc_pixels(w*h);
// //   std::vector< std::vector< struct pix > > cc_boundary(w*h);
// //
// //   for (int j = 0; j < h; j++)
// //      for (int i = 0; i < w; i++)
// //      {
// //         int idx = j*w + i;
// //         int cc = rep[idx];
// //         if (cc>=0) { // if it is a region 
// //            struct pix pp;
// //            pp.x=i;
// //            pp.y=j;
// //            pp.v=img[idx];
// //            cc_pixels[cc].push_back(pp);
// //            extract_cc_boundary_at_pix(w, h, i, j, &(rep.front()), img, &(cc_boundary[cc]));
// //         }
// //      }
// //   printf("%d\n", n_cc);
// //
// //   for (int i=0;i<w*h;i++) {
// //      out[i] = img[i];
// //      outmsk[i] = msk[i];
// //   }
// //
// //   // for each connected component
// //   std::vector< float > data(3*w*h);
// //   for(int t =0;t < n_cc; t++) 
// //   {
// //      int cc = cc_ids[t];
// //      std::vector< struct pix > pixels   =  cc_pixels[cc];
// //      std::vector< struct pix > boundary =  cc_boundary[cc];
// //      int l=boundary.size();
// //
// //      // copy the data at the boundary and compute the affine model
// //      double model[3];
// //      for(int i = 0;i<l;i++) {
// //         data[i*3+0]=(float)boundary[i].x;
// //         data[i*3+1]=(float)boundary[i].y;
// //         data[i*3+2]=(float)boundary[i].v;
// //      }
// //      generate_model ( 3, &(data.front()), l, model);
// //
// //      // compute stuff
// //      float maxerr = 0;
// //      float sum = 0;
// //      float vmax = -INFINITY;
// //      float vmin = INFINITY;
// //      for(int i = 0;i<l;i++){
// //         sum += boundary[i].v;
// //         vmax = fmax(boundary[i].v,vmax);
// //         vmin = fmin(boundary[i].v,vmin);
// //
// //         float err =  eval_model(3, &(data[i*3]), model) - data[i*3+2];
// //         maxerr = fmax(maxerr, fabs(err));
// //      }
// //      sum/=l;
// //
// //      // interpolate the region if necessary
// //      if(LR_REVERSE()*model[0] < 1.0 )      //      if(model[0] > -1.0 )
// //      if(maxerr< MAXERR() && vmax-vmin< MAXDIFF() ) {
// //         for(int i = 0;i<pixels.size();i++) {
// //            int idx = pixels[i].x + pixels[i].y*w;
// //            float vv[3] = {pixels[i].x, pixels[i].y, 0};
// //            out[idx] = eval_model(3, vv, model);
// //            outmsk[idx] = 254;
// //         }
// //         //      printf("id:%d sz:%d per:%d vmax-vmin:%f avg:%f maxerr:%f\n", 
// //         //            cc, (int) pixels.size(), (int)boundary.size(), vmax-vmin, sum, maxerr);
// //      }
// //   }
// 
// 	/* write the results */
//    iio_save_image_float_split(out_file, out,w,h,nch);
//    free(out);
//    free(img);
// }
// 
// 
