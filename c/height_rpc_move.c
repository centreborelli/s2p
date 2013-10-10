#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "vvector.h"
#include "iio.h"
#include "rpc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#define max(a,b) (((a)>(b))?(a):(b))

#include "read_matrix.c"

int writeMatrix_FILE(double H[3][3], FILE *f) {
  fprintf(f, "[");
  fprintf(f, "%20f %20f %20f ;", H[0][0],  H[0][1],  H[0][2]);
  fprintf(f, "%20f %20f %20f ;", H[1][0],  H[1][1],  H[1][2]);
  fprintf(f, "%20f %20f %20f ;", H[2][0],  H[2][1],  H[2][2]);
  fprintf(f, "]\n");
  return 1;
}

int printMatrix(double H[3][3]) {
  FILE *f = stdout;
  writeMatrix_FILE(H, f);
  return 1;
}

int writeMatrix(double H[3][3] , char *file) {
   FILE *f;
   if(f = fopen(file,"w")){
      writeMatrix_FILE(H, f);
      fclose(f);
      return 1;
   }
   else return 0;
}


/// Find endianness of the system
static char * endian() {
   int i=0;
   *((unsigned char*)&i)=1;
   return (i==1? "little_endian": "big_endian");
}


FILE* writePLY_header(FILE* f,int npoints,int ntriangles) {
   int BINARY=0;   //# UNSUPPORTED
   fprintf(f,"ply\n");
   if(BINARY) {
      fprintf(f,"format binary_%s_endian 1.0\n", endian());
   }
   else {
      fprintf(f,"format ascii 1.0\n");
      fprintf(f,"comment created by RPC\n");
      fprintf(f,"element vertex %d \n",npoints);
      fprintf(f,"property float x\n");
      fprintf(f,"property float y\n");
      fprintf(f,"property float z\n");
      fprintf(f,"property uchar red\n");
      fprintf(f,"property uchar green\n");
      fprintf(f,"property uchar blue\n");
      if (ntriangles >0) {
         fprintf(f,"element face %d \n", ntriangles);
         fprintf(f,"property list uchar int vertex_index\n");
      }
      fprintf(f,"end_header\n");
   }
   return f;
}

void applyHom(double outv[3], double M[3][3], double v[3]) {
    MAT_DOT_VEC_3X3(outv,M,v);
    outv[0]/=outv[2]; outv[1]/=outv[2]; outv[2]/=outv[2];
}




static int main_height_rpc_move(int c, char *v[])
{
	if (c < 8) {
		fprintf(stderr, "usage:\n\t"
            "transfer the height map A (and mask A) from the image geometry given by rpca and Ha to the\n\t"
            "geometry given by rpcb and Hb. Optionally indicate the output image size setting NX and NY\n\t"
				"%s rpca  Ha heightA mskA  rpcb Hb  outheightB outmskB  NX  NY"
				//0    1   2   3      4    5      6       7        8     9  10
				"\n", *v);
		return EXIT_FAILURE;
	}



   //read the parameters
   struct rpc rpca[1];
   struct rpc rpcb[1];
   double Ha[3][3], Hb[3][3];
   float *heightA, *mskA;
   int nx,ny,nch;
   char *f_outheightB=NULL;
   char *f_outmskB=NULL;
   int outnx,outny;

   int i = 1;
   read_rpc_file_xml(rpca, v[i]);                              i++; //1
   read_matrix(Ha,v[i]);                                       i++; //2
   heightA = iio_read_image_float_split(v[i], &nx, &ny, &nch); i++; //3
   mskA = iio_read_image_float_split(v[i], &nx, &ny, &nch);    i++; //4
   read_rpc_file_xml(rpcb, v[i]);                              i++; //5
   read_matrix(Hb,v[i]);                                       i++; //6
   f_outheightB = v[i];                                        i++; //7
   f_outmskB    = v[i];                                        i++; //8
   outnx = (c>i)? atoi(v[i]): nx;                              i++; //9
   outny = (c>i)? atoi(v[i]): ny;                              i++; //10



   float *outmskB = calloc(outnx*outny,sizeof(*outmskB));
   float *outheightB = calloc(outnx*outny,sizeof(*outheightB));

   double det;
   double invHa[3][3];
   double invHb[3][3];
   INVERT_3X3(invHa, det, Ha);
   INVERT_3X3(invHb, det, Hb);



   int x,y;
   for(y=0;y<outny*outnx;y++){
      outmskB[y] = 0;   // initialize (non visible)
      outheightB[y] = -1./0.; // initialize (infinity)
   }
   for(y=0;y<ny;y++){
      for(x=0;x<nx;x++){
         int pos = x + nx*y;
         if (  mskA[pos] > 0) {
            double hA = heightA[pos];
            double p0[3] = {x,y,1};
            double p1[3], q0[3], q1[3];

            applyHom(q0,invHa,p0);

	         double tmpGround[2];
	         eval_rpc(tmpGround, rpca, q0[0], q0[1] , hA);
	         eval_rpci(q1, rpcb, tmpGround[0], tmpGround[1], hA);

            q1[2]=1; // fix homogeneous
            applyHom(p1,Hb,q1);

            p1[0] = round(p1[0]);
            p1[1] = round(p1[1]);
            if(p1[0]>=0 && p1[0]<outnx && p1[1]>=0 && p1[1]<outny){
               int pos_p1 = (int)(p1[0] + outnx*p1[1]);
               outheightB[pos_p1] = max(hA , outheightB[pos_p1]);
               outmskB[pos_p1] = 1;
            }

         }
      }
   }


   // save re-projected disparity and mask
   iio_save_image_float_vec(f_outheightB, outheightB, outnx,outny, 1);
   iio_save_image_float_vec(f_outmskB, outmskB, outnx,outny, 1);


	return 0;
}

int main(int c, char *v[])
{
	return main_height_rpc_move(c, v);
}
