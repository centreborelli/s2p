#ifndef IMG_H_
#define IMG_H_
#include <vector>
#include <assert.h>

struct Img
{
   std::vector<float > data;
   union{
      int sz[3];
      struct{
         union { int ncol; int nx; };
         union { int nrow; int ny; };
         int nch;
      };
   };
//   int nx;
//   int ny;
//   int nch;
   int npix;


	Img(int nx, int ny, int nch=1);
	Img(float *copydata, int nx, int ny, int nch=1);

	inline Img() {nx=0;ny=0;nch=0;npix=0;}

   inline float operator[](int i) const { assert(i>=0 && i < npix*nch); return data[i];}
   inline float& operator[](int i) { assert(i>=0 && i < npix*nch); return data[i];}
   inline float operator()(int i) const { assert(i>=0 && i < npix*nch); return data[i];}
   inline float& operator()(int i) { assert(i>=0 && i < npix*nch); return data[i];}
   inline float operator()(int x, int y, int c = 0) const { int i=x+y*nx+c*npix; assert(i>=0 && i < npix*nch); return data[i];}
   inline float& operator()(int x, int y, int c = 0) { int i=x+y*nx+c*npix; assert(i>=0 && i < npix*nch); return data[i];}

   inline float& val(int i, int j, int c) { 
      assert(i >= 0 && i < nx && 
             j >= 0 && j < ny &&
             c >= 0 && c < nch  ); 
      return data[i + j*nx + c*nx*ny];
   }

   //inline float val(int x, int y, int c) const { return data[x+y*nx+c*nx*ny];} 
   inline float val(int i, int j, int c) const { 
      assert(i >= 0 && i < nx && 
             j >= 0 && j < ny &&
             c >= 0 && c < nch  ); 
      return data[i + j*nx + c*nx*ny];
   }


//   private:
//   Img(const Img&);      // disable copy constructor
//   void operator=(const Img&);
//	  Img& operator= (const Img&p);

};

#endif /* IMG_H_ */
