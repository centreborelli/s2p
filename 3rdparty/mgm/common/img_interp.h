#ifndef IMG_H_
#define IMG_H_
#include <vector>
#include <assert.h>

enum Interpolator { nearest, bilinear, bicubic, spline5, spline7, spline11 };

inline float bilinear_interpolation_cell(const float p[2][2], const float x, const float y)
{
   float r = 0;
   r += p[0][0] * (1-x) * (1-y);
   r += p[1][0] * ( x ) * (1-y);
   r += p[0][1] * (1-x) * ( y );
   r += p[1][1] * ( x ) * ( y );
   return r;
}

inline float cubic_interpolation(const float v[4], const float x)
{
   return v[1] + 0.5 * x*(v[2] - v[0]
         + x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
         + x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static float bicubic_interpolation_cell(const float p[4][4], const float x, const float y)
{
   float v[4];
   v[0] = cubic_interpolation(p[0], y);
   v[1] = cubic_interpolation(p[1], y);
   v[2] = cubic_interpolation(p[2], y);
   v[3] = cubic_interpolation(p[3], y);
   return cubic_interpolation(v, x);
}

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

//   std::vector<float > zdata;
//   Interpolator interp_type;
//   int interp_zoom; // TODO zdata may be zoomed to simplify the interpolation
//

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


// Img_subpix: prefilter on request at the first subpixel access. 
//               Img.getpixel(x, y, c), Img.getpixel_int (x,y,c). 
//               pu=interp_prefilter(u, type), 
//               float=intrep_sample (pu, x, y, c, type)

   void interp_prefilter(enum Interpolator); // TODO not implemented

   // default boundary condition: neumann
inline float getpixel_int(int x, int y, int c) const {
   if(x<0)   x=0;
   if(x>=nx) x=nx-1;
   if(y<0)   y=0;
   if(y>=ny) y=ny-1;
   if(c<0)   c=0;
   if(c>=nch)c=nch-1;
   return (*this)(x,y,c);
}
inline float getpixel(int x, int y, int c) const {
   return getpixel_int(x,y,c);
}

inline float getpixel(const float x, const float y, const int c) const {
   if((float)(int) x == x && (float)(int) y == y) 
      getpixel_int((int)x, (int)y, c);
   return interp_bilinear(x,y,c);
   return interp_bicubic(x,y,c);
//   return interp_nearest(x,y,c);
}



inline float interp_bilinear(const float x, const float y, const int ch) const
{
   int ix = floor(x);
   int iy = floor(y);

   float c[2][2];
   for (int j = 0; j < 2; j++)
      for (int i = 0; i < 2; i++)
         c[i][j] = getpixel_int(ix + i, iy + j, ch);
   return  bilinear_interpolation_cell(c, x - ix, y - iy);
}


float interp_bicubic(const float x, const float y, const int ch) const
{
   int ix = floor(x-1);
   int iy = floor(y-1);

   float c[4][4];
   for (int j = 0; j < 4; j++)
      for (int i = 0; i < 4; i++)
         c[i][j] = getpixel_int(ix + i, iy + j, ch);
   return  bicubic_interpolation_cell(c, x -1- ix, y -1- iy);
}

float interp_nearest(const float x, const float y, const int ch) const
{
   int ix = round(x);
   int iy = round(y);
   return getpixel_int(ix, iy, ch);
}


};

#endif /* IMG_H_ */
