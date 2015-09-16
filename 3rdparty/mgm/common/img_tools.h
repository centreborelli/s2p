//// a structure to wrap images 
#ifndef IMG_TOOLS_H_
#define IMG_TOOLS_H_

#include "img.h"
#include <algorithm>
#include "point.h"
extern "C" {
#include "iio.h"
}

/************ IMG IO  **************/

struct Img iio_read_vector_split(char *nm)
{
	struct Img out;
	float *tmpout = iio_read_image_float_split(nm, &out.nx, &out.ny, &out.nch);
	out.data.assign(tmpout,tmpout + out.nx * out.ny * out.nch);
	out.npix = out.nx * out.ny;
	free (tmpout);
	return out;
}



void iio_write_vector_split(char *nm, struct Img &out)
{
	// .front() -> .data() in C++11
	iio_save_image_float_split(nm, &(out.data.front()), out.nx, out.ny, out.nch);
}


void remove_nonfinite_values_Img(struct Img &u, float newval) 
{
   for(int i=0;i<u.npix*u.nch;i++) 
      if (!isfinite(u[i])) u[i] = newval; 
}


/************ IMG ACCESS **************/

inline float val(const struct Img &u, const Point p, const int ch=0) {
	int x = p.x;
	int y = p.y;
	return u.data[x+u.nx*y+ ch*u.npix];
}

inline int check_inside_image(const Point p, const struct Img &u) {
	int nx = u.nx;
	int ny = u.ny;
	float x = p.x;
	float y = p.y;
	if(x>=0 && y>=0 && x<nx && y<ny) return 1;
	else return 0;
}

inline float valnan(const struct Img &u, const Point p, const int ch=0)
{
	return check_inside_image(p, u) ? val(u, p, ch) : NAN;
}

inline float valzero(const struct Img &u, const Point p, const int ch=0)
{
	return check_inside_image(p, u) ? u.val(p.x, p.y, ch) : 0;
}

inline float valzero(const struct Img &u, const int x, const int y, const int ch=0)
{
	return check_inside_image(Point(x,y), u) ? u.val(x,y,ch) : 0;
}

inline float valneumann(const struct Img &u, const int x, const int y, const int ch=0)
{  
   int xx=x, yy=y;
   xx = x >=  0  ? xx : 0;
   xx = x < u.nx ? xx : u.nx - 1;
   yy = y >=  0  ? yy : 0;
   yy = y < u.ny ? yy : u.ny - 1;
	return u.val(xx,yy,ch);
}

/************ IMG PROC **************/

struct Img compute_insensity_image(struct Img &u) {
   int nx = u.nx;
   int ny = u.ny;
   int nch= u.nch;
   struct Img Intensity(nx,ny);
   for(int i=0;i<nx*ny;i++) Intensity[i]=0;

   for (int c=0;c<nch;c++)
   for(int j=0;j<ny;j++)
   for(int i=0;i<nx;i++)
      Intensity[i+j*nx] += u.val(i,j,c);

   for(int i=0;i<nx*ny;i++) Intensity[i]/=nch;

   return Intensity;
}

struct Img apply_filter(struct Img &u, struct Img &filter) {
   struct Img fu(u);
   int hfnx  = filter.nx  / 2;
   int hfny  = filter.ny  / 2;
   int hfnch = filter.nch / 2;

   for(int c = 0; c < u.nch; c++)
   for(int j = 0; j < u.ny ; j++)
   for(int i = 0; i < u.nx ; i++) {
      float v = 0;
      for (int cc = 0; cc < filter.nch; cc++)
      for (int jj = 0; jj < filter.ny ; jj++)
      for (int ii = 0; ii < filter.nx ; ii++) {
         v += valneumann(u, i + ii - hfnx,
                            j + jj - hfny,
                            c + cc - hfnch) *  
              filter.val(ii, jj, cc);
      }
      fu.val(i,j,c) = v;
   }

   return fu;
}

struct Img apply_filter(struct Img &u, float ff[], int fnx, int fny, int fnc) {
   struct Img f(fnx,fny,fnc);
   for(int i=0;i<fnx*fny*fnc;i++) f[i] = ff[i];
   return apply_filter(u,f);
}

//struct Img sobel_x(struct Img &u) {
//   struct Img f(3,3,1);
//   float ff[] = {-1,0,1, -1,0,1, -1,0,1};
//   for(int i=0;i<9;i++) f[i] = ff[i];
//   return apply_filter(u,f);
//}



static float unnormalized_gaussian_function(float sigma, float x) {
	return exp(-x*x/(2*sigma*sigma));
}

#define KWMAX 39
static int gaussian_kernel_width(float sigma) {
	float radius = 3 * fabs(sigma);
	int r = ceil(1 + 2*radius);
	if (r < 1) r = 1;
	if (r > KWMAX) r = KWMAX;
	return r;
}

static void fill_gaussian_kernel(float *k, int w, int h, float s) {
	int cw = (w - 1)/2;
	int ch = (h - 1)/2;
   float m = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
      float v  = unnormalized_gaussian_function(s, hypot(i-cw,j-ch));
		k[j*w+i] = v;
      m       += v;
   }
	for (int i = 0; i < w*h; i++)  {
      k[i] /= m;
   }
}

struct Img gblur_truncated(struct Img &u, float sigma) {
   // determine the size of the kernel
   int rad = gaussian_kernel_width(sigma);
   struct Img fx(rad,1,1), fy(1,rad,1);
   fill_gaussian_kernel(&(fx[0]), rad, 1, sigma);
   fill_gaussian_kernel(&(fy[0]), 1, rad, sigma);
   struct Img tmp = apply_filter(u, fx);
   return apply_filter(tmp, fy);
}


std::pair<float, float> image_minmax(struct Img &u){
   // returns global (finite) min and max of an image
   int nx = u.nx;
   int ny = u.ny;
   int nch= u.nch;
   float gmin = INFINITY; float gmax = -INFINITY;
   for (int c=0;c<nch;c++)
   for (int j=0;j<ny;j++)
   for (int i=0;i<nx;i++) {
      float v = val(u,Point(i,j), c);
      if (isfinite(v)) { 
         if (v < gmin) gmin = v;   
         if (v > gmax) gmax = v;
      }
   }

   return std::pair<float, float> (gmin, gmax);
}

/// Median filter
struct Img median_filter(struct Img &u, int radius) {
    struct Img M(u);
    int size=2*radius+1;
    size *= size;
    std::vector<float> v(size);
    for(int k=0; k<M.nch; k++) 
    for(int y=0; y<M.ny; y++)
    for(int x=0; x<M.nx; x++)
    {
       int n=0;
       for(int j=-radius; j<=radius; j++)
          if(0<=j+y && j+y<u.ny)
             for(int i=-radius; i<=radius; i++)
                if(0<=i+x && i+x<u.nx)
                   v[n++] = M.val(i+x,j+y,k);
       std::nth_element(v.begin(), v.begin()+n/2, v.end());
       M.val(x,y,k) = v[n/2];
    }
    return M;
}

#endif // IMG_TOOLS_H_
