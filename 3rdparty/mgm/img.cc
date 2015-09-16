#include "img.h"
#include <vector>
#include <stdio.h>


Img::Img(int nx, int ny, int nch)
{
	this->nx = nx;
	this->ny = ny;
   this->nch = nch;
   this->npix = nx*ny;
   this->data = std::vector<float >(nx*ny*nch,0);
}

// copy row major data and create a new image
// the data is assumed to be organized in color planes
Img::Img(float *copydata, int nx, int ny, int nch)
{
	this->nx = nx;
	this->ny = ny;
   this->nch = nch;
   this->npix = nx*ny;
   this->data = std::vector<float >(nx*ny*nch,0);
   for(int i=0;i<nx*ny*nch;i++)
      this->data[i] = copydata[i];
}

//int main(){
//   struct Img c,d;
//   d.data.push_back(-1.1);
//   d.npix++;
//   struct Img a = Img(1002,1003,3);
//   a[0]=10;
//   struct Img b(1000,1000,3);
//   b[0]=20;
//   c = b;
//   c[0]++;
//   c[1]=-10;
//   printf("%f %f %f %f %f %f %d\n", a[0],b[0],c[0],c[1],a.val(10,10,1),d[0],d.nx);
//   printf("%d %d %d %d %d\n", a.nx,a.ny,a.sz[0],a.ncol, c.sz[0]);
//}
