// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "flimage.h"



//////////////////////////////////////////////// Class flimage
//// Construction
flimage::flimage() : width(0), height(0), p(0) 
{
}	 

flimage::flimage(int w, int h) : width(w), height(h), p(new float[w*h]) 
{
	for (int j=width*height-1; j>=0 ; j--) p[j] = 0.0;
}	


flimage::flimage(int w, int h, float v) : width(w), height(h), p(new float[w*h]) 
{
	for (int j=width*height-1; j>=0 ; j--) p[j] = v;
}


flimage::flimage(int w, int h, float* v) : width(w), height(h), p(new float[w*h]) 
{
	for (int j=width*height-1; j>=0 ; j--) p[j] = v[j];
}


void flimage::create(int w, int h)
{ 
 	erase();
	width = w; height = h; 
	p = new float[w*h];	
	for (int j=width*height-1; j>=0 ; j--) p[j] = 0.0;
}

void flimage::create(int w, int h, float* v)
{
 	erase();
	width = w; height = h;  p = new float[w*h]; 
	for (int j=width*height-1; j>=0 ; j--) p[j] = v[j];
}


flimage::flimage(const flimage& im) : width(im.width), height(im.height), p(new float[im.width*im.height]) 
{
	for (int j=width*height-1; j>=0 ; j--) p[j] = im.p[j];
}

flimage&  flimage::operator= (const flimage& im)
{	
	if (&im == this) {
		return *this;
	}
	
	if (width != im.width || height != im.height)
	{  			
	  	erase();
		width = im.width; height=im.height; p = new float[width*height];
	}
	
	for (int j=width*height-1; j>=0 ; j--) p[j] = im.p[j];
	return *this;	
}


//// Destruction
void flimage::erase() 
{
	width = height = 0; 
	if (p) delete[] p; 
	p=0;
} 

flimage::~flimage()
{
	erase();
}




