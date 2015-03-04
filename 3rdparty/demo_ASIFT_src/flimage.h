// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef _FLIMAGE_H_
#define _FLIMAGE_H_

#include <iostream>
#include <string>

class flimage {
	
private:
	
	int	width, height;	// image size
	float*	p;		// array of color levels: level of pixel (x,y) is p[y*width+x]
	
public:
	
	
	//// Construction
	flimage();
	flimage(int w, int h);
	flimage(int w, int h, float v);
	flimage(int w, int h, float* v);
	flimage(const flimage& im);
	flimage& operator= (const flimage& im);
	
	
	void create(int w, int h);
	void create(int w, int h, float *v);
	
	//// Destruction
	void erase();
	~flimage();
	
	//// Get Basic Data	
	int nwidth() const {return width;} 	// image size
	int nheight() const {return height;} 
	
	/// Access values
	float* getPlane() {return p;}	// return the adress of the array of values 
	
	float operator()(int x, int y) const {return p[ y*width + x ];} 	// acces to the (x,y) value
	float& operator()(int x, int y) {return p[ y*width + x ];}	// by value (for const images) and by reference
	
	
};


#endif

