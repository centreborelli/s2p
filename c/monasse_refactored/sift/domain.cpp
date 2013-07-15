// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// David Lowe  "Method and apparatus for identifying scale invariant 
// features in an image and use of same for locating an object in an 
// image",  U.S. Patent 6,711,293.
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.

#include "domain.h"
#include "library.h"
#include "splines.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* extract image value (even outside image domain) */
static float v(float *in,int x,int y,float bg, int width, int height)
{
  if (x<0 || x>=width || y<0 || y>=height)
    return(bg); else return(in[y*width+x]);
}

void apply_zoom(float *input, float *out, float zoom, int order, int width, int height)
{

	int nwidth = (int)( zoom * (float) width);
	int nheight = (int)( zoom * (float) height);

	float *coeffs;
	float *ref;

	float  cx[12],cy[12],ak[13];

	if (order!=0 && order!=1 && order!=-3 && 
      order!=3 && order!=5 && order!=7 && order!=9 && order!=11)
    	{	
		printf("unrecognized interpolation order.\n");
		exit(-1);
	}

        if (order>=3) {

		coeffs = new float[width*height];
    		finvspline(input,order,coeffs,width,height);
    		ref = coeffs;
    		if (order>3) init_splinen(ak,order);

	} else 
	{
    		coeffs = NULL;
    		ref = input;
  	}

	int xi,yi;
	float xp,yp;
	float res;
	int n1,n2;
	float bg = 0.0f;
	float p=-0.5;
	for(int i=0; i < nwidth; i++)
		for(int j=0; j < nheight; j++)
		{

			xp =  (float) i / zoom;
			yp =  (float) j / zoom;

			if (order == 0) { 
	
				xi = (int)floor((double)xp); 
				yi = (int)floor((double)yp);
		
				if (xi<0 || xi>=width || yi<0 || yi>=height)
		 			 res = bg; 
				else res = input[yi*width+xi];
	
     			 } else { 
	
		
				if (xp<0. || xp>=(float)width || yp<0. || yp>=(float)height) res=bg; 
				else {
					xp -= 0.5; yp -= 0.5;
					int xi = (int)floor((double)xp); 
					int yi = (int)floor((double)yp);
					float ux = xp-(float)xi;
	  				float uy = yp-(float)yi;

					switch (order) 
	   				{
					    	case 1: /* first order interpolation (bilinear) */
	      						n2 = 1;
							cx[0]=ux;	cx[1]=1.f-ux;
							cy[0]=uy; cy[1]=1.f-uy;
							break;
						
						case -3: /* third order interpolation (bicubic Keys' function) */
							n2 = 2;
							keys(cx,ux,p);
							keys(cy,uy,p);
							break;

						case 3: /* spline of order 3 */
							n2 = 2;
							spline3(cx,ux);
							spline3(cy,uy);
							break;

						default: /* spline of order >3 */
							n2 = (1+order)/2;
							splinen(cx,ux,ak,order);
							splinen(cy,uy,ak,order);
							break;
					}
	  
	  				res = 0.; n1 = 1-n2;
	 				if (xi+n1>=0 && xi+n2<width && yi+n1>=0 && yi+n2<height) {
	    				
						int adr = yi*width+xi; 
	    					for (int dy=n1;dy<=n2;dy++) 
	    					  for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*ref[adr+width*dy+dx];
	  				} else 
	
					   	for (int dy=n1;dy<=n2;dy++)
						      for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*v(ref,xi+dx,yi+dy,bg,width,height);
					
      			}

		}	

		out[j*nwidth+i] = res;
    	
	}
    if(ref != input)
        delete [] ref;
}
