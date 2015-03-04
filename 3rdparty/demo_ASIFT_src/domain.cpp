// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "domain.h"


#define DEBUG 0




void apply_zoom(float *input, float *out, float zoom, int order, int width, int height)
{

	int nwidth = (int)( zoom * (float) width);
	int nheight = (int)( zoom * (float) height);

	float *coeffs;
	float *ref;

	float  cx[12],cy[12],ak[13];

	// Guoshen Yu, 2010.09.22, Windows versions
	vector<float> input_vec, coeffs_vec, ref_vec;
	input_vec = vector<float>(width*height);
	coeffs_vec = vector<float>(width*height);
	ref_vec = vector<float>(width*height);
	for (int i = 0; i < width*height; i++)
		input_vec[i] = input[i];

	if (order!=0 && order!=1 && order!=-3 && 
      order!=3 && order!=5 && order!=7 && order!=9 && order!=11)
    	{	
		printf("unrecognized interpolation order.\n");
		exit(-1);
	}

        if (order>=3) {

			coeffs = new float[width*height];

			// Guoshen Yu, 2010.09.21, Windows version
    		//finvspline(input,order,coeffs,width,height);
			finvspline(input_vec,order,coeffs_vec,width,height);
			for (int i = 0; i < width*height; i++)
				coeffs[i] = coeffs_vec[i];

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
							cx[0]=ux;	cx[1]=1.-ux;
							cy[0]=uy; cy[1]=1.-uy;
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
	
						// Guoshen Yu, 2010.09.21, Windows
						for (int i = 0; i < width*height; i++)
							ref_vec[i] = ref[i];

					   	for (int dy=n1;dy<=n2;dy++)
						      for (int dx=n1;dx<=n2;dx++) 
							// Guoshen Yu, 2010.09.21, Windows
							// res += cy[n2-dy]*cx[n2-dx]*v(ref,xi+dx,yi+dy,bg,width,height);
							res += cy[n2-dy]*cx[n2-dx]*v(ref_vec,xi+dx,yi+dy,bg,width,height);
      			}

		}	

		out[j*nwidth+i] = res;
    	
	}

}

