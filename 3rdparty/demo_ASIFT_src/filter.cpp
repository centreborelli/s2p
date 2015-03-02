// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "filter.h"


/////////////////////////////////////////////////////////////// Build Gaussian filters
float * directional_gauss_filter(float xsigma, float ysigma, float angle, int *kwidth, int *kheight)
{
	
	
	int ksize = (int)(2.0 * 2.0 * MAX(xsigma, ysigma) + 1.0);
	float *kernel = new float[ksize*ksize];
	
	float xsigma2 = xsigma*xsigma;
	float ysigma2 = ysigma*ysigma;
	
	int l2 = ksize/2;
	for(int y = -l2; y <= l2; y++)
		for(int x = -l2; x <= l2; x++)
		{
			
			float a = (float) angle * PI / 180.0f;
			float sina = sin(a);
			float cosa = cos(a);
			
			float ax = (float) x * cosa + (float) y * sina;
			float ay = -(float) x * sina + (float) y * cosa;
			kernel[(y+l2) * ksize + x + l2] =  exp(-(ax*ax)/(2.0f*xsigma2)  - (ay*ay)/(2.0f*ysigma2) );  
			
		}
	
	
	float sum=0.0;
	for(int i=0; i < ksize*ksize; i++) sum += kernel[i];
	for(int i=0; i < ksize*ksize; i++) kernel[i] /= sum;
	
	*kwidth = ksize;
	*kheight = ksize;
	
	return kernel;
}






/* Convolution with a kernel */
/* No padding applied to the image */
void convol(float *u,float *v,int width,int height,float *kernel,int kwidth,int kheight)
{
  
//  float S;
//  int K2,L2,m,n,kmin,kmax,lmin,lmax,l,k;
 
  int K2 = kwidth / 2;
  int L2 = kheight / 2;

  for(int y=0 ; y < height; y++) 
    for (int x=0 ; x < width; x++) {

	float S = 0.0;
//      kmax = MIN(kwidth-1,n+K2);
//      kmin = MAX(0,1+n+K2-width);
//      lmax = MIN(kheight-1,m+L2);
//      lmin = MAX(0,1+m+L2-height);
	
	for (int l = -L2; l <= L2; l++) 
		for (int k = -K2 ; k<= K2; k++)
		{ 
			int px=x+k;
			int py=y+l;

			if (px>=0 && px < width && py>=0 && py<height)
			  S += u[width*py + px] * kernel[kwidth*(l+L2) + k+K2];
		}

      v[y*width+x] = (float) S;

    }
}



void median(float *u,float *v, float radius, int niter, int width,int height)
{
    

	int iradius = (int)(radius+1.0);
	int rsize=(2*iradius+1)*(2*iradius+1);
  
	float * vector = new float[rsize];
	float * index = new float[rsize];
  
	for(int n=0; n< niter;n++){
		
		for(int x=0;x<width;x++)
			for(int y=0;y<height;y++){
      
				int count=0;
				for(int i=-iradius;i<=iradius;i++)
					for(int j=-iradius;j<=iradius;j++)
						if ((float) (i*i + j*j) <= iradius*iradius){
		
							int x0=x+i;
							int y0=y+j;

							if (x0>=0 && y0>=0 && x0 < width && y0 < height) { 
								vector[count] = u[y0*width+x0];
								index[count] = count;
								count++; 
							}
				}

				quick_sort(vector,index,count);
				v[y*width+x] = vector[count/2];

		}

		copy(v,u,width*height);
	}

	delete[] vector;
	delete[] index;
  
}



void remove_outliers(float *igray,float *ogray,int width, int height)
{
  

	int bloc=1;
	int bsize = (2*bloc+1)*(2*bloc+1)-1;
	for(int x=bloc;x<width-bloc;x++)
		for(int y=bloc;y<height-bloc;y++) {
	  
			int l = y*width+x;

			int countmax=0;
			int countmin=0;
			float valueg0 = igray[l];
    
	//		float distmin = MAXFLOAT;
			float distmin = FLT_MAX; // Guoshen Yu
			float green = igray[l];

			for(int i=-bloc;i<=bloc;i++)
				for(int j=-bloc;j<=bloc;j++)
					if ((i!=0 || j!=0)){
		
					int l0 = (y+j)*width+x+i;
						
					int valueg = (int) igray[l0];
					
					if (valueg0>valueg) countmax++;
					if (valueg0<valueg) countmin++;
						
						
					float dist = fabsf(valueg - valueg0);
			
					if (dist < distmin) {distmin=dist;green=valueg;}
							
			}
	  
	  
			if (countmin == bsize || countmax == bsize ) 	ogray[l]=green;
			else ogray[l] = igray[l];
	
	  
		}

 
}



/* Convolution with a  separable kernel */
/* boundary condition: 0=zero, 1=symmetry */
void separable_convolution(float *u, float *v, int width, int height,float * xkernel, int xsize,float *ykernel,int ysize,int boundary)
{

  int width2 = 2*width;
  int height2 = 2*height;

  float *tmp = (float *) malloc(width*height*sizeof(float));
 

  /* convolution along x axis */
  float sum = 0.0;
  int org = xsize / 2;
  for (int y=height;y--;) 
    for (int x=width;x--;) {
	
      sum = 0.0;	
      for (int i=xsize;i--;) {
	int s = x-i+org;
	switch(boundary) {
	  
	case 0: 
	  if (s>=0 && s<width) sum += xkernel[i]*u[y*width+s];
	  break;
      
	case 1: 

	  while (s<0) s+=width2;
	  while (s>=width2) s-=width2;
	  if (s>=width) s = width2-1-s;
	  sum += xkernel[i]*u[y*width+s];
	  break;

	}
      }
      tmp[y*width+x] = sum;
    }
  
  /* convolution along y axis */
  org = ysize / 2;
  for (int y=height;y--;) 
    for (int x=width;x--;) {

      sum=0.0;
      for (int i=ysize;i--;) {
	int s = y-i+org;
	switch(boundary) {
	case 0: 
	  if (s>=0 && s<height) sum += ykernel[i]*tmp[s*width+x];
	  break;
	case 1: 
	  while (s<0) s+=height2;
	  while (s>=height2) s-=height2;
	  if (s>=height) s = height2-1-s;
	  sum += ykernel[i]*tmp[s*width+x];
	  break;
	}
      }
      v[y*width+x] = sum;
    }
  
  free(tmp);
}


void gaussian_convolution(float *u, float *v, int width, int height, float sigma)
{

	int ksize;	
	float * kernel;

	ksize = (int)(2.0 * 4.0 * sigma + 1.0);
	kernel = gauss(1,sigma,&ksize);

	int boundary = 1;

	copy(u,v,width*height);
	horizontal_convolution(v, v, width, height, kernel, ksize, boundary);
    vertical_convolution(v, v, width, height,  kernel,  ksize, boundary);
	delete[] kernel; /*memcheck*/
}


void gaussian_convolution(float *u, float *v, int width, int height, float sigma, int ksize)
{
	float * kernel;
	kernel = gauss(1,sigma,&ksize);

	int boundary = 1;

	copy(u,v,width*height);
	horizontal_convolution(v, v, width, height, kernel, ksize, boundary);
    	vertical_convolution(v, v, width, height,  kernel,  ksize, boundary);
}


void fast_separable_convolution(float *u, float *v, int width, int height,float * xkernel, int xsize,float *ykernel,int ysize,int boundary)
{
    copy(u,v,width*height);

    horizontal_convolution(v, v, width, height, xkernel, xsize, boundary);
    vertical_convolution(v, v, width, height,  ykernel,  ysize, boundary);

}

/* Loop unrolling simply sums 5 multiplications 
   at a time to allow the compiler to schedule
   operations better and avoid loop overhead.
*/
void buffer_convolution(float *buffer,float *kernel,int size,int ksize)
{

    for (int i = 0; i < size; i++) {

      float sum = 0.0;
      float *bp = &buffer[i];
      float *kp = &kernel[0];
      
		
      /* Loop unrolling: do 5 multiplications at a time. */
//      int k=0;
		
	  for(int k = 0; k < ksize; k++) 
		 sum += *bp++ * *kp++;
		
  //    for(;k + 4 < ksize;  bp += 5, kp += 5, k += 5) 
//	        sum += bp[0] * kp[0] +  bp[1] * kp[1] + bp[2] * kp[2] +
  //      	       bp[3] * kp[3] +  bp[4] * kp[4];

      /* Do multiplications at a time on remaining items. */
//      for(; k < ksize; bp++ , kp++, k++)  sum += *bp * (*kp);

      buffer[i] = sum;
    }
}



/* Convolve image with the 1-D kernel vector along image rows.  This
   is designed to be as efficient as possible.  
*/
void horizontal_convolution(float *u, float *v, int width, int height, float *kernel, int ksize, int boundary)
{

    int halfsize = ksize / 2;
    int buffersize = width + ksize;
    float *buffer = new float[buffersize];

    for (int r = 0; r < height; r++) {

	/// symmetry
	int l = r*width;
	if (boundary == 1)
        	for (int i = 0; i < halfsize; i++)
            		buffer[i] = u[l + halfsize - 1 - i ];
	else
		for (int i = 0; i < halfsize; i++)
            		buffer[i] = 0.0;


        for (int i = 0; i < width; i++)
            buffer[halfsize + i] = u[l + i];


	if (boundary == 1)
	        for (int i = 0; i <  halfsize; i++)
        	    buffer[i + width + halfsize] = u[l + width - 1 - i];
	else 
		for (int i = 0; i <  halfsize; i++)
        	    buffer[i + width + halfsize] = 0.0;

        buffer_convolution(buffer, kernel, width, ksize);
        for (int c = 0; c < width; c++)
          v[r*width+c] = buffer[c];
    }
	delete[] buffer; /*memcheck*/
}



void vertical_convolution(float *u, float *v, int width, int height, float *kernel,int ksize, int boundary)
{
    int halfsize = ksize / 2;
    int buffersize = height + ksize;
    float *buffer = new float[buffersize];

    for (int c = 0; c < width; c++) {

	if (boundary == 1)
		for (int i = 0; i < halfsize; i++)
				buffer[i] = u[(halfsize-i-1)*width + c];
	else 
		for (int i = 0; i < halfsize; i++)
				buffer[i] = 0.0f;

      for (int i = 0; i < height; i++)
		buffer[halfsize + i] = u[i*width + c];

	if (boundary == 1)
		for (int i = 0; i < halfsize; i++)
				buffer[halfsize + height + i] = u[(height - i - 1)*width+c];
	else
		for (int i = 0; i < halfsize; i++)
			buffer[halfsize + height + i] = 0.0f;

      buffer_convolution(buffer, kernel, height, ksize);

      for (int r = 0; r < height; r++)
          v[r*width+c] = buffer[r];

	}
    delete[] buffer; /*memcheck*/
}



void heat(float *input, float *out, float step, int niter, float sigma, int width, int height)
{

	int i,j,n,ksize,size,im,i1,j1,jm;
	float *kernel = NULL, *laplacian = NULL, *convolved = NULL;
	
	
	size = width*height;
	
	if (sigma > 0.0) kernel = gauss(0,sigma,&ksize);		  		
		
	laplacian = (float *) malloc(size*sizeof(float));
	convolved = (float *) malloc(size*sizeof(float));
	
	
	
	for(n=0; n < niter; n++)
	{
	
		
		if (sigma > 0.0)
		{
			
			separable_convolution(input,convolved,width,height, kernel, ksize,kernel,ksize,1);
			
			for(i=0; i< size; i++) laplacian[i] = convolved[i] - input[i];
		
		} else
		{
		
		
			for (i=0; i < width;i++)
			  for (j=0; j< height ;j++) 
				{
				
				  if (j==0) jm=1; else jm=j-1;
				  if (j==height-1) j1=height-2; else j1=j+1;
				
				  if (i==0) im=1; else im=i-1;
				  if (i==width-1) i1=width-2; else i1=i+1;
				
				  laplacian[j*width + i] =  - 4.0  * input[width*j+i] + input[width*j+im]+ input[width*j+i1]+input[width*jm + i] + input[width*j1 + i];
			     }
		}
		
		
		
		for(i=0; i < size; i++) out[i] = input[i] + step * laplacian[i];
		
		copy(out,input,size);
		
	}
	
	
	free(laplacian);
	free(convolved);
	if (kernel) free(kernel);
	
}


		

