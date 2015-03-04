// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#include "library.h"



void wxwarning(const char * message, const char *function,const char *file)
{
	
	printf("warning :: %s :: %s :: %s\n", file, function, message);
	
}


void wxerror(const char * message, const char *function, const char *file)
{
	
	printf("error :: %s :: %s :: %s\n", file, function, message);
	exit(-1);
	
}




double fsqr(double a) {	return a*a; }


void  fill_exp_lut(float *lut, int size)
{
  for(int i=0; i< size;i++)   lut[i]=expf( - (float) i / LUTPRECISION);
}

  
  
    
/* Looks for f(dif) in the lut */
float slut(float dif,float *lut)
{
  if (dif >= (float) LUTMAX) return 0.0;
  	
  int  x= (int) floor(dif*LUTPRECISION);
 	
  float y1=lut[x];
  float y2=lut[x+1];
 
  return y1 + (y2-y1)*(dif*LUTPRECISION - (float) x); 
}




float max(float *u,int *pos, int size)
{  
	float max=u[0];
	if (pos) *pos=0;  
	for(int i=1; i<size; i++)  if(u[i]>max){ max=u[i]; if (pos) *pos=i; } 
	return max;
}


float min(float *u,int *pos,int size) 
{
	float min=u[0];  
	if (pos) *pos=0;
	for(int i=1;i<size;i++)	if(u[i]<min){ min=u[i]; if (pos) *pos=i; }
	return min;
}


void max_u_v(float *u,float *v,int size) { for(int i=0;i<size;i++)  u[i]=MAX(u[i],v[i]);}

void max_u_k(float *u, float k, int size) { for(int i=0;i<size;i++)  u[i]=MAX(u[i],k);}

void min_u_v(float *u, float *v,int size) { for(int i=0;i<size;i++)  u[i]=MIN(u[i],v[i]);}

void min_u_k(float *u,float k,int size) { int i=0; for(i=0;i<size;i++)  u[i]=MIN(u[i],k);}

void abs(float *u,float *v,int size){ for(int i=0;i<size ;i++)  v[i] = fabsf(u[i]); }

void copy(float *u,float *v,int size) {  for(int i=0; i<size ;i++)  v[i]=u[i]; }

void clear(float *u,float value,int size) { for(int i=0; i<size; i++) u[i]=value;  }

void combine(float *u,float a,float *v,float b, float *w,  int size)  { for(int i=0;i<size ;i++)   w[i]= a*u[i] + b*v[i];  }

void multiple(float *u,float multiplier,int size) { for(int i=0;i<size;i++)  u[i]=multiplier*u[i]; }

float scalar_product(float *u, float *v, int n)
{
	
	float aux = 0.0f;
	for(int i=0; i < n; i++)
		aux += u[i] * v[i];
	
	return aux;
	
}



float  var(float *u,int size)
{
	
	float *ptru=&u[0];
	float mean=0.0;
	float mean2 = 0.0;
	for(int i=0;i<size;i++,ptru++) { mean +=  *ptru; mean2 += *ptru * (*ptru);}

	mean/=(float) size; mean2/=(float) size;
	float var = mean2- mean*mean;

	return var;
}



float  mean(float *u,int size)
{

	float *ptru=&u[0];
	float mean=0.0;
	for(int i=0; i<size; i++,ptru++)  mean +=  *ptru;
	mean/=(float) size;
	return mean;
}



float  median(float *u,int size)
{
	
	float *vector = new float[size];
	float *index = new float[size];
	

	float *ptru=&u[0];
	for(int i=0; i<size; i++,ptru++)  { vector[i] = *ptru; index[i] = (float) i;}
	quick_sort(vector,index,size);	
	
	float median;
	if (size%2==1)	median = vector[size/2];
	else	median = (vector[size/2] + vector[(size/2) - 1])/2;
	
	delete[] vector; delete[] index;
	return median;
}




int normalize(float *u,int size) 
{
  
  float some = 0.0;  
  for(int i=0;i<size;i++) some+=u[i];
  if (some != 0.0) {	for(int i=0;i<size;i++) u[i]/=some;  return 1;}  
  else return 0;	
}
	
	

float  nearest(float *u,float value,int *pos,int size)
{
	float mindif =  fabsf(value - u[0]);
	float minvalue = u[0];
	if (pos) *pos = 0;

	for(int i=0;i<size;i++){

		float dif = fabsf(value - u[i]);
		if (dif < mindif) { mindif=dif; minvalue=u[i]; if (pos!=NULL) *pos=i;} 

	}

	return minvalue;
}


void binarize(float *u, float *v,float value, int inverse, int size)
{ 
	for(int i=0;i<size;i++){
 
		if (u[i] >= value && !inverse) 	v[i]= 255.0;
		else if (u[i] <= value && inverse)  v[i]= 255.0;
		else v[i]= 0.0;
	}
}



float *  gauss(int sflag,float std,int *size)
{

   float *u,prec = 4.0,shift;
   double v;
//   int n,i,flag;
   int n,i;
   int flag = 1; //Guoshen Yu


   if (sflag) n=*size;
   else
	n = 1+2*(int)ceil((double)std*sqrt(prec*2.*log(10.)));   
   
   u =new float[n];

   if (n==1) 
    u[0]=1.0;
   else{

      shift = 0.5*(float)(n-1);

      for (i=(n+1)/2;i--;) {
         v = ((double)i - (double) shift)/(double)std;
         u[i] = u[n-1-i] = (float) exp(-0.5*v*v); 
      }
   }	

//   if (flag = normalize(u,n)) {
	if (flag == normalize(u,n)) 
	{
		*size=n;
        return u;
    } 
	else 
	{
		printf("ERROR: _gauss: _normalize: normalization equals zero.\n");
		delete[] u; /*memcheck*/
		return 0; // Guoshen Yu
    }
}



/* Quicksort,  values in arr are set in increasing order and brr elements are switched at the same time*/
void FSWAP(float *x,float *y)
{
  float aux;
  aux=*x;
  *x=*y;
  *y=aux;
}



void quick_sort(float *arr,float *brr,int n)
{
  int M=7,NSTACK=50;
  int i,ir,j,k,jstack=-1,l=0;
  float a,b;
  int istack[50];
  
  ir=n-1;


  for(;;){
    if(ir-l<M){
      for(j=l+1;j<=ir;j++){
	a=arr[j];
	b=brr[j];
	for(i=j-1;i>=l;i--){
	  if (arr[i]<=a) break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;

      }

      if (jstack<0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {

      k=(l+ir) >> 1;
      FSWAP(&arr[k],&arr[l+1]);
      FSWAP(&brr[k],&brr[l+1]);
      if (arr[l]>arr[ir]){
	FSWAP(&arr[l],&arr[ir]);
	FSWAP(&brr[l],&brr[ir]);
      }
      if (arr[l+1]>arr[ir]){
	FSWAP(&arr[l+1],&arr[ir]);
	FSWAP(&brr[l+1],&brr[ir]);
      }
      if (arr[l]>arr[l+1]){
	FSWAP(&arr[l],&arr[l+1]);
	FSWAP(&brr[l],&brr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for(;;){
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i) break;
	FSWAP(&arr[i],&arr[j]);
	FSWAP(&brr[i],&brr[j]);
      }

      arr[l+1]=arr[j];
      arr[j]=a;
      brr[l+1]=brr[j];
      brr[j]=b;
      jstack+=2;

      if (jstack>=NSTACK) { printf("Stack too small\n"); exit(-1);}
      if (ir-i+1>=j-l){
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }

    }
  }

  
}



/// histogram of values. 'n' (number of bins) or  's' (step) must be selected in flag while the other value is filled 
float * histo(float* input, float *iminim, float *imaxim, int *n, float *s, int size, char flag)
{
		
	if (flag != 's' && flag != 'n') 	{ printf("Warning (histo): Please select s or n as flag\n");  return NULL;	}
			
	float minim;
	if (iminim) minim = *iminim;
	else minim = min(input, NULL, size);

	float maxim;
	if (imaxim) maxim = *imaxim;
	else maxim = max(input, NULL, size);

	int num;
	float step;
	if (flag == 'n')
	{
		 num = *n;
		 step = (maxim-minim)/ (float)num;
		 *s = step;
	} else
	{
		step = *s;
		num = (int)(0.5+(maxim-minim)/step);
		*n = num;
	}
	
	float *histo = new float[num];
	clear(histo,0.0,num);
		
	
	for(int i=0; i < size; i++)
	{
	
		int cell = (int) floorf((input[i]-minim) / step);
	
		if (cell < 0) cell = 0;
		if (cell >= num) cell = num - 1;
		
		histo[cell]++;
	}
	
	
	return histo;
	
}





/////////////////////////////////////////////////////////////////////////////////////////////////////   Used and checked image functions

void compute_gradient_orientation(float* igray,float *grad, float *ori, int width, int height)
{

    float xgrad, ygrad;
    int rows, cols, r, c;

    rows = height;
    cols = width;
    

    for (r = 0; r < rows; r++)
      for (c = 0; c < cols; c++) {
        if (c == 0)
          xgrad = 2.0 * (igray[r*cols+c+1] - igray[r*cols+c]);
        else if (c == cols-1)
          xgrad = 2.0 * (igray[r*cols+c] - igray[r*cols+c-1]);
        else
          xgrad = igray[r*cols+c+1] - igray[r*cols+c-1];
        if (r == 0)
          ygrad = 2.0 * (igray[r*cols+c] - igray[(r+1)*cols+c]);
        else if (r == rows-1)
          ygrad = 2.0 * (igray[(r-1)*cols+c] - igray[r*cols+c]);
        else
          ygrad = igray[(r-1)*cols+c] - igray[(r+1)*cols+c];
        

        if (grad) grad[r*cols+c] = (float)sqrt((double)(xgrad * xgrad + ygrad * ygrad));
        if (ori) ori[r*cols+c] = (float)atan2 (-(double)ygrad,(double)xgrad);
      
      }
}




void sample(float *igray,float *ogray, float factor, int width, int height)
{

	int swidth = (int)((float) width / factor);
	int sheight = (int)((float) height / factor);
	
	for(int j=0; j < sheight; j++)
	 for(int i=0; i < swidth; i++)
		ogray[j*swidth + i] = igray[(int)((float) j * factor) * width + (int) ((float) i*factor)];

}


void sample_aglomeration(float *igray,float *ogray, float factor, int width, int height)
{

	int swidth = (int)((float) width / factor);
	int sheight = (int)((float) height / factor);
	int ssize = swidth * sheight;

	clear(ogray,0.0,swidth*sheight);

	for(int j=0; j < height; j++)
	 for(int i=0; i < width; i++){
	 	int index = (int)((float) j / factor) * swidth + (int) ((float) i / factor);
		if (index < ssize) ogray[index]  += igray[j*width+i];
	}

	factor *= factor;
	for(int i = 0; i < swidth*sheight; i++)
		ogray[i] /= factor; 


}


/*
void extract(float *igray,float *ogray, int ax, int ay,int cwidth, int cheight,int width, int height)
{
	for(int j=0; j < cheight; j++)
		for(int i=0; i < cwidth; i++)	
			ogray[j*cwidth + i] = igray[(ay+j)*width + ax+i];	
}
 */


void gray(float *red, float *green,float *blue, float *out, int width, int height)
{
	for(int i=width*height-1; i>0; i--)	out[i] = (red[i] + green[i] + blue[i]) /3.0;
}


/*
float l2_distance(float *u0,float *u1,int i0,int j0,int i1,int j1,int radius,int width,int height)
{

	int wsize=(2*radius+1)*(2*radius+1);

	float dist=0.0;       
	for (int s=-radius; s<= radius; s++){
	
		int l = (j0+s)*width + (i0-radius);
		float *ptr0 = &u0[l];
	
		l = (j1+s)*width + (i1-radius);
		float *ptr1 = &u1[l];
	
		for(int r=-radius;r<=radius;r++,ptr0++,ptr1++){	float dif = (*ptr0 - *ptr1); dist += (dif*dif); }

	}

	dist/=(float) wsize;
	return dist;
}

 
float l2_distance_non_normalized(float *u0,float *u1,int i0,int j0,int i1,int j1,int radius,int width,int height)
{


	float dist=0.0;       
	for (int s=-radius; s<= radius; s++){
	
		int l = (j0+s)*width + (i0-radius);
		float *ptr0 = &u0[l];
	
		l = (j1+s)*width + (i1-radius);
		float *ptr1 = &u1[l];
	
		for(int r=-radius;r<=radius;r++,ptr0++,ptr1++){	float dif = (*ptr0 - *ptr1); dist += (dif*dif); }

	}

	return dist;
}


float l2_distance_nsq(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int yradius,int width,int height)
{

	int wsize=(2*xradius+1)*(2*yradius+1);

	float dist=0.0;       
	for (int s=-yradius; s <= yradius; s++){
	
		int l = (j0+s)*width + (i0-xradius);
		float *ptr0 = &u0[l];
	
		l = (j1+s)*width + (i1-xradius);
		float *ptr1 = &u1[l];
	
		for(int r=-xradius;r<=xradius;r++,ptr0++,ptr1++)
		{
			float dif = (*ptr0 - *ptr1); dist += (dif*dif); 
		}

	}

	dist/=(float) wsize;
	return dist;
}




float weighted_l2_distance(float *u0,float *u1,int i0,int j0,int i1,int j1,int width,int height,float *kernel,int radius)
{


	float *ptrk=&kernel[0];
	float dist=0.0;       
	for (int s=-radius; s<= radius; s++){
	
		int l = (j0+s)*width + (i0-radius);
		float *ptr0 = &u0[l];
	
		l = (j1+s)*width + (i1-radius);
		float *ptr1 = &u1[l];
	
	
		for(int r=-radius;r<=radius;r++,ptr0++,ptr1++,ptrk++){ float dif = (*ptr0 - *ptr1); dist += *ptrk*(dif*dif); }
	
	}

	return dist;
}



float weighted_l2_distance_nsq(float *u0,float *u1,int i0,int j0,int i1,int j1,int width,int height,float *kernel,int xradius, int yradius)
{


	float *ptrk=&kernel[0];
	float dist=0.0;       
	for (int s=-yradius; s<= yradius; s++){
	
		int l = (j0+s)*width + (i0-xradius);
		float *ptr0 = &u0[l];
	
		l = (j1+s)*width + (i1-xradius);
		float *ptr1 = &u1[l];
	
	
		for(int r=-xradius;r<=xradius;r++,ptr0++,ptr1++,ptrk++){ float dif = (*ptr0 - *ptr1); dist += *ptrk*(dif*dif); }
	
	}

	return dist;
}
*/



/// RGV to YUV conversion
void rgb2yuv(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height) 
{
	int size=height*width;
  
	for(int i=0;i<size;i++){
		y[i] = COEFF_YR * r[i] + COEFF_YG * g[i] + COEFF_YB * b[i];
		u[i] = r[i] - y[i];
		v[i] = b[i] - y[i];
	}
}



void rgb2yuv(float *r,float *g,float *b,float *y,float *u,float *v,float yR, float yG, float yB, int width,int height) 
{
	int size=height*width;
	

	for(int i=0;i<size;i++){
		y[i] = yR * r[i] + yG * g[i] + yB * b[i];
		u[i] = r[i] - y[i];
		v[i] = b[i] - y[i];
	}
}



/// YUV to RGB conversion
void yuv2rgb(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height)  
{
	int size=height*width;
  
	for(int i=0;i<size;i++){
		g[i] = (y[i] - COEFF_YR * (u[i] + y[i]) - COEFF_YB * (v[i] + y[i])) / COEFF_YG;
		r[i] = u[i] + y[i];
		b[i] = v[i] + y[i];
	}
}


void yuv2rgb(float *r,float *g,float *b,float *y,float *u,float *v,float yR, float yG, float yB, int width,int height)  
{
	int size=height*width;
	
	for(int i=0;i<size;i++){
		g[i] = (y[i] - yR * (u[i] + y[i]) - yB * (v[i] + y[i])) / yG;
		r[i] = u[i] + y[i];
		b[i] = v[i] + y[i];
	}
}


void draw_line(float *igray, int a0, int b0, int a1, int b1, float value, int width, int height)
{

  int bdx,bdy;
  int sx,sy,dx,dy,x,y,z,l;
  
  bdx = width;
  bdy = height;

  if (a0 < 0) a0=0; 
  else if (a0>=bdx) a0=bdx-1;
   
  if (a1<0)  a1=0; 
  else  if (a1>=bdx)   a1=bdx-1;
	  
  if (b0<0) b0=0; 
  else if (b0>=bdy) b0=bdy-1;
   
  if (b1<0) 	b1=0; 
  else if (b1>=bdy) b1=bdy-1; 

  if (a0<a1) { sx = 1; dx = a1-a0; } else { sx = -1; dx = a0-a1; }
  if (b0<b1) { sy = 1; dy = b1-b0; } else { sy = -1; dy = b0-b1; }
  x=0; y=0;
  
  if (dx>=dy) 
    {
      z = (-dx) / 2;
      while (abs(x) <= dx) 
	{

	  l =  (y+b0)*bdx+x+a0; 
	
	  igray[l] = value;
	  
	  x+=sx;
	  z+=dy;
	  if (z>0) { y+=sy; z-=dx; }

	} 

    }
  else 
    {
      z = (-dy) / 2;
      while (abs(y) <= dy) {

	l = (y+b0)*bdx+x+a0;
  	igray[l] = value;
 
	y+=sy;
	z+=dx;
	if (z>0) { x+=sx; z-=dy; }
      }
    }

}



void draw_square(float *igray, int a0, int b0, int w0, int h0, float value, int width, int height)
{

	draw_line(igray,a0,b0,a0+w0,b0,value,width,height);
	draw_line(igray,a0,b0,a0,b0+h0,value,width,height);
	draw_line(igray,a0+w0,b0,a0+w0,b0+h0,value,width,height);
	draw_line(igray,a0,b0+h0,a0+w0,b0+h0,value,width,height);

}



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/*


*/
/*
*/





/*

void _sign(float *u,float *v, int size)
{

	int i=0;
 
	for(i=0;i<size;i++){
		 
		if (u[i]>0) v[i] = 1.0;
		else if (u[i]<0) v[i]=-1.0;
		else v[i]=0.0;
	}
	

}

	
	
	
	





void _multiple(float *u,float multiplier,int size) 
{
   int i=0;
   float *ptru;

   ptru=&u[0];
   for(i=0;i<size;i++,ptru++)  *ptru=multiplier*(*ptru);

}


void _product(float *u,float *v,int size) 
{
  int i;
  float *ptru,*ptrv;
  
  ptru=&u[0];
  ptrv=&v[0];
  for(i=0;i<size;i++,ptru++,ptrv++)   *ptru= *ptru*(*ptrv);
}



int _is_increasing(float *u,float tolerance, int size)
{

	int i=1;
	while (i < size && u[i] > (u[i-1] - tolerance))
	{
		i++;
	}
	
	if (i==size) return 1;
	else return 0; 
	
}








void _offset(float *u,float offset,int size)  
{
  int i=0;
  float *ptru;

  ptru=&u[0];
  for(i=0;i<size;i++,ptru++)  *ptru=*ptru + offset;

}












		
		
void _threshold(float *u, float *v,float valuem,float valueM, int size)
{

	int i;
 
	for(i=0;i<size;i++){
 
		if (u[i] >= valueM) 	v[i]= valueM;
		else if (u[i] <= valuem)  v[i]= valuem;
		else v[i] = u[i];
			
	}
	
}
		



void _absdif(float *u, float *v,int size)  
{
	int i=0;
 
	for(i=0;i<size;i++)  u[i] = (float) fabsf( u[i] -  v[i]);
}






	
	
	
float *  _diag_gauss(int sflag,float std,int *size) //Create a 1d gauss kernel of standard deviation std  (megawave2)
{

	float *u,prec = 4.0,shift;
	double v;
	int n,i,flag;

	if (sflag) n=*size;
	else
		n = 1+2*(int)ceil((double)std*sqrt(prec*2.*log(10.)));   
   
	u =(float *) malloc(n*sizeof(float));

	if (n==1) 
		u[0]=1.0;
	else{

		shift = 0.5*(float)(n-1);

		for (i=(n+1)/2;i--;) {
			
			v = ((double)i - (double) shift)/(double)std;
			
			u[i] = u[n-1-i] = (float) exp(-2.0*0.5*v*v);  // 2.0 because distances are in the diagonal 
		
		}
	}	

	if (flag = _normalize(u,n)) {
		*size=n;
		return u;
	} else {
		printf("ERROR: mdSigGaussKernel: mdSigNormalize: normalization equals zero.\n Try to reduce std.\n");
	}
}
	








void  _quant(float *u,float *v,float lambda,int size)
{

	int i,n;
	float a=lambda/2;

	for(i=0;i<size;i++){
		
		n = (int) floorf(u[i] / a);
		if (n%2==0)
			v[i] = (float) n * a;
		else 
			v[i] = (float) (n+1) * a;
	}
	
}



void  _projectquant(float *u,float *v,float lambda,int size)
{
	
	int i;
	float a=lambda/2;

	for(i=0;i<size;i++){
		
		if (v[i] < u[i] - a) v[i]=u[i] - a;
		else if (v[i] > u[i] + a) v[i]=u[i] + a;
	
	}
	
}

*/






