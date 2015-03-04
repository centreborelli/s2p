// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef _LIBRARY_H_
#define _LIBRARY_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// #include <unistd.h>
#include <float.h>


#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

#define LUTMAX 30
#define LUTPRECISION 1000.0

#define TINY 1.0e-10
//#define MAXFLOAT 10000000.0


#define IRAC8    0.35355339  /* 1/sqrt(8)   */
#define IRAC2P2  0.29289322  /* 1/(sqrt(2)+2) */
#define IRAC2    0.70710678  /* 1/sqrt(2) */
#define RAC8P4   6.8284271   /* sqrt(8)+4 */
#define RADIANS_TO_DEGREES (180.0/M_PI)

#define PI 3.14159


#define COEFF_YR 0.299
#define COEFF_YG 0.587
#define COEFF_YB 0.114




///////////////////////////////////////////////////////////////// Errors
void wxwarning(const char * message,const char *function, const char *file);
void wxerror(const char * message, const char *function, const char *file);


/////////////////////////////////////////////////////////////////////////////////////////////////////   Used and checked mathematical functions
double fsqr(double a);

void  fill_exp_lut(float *lut,int size); /* Fills exp(x) for x great or equal than zero*/
float slut(float dif,float *lut); /* We look for f(dif) in the lut*/




/////////////////////////////////////////////////////////////////////////////////////////////////////   Used and checked mxsignal functions
float max(float *u,int *pos, int size);    /// Max(u), pos contains the index of the maximum 
float min(float *u,int *pos, int size);    /// Min(u), pos contains the index of the minimum 	
void max_u_v(float *u,float *v,int size); 
void max_u_k(float *u,float k,int size);
void min_u_v(float *u,float *v,int size);
void min_u_k(float *u,float k,int size);



void abs(float *u,float *v,int size);	 	/// v = abs(u)

void  copy(float *u,float *v,int size);			/// v = u
void  clear(float *u, float value ,int size);        	/// u = k 
void  combine(float *u,float a,float *v,float b, float *w,  int size);    /// w = a*u + b*v 
void  multiple(float *u,float multiplier,int size);         	/// u = K * u  


float scalar_product(float *u, float *v, int n);

// float lpdist(float *u,float *v,float *mask,int pow,int size);   
// float lpnorm(float *u,int fpow,int size);   

float mean(float *u,int size);
float var(float *u,int size);
float median(float *u,int size);

float nearest(float *u,float value,int *pos,int size);  /// Returns the nearest value in the vector u and its position if selected
void binarize(float *u, float *v,float value, int inverse, int size);    /// v = 255 if u > value 0 else  
int normalize(float *u,int size);    			/// u = u / sum_i u(i).   Returns 0 if mxsignal sum equals zero 



float *  gauss(int sflag,float std,int *size);  	/// Create a 1d gauss kernel of standard deviation std

// void addnoise(float *u,float *v,float std,long int randinit, int size);    
// void addnoise_var_afine(float *u,float *v,float a,float b,long int randinit, int size);    


void quick_sort(float *arr,float *brr,int n);  		/// Quicksort


/// histogram of values. 'n' (number of bins) or  's' (step) must be selected in flag while the other value is filled 
float * histo ( float* input, float *iminim, float *imaxim, int *n, float *s, int size, char flag );	





/////////////////////////////////////////////////////////////////////////////////////////////////////   Used and checked image functions

void compute_gradient_orientation(float* igray,float *grad, float *ori, int width, int height);

// void extract ( float *igray,float *ogray, int ax, int ay,int cwidth, int cweight,int width, int height );

void sample ( float *igray,float *ogray, float factor ,int width, int height);
void sample_aglomeration(float *igray,float *ogray, float factor, int width, int height);
 
void gray ( float *red, float *green,float *blue, float *out, int width, int height );

/*
float l2_distance ( float * u0,float *u1,int i0,int j0,int i1,int j1,int radius,int width,int height);
float l2_distance_non_normalized(float *u0,float *u1,int i0,int j0,int i1,int j1,int radius,int width,int height);
float weighted_l2_distance ( float *u0,float *u1,int i0,int j0,int i1,int j1,int width,int height,float * kernel,int radius );
float l2_distance_nsq(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int yradius,int width,int height);
float weighted_l2_distance_nsq(float *u0,float *u1,int i0,int j0,int i1,int j1,int width,int height,float *kernel,int xradius, int yradius);
*/

void rgb2yuv ( float *r,float *g,float *b,float *y,float *u,float *v,int width,int height );  
void yuv2rgb ( float *r,float *g,float *b,float *y,float *u,float *v,int width,int height ); 

void rgb2yuv(float *r,float *g,float *b,float *y,float *u,float *v,float yR, float yG, float yB, int width,int height);
void yuv2rgb(float *r,float *g,float *b,float *y,float *u,float *v,float yR, float yG, float yB, int width,int height);


void draw_line(float *igray, int a0, int b0, int a1, int b1, float value, int width, int height);
// void draw_circle(float *igray, int pi,int pj,float radius, float value, int width, int height);
void draw_square(float *igray, int a0, int b0, int w0, int h0, float value, int width, int height);



#endif // _LIBRARY_H_


/////////////////////////////////////// Not often used and not checked

//

//





/*		
*/




//void md_fsig_absdif(float *u,float *v,int size);  	/// v = abs(v-u)



//void md_fsig_sign(float *u,float *v, int size);    	//// ¿¿¿¿¿ ?????
//int md_fsig_is_increasing(float *u,float tolerance, int size);    ///// ¿¿¿¿¿ ?????	








//void md_fsig_multiple(float *u,float multiplier,int size);         	// u = K * u  
//void md_fsig_product(float *u,float *v,int size);        		// u = u * v 
//void md_fsig_offset(float *u,float offset, int size);    		// u = u - K  



//


//
//void md_fsig_threshold(float *u, float *v,float valuem,float valueM, int size);  // threshold into (m,M) 



//--- Conversion ---
//void md_fsig_float2char(float min, float max, float *u, float *v, int size);    // Linear Conversion between  (min,max) and (0,255)  



  
//void md_fsig_addnoise(float *u,float *v,float std,long int randinit, int size);    // Add gaussian noise of standard deviation sigma  



// v quantified mxmximage u with interval length lambda
// if u \in ( (n - 1/2) l, (n + 1/2 ) l )  - >  v = n l  
// n= 0 -> (-1/2 * l, 1/2 * l) //
//void  md_fsig_quant(float *u, float *v, float lambda,  int size);


// v is projected to the space of quantizations of length lambda that gives u
//void  md_fsig_projectquant(float *u,float *v,float lambda, int size);
