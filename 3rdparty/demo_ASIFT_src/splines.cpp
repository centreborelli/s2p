// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "splines.h"

double initcausal(double *c,int n,double z)
{
  double zk,z2k,iz,sum;
  int k;

  zk = z; iz = 1./z;
  z2k = pow(z,(double)n-1.);
  sum = c[0] + z2k * c[n-1];
  z2k = z2k*z2k*iz;
  for (k=1;k<=n-2;k++) {
    sum += (zk+z2k)*c[k];
    zk *= z;
    z2k *= iz;
  }
  return (sum/(1.-zk*zk));
}

double initanticausal(double *c,int n,double z)
{
  return((z/(z*z-1.))*(z*c[n-2]+c[n-1]));
}


void invspline1D(double *c,int size,double *z,int npoles)
{
  double lambda;
  int n,k;

  /* normalization */
  for (k=npoles,lambda=1.;k--;) lambda *= (1.-z[k])*(1.-1./z[k]);
  for (n=size;n--;) c[n] *= lambda;

  /*----- Loop on poles -----*/
  for (k=0;k<npoles;k++) {

    /* forward recursion */
    c[0] = initcausal(c,size,z[k]);
    for (n=1;n<size;n++) 
      c[n] += z[k]*c[n-1];

    /* backwards recursion */
    c[size-1] = initanticausal(c,size,z[k]);
    for (n=size-1;n--;) 
      c[n] = z[k]*(c[n+1]-c[n]);
    
  }
}


/*------------------------------ MAIN MODULE ------------------------------*/

// void finvspline(float *in,int order,float *out, int width, int height)

// Guoshen Yu, 2010.09.21, Windows version
void finvspline(vector<float> &in,int order,vector<float>& out, int width, int height)
// void finvspline(float *in,int order,float *out, int width, int height)
{
  double *c,*d,z[5];
  int npoles,nx,ny,x,y;
 
  ny = height; nx = width;

  /* initialize poles of associated z-filter */
  switch (order) 
    {
    case 2: z[0]=-0.17157288;  /* sqrt(8)-3 */
      break;

    case 3: z[0]=-0.26794919;  /* sqrt(3)-2 */ 
      break;

    case 4: z[0]=-0.361341; z[1]=-0.0137254;
      break;

    case 5: z[0]=-0.430575; z[1]=-0.0430963;
      break;
      
    case 6: z[0]=-0.488295; z[1]=-0.0816793; z[2]=-0.00141415;
      break;

    case 7: z[0]=-0.53528; z[1]=-0.122555; z[2]=-0.00914869;
      break;
      
    case 8: z[0]=-0.574687; z[1]=-0.163035; z[2]=-0.0236323; z[3]=-0.000153821;
      break;

    case 9: z[0]=-0.607997; z[1]=-0.201751; z[2]=-0.0432226; z[3]=-0.00212131;
      break;
      
    case 10: z[0]=-0.636551; z[1]=-0.238183; z[2]=-0.065727; z[3]=-0.00752819;
      z[4]=-0.0000169828;
      break;
      
    case 11: z[0]=-0.661266; z[1]=-0.27218; z[2]=-0.0897596; z[3]=-0.0166696; 
      z[4]=-0.000510558;
      break;
      
     default:
      printf("finvspline: order should be in 2..11.\n");
      exit(-1);
    }

  npoles = order/2;

  /* initialize double array containing image */
  c = (double *)malloc(nx*ny*sizeof(double));
  d = (double *)malloc(nx*ny*sizeof(double));
  for (x=nx*ny;x--;) 
    c[x] = (double)in[x];

  /* apply filter on lines */
  for (y=0;y<ny;y++) 
    invspline1D(c+y*nx,nx,z,npoles);

  /* transpose */
  for (x=0;x<nx;x++)
    for (y=0;y<ny;y++) 
      d[x*ny+y] = c[y*nx+x];
      
  /* apply filter on columns */
  for (x=0;x<nx;x++) 
    invspline1D(d+x*ny,ny,z,npoles);

  /* transpose directy into image */
  for (x=0;x<nx;x++)
    for (y=0;y<ny;y++) 
      out[y*nx+x] = (float)(d[x*ny+y]);

  /* free array */
  free(d);
  free(c);
}



/* extract image value (even outside image domain) */
//float v(float *in,int x,int y,float bg, int width, int height)
// Guoshen Yu, 2010.09.21, Windows version
float v(vector<float>& in,int x,int y,float bg, int width, int height)
// float v(float *in, int x,int y,float bg, int width, int height)
{
  if (x<0 || x>=width || y<0 || y>=height)
    return(bg); else return(in[y*width+x]);
}


/* c[] = values of interpolation function at ...,t-2,t-1,t,t+1,... */

/* coefficients for cubic interpolant (Keys' function) */
void keys(float *c,float t,float a)
{
  float t2,at;

  t2 = t*t;
  at = a*t;
  c[0] = a*t2*(1.0-t);
  c[1] = (2.0*a+3.0 - (a+2.0)*t)*t2 - at;
  c[2] = ((a+2.0)*t - a-3.0)*t2 + 1.0;
  c[3] = a*(t-2.0)*t2 + at;
}

/* coefficients for cubic spline */
void spline3(float *c,float t)
{
  float tmp;

  tmp = 1.-t;
  c[0] = 0.1666666666*t*t*t;
  c[1] = 0.6666666666-0.5*tmp*tmp*(1.+t);
  c[2] = 0.6666666666-0.5*t*t*(2.-t);
  c[3] = 0.1666666666*tmp*tmp*tmp;
}

/* pre-computation for spline of order >3 */
void init_splinen(float *a,int n)
{
  int k;

  a[0] = 1.;
  for (k=2;k<=n;k++) a[0]/=(float)k;
  for (k=1;k<=n+1;k++)
    a[k] = - a[k-1] *(float)(n+2-k)/(float)k;
}

/* fast integral power function */
float ipow(float x,int n)
{
  float res;

  for (res=1.;n;n>>=1) {
    if (n&1) res*=x;
    x*=x;
  }
  return(res);
}

/* coefficients for spline of order >3 */
void splinen(float *c,float t,float *a,int n)
{
  int i,k;
  float xn;
  
  memset((void *)c,0,(n+1)*sizeof(float));
  for (k=0;k<=n+1;k++) { 
    xn = ipow(t+(float)k,n);
    for (i=k;i<=n;i++) 
      c[i] += a[i-k]*xn;
  }
}


