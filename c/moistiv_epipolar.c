// somewhat de-uglyfied code from moistiv

#include <math.h>

#include "fail.c"

#ifndef M_PI
#define M_PI 3.1416
#endif


/*--------------------------- MegaWave2 module  -----------------------------*/
/* mwcommand
 name = {stereomatch};
 author = {"Lionel Moisan"};
 version = {"1.1"};
 function = {"Detect and rate the best stereo correspondence between 2D point matches"};
 usage = {
 'v'->verb          "verbose",
 's'->stop          "stop as soon as the first meaningful correspondence is found",
 't':[t=10000]->t   "maximum number of ransac trials (default: 10000)",
 'n'->n_flag        "in order NOT to reinitialize the random seed",
 'm':[mode=3]->mode "mode: 0=deterministic 1=ransac 2=optimized ransac (ORSA) 3=automatic",
 u1->u1             "input: image 1 (used for dimensions)",
 p1->p1             "input: 1st Flist of 2D points",
 p2->p2             "input: 2nd Flist of 2D points",
 f<-f               "output: fundamental matrix (3x3 Flist)",
 index<-index       "output: indexes of matching pairs (1-Flist)",
 lnfa<-stereomatch  "output: meaningfulness of the matching (-log10(nfa))"
};
*/

/*----------------------------------------------------------------------
 v1.0 (10/2007): initial version from private file stereomatch3.c (LM)
 v1.1 (11/2008): useless lines removed
----------------------------------------------------------------------*/

//#include <stdio.h>
//#include <math.h>
//#include <time.h>
//#include <unistd.h>
//#include  "mw.h"


/*-------------------- GENERAL PURPOSE ROUTINES --------------------*/

/* routines for vectors and matrices */

static
float *vector(int nl, int nh)
{
  float *v;

  v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
  if (!v) fail("allocation failure in vector()");
  return v-nl;
}

static
float **matrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  float **m;

  m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
  if (!m) fail("allocation failure 1 in matrix()");
  m -= nrl;
  for(i=nrl;i<=nrh;i++) {
    m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
    if (!m[i]) fail("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

static
void free_vector(float *v, int nl, int nh)
{
  free((char*) (v+nl));
}

static
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

/* Singular Value Decomposition routine */

//static float at,bt,ct;
#define PYTHAG(a,b) hypot(a,b)
//#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ?
//(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static
void svdcmp(float **a, int m, int n, float *w, float **v)
{
  int flag,i,its,j,jj,k,l,nm;
  float c,f,h,s,x,y,z;
  float anorm=0.0,g=0.0,scale=0.0;
  float *rv1;

  if (m<n) fail("SVDCMP: You must augment A with extra zero rows");
  rv1=vector(1,n);
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	if (i != n) {
	  for (j=l;j<=n;j++) {
	    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	    f=s/h;
	    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	  }
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	if (i != m) {
	  for (j=l;j<=m;j++) {
	    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	  }
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=n;i>=1;i--) {
    l=i+1;
    g=w[i];
    if (i < n)
      for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      if (i != n) {
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	  f=(s/a[i][i])*g;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else {
      for (j=i;j<=m;j++) a[j][i]=0.0;
    }
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if (fabs(rv1[l])+anorm == anorm) {
	  flag=0;
	  break;
	}
	if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  if (fabs(f)+anorm != anorm) {
	    g=w[i];
	    h=PYTHAG(f,g);
	    w[i]=h;
	    h=1.0/h;
	    c=g*h;
	    s=(-f*h);
	    for (j=1;j<=m;j++) {
	      y=a[j][nm];
	      z=a[j][i];
	      a[j][nm]=y*c+z*s;
	      a[j][i]=z*c-y*s;
	    }
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
	}
	break;
      }
      if (its == 30) fail("No convergence in 30 SVDCMP iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=PYTHAG(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=PYTHAG(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y=y*c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=PYTHAG(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=(c*g)+(s*y);
	x=(c*y)-(s*g);
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_vector(rv1,1,n);
}

#undef SIGN
#undef MAX
#undef PYTHAG


/* Compute the real roots of a third order polynomial */
/* returns 1 or 3, the number of roots found */

static
int FindCubicRoots(float coeff[4], float x[3])
{
  float a1 = coeff[2] / coeff[3];
  float a2 = coeff[1] / coeff[3];
  float a3 = coeff[0] / coeff[3];

  double Q = (a1 * a1 - 3 * a2) / 9;
  double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
  double Qcubed = Q * Q * Q;
  double d = Qcubed - R * R;

  /* Three real roots */
  if (d >= 0) {
    double theta = acos(R / sqrt(Qcubed));
    double sqrtQ = sqrt(Q);
    x[0] = -2 * sqrtQ * cos( theta             / 3) - a1 / 3;
    x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
    x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
    return (3);
  }

  /* One real root */
  else {
    double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
    if (R > 0)
      e = -e;
    x[0] = (e + Q / e) - a1 / 3.;
    return (1);
  }
}


///* logarithm (base 10) of binomial coefficient */
//float logcombi(k,n)
//     int k,n;
//{
//  double r;
//  int i;
// 
//  if (k>=n || k<=0) return(0.);
//  if (n-k<k) k=n-k;
//  r = 0.;
//  for (i=1;i<=k;i++) 
//    r += log10((double)(n-i+1))-log10((double)i);
//
//  return((float)r);
//}
//
///* tabulate logcombi(.,n) */
//float *makelogcombi_n(n)
//     int n;
//{
//  float *l;
//  int k;
//
//  l = (float *)malloc((n+1)*sizeof(float));
//  for (k=0;k<=n;k++) l[k]=logcombi(k,n);
//
//  return(l);
//}
//
///* tabulate logcombi(k,.) */
//float *makelogcombi_k(k,nmax)
//     int k,nmax;
//{
//  float *l;
//  int n;
//
//  l = (float *)malloc((nmax+1)*sizeof(float));
//  for (n=0;n<=nmax;n++) l[n]=logcombi(k,n);
//
//  return(l);
//}
//
//
///* get a (sorted) random 7-uple of 0..n-1 */
//void random_p7(k,n)
//     int *k,n;
//{
//  int i,j,j0,r;
//
//  for (i=0;i<7;i++) {
//    r = (rand()>>3)%(n-i);
//    for (j=0;j<i && r>=k[j];j++) r++;
//    j0 = j;
//    for (j=i;j>j0;j--) k[j]=k[j-1];
//    k[j0]=r;
//  }
//}

/*-------------------- END OF GENERAL PURPOSE ROUTINES --------------------*/


///* float comparison for qsort() */
//static int compf(i,j)
//     void *i,*j;
//{ 
//  float a,b;
//
//  a = *((float *)i);
//  b = *((float *)j);
//  return(a<b?-1:(a>b?1:0));
//}
//
///* find the increasing sequence of squared distances to epipolar lines */
///* e[n*2] = distances, e[n*2+1] = indexes (to cast into an int) */
//
//void matcherrorn(F,p1,p2,e) 
//     float **F;
//     Flist p1,p2;
//     float *e;
//{
//  int i;
//  float x,y,a,b,c,d;
//
//  for (i=0;i<p1->size;i++) {
//    x = p1->values[i*2];
//    y = p1->values[i*2+1];
//    a = F[1][1]*x+F[1][2]*y+F[1][3];
//    b = F[2][1]*x+F[2][2]*y+F[2][3];
//    c = F[3][1]*x+F[3][2]*y+F[3][3];
//    d = (a*p2->values[i*2]+b*p2->values[i*2+1]+c);
//    e[i*2] = (d*d)/(a*a+b*b);
//    e[i*2+1] = (float)i;
//  }
//  qsort(e,p1->size,2*sizeof(float),compf);
//}


/*---------- compute the epipolar geometry associated to 7 pairs ----------*/
/*                                                                         */
/*  INPUT: the points are (m1[k[i]*2],m1[k[i]*2+1]), m2... 0<i<7           */
/*                                                                         */
/*  OUTPUT:                                                                */
/*             return the number of roots found, stored in z[]             */
/*   the epipolar matrices are F1[i][j]+z[k]*F2[i][j], 1<=i,j<=3, 0<=k<m   */


/* global (intermediate) variables for epipolar() */
//float *w,**c,**v,a[4];

int moistiv_epipolar(float *m1, float *m2, int *k, float *z,
		float F1[4][4], float F2[4][4])
{
  int i,j,i2,i3,imin1,imin2;
  float wmin1,wmin2;
  float **c  = matrix(1,9,1,9);
  float *w  = vector(1,9);
  float **v  = matrix(1,9,1,9);
  float a[4];

  /* build 9xn matrix from point matches */
  for (i=0;i<7;i++) {
    c[i+1][1] = m1[k[i]*2  ]*m2[k[i]*2  ];
    c[i+1][2] = m1[k[i]*2+1]*m2[k[i]*2  ];
    c[i+1][3] =                      m2[k[i]*2  ];
    c[i+1][4] = m1[k[i]*2  ]*m2[k[i]*2+1];
    c[i+1][5] = m1[k[i]*2+1]*m2[k[i]*2+1];
    c[i+1][6] =                      m2[k[i]*2+1];
    c[i+1][7] = m1[k[i]*2  ];
    c[i+1][8] = m1[k[i]*2+1];
    c[i+1][9] = 1.;
  }
  for (i=1;i<=9;i++) c[8][i] = c[9][i] = 0.;
  
  /* SVD */
  svdcmp(c,9,9,w,v);
  
  /* look for the two smallest eigenvalue of c'c */
  if (w[1]<w[2]) {
    imin1 = 1;
    imin2 = 2;
  } else {
    imin2 = 1;
    imin1 = 2;
  }
  wmin1 = w[imin1];
  wmin2 = w[imin2];
  for (i=3;i<=9;i++) {
    if (w[i]<wmin1) {imin2=imin1; wmin2=wmin1; wmin1=w[i]; imin1=i;}
    else if (w[i]<wmin2) {wmin2=w[i]; imin2=i;}
  }
  
  /* build basis of solutions */
  for (i=1;i<=3;i++) 
    for (j=1;j<=3;j++) {
      F1[i][j] = v[(i-1)*3+j][imin1];
      F2[i][j] = v[(i-1)*3+j][imin2]-F1[i][j];
    }
  
  /* build cubic polynomial P(x)=det(F1+xF2) */
  a[0] = a[1] = a[2] = a[3] = 0.;
  for (i=1;i<=3;i++) {
    i2 = i%3+1;
    i3 = i2%3+1;
    a[0] += F1[i][1]*F1[i2][2]*F1[i3][3];
    a[1] += 
      F2[i][1]*F1[i2][2]*F1[i3][3]+
      F1[i][1]*F2[i2][2]*F1[i3][3]+
      F1[i][1]*F1[i2][2]*F2[i3][3];
    a[2] += 
      F1[i][1]*F2[i2][2]*F2[i3][3]+
      F2[i][1]*F1[i2][2]*F2[i3][3]+
      F2[i][1]*F2[i2][2]*F1[i3][3];
    a[3] += F2[i][1]*F2[i2][2]*F2[i3][3];
  }
  for (i=1;i<=3;i++) {
    i2 = (i+1)%3+1;
    i3 = (i2+1)%3+1;
    a[0] -= F1[i][1]*F1[i2][2]*F1[i3][3];
    a[1] -= 
      F2[i][1]*F1[i2][2]*F1[i3][3]+
      F1[i][1]*F2[i2][2]*F1[i3][3]+
      F1[i][1]*F1[i2][2]*F2[i3][3];
    a[2] -= 
      F1[i][1]*F2[i2][2]*F2[i3][3]+
      F2[i][1]*F1[i2][2]*F2[i3][3]+
      F2[i][1]*F2[i2][2]*F1[i3][3];
    a[3] -= F2[i][1]*F2[i2][2]*F2[i3][3];
  }


  free_matrix(c,1,9,1,9);
  free_matrix(v,1,9,1,9);
  free_vector(w,1,9);

  return(FindCubicRoots(a,z));
}


///*------------------------------ MAIN MODULE ------------------------------*/
//
//
///* NOTE: if f=NULL, the fundamental matrix is not returned */
///* idem if index=NULL */
///* if *mode=3, the mode chosen (0 or 2) is returned in *mode */
//
//float stereomatch(u1,p1,p2,f,index,t,verb,n_flag,mode,stop)
//     Cimage u1;
//     Flist p1,p2,f,index;
//     int *t,*verb,*n_flag,*mode,*stop;
//{
//  int i,j,i0,k[8],idk[8],*id,m,n,l,minicur,miniall,delete_index,nid;
//  int niter,maxniter,better,cont,optimization;
//  float **F1,**F2,**F,nx,ny,z[3],minepscur,minepsall,nfa;
//  float norm,minlogalphacur,minlogalphaall,logalpha,logalpha0;
//  float *e,*logcn,*logc7,loge0;
//
//  /* initialize random seed if necessary */
//  if (!n_flag) srand48( (long int) time (NULL) + (long int) getpid() );
//  
//  /* check sizes */
//  if (p1->size!=p2->size || p1->size<7) 
//    mwerror(FATAL,1,"Inconsistent sizes. ");
//  n = p1->size;
//  
//  /* tabulate logcombi */
//  loge0 = (float)log10(3.*(double)(n-7));
//  logcn = makelogcombi_n(n);
//  logc7 = makelogcombi_k(7,n); 
//  
//  /* choose mode */
//  if (*mode==3) {
//    if (logcn[7]<=(float)log10((double)(*t)))
//      *mode=0; 
//    else *mode=2;
//  }
//  if (verb) 
//    switch(*mode) {
//    case 0: 
//      i = (int)(0.5+pow(10.,logc7[n]));
//      printf("I will use deterministic mode (systematic search).\n");
//      printf("I have to test %d bases\n",i);
//      break;
//    case 1:
//      printf("I will use pure stochastic mode with no optimization.\n");
//      break;
//    case 2:
//      printf("I will use optimized stochastic mode (ORSA).\n");
//    }
//
//  /* normalize coordinates */
//  nx = (float)u1->ncol;
//  ny = (float)u1->nrow;
//  norm = 1./(float)sqrt((double)(nx*ny));
//  logalpha0 = (float)(log10(2.)+0.5*log10((double)((nx*nx+ny*ny)*norm*norm)));
//  for (i=0;i<n;i++) {
//    p1->values[i*2  ] =  (p1->values[i*2  ]-0.5*nx)*norm;
//    p1->values[i*2+1] =  (p1->values[i*2+1]-0.5*ny)*norm;
//    p2->values[i*2  ] =  (p2->values[i*2  ]-0.5*nx)*norm;
//    p2->values[i*2+1] =  (p2->values[i*2+1]-0.5*ny)*norm;
//  }
//
//  /* allocate and initialize memory */
//  c  = matrix(1,9,1,9);
//  w  = vector(1,9);
//  v  = matrix(1,9,1,9);
//  F  = matrix(1,3,1,3);
//  F1  = matrix(1,3,1,3);
//  F2  = matrix(1,3,1,3);
//  if(f) {
//    mw_change_flist(f,3,3,3);
//    mw_clear_flist(f,0.);
//  }
//  delete_index = (index?0:1);
//  index = mw_change_flist(index,n,0,1);
//  e = (float *)malloc(2*n*sizeof(float));
//  id = (int *)malloc(n*sizeof(int));
//  for (i=0;i<n;i++) id[i]=i;
//  nid = n;
//
//  maxniter = (*mode==0?*t:*t-(*t)/10);
//  minlogalphaall = minepsall = 10000.;
//  niter = optimization = 0;
//  i0=0; k[0]=-1; k[7]=n;
//
//  /********** MAIN LOOP **********/
//  do {
//
//    niter++;
//
//     /* build next list of points */
//    if (*mode) random_p7(k,nid);
//    else {
//      k[i0]++;
//      for (i=i0+1;i<=6;i++) k[i]=k[i-1]+1;
//    }
//    for (i=0;i<7;i++) idk[i]=id[k[i]];
//
//    /* find epipolar transform */
//    m = epipolar(p1->values,p2->values,idk,z,F1,F2);
//
//   /* loop on roots */
//    for (;m--;) {
//
//      for (i=1;i<=3;i++) 
//	for (j=1;j<=3;j++) 
//	  F[i][j] = F1[i][j]+z[m]*F2[i][j];
//      
//      /* sort errors */
//      matcherrorn(F,p1,p2,e);
//      
//      /* find most meaningful subset */
//      minepscur = minlogalphacur = 10000.;
//      for (i=7;i<n;i++) {
//	logalpha = logalpha0+0.5*(float)log10((double)e[i*2]);
//	nfa = loge0+logalpha*(float)(i-6)+logcn[i+1]+logc7[i+1];
//	if (nfa<minepscur) {
//	  minepscur = nfa;
//	  minicur = i;
//	  minlogalphacur = logalpha; 
//	}
//      }
//      
//      if (minepscur<minepsall) {
//	/* store best result so far */
//	better = 1;
//	minepsall = minepscur;
//	minlogalphaall = minlogalphacur; 
//	miniall = minicur;
//	if (f) 
//	  for (l=1;l<=3;l++) 
//	    for (j=1;j<=3;j++) 
//	      f->values[(l-1)*f->dim+j-1] = F[l][j];
//	for (i=0;i<=minicur;i++) 
//	  index->values[i] = e[i*2+1];
//	index->size = minicur+1;
//      } else better=0;
// 
//      if (*mode==2 && ((better && minepsall<0.) || 
//		       (niter==maxniter && !optimization))) {
//	if (!optimization) maxniter = niter + (*t)/10;
//	optimization = 1;
//	/* final optimization */
//	if (verb) {
//	  printf("   nfa=%f size=%d (niter=%d)\n",minepsall,miniall+1,niter);
//	  printf("optimization...\n");
//	}
//	nid = miniall+1;
//	for (j=0;j<=miniall;j++)
//	  id[j] = (int)index->values[j];
//      }
//
//    }
//
//   /* prepare next list of points */
//    if (*mode==0) 
//      for (i0=6;i0>=0 && k[i0]==k[i0+1]-1;i0--);
//
//    if (stop && minepsall<0.) cont=0;
//    else if (*mode==0) cont=(i0>=0?1:0); 
//    else cont=(niter<maxniter?1:0);
//
//  } while (cont);
//
//  if (verb) 
//    printf("best matching found:  %d points  log(alpha)=%f  (%d iterations)\n",
//	   miniall+1,minlogalphaall,niter);
//
//  /* free memory */
//  free(id);
//  free(e);
//  if (delete_index) mw_delete_flist(index);
//  free_matrix(F2,1,3,1,3);
//  free_matrix(F1,1,3,1,3);
//  free_matrix(F,1,3,1,3);
//  free_matrix(v,1,9,1,9);
//  free_vector(w,1,9);
//  free_matrix(c,1,9,1,9);
//  free(logc7);
//  free(logcn);
//
//  return(minepsall);
//}
