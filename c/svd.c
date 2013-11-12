/*  CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */

/* Note: this file contains code from the following original ccmath's files:
 * 	svduv.c
 * 	sv2uv.c
 * 	qrbdv.c
 * 	ldumat.c
 * 	ldvmat.c
 */
#include <math.h>
#include <stdlib.h>


static int qrbdv(double *dm,double *em,double *um,int mm,double *vm,int m)
{ int i,j,k,n,jj,nm;
  double u,x,y,a,b,c,s,t,w,*p,*q;
  for (j=1,t=fabs(dm[0]); j<m ;++j)
    if((s=fabs(dm[j])+fabs(em[j-1]))>t) t=s;
  t*=1.e-15; n=100*m; nm=m;
  for(j=0; m>1 && j<n ;++j){
    for(k=m-1; k>0 ;--k){
      if(fabs(em[k-1])<t) break;
      if(fabs(dm[k-1])<t){
        for(i=k,s=1.,c=0.; i<m ;++i){
          a=s*em[i-1]; b=dm[i]; em[i-1]*=c;
          dm[i]=u=sqrt(a*a+b*b); s= -a/u; c=b/u;
          for(jj=0,p=um+k-1; jj<mm ;++jj,p+=mm){
            q=p+i-k+1;
            w=c* *p+s* *q; *q=c* *q-s* *p; *p=w;
           }
	     }
        break;
       }
     }
    y=dm[k]; x=dm[m-1]; u=em[m-2];
    a=(y+x)*(y-x)-u*u; s=y*em[k]; b=s+s;
    u=sqrt(a*a+b*b);
    if(u!=0.){
      c=sqrt((u+a)/(u+u));
	  if(c!=0.) s/=(c*u); else s=1.;
	  for(i=k; i<m-1 ;++i){
        b=em[i];
        if(i>k){
          a=s*em[i]; b*=c;
	      em[i-1]=u=sqrt(x*x+a*a);
    	  c=x/u; s=a/u;
         }
        a=c*y+s*b; b=c*b-s*y;
        for(jj=0,p=vm+i; jj<nm ;++jj,p+=nm){
          w=c* *p+s* *(p+1); *(p+1)=c* *(p+1)-s* *p; *p=w;
         }
        s*=dm[i+1]; dm[i]=u=sqrt(a*a+s*s);
        y=c*dm[i+1]; c=a/u; s/=u;
        x=c*b+s*y; y=c*y-s*b;
        for(jj=0,p=um+i; jj<mm ;++jj,p+=mm){
          w=c* *p+s* *(p+1); *(p+1)=c* *(p+1)-s* *p; *p=w;
         }
	   }
     }
    em[m-2]=x; dm[m-1]=y;
    if(fabs(x)<t) --m;
    if(m==k+1) --m; 
   }
  return j;
}

static void ldumat(double *a,double *u,int m,int n)
{ double *p0,*q0,*p,*q,*w;
  int i,j,k,mm;
  double s,h;
  w=(double *)calloc(m,sizeof(double));
  for(i=0,mm=m*m,q=u; i<mm ;++i) *q++ =0.;
  p0=a+n*n-1; q0=u+m*m-1; mm=m-n; i=n-1;
  for(j=0; j<mm ;++j,q0-=m+1) *q0=1.;
  if(mm==0){ p0-=n+1; *q0=1.; q0-=m+1; --i; ++mm;}
  for(; i>=0 ;--i,++mm,p0-=n+1,q0-=m+1){
    if(*p0!=0.){
      for(j=0,p=p0+n,h=1.; j<mm ;p+=n) w[j++]= *p;
      h= *p0; *q0=1.-h;
      for(j=0,q=q0+m; j<mm ;q+=m) *q= -h*w[j++];
      for(k=i+1,q=q0+1; k<m ;++k){
	for(j=0,p=q+m,s=0.; j<mm ;p+=m) s+=w[j++]* *p;
	s*=h;
	for(j=0,p=q+m; j<mm ;p+=m) *p-=s*w[j++];
        *q++ = -s;
       }
     }
    else{
      *q0=1.;
      for(j=0,p=q0+1,q=q0+m; j<mm ;++j,q+=m) *q= *p++ =0.;
     }
   }
  free(w);
}

static void ldvmat(double *a,double *v,int n)
{ double *p0,*q0,*p,*q,*qq;
  double h,s;
  int i,j,k,mm;
  for(i=0,mm=n*n,q=v; i<mm ;++i) *q++ =0.;
  *v=1.; q0=v+n*n-1; *q0=1.; q0-=n+1;
  p0=a+n*n-n-n-1;
  for(i=n-2,mm=1; i>0 ;--i,p0-=n+1,q0-=n+1,++mm){
    if(*(p0-1)!=0.){
      for(j=0,p=p0,h=1.; j<mm ;++j,++p) h+= *p* *p;
      h= *(p0-1); *q0=1.-h;
      for(j=0,q=q0+n,p=p0; j<mm ;++j,q+=n) *q= -h* *p++; 
      for(k=i+1,q=q0+1; k<n ;++k){
        for(j=0,qq=q+n,p=p0,s=0.; j<mm ;++j,qq+=n) s+= *qq* *p++;
        s*=h;
        for(j=0,qq=q+n,p=p0; j<mm ;++j,qq+=n) *qq-=s* *p++;
        *q++ = -s;
       }
     }
    else{
      *q0=1.;
      for(j=0,p=q0+1,q=q0+n; j<mm ;++j,q+=n) *q= *p++ =0.;
     }
   }
}

static int svduv(double *d,double *a,double *u,int m,double *v,int n)
{ double *p,*p1,*q,*pp,*w,*e;
  double s,h,r,t,sv;
  int i,j,k,mm,nm,ms;
  if(m<n) return -1;
  w=(double *)calloc(m+n,sizeof(double)); e=w+m;
  for(i=0,mm=m,nm=n-1,p=a; i<n ;++i,--mm,--nm,p+=n+1){
    if(mm>1){ sv=h=0.;
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; t=1./(w[0]+=h);
        sv=1.+fabs(*p/h);
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  r*=s;
	  for(j=0,q=p+k; j<mm ;q+=n) *q-=r*w[j++];
	 }
        for(j=1,q=p; j<mm ;) *(q+=n)=t*w[j++];
       }
      *p=sv; d[i]= -h;
     }
    if(mm==1) d[i]= *p;
    p1=p+1; sv=h=0.;
    if(nm>1){
      for(j=0,q=p1,s=0.; j<nm ;++j,++q) s+= *q* *q;
      if(s>0.){
	h=sqrt(s); if(*p1<0.) h= -h;
        sv=1.+fabs(*p1/h);
	s+= *p1*h; s=1./s; t=1./(*p1+=h);
	for(k=n,ms=n*(m-i); k<ms ;k+=n){
	  for(j=0,q=p1,pp=p1+k,r=0.; j<nm ;++j) r+= *q++ * *pp++;
	  r*=s;
	  for(j=0,q=p1,pp=p1+k; j<nm ;++j) *pp++ -=r* *q++;
	 }
        for(j=1,q=p1+1; j<nm ;++j) *q++ *=t;
       }
      *p1=sv; e[i]= -h;
     }
    if(nm==1) e[i]= *p1;
   }
  ldvmat(a,v,n); ldumat(a,u,m,n);
  qrbdv(d,e,u,m,v,n);
  for(i=0; i<n ;++i){
    if(d[i]<0.){ d[i]= -d[i];
      for(j=0,p=v+i; j<n ;++j,p+=n) *p= - *p;
     }
   }
  free(w);
  return 0;
}

static int sv2uv(double *d,double *a,double *u,int m,double *v,int n)
{ double *p,*p1,*q,*pp,*w,*e;
  double s,t,h,r,sv;
  int i,j,k,mm,nm,ms;
  if(m<n) return -1;
  w=(double *)calloc(m+n,sizeof(double)); e=w+m;
  for(i=0,mm=m,p=a; i<n ;++i,--mm,p+=n+1){
    if(mm>1){ sv=h=0.;
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; t=1./(w[0]+=h);
        sv=1.+fabs(*p/h);
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  r=r*s;
	  for(j=0,q=p+k; j<mm ;q+=n) *q-=r*w[j++];
	 }
	for(j=1,q=p; j<mm ;) *(q+=n)=w[j++]*t;
       }
      *p=sv; d[i]= -h;
     }
    if(mm==1) d[i]= *p;
   }
  ldumat(a,u,m,n);
  for(i=0,q=a; i<n ;++i){
    for(j=0; j<n ;++j,++q){
      if(j<i) *q=0.;
      else if(j==i) *q=d[i];
     }
   }
  for(i=0,mm=n,nm=n-1,p=a; i<n ;++i,--mm,--nm,p+=n+1){
    if(i && mm>1){ sv=h=0.;
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; t=1./(w[0]+=h);
        sv=1.+fabs(*p/h);
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  for(j=0,q=p+k,r*=s; j<mm ;q+=n) *q-=r*w[j++];
	 }
        for(k=0,p1=u+i; k<m ;++k,p1+=m){
          for(j=0,q=p1,r=0.; j<mm ;) r+=w[j++]* *q++;
	  for(j=0,q=p1,r*=s; j<mm ;) *q++ -=r*w[j++];
	 }
       }
      *p=sv; d[i]= -h;
     }
    if(mm==1) d[i]= *p;
    p1=p+1;
    if(nm>1){ sv=h=0.;
      for(j=0,q=p1,s=0.; j<nm ;++j,++q) s+= *q* *q;
      if(s>0.){
	h=sqrt(s); if(*p1<0.) h= -h;
        sv=1.+fabs(*p1/h);
	s+= *p1*h; s=1./s; t=1./(*p1+=h);
	for(k=n,ms=n*(n-i); k<ms ;k+=n){
	  for(j=0,q=p1,pp=p1+k,r=0.; j<nm ;++j) r+= *q++ * *pp++;
	  for(j=0,q=p1,pp=p1+k,r*=s; j<nm ;++j) *pp++ -=r* *q++;
	 }
	for(j=1,q=p1+1; j<nm ;++j) *q++ *=t;
       }
      *p1=sv; e[i]= -h;
     }
    if(nm==1) e[i]= *p1;
   }
  ldvmat(a,v,n);
  qrbdv(d,e,u,m,v,n);
  for(i=0; i<n ;++i){
    if(d[i]<0.){ d[i]= -d[i];
      for(j=0,p=v+i; j<n ;++j,p+=n) *p= - *p;
     }
   }
  free(w);
  return 0;
}

// d: output size n (singular values)
// a: input size n*m (m>=n, will get trashed)
// u: output size m*m
// v: output size n*n
int svd(double *d, double *a, double *u, int m, double *v, int n)
{
	int (*f)(double*, double*, double*, int, double*, int);
	f = 2*m > 3*n ? sv2uv : svduv;
	return f(d, a, u, m, v, n);
}

#include <assert.h>
#include <stdio.h>

// d: output size n (singular values)
// a: input size n*m (m>=n)
// u: output size m*m
// v: output size n*n
int svd_float(float *d, float *a, float *u, int m, float *v, int n)
{
	assert(m >= n);
	double *dd = malloc(n * sizeof(double));
	double *da = malloc(n * m * sizeof(double));
	double *du = malloc(m * m * sizeof(double));
	double *dv = malloc(n * n * sizeof(double));
	if (!(dd && da && du && dv)) exit(fprintf(stderr,"out of mem"));
	for (int i = 0; i < n*m; i++) da[i] = a[i];
	int r = svd(dd, da, du, m, dv, n);
	if (d) for (int i = 0; i < n; i++) d[i] = dd[i];
	if (u) for (int i = 0; i < m*m; i++) u[i] = du[i];
	if (v) for (int i = 0; i < n*n; i++) v[i] = dv[i];
	free(du);
	free(dv);
	free(da);
	free(dd);
	return r;
}

// d: output size n (singular values)
// a: input size n*m (m>=n)
// u: output size m*m
// v: output size n*n
int svd_double(double *d, double *a, double *u, int m, double *v, int n)
{
	assert(m >= n);
	double *dd = malloc(n * sizeof(double));
	double *da = malloc(n * m * sizeof(double));
	double *du = malloc(m * m * sizeof(double));
	double *dv = malloc(n * n * sizeof(double));
	if (!(dd && da && du && dv)) exit(fprintf(stderr,"out of mem"));
	for (int i = 0; i < n*m; i++) da[i] = a[i];
	int r = svd(dd, da, du, m, dv, n);
	if (d) for (int i = 0; i < n; i++) d[i] = dd[i];
	if (u) for (int i = 0; i < m*m; i++) u[i] = du[i];
	if (v) for (int i = 0; i < n*n; i++) v[i] = dv[i];
	free(du);
	free(dv);
	free(da);
	free(dd);
	return r;
}
