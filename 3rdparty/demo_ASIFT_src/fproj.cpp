// Copyright (c) 2007 Lionel Moisan <Lionel.Moisan@parisdescartes.fr>

#include <stdio.h>
#include <math.h>
#include "splines.h"
#include "fproj.h"




/*------------------------ MAIN MODULE ---------------------------------*/

//void fproj(float *in, float *out, int nx, int ny, int *sx, int *sy, float *bg, int *o, float *p, char *i, float X1, float Y1, float X2, float Y2, float X3, float Y3, float *x4, float *y4)
void fproj(vector<float>& in, vector<float>& out, int nx, int ny, int *sx, int *sy, float *bg, int *o, float *p, char *i, float X1, float Y1, float X2, float Y2, float X3, float Y3, float *x4, float *y4)
/*     Fimage in,out;
     int    *sx,*sy,*o;
     char   *i;
     float  *bg,*p,X1,Y1,X2,Y2,X3,Y3,*x4,*y4; */
{
/*  int    n1,n2,nx,ny,x,y,xi,yi,adr,dx,dy;*/
  int    n1,n2,x,y,xi,yi,adr,dx,dy;
  float  res,xx,yy,xp,yp,ux,uy,a,b,d,fx,fy,x12,x13,y12,y13;
  float  cx[12],cy[12],ak[13];
 /* Fimage ref,coeffs; */
//  float *ref, *coeffs;
	vector<float> ref, coeffs;


  /* CHECK ORDER */
  if (*o!=0 && *o!=1 && *o!=-3 && 
      *o!=3 && *o!=5 && *o!=7 && *o!=9 && *o!=11)
  /*  mwerror(FATAL,1,"unrecognized interpolation order.\n"); */
  {	
	  printf("unrecognized interpolation order.\n");
	  exit(-1);
  }

  /* ALLOCATE NEW IMAGE */
/*  nx = in->ncol; ny = in->nrow; */
/*  out = mw_change_fimage(out,*sy,*sx); 
  if (!out) mwerror(FATAL,1,"not enough memory\n"); */


  if (*o>=3) {
/*    coeffs = mw_new_fimage();
    finvspline(in,*o,coeffs); */
	  
//	  coeffs = new float[nx*ny];
	  
	  coeffs = vector<float>(nx*ny);
	  
	finvspline(in,*o,coeffs,nx,ny);
	  
    ref = coeffs;
    if (*o>3) init_splinen(ak,*o);
  } else {
//    coeffs = NULL;
    ref = in;
  }


  /* COMPUTE NEW BASIS */
  if (i) {
    x12 = (X2-X1)/(float)nx;
    y12 = (Y2-Y1)/(float)nx;
    x13 = (X3-X1)/(float)ny;
    y13 = (Y3-Y1)/(float)ny;
  } else {
    x12 = (X2-X1)/(float)(*sx);
    y12 = (Y2-Y1)/(float)(*sx);
    x13 = (X3-X1)/(float)(*sy);
    y13 = (Y3-Y1)/(float)(*sy);
  }



  if (y4) { 
    xx=((*x4-X1)*(Y3-Y1)-(*y4-Y1)*(X3-X1))/((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1));
    yy=((*x4-X1)*(Y2-Y1)-(*y4-Y1)*(X2-X1))/((X3-X1)*(Y2-Y1)-(Y3-Y1)*(X2-X1));
    a = (yy-1.0)/(1.0-xx-yy);
    b = (xx-1.0)/(1.0-xx-yy);
  } 
  else 
    {
      a=b=0.0;
    }

     


  /********** MAIN LOOP **********/

  for (x=0;x<*sx;x++) 
    for (y=0;y<*sy;y++) {
      
      /* COMPUTE LOCATION IN INPUT IMAGE */
      if (i) {
	xx = 0.5+(((float)x-X1)*y13-((float)y-Y1)*x13)/(x12*y13-y12*x13);
	yy = 0.5-(((float)x-X1)*y12-((float)y-Y1)*x12)/(x12*y13-y12*x13);
	d = 1.0-(a/(a+1.0))*xx/(float)nx-(b/(b+1.0))*yy/(float)ny;
	xp = xx/((a+1.0)*d);
	yp = yy/((b+1.0)*d);
      } else {
	fx = (float)x + 0.5;
	fy = (float)y + 0.5;
	d = a*fx/(float)(*sx)+b*fy/(float)(*sy)+1.0;
	xx = (a+1.0)*fx/d;
	yy = (b+1.0)*fy/d;
	xp = X1 + xx*x12 + yy*x13;
	yp = Y1 + xx*y12 + yy*y13;
      }


      /* INTERPOLATION */
      
      if (*o==0) { 
	
	/* zero order interpolation (pixel replication) */
	xi = (int)floor((double)xp); 
	yi = (int)floor((double)yp);
/*	if (xi<0 || xi>=in->ncol || yi<0 || yi>=in->nrow)*/
	if (xi<0 || xi>=nx || yi<0 || yi>=ny)
	  res = *bg; 
	else 
		/* res = in->gray[yi*in->ncol+xi]; */
		res = in[yi*nx+xi];
      } else { 
	
	/* higher order interpolations */
	if (xp<0. || xp>(float)nx || yp<0. || yp>(float)ny) res=*bg; 
	else {
	  xp -= 0.5; yp -= 0.5;
	  xi = (int)floor((double)xp); 
	  yi = (int)floor((double)yp);
	  ux = xp-(float)xi;
	  uy = yp-(float)yi;
	  switch (*o) 
	    {
	    case 1: /* first order interpolation (bilinear) */
	      n2 = 1;
	      cx[0]=ux;	cx[1]=1.-ux;
	      cy[0]=uy; cy[1]=1.-uy;
	      break;
	      
	    case -3: /* third order interpolation (bicubic Keys' function) */
	      n2 = 2;
	      keys(cx,ux,*p);
	      keys(cy,uy,*p);
	      break;

	    case 3: /* spline of order 3 */
	      n2 = 2;
	      spline3(cx,ux);
	      spline3(cy,uy);
	      break;

	    default: /* spline of order >3 */
	      n2 = (1+*o)/2;
	      splinen(cx,ux,ak,*o);
	      splinen(cy,uy,ak,*o);
	      break;
	    }
	  
	  res = 0.; n1 = 1-n2;
	  /* this test saves computation time */
	  if (xi+n1>=0 && xi+n2<nx && yi+n1>=0 && yi+n2<ny) {
	    adr = yi*nx+xi; 
	    for (dy=n1;dy<=n2;dy++) 
	      for (dx=n1;dx<=n2;dx++) 
/*		res += cy[n2-dy]*cx[n2-dx]*ref->gray[adr+nx*dy+dx];*/
			  res += cy[n2-dy]*cx[n2-dx]*ref[adr+nx*dy+dx];
	  } else 
	    for (dy=n1;dy<=n2;dy++)
	      for (dx=n1;dx<=n2;dx++) 
/*		res += cy[n2-dy]*cx[n2-dx]*v(ref,xi+dx,yi+dy,*bg); */
	  res += cy[n2-dy]*cx[n2-dx]*v(ref,xi+dx,yi+dy,*bg,nx,ny); 
	}
      }		
      /* out->gray[y*(*sx)+x] = res; */
		out[y*(*sx)+x] = res;
    }
  //if (coeffs) 
	  /* mw_delete_fimage(coeffs); */
	//  delete[] coeffs; 
}

