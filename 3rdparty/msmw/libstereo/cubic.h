#include <math.h>
//double cubicInterpolate (double p[4], double x) {
//   assert(x<=1);
//   assert(x>=0);
//   return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
//}
//
//// find the minimum in the iterval [0,1] of a cubic interpolated function
//double cubicMinimum(double p[4]) {
//      // trivial minimums
//      double cmin = fmin(p[1],p[2]);
//
//      double a,b,c,z1,z2,discr;
//      // coefficients of:  ax^2 +bx +c = 0
//      a = 0.5 * 3.0 * (3.0*(p[1] - p[2]) + p[3] - p[0]);
//      b = 2.0 * p[0] - 5.0*p[1] + 4.0*p[2] - p[3];
//      c = 0.5 * (p[2] - p[0]);
//
//      // discriminant
//      discr=b*b -4.0 *a *c;
//      if (discr >=0) {
//         z1 = (-b+ sqrt(discr))/(2.0*a);
//         z2 = (-b- sqrt(discr))/(2.0*a);
//         if(z1>0.0 && z1<1.0) cmin = fmin(z1,cmin);
//         if(z2>0.0 && z2<1.0) cmin = fmin(z2,cmin);
//      }
//      return cmin;
//}
      


inline float cubicInterpolate (const float p[4], const float x) {
//   assert(x<=1);
//   assert(x>=0);
   return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}


// find the minimum in the iterval [0,1] of a cubic interpolated function
float cubicMinimum(float p[4], float *out_pmin, float *out_xmin) {
      // trivial minimums
      float pmin,xmin;
      if(p[1] < p[2]) {
         pmin=p[1];
         xmin=0.0;
      } else {
         pmin=p[2];
         xmin=1.0;
      }

      double a,b,c,z1,z2,discr;
      // coefficients of:  ax^2 +bx +c = 0
      a = 0.5 * 3.0 * (3.0*(p[1] - p[2]) + p[3] - p[0]);
      b = 2.0 * p[0] - 5.0*p[1] + 4.0*p[2] - p[3];
      c = 0.5 * (p[2] - p[0]);

      // discriminant
      discr=b*b -4.0 *a *c;
      if (discr >=0) {
         z1 = (-b+ sqrt(discr))/(2.0*a);
         z2 = (-b- sqrt(discr))/(2.0*a);
         if(z1>0.0 && z1<1.0) {
            float tmp = cubicInterpolate(p,z1);
            if(tmp < pmin) {
               pmin=tmp;
               xmin=z1;
            }
         }
         if(z2>0.0 && z2<1.0) {
            float tmp = cubicInterpolate(p,z2);
            if(tmp < pmin) {
               pmin=tmp;
               xmin=z2;
            }
         }
      }

//      if( *out_pmin < pmin) {
         *out_pmin = pmin;
         *out_xmin = xmin;
//      }
      return pmin;
}
