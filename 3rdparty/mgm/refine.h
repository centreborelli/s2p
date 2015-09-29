#include <math.h>

float ParabolafitMinimumOpenCV(const float v[3], float *v_min, float *x_min) {
   // P(-1) = v[0]
   // P(0)  = v[1]
   // P(1)  = v[2]

   // if we can't fit a parabola in the range [-1,1] then we leave the center
   if(v[1]>v[0] && v[1]>v[2]) {
      *x_min=0;
      *v_min=v[1];
      return 0;
   }
   //ax^2 + bx + c = v(x)
   //--------------------
   //c = v[1]
   //a = (v[2] - 2*v[1] + v[0])/2
   //b = (v[2] -   v[0])/2
   float c = v[1];
   float b = (v[2]-v[0])/2;
   float a = (v[2]-2*v[1]+v[0])/2;

   // minimum at: x := -b/2a
   // /// THE FOLLOWING 3 LINES DON'T MAKE ANY SENSE!
   a*=2; b*=2;
   a = a > 1.0 ? a : 1.0;
   float      x = (-b+a)/(2*a);

   if(x >  1) x =  1;
   if(x < -1) x = -1;

   *v_min = (a*x + b)*x + c;
   *x_min = x;
   return x;
}

float ParabolafitMinimum(const float v[3], float *v_min, float *x_min) {
   // P(-1) = v[0]
   // P(0)  = v[1]
   // P(1)  = v[2]

   // if we can't fit a parabola in the range [-1,1] then we leave the center
   if(v[1]>v[0] && v[1]>v[2]) {
      *x_min=0;
      *v_min=v[1];
      return 0;
   }
   //ax^2 + bx + c = v(x)
   //--------------------
   //c = v[1]
   //a = (v[2] - 2*v[1] + v[0])/2
   //b = (v[2] -   v[0])/2
   float c = v[1];
   float b = (v[2]-v[0])/2;
   float a = (v[2]-2*v[1]+v[0])/2;

   // minimum at: x := -b/2a
   float      x = -b/(2*a);
   if(x >  1) x =  1;
   if(x < -1) x = -1;

   *v_min = (a*x + b)*x + c;
   *x_min = x;
   return x;
}

float VfitMinimum(const float v[3], float *v_min, float *x_min) {
   // P(-1) = v[0]
   // P(0)  = v[1]
   // P(1)  = v[2]

   // if we can't fit a V in the range [-1,1] then we leave the center
   if( (v[1] > v[0])  && (v[1] > v[2]) ) {
      *v_min = v[1];
      *x_min = 0;
      return 0;
   }

   // y = P(1) + (x - 1) * slope
   // y = P(-1) + (x - (-1)) * (-slope)
   // x = (P(-1) - P(1))  / (2*slope)
   float slope = v[2] - v[1];
   if ( (v[2] - v[1]) < (v[0] - v[1]) ) 
      slope = v[0] - v[1];

   *x_min = (v[0] - v[2]) / (2*slope);
   *v_min = v[2] + (*x_min - 1) * slope;
   return *x_min;
}

float cubicInterpolate (const float p[4], const float x) {
//   assert(x<=1);
//   assert(x>=0);
   return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}


// find the minimum in the iterval [0,1] of a cubic interpolated function
float CubicfitMinimum(const float p[4], float *out_pmin, float *out_xmin) {
      // trivial minima
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
