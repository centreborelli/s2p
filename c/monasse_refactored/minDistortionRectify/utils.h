#ifndef RECTIFY_MD_UTILS_H
#define RECTIFY_MD_UTILS_H
#include <vector>
#include <cstdio>
#include <iostream>
#include <assert.h>
#include "libNumerics/numerics.h"
#include "libNumerics/matrix.h"


#define DEBUG(x) (std::cout<< #x << " = " << x <<std::endl)
#define EPSILON 1e-9

using namespace libNumerics;

namespace rect{

/// Simple class to compute polynomial expressions.
class poly : public std::vector<flnum> {
  public:

    poly(): std::vector<flnum>(){
    }

    poly(int n) : std::vector<flnum>(n){

    }

    poly operator+(const poly& a) const{
      unsigned int tsize = std::max(a.size(),size());
      poly tmp(tsize);
      for (unsigned int i = 0; i < tsize; i++) {
        tmp[i]=0.;
        if(i<a.size()) tmp[i]+=a[i];
        if(i<size()) tmp[i]+=at(i);
      } // i
      return tmp;
    }

    poly operator-(const poly& a) const{
      unsigned int tsize = std::max(a.size(),size());
      poly tmp(tsize);
      for (unsigned int i = 0; i < tsize; i++) {
        tmp[i]=0.;
        if(i<a.size()) tmp[i]-=a[i];
        if(i<size()) tmp[i]+=at(i);
      } // i
      return tmp;
    }
    
    poly& operator*=(const flnum scalar){
      for (unsigned int i = 0; i < size(); i++) {
        at(i)*=scalar;
      } // i
      return *this;
    }

    poly operator*(const flnum scalar) const{
      poly tmp(size());
      for (unsigned int i = 0; i < tmp.size(); i++) {
        tmp[i]=scalar*at(i);
      } // i
      return tmp;
    }

    poly ruffini(double r, double & s) const{
      int n = size()-1;
      poly ret(n);
      ret[n-1] = at(n);
      for (int i = n-2; i >= 0 ; i--) {
         ret[i] = ret[i+1]*r + at(i+1);
      } // i
      s = at(0)+ret[0] *r;
      return ret;
    }


    poly operator*(const poly& a) const{
      unsigned int tsize = a.size()+size()-1;
      poly tmp(tsize);

      for (unsigned int i = 0; i < tmp.size(); i++) {
        tmp[i]=0.;
      } // i

      for (unsigned int i = 0; i < a.size(); i++) {
        for (unsigned int j = 0; j < size(); j++) {
          tmp[i+j]+=a[i]*at(j);
        } // j
      } // i
      return tmp;

    }

    /// Computes a polynomial from a 2X2 matrix such that
    /// ret = [x 1]^T A [x 1]
    static poly fromMatrix(const matrix<flnum> & A){
      assert(A.nrow()==2 && A.ncol()==2);
      poly tmp(3);
      tmp[0]=A(1,1);
      tmp[1]=A(0,1)+A(1,0);
      tmp[2]=A(0,0);
      return tmp;
    }
    
    /// Computes a polynomial from a 2X2 matrix such that
    /// ret = [1 x]^T A [1 x]
    static poly fromMatrix2(const matrix<flnum> & A){
      assert(A.nrow()==2 && A.ncol()==2);
      poly tmp(3);
      tmp[0]=A(0,0);
      tmp[1]=A(0,1)+A(1,0);
      tmp[2]=A(1,1);
      return tmp;
    }

    /// Returns the derivative of a polynomial;
    poly deriv(){
      poly tmp(1);
      tmp[0]=0;
      if (size()==1) 
        return tmp;
      tmp.resize(size()-1);
      for (unsigned int i = 0; i < size()-1; i++) {
        tmp[i]=(i+1)*at(i+1);
      } // i
      return tmp;
    }

    /// Computes the polynomial evaluation at x.
    flnum eval(double x){
      int i, n = size()-1;
      double p;
      p = at(n)*x+at(n-1);
      for(i=n-2;i>=0;i--){
        p=at(i)+p*x;
      }
      return p;
    }

    /// Auxiliary method for debugging purposes
    void print(const char name[]){
      printf("%s = ",name);
      for (int i = size()-1; i >= 0; i--){
        printf("%.1lfx^%d ", at(i),i);
          if (i>0) printf("+ ");
      } // i
      printf("\n");
    }

};

/// Computes the area in the triangle defined by vectors a and b.
flnum getAreaP(const vector<flnum> & a,const vector<flnum> & b);

/// Computes vector e such that F*e=0
vector<flnum> computeEpipole(const matrix<flnum> & F);

/// Cross Product matrix, such that [e]_x f = cross(e,f)
matrix<flnum> getCrossProdMatrix(const vector<flnum> & e);

/// cholesky decomposition of the 2X2 matrix.
matrix<flnum> cholesky2(const matrix<flnum> & A);

//! Finds the Z that minimises (z^T*A1*z)/(z^T*B1*z) + (z^T*A2*z)/(z^T*B2*z)
void minimize_sum_polyratio(matrix<flnum> A1, matrix<flnum> B1, matrix<flnum> A2, matrix<flnum> B2, vector<flnum> & ret, bool verbose=true);

}

#endif 

