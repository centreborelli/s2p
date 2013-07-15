#include "minDistortionRectify/utils.h"
#include <assert.h>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <cfloat>

namespace rect{

/// Finds polynomial root near x0.
/// Newton-raphson method
double solve_poly(const poly &c, double x0, bool verbose){
  int N= 4000000;
  int i,n = c.size()-1;
  double p=DBL_MAX,pd=DBL_MAX,x = x0;

  double x1=DBL_MAX ,x2=DBL_MAX;
  double p1=DBL_MAX ,p2=DBL_MAX;
  double pd1=DBL_MAX,pd2=DBL_MAX;

  while(N--){
    pd2=pd1;
    pd1=pd;
    p2=p1;
    p1=p;

    p = c[n]*x+c[n-1];
    pd=c[n];
    for(i=n-2;i>=0;i--){
      pd=p+pd*x;
      p=c[i]+p*x;
    }

    if(std::abs(pd)<1e-6 || x==x1){
      break;
    }

    x2=x1;
    x1=x;
    x-=p/pd;
  }

  if (verbose){
    std::cerr << "p = "<<p << " , " <<p1 << " , " <<p2 << " ... " <<std::endl;
    std::cerr << "p'= "<<pd<< " , " <<pd1<< " , " <<pd2<< " ... " <<std::endl;
    std::cerr << "x = "<<x << " , " <<x1 << " , " <<x2 << " ... " <<std::endl;
  }
  return x;
}

vector<flnum> computeEpipole(const matrix<flnum> & F){
  SVD svd(F);
  // picks the vector corresponding to the lowest singular value of F.
  vector<flnum> y(svd.V().col(2));

  // Try to return e such that w component is equal to 1.
  // Otherwise, returns one such that w-component is zero and
  // the maximum element is 1.
  double u = std::abs(y(0));
  double v = std::abs(y(1));
  double w = std::abs(y(2));
  double amax = std::max(u,std::max(v,w));
  if( w == 0. || amax/w > 1e+7){
    y(2)=0.;
    return y*(1./amax);
  }
  return y/y(2);
}



/// Finds partial solution for z. 
/// Maximises expression z^TBz/z^TAz 
void minimize_polyratio2( matrix<flnum> A, matrix<flnum> B, vector<flnum> & ret ){
  matrix<flnum> D ( cholesky2(A) );
  matrix<flnum> M ( D.inv().t() * B * D.inv());
  //Extracting Eigenvector
  SVD svd(M);
  //Picks vector corresponding to the highest singular value
  vector<flnum> y(svd.V().col(0)); 
  ret = D.inv()*y;
}


// Minimizes the rational expression p/q given an initial guess.
// Returns the minimum value found.
double solve_r4(poly p, poly q, const double x0  ,double& xmin, bool verbose){

  poly deriv = p*q.deriv()*(-1) + q*p.deriv();

  // Try to solve derivative=0 such that it does not solves the 
  // denominator q at the same time.
  std::vector<double> candidates;
  while(1){
    double candidate = solve_poly(deriv, x0, verbose);
    candidates.push_back(candidate);
    double residue=0;
    deriv = deriv.ruffini(candidate,residue);
    if (deriv.size() == 1){
      break;
    }
  }

  xmin = DBL_MAX;
  double best_val = DBL_MAX;
  for (unsigned int i = 0; i < candidates.size(); i++) {
    double candidate = candidates[i];
    double val = std::abs(p.eval(candidate)/q.eval(candidate));
    if (val < best_val){
      best_val = val;
      xmin = candidate;
    }
  } // i
  return p.eval(xmin)/q.eval(xmin);
}


void minimize_sum_polyratio(matrix<flnum> A1, matrix<flnum> B1, matrix<flnum> A2, matrix<flnum> B2, vector<flnum> & ret, bool verbose){
  // Computing initial guess of lambda
  vector<flnum> zs1(2); 
  vector<flnum> zs2(2); 
  minimize_polyratio2(A1, B1, zs1);
  minimize_polyratio2(A2, B2, zs2);  
  vector<flnum> z0 =  0.5*(zs1/sqrt(zs1.qnorm()) + zs2/sqrt(zs2.qnorm()));

  //first case: mu=1
  double lambda0 = z0(0)/z0(1);
  poly p1 = poly::fromMatrix(A1);
  poly q1 = poly::fromMatrix(B1);
  poly p2 = poly::fromMatrix(A2);
  poly q2 = poly::fromMatrix(B2);
  poly pnum=(p1*q2 +q1*p2);
  poly pden=(q2*q1);
  double lambda = DBL_MAX;
  double lambda_distortion = solve_r4(pnum,pden,lambda0, lambda, verbose);


  //second case: lambda=1
  double mu0 = z0(1)/z0(0);
  p1 = poly::fromMatrix2(A1);
  q1 = poly::fromMatrix2(B1);
  p2 = poly::fromMatrix2(A2);
  q2 = poly::fromMatrix2(B2);
  pnum=(p1*q2 +q1*p2);
  pden=(q2*q1);
  double mu = DBL_MAX;
  double mu_distortion =solve_r4(pnum,pden,mu0, mu, verbose);

  if (verbose){
    std::cerr<< "lambda0 -> lambda* :  " << lambda0 <<" -> " << lambda << std::endl;
    std::cerr<< "Distortion(lambda) = " << lambda_distortion << std::endl;
    std::cerr<< std::endl;
    std::cerr<< "mu0 -> mu* :  " << mu0 <<" -> " << mu << std::endl;
    std::cerr<< "Distortion(mu) = " << mu_distortion << std::endl;
    std::cerr<< std::endl;
  }

  if (lambda_distortion < mu_distortion)
  {
    ret(0) = lambda;
    ret(1) = 1.;
  }else
  {
    ret(0) = 1;
    ret(1) = mu;
  }
}

flnum getAreaP(const vector<flnum> & a,const vector<flnum> & b){
  return .5*std::abs(a(0)*b(1)- a(1)*b(0));
}


matrix<flnum> getCrossProdMatrix(const vector<flnum> & e){
  matrix<flnum> A(3,3);
  A=0;
  A(0,1) = -e(2);
  A(1,0) = e(2);
  A(0,2) = e(1);
  A(2,0) = -e(1);
  A(1,2) = -e(0);
  A(2,1) = e(0);
  return A;
}

matrix<flnum> cholesky2(const matrix<flnum> & A){
  assert(A.nrow()==2&&A.ncol()==2);
  matrix<flnum> ret(2,2);
  ret(0,0) = sqrt(A(0,0));
  ret(0,1) = A(0,1)/ret(0,0);
  ret(1,0) = 0;
  ret(1,1) = sqrt(A(1,1) - ret(0,1)*ret(0,1));
  return ret;
}

}


