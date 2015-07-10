//
// C++ Implementation: stereomatch
//
// Description: eliminate the false matches with epipolar geometry constraint. 
//		See http://www.math-info.univ-paris5.fr/~moisan/epipolar/
//
// Copyright (c) 2007 Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
// Changelog : 2011 Use Eigen SVD <Pierre Moulon>
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "orsa.h"
#include <third_party/Eigen/Cholesky>
#include <third_party/Eigen/Core>
#include <third_party/Eigen/Eigenvalues>
#include <third_party/Eigen/LU>
#include <third_party/Eigen/QR>
#include <third_party/Eigen/SVD>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*-------------------- GENERAL PURPOSE ROUTINES --------------------*/

/* routines for vectors and matrices */

float *vector(int nl, int nh)
{
  float *v;
  
  v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
  if (!v) {
    // mwerror(FATAL,1,"allocation failure in vector()");
    fprintf(stderr, "allocation failure in vector()\n");
    exit(EXIT_FAILURE); /* indicate failure.*/
  }
  return v-nl;
}

float **matrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  float **m;
  
  m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
  if (!m) {
    // mwerror(FATAL,1,"allocation failure 1 in matrix()");
    fprintf(stderr, "allocation failure 1 in matrix()\n");
    exit(EXIT_FAILURE); /* indicate failure.*/
  }
  m -= nrl;
  for(i=nrl;i<=nrh;i++) {
    m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
    if (!m[i]) {
      // mwerror(FATAL,1,"allocation failure 2 in matrix()");
      fprintf(stderr, "allocation failure 2 in matrix()\n");
      exit(EXIT_FAILURE); /* indicate failure.*/
    }
    m[i] -= ncl;
  }
  return m;
}

void free_vector(float *v, int nl, int nh)
{
  free((char*) (v+nl));
	
	nh = nh; // to remove the warning "unused parameter ‘nh’"
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
	
  nch = nch; // to remove the warning "unused parameter ‘nh’"
}

/* Compute the real roots of a third order polynomial */
/* returns 1 or 3, the number of roots found */

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


/* logarithm (base 10) of binomial coefficient */
float logcombi(int k, int n)
{
  double r;
  int i;
 
  if (k>=n || k<=0) return(0.);
  if (n-k<k) k=n-k;
  r = 0.;
  for (i=1;i<=k;i++) 
    r += log10((double)(n-i+1))-log10((double)i);

  return((float)r);
}

/* tabulate logcombi(.,n) */
float *makelogcombi_n(int n)
{
  float *l;
  int k;

  l = (float *)malloc((n+1)*sizeof(float));
  for (k=0;k<=n;k++) l[k]=logcombi(k,n);

  return(l);
}

/* tabulate logcombi(k,.) */
float *makelogcombi_k(int k, int nmax)
{
  float *l;
  int n;

  l = (float *)malloc((nmax+1)*sizeof(float));
  for (n=0;n<=nmax;n++) l[n]=logcombi(k,n);

  return(l);
}


/* get a (sorted) random 7-uple of 0..n-1 */
void random_p7(int *k, int n)
{
  int i,j,j0,r;

  for (i=0;i<7;i++) {
    r = (rand()>>3)%(n-i);
    for (j=0;j<i && r>=k[j];j++) r++;
    j0 = j;
    for (j=i;j>j0;j--) k[j]=k[j-1];
    k[j0]=r;
  }
}

/*-------------------- END OF GENERAL PURPOSE ROUTINES --------------------*/


/* float comparison for qsort() */
//According to http://www.cplusplus.com/reference/clibrary/cstdlib/qsort/, 
//we should have: void qsort ( void * base, size_t num, size_t size, int ( * comparator ) ( const void *, const void * ) ); that means, for "qsort", the "comparator" has two constant void* type input parameters
// static int compf(void *i, void *j)
int compf(const void *i, const void *j)
{ 
  float a,b;

  a = *((float *)i);
  b = *((float *)j);
  return(a<b?-1:(a>b?1:0));
}



/* find the increasing sequence of squared distances to epipolar lines */
/* e[n*2] = distances, e[n*2+1] = indexes (to cast into an int) */

//void matcherrorn(float **F, Flist p1, Flist p2, float *e)
void matcherrorn(float **F, const std::vector<float>& p1, const std::vector<float>& p2, float *e)
{
	int i;
	double x,y,a,b,c,d; // Guoshen Yu, double precision is needed. When the two images are identical, the error under float precision is 0 => log(error)=-inf. 
	
	int pt_size = (p1.size())/2;
	
	for (i = 0; i < pt_size; i++) {
		x = (double) p1[i*2];
		y = (double) p1[i*2+1];
		a = (double) F[1][1]*x+(double) F[1][2]*y+(double) F[1][3];
		b = (double) F[2][1]*x+(double) F[2][2]*y+(double) F[2][3];
		c = (double) F[3][1]*x+(double) F[3][2]*y+(double) F[3][3];
		d = (a*(double) p2[i*2]+b*(double) p2[i*2+1]+c);
		e[i*2] = (d*d)/(a*a+b*b); 
		  	  
		e[i*2+1] = (float)i;
	}
	qsort(e, pt_size, 2*sizeof(float), compf);
}


/*---------- compute the epipolar geometry associated to 7 pairs ----------*/
/*                                                                         */
/*  INPUT: the points are (m1[k[i]*2],m1[k[i]*2+1]), m2... 0<i<7           */
/*                                                                         */
/*  OUTPUT:                                                                */
/*             return the number of roots found, stored in z[]             */
/*   the epipolar matrices are F1[i][j]+z[k]*F2[i][j], 1<=i,j<=3, 0<=k<m   */

// int epipolar(float *m1, float *m2, int *k, float *z, float **F1, float **F2)
int epipolar(std::vector<float>& m1, std::vector<float>& m2, int *k, float *z, float **F1, float **F2)
{
  float a[4];
  int i,j,i2,i3;

  typedef Eigen::MatrixXf Mat;
  Mat c(7,9);
  /* build 9xn matrix from point matches */
  for (i=0;i<7;i++) {
    c(i,0) = m1[k[i]*2  ]*m2[k[i]*2  ];
    c(i,1) = m1[k[i]*2+1]*m2[k[i]*2  ];
    c(i,2) =                      m2[k[i]*2  ];
    c(i,3) = m1[k[i]*2  ]*m2[k[i]*2+1];
    c(i,4) = m1[k[i]*2+1]*m2[k[i]*2+1];
    c(i,5) =                      m2[k[i]*2+1];
    c(i,6) = m1[k[i]*2  ];
    c(i,7) = m1[k[i]*2+1];
    c(i,8) = 1.;
  }
  
  // SVD  
  Eigen::JacobiSVD<Mat> svd(c, Eigen::ComputeFullV);
  // look for the two smallest eigenvalue of c'c 
  typedef Eigen::Matrix<float, 9, 1> Vec9;
  Vec9 F1Vec = svd.matrixV().col(c.cols()-1);
  Vec9 F2Vec = svd.matrixV().col(c.cols()-2);
  
  /* build basis of solutions */
  int cpt = 0;
  for (i=1;i<=3;i++) 
    for (j=1;j<=3;j++)
    {
      F1[i][j] = F1Vec(cpt);
      F2[i][j] = F2Vec(cpt);
      cpt++;
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
  
    return(FindCubicRoots(a,z));
}

void divide_match(const std::vector<Match>& matches, std::vector<float>& p1, std::vector<float>& p2)
{
  float x1, y1, x2, y2;

  p1.clear();
  p2.clear();
  p1.reserve(2 * matches.size());
  p2.reserve(2 * matches.size());
  std::vector<Match>::const_iterator it=matches.begin();
  for(; it != matches.end(); ++it) {
    x1 = (*it).x1; y1 = (*it).y1;
    x2 = (*it).x2; y2 = (*it).y2;
    p1.push_back(x1); p1.push_back(y1);
    p2.push_back(x2); p2.push_back(y2);
  }
}


// float stereomatch(int img_x, int img_y, int size_pt, float* p1, float* p2, float** f, float* index, int* t, int* verb, int* n_flag, int* mode, int* stop)
// float stereomatch(const wxImage& u1, std::vector<float>& p1, std::vector<float>& p2, std::vector<SmallVector<float,3> >& f, std::vector<float>& index, int* t, int* verb, int* n_flag, int* mode, int* stop)
//int main(int argc, char** argv)
float orsa(int width, int height, std::vector<Match>& match, std::vector<float>& index, int t_value, int verb_value, int n_flag_value, int mode_value, int stop_value)
{
 //   int width = 0, height = 0;
  //  int t_value, verb_value, n_flag_value, mode_value, stop_value;
    int *t, *verb, *n_flag, *mode, *stop;
    
    t = (int*)malloc(sizeof(int)); // maximum number of ransac trials
    verb = (int*)malloc(sizeof(int)); //verbose
    n_flag = (int*)malloc(sizeof(int)); // in order NOT to reinitialize the random seed
    mode = (int*)malloc(sizeof(int)); // mode: 0=deterministic 1=ransac 2=optimized ransac (ORSA) 3=automatic
    stop = (int*)malloc(sizeof(int)); // stop as soon as the first meaningful correspondence is found   

    if(width <=0 || height <= 0) {
        std::cerr << "Wrong dimensions of image" << std::endl;
        return 1;
    }
    
    std::vector<float> p1(2*match.size()), p2(2*match.size()), p1_backup(2*match.size()), p2_backup(2*match.size());

    divide_match(match, p1, p2);
    p1_backup = p1;
    p2_backup = p2;
    
    libNumerics::matrix<libNumerics::flnum> f(3, 3);
    f = 0;
	index = std::vector<float>(match.size());
	// Guoshen Yu, 2010.09.23
	// index.clear();

    if(t_value <= 0) {
      std::cerr << "t should be greater than 0" << std::endl;
      return 1;
    }
    *t = t_value;

    if(verb_value == 0) {
      free(verb);
      verb = NULL;
    }
    else
      *verb = verb_value;
    if(verb_value != 1 && verb_value != 0) {
      std::cerr << "verb can only be 0 or 1" << std::endl;
      return 1;
    }
    
    if(n_flag_value == 0) {
      free(n_flag);
      n_flag = NULL;
    }
    else
      *n_flag = n_flag_value;
    if(n_flag_value != 1 && n_flag_value != 0) {
      std::cerr << "n_flag can only be 0 or 1" << std::endl;
      return 1;
    }

    if(mode_value != 0 && mode_value != 1 && mode_value != 2 && mode_value != 3) {
      std::cerr << "mode can only be 0 or 1 or 2 or 3" << std::endl;
      return 1;
    }
    *mode = mode_value;
    
    if(stop_value == 0) {
      free(stop);
      stop = NULL;
    }
    else
      *stop = stop_value;
    if(stop_value != 1 && stop_value != 0) {
      std::cerr << "stop can only be 0 or 1" << std::endl;
      return 1;
    }


  int i,j,i0,k[8],idk[8],*id,m,n,l,minicur=0,miniall=0,delete_index,nid;
  int niter,maxniter,better,cont,optimization;
  float **F1,**F2,**F,nx,ny,z[3],minepscur,minepsall,nfa;
  float norm,minlogalphacur,minlogalphaall,logalpha,logalpha0;
  float *e,*logcn,*logc7,loge0;

  /* initialize random seed if necessary */
  // if (!n_flag) srand48( (long int) time (NULL) + (long int) getpid() );
 // if (!n_flag) srand( (long int) time (NULL) + (long int) getpid() );

  // Guoshen Yu, 2010.09.21: remove getpid which does not exist under Windows
   if (!n_flag) srand( (long int) time (NULL) );
  
  /* check sizes */
  if (p1.size() != p2.size() || p1.size() < 14) {
    fprintf(stderr, "Inconsistent sizes.\n");
    exit(EXIT_FAILURE); /* indicate failure.*/
  } 
  n = p1.size()/2;
  
  /* tabulate logcombi */
  loge0 = (float)log10(3.*(double)(n-7));
  logcn = makelogcombi_n(n);
  logc7 = makelogcombi_k(7,n); 
  
  /* choose mode */
  if (*mode==3) {
    if (logcn[7]<=(float)log10((double)(*t)))
      *mode=0; 
    else *mode=2;
  }
  if (verb) 
    switch(*mode) {
      case 0: 
//	i = (int)(0.5+pow(10.,logc7[n]));
    // Guoshen Yu, 2010.09.22, Windows version
    i = (int)(0.5+pow(10., (double)(logc7[n])));
	printf("I will use deterministic mode (systematic search).\n");
	printf("I have to test %d bases\n",i);
	break;
      case 1:
	printf("I will use pure stochastic mode with no optimization.\n");
	break;
      case 2:
	printf("I will use optimized stochastic mode (ORSA).\n");
    }

    /* normalize coordinates */
    nx = (float)width;
    ny = (float)height;
    norm = 1./(float)sqrt((double)(nx*ny));
    logalpha0 = (float)(log10(2.)+0.5*log10((double)((nx*nx+ny*ny)*norm*norm)));
    for (i=0;i<n;i++) {
	 p1[i*2  ] =  (p1[i*2  ]-0.5*nx)*norm;
	 p1[i*2+1] =  (p1[i*2+1]-0.5*ny)*norm;
	 p2[i*2  ] =  (p2[i*2  ]-0.5*nx)*norm;
	 p2[i*2+1] =  (p2[i*2+1]-0.5*ny)*norm;
    }

    /* allocate and initialize memory */
    F  = matrix(1,3,1,3);
    F1  = matrix(1,3,1,3);
    F2  = matrix(1,3,1,3);
	
      
    //delete_index = (index?0:1);
    delete_index = 0;
    
    e = (float *)malloc(2*n*sizeof(float));
    id = (int *)malloc(n*sizeof(int));
    for (i=0;i<n;i++) id[i]=i;
    nid = n;

    maxniter = (*mode==0?*t:*t-(*t)/10);
    minlogalphaall = minepsall = 10000.;
    niter = optimization = 0;
    i0=0; k[0]=-1; k[7]=n;
	
    /********** MAIN LOOP **********/
    do {

      niter++;

      /* build next list of points */
      if (*mode) random_p7(k,nid);
      else {
	k[i0]++;
	for (i=i0+1;i<=6;i++) k[i]=k[i-1]+1;
      }
      for (i=0;i<7;i++) idk[i]=id[k[i]];

      /* find epipolar transform */
      m = epipolar(p1,p2,idk,z,F1,F2);
      
      /* loop on roots */
      for (;m--;) {

	for (i=1;i<=3;i++) 
	  for (j=1;j<=3;j++) 
	    F[i][j] = F1[i][j]+z[m]*F2[i][j];
      
	/* sort errors */
	matcherrorn(F,p1,p2,e);
      
		
	

	/* find most meaningful subset */
	minepscur = minlogalphacur = 10000.;
	for (i=7;i<n;i++) {
	  logalpha = logalpha0+0.5*(float)log10((double)e[i*2]);
	  nfa = loge0+logalpha*(float)(i-6)+logcn[i+1]+logc7[i+1];
	  if (nfa<minepscur) {
	    minepscur = nfa;
	    minicur = i;
	    minlogalphacur = logalpha; 
	  }
	}
	if (minepscur<minepsall) {
	  /* store best result so far */
	  better = 1;
	  minepsall = minepscur;
	  minlogalphaall = minlogalphacur; 
	  miniall = minicur;
	  // if (f) 
	    for (l=1;l<=3;l++) 
	      for (j=1;j<=3;j++) 
		f(l-1, j-1) = F[l][j];

    // Guoshen Yu, 2010.09.22
	//  for (i=0;i<=minicur;i++) 
	for (i=0;i<minicur;i++) 
	 {		
	    index[i] = e[i*2+1];		
	 }
	} else better=0;


	if (*mode==2 && ((better && minepsall<0.) || 
		    (niter==maxniter && !optimization))) {
	  if (!optimization) maxniter = niter + (*t)/10;
	  optimization = 1;
	  /* final optimization */
	  if (verb) {
	    printf("   nfa=%f size=%d (niter=%d)\n",minepsall,miniall+1,niter);
	    printf("optimization...\n");
	  }
	  nid = miniall+1;

	  // Guoshen Yu, 2010.09.22
	  // for (j=0;j<=miniall;j++)
	 for (j=0;j<miniall;j++)
	    id[j] = (int)(index[j]);
		    }
      }

      /* prepare next list of points */
      if (*mode==0) 
		  for(i0=6;i0>=0 && k[i0]==k[i0+1]-1;i0--){};

      if (stop && minepsall<0.) cont=0;
      else if (*mode==0) cont=(i0>=0?1:0); 
      else cont=(niter<maxniter?1:0);

    } while (cont);

	
    //erase "index", only get the index of the meaningful matchings
    index.erase(index.begin()+miniall+1, index.end());
    if (verb) 
      printf("best matching found:  %d points  log(alpha)=%f  (%d iterations)\n",
	     miniall+1,minlogalphaall,niter);
   



    /* free memory */
    free(id);
    free(e);
    // if (delete_index) mw_delete_flist(index);
    free_matrix(F2,1,3,1,3);
    free_matrix(F1,1,3,1,3);
    free_matrix(F,1,3,1,3);
    free(logc7);
    free(logcn);

    if(t) free(t); 
    if(verb) free(verb); 
    if(n_flag) free(n_flag); 
    if(mode) free(mode); 
    if(stop) free(stop);

//    return 0;
	return(minepsall);
}
