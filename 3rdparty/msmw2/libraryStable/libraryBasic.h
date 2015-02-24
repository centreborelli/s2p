#ifndef _LIBRARY_BASIC_H_
#define _LIBRARY_BASIC_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cmath>
#include <cassert>
#include <vector>
#include <getopt.h>


extern "C" {
#include "iio.h"
}


#define USE_FFTW
#ifdef USE_FFTW
#include <fftw3.h>
#endif



namespace libIIPStable {
	
    
    
#ifndef PI	
#define PI 3.14159265358979323846264338327
#endif
    
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
    
#define dTiny 1e-10
#define fTiny 0.00000001f
//#define fLarge 100000000.0f
#define fLarge INFINITY
#define dLarge 1e+10
    
#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )
    
#define iipHorizontal 0
#define iipVertical 1
    

    
    
    
    //
	//! Class definitions
	//
	
	class cflimage;
	class flimage;
	class laMatrix;
    
    
    
	//
	//! Value operations
	//
	
	void fpClear(float *fpI,float fValue, int iLength);
	void fpCopy(float *fpI,float *fpO, int iLength);
	
	
	void dpCopy(float *dpI,float *dpO, int iLength);
	void dpClear(double *dpI,double dValue, int iLength);
	
	
	float fpMax(float *u,int *pos, int size);
	float fpMin(float *u,int *pos,int size);
	
	double dpMax(double *u,int *pos, int size);
	double dpMin(double *u,int *pos,int size);
	
	float fpVar(float *u,int size);	
	float fpMean(float *u,int size);
    float fpMedian(float *u,int size);
    
	double dpVar(double *u,int size);	
	double dpMean(double *u,int size);
	
	
    void fpCombine(float *u,float a,float *v,float b, float *w,  int size);
    void dpCombine(double *u,double a,double *v,double b, double *w,  int size);
	
    
    
    
    void fiImageDrawCircle(float *igray, int pi,int pj, double radius, float value, int width, int height);
    void fiImageDrawLine(float *igray, int a0, int b0, int a1, int b1, float value, int width, int height);
    
    

    void fpBinarize(float *u, float *v, float value, int inverse, int size);
    
    
    float fpDistLp(float *u, float *v, int p, int size);
    float fpDistLp(float *u, float *v, float *m, int p, int size);

    
    
	//
	//! Float pointer ordering
	//
	
	void fpQuickSort(float *fpI, int iLength, int inverse = 0 );
	void fpQuickSort(float *fpI, float *fpO, int iLength, int inverse = 0);
	
	
	
	//
	//! Gradient based
	//
	
    void fiComputeImageGradient(float * fpI,float *fpXgrad, float *fpYgrad, float *fpGrad, float * fpOri, int iWidth, int iHeight, char cType = 'f');
    
	void fiComputeImageGradient(float * fpI, float *fpGrad, float * fpOri, int iWidth, int iHeight, char cType = 'f');
    
	void fiComputeImageGradient(float * fpI, float *fpGrad, int iWidth, int iHeight, char cType = 'f');
	
    
	
	//
	//! Noise
	//
	
	void fpAddNoiseGaussian(float *u, float *v, float std, long int randinit, int size); 
	void fpAddNoiseGaussianAfine(float *u,float *v, float a,float b,long int randinit, int size);
	
	
	
	//
	//! Histogram
	//
	
	float* fpHisto(float* input, float *iminim, float *imaxim, int *n, float *s, int size, char flag);
    
#define HSP_NB_CASES 500    	
	void fk_histogram_midway(float *in1, float *in2, float *out1, float *out2, int width1, int height1, int width2, int height2);
	void fk_histogram_midway_sequence(float **in, float **out, int nImages, int width, int height);
	void fk_histogram_specification(float *in1, float *in2, float *out, int width1, int height1, int width2, int height2);
	
	
	
    
	
    //
	//! Image Convolution
	//
    
	float* fiFloatGaussKernel(float std, int & size);
	float* fiFloatDirectionalGaussKernel(float xsigma, float ysigma, float angle, float *kernel, int kwidth, int kheight);
	
    void fiFloatHorizontalConvolution(float *u, float *v, int width, int height, float *kernel, int ksize, int boundary);
    void fiFloatVerticalConvolution(float *u, float *v, int width, int height, float *kernel,int ksize, int boundary);
	void fiGaussianConvol(float *u, float *v, int width, int height, float sigma);
	void fiConvol(float *u,float *v,int width,int height,float *kernel,int kwidth,int kheight);
	void fiSepConvol(float *u,float *v,int width,int height,float *xkernel, int xksize, float *ykernel, int yksize);
	
	
    
    
    
    
	//
    //! Sampling functions
    //
    void fiImageSample(float *igray,float *ogray, int factor, int width, int height);
    void fiImageSampleCenter(float *igray,float *ogray, int factor, int width, int height);
    void fiImageSampleAglomeration(float *igray,float *ogray, int factor, int width, int height);
    


    
	//
    //! Patch Processing
    //
    
    void fiPatchStatistics(float *fpIn, float *fpMinV, float *fpMaxV, float *fpMeanV, float *fpVarV, float *fpMedianV, float fRadius, int iWidth, int iHeight);
    void fiComputeIntegralImage(float *in, float *out, int width, int height);
    
    void fiPatchMin(float *fpIn, float *fpMinV, float fRadius, int iWidth, int iHeight);
    void fiPatchMax(float *fpIn, float *fpMaxV, float fRadius, int iWidth, int iHeight);
    void fiPatchMean(float *fpIn, float *fpMeanV, float fRadius, int iWidth, int iHeight);
    void fiPatchVar(float *fpIn, float *fpVarV, float fRadius, int iWidth, int iHeight);
    void fiPatchMedian(float *fpIn, float *fpMedianV, float fRadius, int iWidth, int iHeight);
    
    
    
    //
	//! Image Conversion
	//
	
#define COEFF_YR 0.299
#define COEFF_YG 0.587
#define COEFF_YB 0.114
	void fiRgb2Yuv(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height);
	void fiYuv2Rgb(float *r,float *g,float *b,float *y,float *u,float *v, int width,int height);
	
	
	void fiRgb2YuvO(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height);
	void fiYuvO2Rgb(float *r,float *g,float *b,float *y,float *u,float *v, int width,int height);
	
	
	
	
    
	
	//
	//! Patch Distances
	//
    
	float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int width0, int width1);
    float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int width0, int width1);
    float fiL2FloatDist(float **u0,float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int channels, int width0, int width1);
	
    
	float fiL2FloatWDist ( float * u0, float *u1, int i0, int j0,int i1,int j1,int xradius, int yradius, float *kernel, int width0, int width1);	
    float fiL2FloatWDist ( float ** u0, float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, float *kernel, int channels, int width0, int width1);
	
 
    
    //
    //! FFT Stuff
    //
    void fft1d(float *Xr,float *Xi,float *Yr,float *Yi, int inverse, int n);
	void fft2d(float *in_re, float *in_im,float *out_re, float *out_im, int i_flag, int width, int height);
	void fiFFTShearPascal(double dAmountOfShear, float *pFloatImageInput,float *pFloatImageOutput, int iAxis, float fDelta, int PleaseDoNotCenter, int width, int height);
	void fiFFTZoom(float *in, float *out, float zoomx, float zoomy, int width, int height);
	
    
    
    
    
    ///////////////////////
	//! Miscellaneous
	///////////////////////
#define LUTMAX 30.0
#define LUTMAXM1 29.0
#define LUTPRECISION 1000.0
	
	
	void  wxFillExpLut(float *lut, int size);
	float wxSLUT(float dif, float *lut);
    
    
    
    
    
    
    //
    //! Splines Stuff
    //
    float vMW(float *in,int x,int y,float bg, int width, int height);
    void keysMW(float *c,float t,float a);
    void spline3MW(float *c,float t);
    void init_splinenMW(float *a,int n);
    float ipowMW(float x,int n);
    void splinenMW(float *c,float t,float *a,int n);
    double initcausalMW(double *c,int n,double z);
    double initanticausalMW(double *c,int n,double z);
    void invspline1DMW(double *c,int size,double *z,int npoles);
	void finvsplineMW(float *in,int order,float *out, int width, int height);
    float evaluate_splineMW(float *input, float *ref, float xp, float yp, float *ak, int order, float bg, int width, int height);
    
    
    
    
    //
    //! Geometrical Transformations
    //
    void compute_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, laMatrix &H);
    int compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, laMatrix &H, int *accorded);
    int compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, laMatrix &H);
    void apply_planar_homography(float *input, int width, int height, laMatrix &H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight);
    void apply_planar_homography_zoom(float *input, int width, int height, laMatrix &H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight, float fZoom);
    void apply_zoom(float *input, float *out, float zoom, int order, int width, int height);    
    
    
	
    
    
    
    //
    //! Image Classes
    //
    
	
	class cflimage {
		
		
	protected:
		
		int d_c, d_w, d_h, d_wh, d_whc;		// channels, width, height, width*height, width*height*channels			
		float *d_v;							// pointer		
		
		
		char name[128];
		float visuMin, visuMax;				// useful only for visualization with mxview
		
		
	public:
		
		
		//
		//! Construction / Memory
		//
		
		cflimage();
		cflimage(int w, int h, int c);
		cflimage(int w, int h, float *igray);
		cflimage(int w, int h, float *ired, float *igreen, float *iblue);
		cflimage(flimage &red, flimage &green, flimage &blue);
		
		
		void create(int w, int h, int c);
		
		cflimage(const cflimage& im);
		
		void erase();
		virtual ~cflimage();
		
		
		
		//
		//! Operators
		// 
		
		cflimage& operator= (const cflimage& im);
		
		
		
		//
		//! Load / Save
		//
		void load(const char *filename);
		void save(const char *filename);
		
		
		
		//
		//! Get Basic Data	
		//
		
		int c() const {return d_c;} 
		int w() const {return d_w;} 	
		int h() const {return d_h;} 
		int wh() const {return d_wh;} 
		int whc() const {return d_whc;} 
		
		
		float * v() {return d_v;}
		float * v(int i) {/*assert(i>= 0 && i < d_c);*/ return &d_v[i * d_wh];}
		
		inline float operator[](int i) const {/*assert(i>=0 && i < d_whc);*/ return d_v[i];} 	
		inline float& operator[](int i) {/*assert(i>=0 && i < d_whc);*/ return d_v[i];}	
		
		
		flimage getChannel(int i);
		operator  flimage();
		
		
		
		int isSameSize(cflimage &inIm);
		
		char * getName(){return name;}
		void  setName(const char *iname){strcpy(name, iname);}
		
		
		float getVisuMin(){return visuMin;}
		float getVisuMax(){return visuMax;}
		
		void setVisuMin(float a){visuMin = a;}
		void setVisuMax(float a){visuMax = a;}

        
		//
		//! Math
		//
		
		cflimage& operator= (float a);
		void operator-=(float a);
		void operator+=(float a);
		void operator*=(float a);
		
		float  min();
		float  max();
		
		float  min_channel (int i);
		float  max_channel (int i);
		
		void  max (float M);
		void  min (float m);	
		
		void normalizeL1();
		
		void rint();
		void abs();
		
		void thre(float m, float M);
		
        
		//
		//! Value operations
		//
		void addGaussianNoise(float std);
        void addGaussianNoiseSD(float std, float sdstd);
        
		
		cflimage gradient(char cType);
		cflimage xgradient(char cType);
		cflimage ygradient(char cType);
		cflimage gradient(cflimage &orientation, char cType);
		cflimage gradient(cflimage &xgrad, cflimage &ygrad, cflimage &orientation, char cType);
        
		
		
		//
		//! Color Conversion
		//
		
		flimage getGray();
		flimage getGray(float wr, float wg, float wb);
		
		void Rgb2Yuv(int iflagOrto);
		void Yuv2Rgb(int iflagOrto);
		
		
		
        cflimage binarize(float value, int inverse);

        
		//
		//! Block operations
		//
		
		cflimage copy(int ipx, int ipy, int iw, int ih);
		void paste(const cflimage &im, int x, int y);
		
		cflimage padding(int w, int h, float fValue);
		cflimage append(cflimage &imIn, int extension);
        
		
        
        
        //
		//! Subsampling & zooming & convolution
		//
		
        cflimage convolveGauss(float fSigma);
		cflimage convolve(flimage &kernel);
		
        
        cflimage subSample(int fFactor, int flagCenter);
		cflimage subSample(int fFactor, float fSigma, int flagCenter);

		cflimage subSampleAglomeration(int iFactor);
 
        
        
        //
        //! Block operations 
        //
        
		cflimage  patchMean(float fRadius);
		cflimage  patchVar(float fRadius);
        cflimage  patchMin(float fRadius);
		cflimage  patchMax(float fRadius);
		cflimage  patchMedian(float fRadius);
		
        
        cflimage  patchMean(flimage &kernel);
        cflimage  patchVar(flimage &kernel);
 
         // with list kernels
        cflimage  patchListMean(flimage &kernel);
       
        //
		//! Geometrical transforms
		//
		
       
        
		cflimage mirror(int Orientation);
		void fftRot(float angle, float xtrans , float ytrans, int flagNoCenter, int flagSymmetric);
		void fftTrans(float xtrans,float ytrans, int flagNoCenter, int flagSymmetric);
		
		void fftShear(float dAmountOfShear, int iOrientation, float fDelta, int flagNoCenter, int flagSymmetric);
		
		
		
		// 
		//! Zooming
		//
		
		cflimage fftUpSample(float fFactorx);
		cflimage fftUpSample(float fFactorx, float fFactory);
		
		
		cflimage upSampleSplines(float fFactor, int order = 0);
		
		
        //
		//! Patch distances
		//
		
		friend float distanceL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy);
		
        
        friend float distancePatchL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h);
		friend float distancePatchWL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel);
		
		
		
		friend float distancePatchL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h, cflimage &mean1, cflimage &mean2);
		
		friend float distancePatchWL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2);
		
		
		friend float distancePatchL1(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h);
		friend float distancePatchWL1(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel);
		

		
		
		friend float distancePatchListWL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &Lkernel);
		friend float distancePatchListWL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &Lkernel,  cflimage &mean1, cflimage &mean2);
        
        
	};
	
	
	
	class flimage : public cflimage
	{
		
	public:
		
      // list structures for representing List kernels
      int   list_len;
		int   *offx;				// pointer		
		int   *offy;				// pointer		
		float *offval;				// pointer		
		
		
		flimage();
		flimage(int w, int h);
		flimage(int w, int h, float *ptr);
		
		flimage(const flimage& im);

		
		void create(int w, int h);
		
		
		using cflimage::operator=;
		flimage& operator=(const flimage& im);
		
		
		float operator()(int i, int j) const {return this->d_v[j * this->d_w + i];} 	
		float& operator()(int i, int j) { return this->d_v[j * this->d_w + i];}	
		
		
		
		virtual ~flimage() {};
		
	};
	

    
    
    
    class cflmovie
    {
        
    public:
        
        
        // Constructor && Destructor 
        cflmovie();
        cflmovie(const char * ifilename);							/// Reads a cflmovie file
        cflmovie(const char * ifilename,  int inframes);   			/// Opens a cflmovie file in writable format
        
        
        cflmovie(const cflmovie & imovie);
        
        ~cflmovie();
        
        
        cflmovie& operator= (const cflmovie& im);
        
        
        // Load
        void load(const char * ifilename);											/// Reads a cflmovie file
        
        
        // Main variables
        int n() {return d_n;}	  
        
        
        
        // Getting image at position fpos 
        cflimage getframe(int fpos);
        
        
        // Write image into film 
        void write(cflimage &frame); 
        
        
    private:
        
        int d_n;
        
        char* filename;
        char* fileadress;				/// Used for reading/writting image names in format .mov
        
        bool  writable;				/// Is current cflmovie writable
        char ** strings;   			/// In case we read each time the image from the file
        
        int pos;        				/// Current position for writting
        std::ofstream outfile;     		/// Output mov file for writing
        
    };
    
    
    
    
    /// List of 3d points with a color
    class flData3D
    {
        
    public:
        
        
        ///////////////// Constructor && Destructor 
        flData3D();
        flData3D(flData3D &inData);
        flData3D(int inpoints, float *idx, float *idy, float *idz, float *idr, float *idg, float *idb);
        flData3D(const char * filename);
        
        ~flData3D();
        
        
        ///////////////// Operations
        flData3D& operator= (const flData3D &in);
        
        void loadFile(const char * filename);   
        void SaveFile(const char * dataname);
        
        
        /////////////// Allocate Desallocate 	
        void allocate_coordinates(int n);
        void allocate_color(int n);
        void desallocate_coordinates();
        void desallocate_color();
        
        
        ////////////// Get pointers
        float *getX(){return dx;}
        float *getY(){return dy;}   
        float *getZ(){return dz;}
        
        float *getR(){return dr;}
        float *getG(){return dg;}
        float *getB(){return db;}
        
        
        ///////////// Get number of points
        int getSize(){return npoints;}
        
        
        bool Ok(){return ok;}
        bool Color(){return color;}
        
        
        
        void normalize();
        //void normalize(float bx, float by, float bz, float amp);
        //void getNormParam(float *bx, float *by,float *bz, float *amp);
        
        
        
    public:
        
        int npoints;       
        float *dx, *dy, *dz;
        float *dr, *dg, *db;
        
        bool ok;
        bool color;
        
    };
    
    
    
    
    
    
    
    
    
    
    //
    //! Numerics Classes
    //
    
    
    inline const  double SQR(const double a) {return a*a;}

	
    inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
    
    
    
    inline void SWAP(double &a, double &b)
	{double dum=a; a=b; b=dum;}
    
    
    
    class laMatrix;
    class laVector;
    
    
    
    class laVector {
    protected:
        int d_n;	// size of array. upper index is nn-1
        double *d_v;
    public:
        
        laVector();
        explicit laVector(int n);								// Zero-based array
        laVector(const double &a, int n);						// initialize to constant value
        laVector(const double *a, int n);						// Initialize to array
        
        laVector(const laVector &rhs);							// Copy constructor
        laVector & operator=(const laVector &rhs);				// assignment
        laVector & operator=(const double &a);					// assign a to every element
        
		double & operator[](const int i);				// i'th element
		const double & operator[](const int i) const;
        
		
		
		double * v();
		
		int size() const;
        ~laVector();
        
        friend class laMatrix;
        
    };
    
    
    
    
	
    class laMatrix 
    {
        
    protected:
        int d_n;			// nrows
        int d_m;			// ncols
        double **d_v;		// pointer
        
    public:
        
        
        //! Construction / Destruction
        
		laMatrix();
        laMatrix(int n, int m);			
        laMatrix(const double &a, int n, int m);	
        laMatrix(const double *a, int n, int m);
        laMatrix(const laMatrix &rhs);		
        
        laMatrix & operator=(const laMatrix &rhs);	
        laMatrix & operator=(const double &a);		
        
        ~laMatrix();
        
        
		
        //! Basic operators
        
        double * operator[](const int i);	//subscripting: pointer to row i
        const double * operator[](const int i) const;
        
        double ** v();
        
        int nrows() const;
        int ncols() const;
        
        
        void create(int n, int m);
        
        
        
        //! Non member Arithmetic Operations
		
		friend laMatrix operator*  (double a, const laMatrix& rhs);                              // scalar matrix product
		friend laMatrix operator/  (const laMatrix& lhs, double a);                              // matrix scalar division
		friend laMatrix operator+  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix sum
		friend laMatrix operator-  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix subtraction
		friend laMatrix operator*  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix product
		
		friend laVector operator*  (const laMatrix& lhs, const laVector & rhs);                  // matrix vector product
		
        
		
		//! Other
		laMatrix transposed();
        
        laMatrix copyBlock(int i0, int j0, int rowb, int colb);
        
        friend class laVector;
        
    };
    
    
    
	void luinv(laMatrix &a, laMatrix &inv);
    void lusolve(laMatrix &a, laVector &x, laVector &b);

    
    void invCC(laMatrix &inv);          //! Inversion of definite positive matrices
    
    
    void compute_svd(laMatrix &A, laMatrix &m_U, laMatrix &m_V, laVector &m_W);
    
    void compute_pca_svd(laMatrix &X, laVector &S, laMatrix &V, laMatrix &U);
    
    
    
    
}
    
    
    







//
//! Parser
//
typedef struct optstruct
{	
    const char *gp;
    int flag;
    const char *defvalue;
    const char *value;
    const char *comment;
    
} OptStruct;



typedef struct parstruct
{
    const char * name;
    char * value;
    const char * comment;
} ParStruct;

typedef struct parstruct2
{
    const char * name;
    char * value;
    const char * comment;
    int isoptional;
} ParStruct2;



int parsecmdline(char *pname,
                 char *function,
                 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
                 std::vector <ParStruct*> & par);



int parsecmdline2(char *pname,
                 char *function,
                 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
                 std::vector <ParStruct2*> & par);



#endif
