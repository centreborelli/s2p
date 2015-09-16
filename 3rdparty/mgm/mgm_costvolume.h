#ifndef COSTVOLUME_H_
#define COSTVOLUME_H_

#include "point.h"

#include "img_interp.h"

#include "img_tools.h"

#include "census_tools.cc"

#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))

// the type of a cost function
typedef float (*cost_t)(Point,Point,const struct Img&,const struct Img&);


inline float computeC_AD_sub( Point p, Point q, const struct Img &u, const struct Img &v) {
	if( !check_inside_image(p,u) ) return INFINITY;
	if( !check_inside_image(q,v) ) return INFINITY;
	float tmp = 0;
	for(int t=0;t<u.nch;t++) {
      float x = u.getpixel(p.x,p.y,t) - v.getpixel(q.x,q.y,t);
      x = __max(x,-x);
      tmp += x ;
   }
	return tmp;
}

inline float computeC_AD( Point p, Point q, const struct Img &u, const struct Img &v) {
	if( !check_inside_image(p,u) ) return INFINITY;
	if( !check_inside_image(q,v) ) return INFINITY;
	float tmp = 0;
	for(int t=0;t<u.nch;t++) {
      float x = val(u,p,t) - val(v,q,t);
      x = __max(x,-x);
      tmp += x ;
   }
	return tmp;
}
inline float computeC_SD( Point p, Point q, const struct Img &u, const struct Img &v) {
   if( !check_inside_image(p,u) ) return INFINITY;
   if( !check_inside_image(q,v) ) return INFINITY;
   float tmp = 0;
   for(int t=0;t<u.nch;t++) {
      float x = val(u,p,t) - val(v,q,t);
      x = __max(x,-x);
      tmp += x * x;
   }
   return tmp;
}


// SMART_PARAMETER(LAMBDA,0.9);

// inline float computeC_COMBI(Point p, Point q, const struct Img &u, const struct Img &v)
// {
// 	float a = computeC_AD(p, q, u, v);
// 	float b = computeC_SC(p, q, u, v);
// 	float l = LAMBDA();
// 	float r = (1 - l) * a + l * b;
// 	return r;
// }


SMART_PARAMETER(CENSUS_NCC_WIN,3)


// fast census (input images must pre-processed by census transform)
inline float computeC_census_on_preprocessed_images( Point p, Point q, const struct Img &u, const struct Img &v)
{
	float r = 0;
	for (int t = 0; t < u.nch; t++)
	{
		float vp = valzero(u, p, t);
		float vq = valzero(v, q, t);

	   r += compute_census_distance_array((uint8_t*)(&vp), (uint8_t*)(&vq), sizeof(float));
	}
   // magic factor each channel contrinutes to r  4 bytes (32bits) 
   // to guarantee that the total cost is below 256 we take r*8 / nch
	return r * 1.0 / u.nch; // magic factor
}


// birchfield and tomasi absolute differences
inline float BTAD( Point p, Point q, int channel, const struct Img &u, const struct Img &v) {
      if( !check_inside_image(p,u) ) return INFINITY;
      if( !check_inside_image(q,v) ) return INFINITY;
      #define min3(a,b,c) (((a)<(b))? (((a)<(c))?(a):(c)) : (((c)<(b))?(c):(b)) )
      #define max3(a,b,c) (((a)>(b))? (((a)>(c))?(a):(c)) : (((c)>(b))?(c):(b)) )

      int t=channel;
      float IL = val(u,p,t);
      float ILp = IL, ILm = IL;
      if (p.x<u.nx-1) ILp = (IL + val(u,p+Point(1,0),t))/2.0;
      if (p.x >= 1  ) ILm = (IL + val(u,p+Point(-1,0),t))/2.0;

      float IR = val(v,q,t);
      float IRp = IR, IRm = IR;
      if (q.x<v.nx-1) IRp = (IR + val(v,q+Point(1,0),t))/2.0;
      if (q.x >= 1  ) IRm = (IR + val(v,q+Point(-1,0),t))/2.0;

      float IminR = min3(IRm,IRp,IR);
      float ImaxR = max3(IRm,IRp,IR);

      float IminL = min3(ILm,ILp,IL);
      float ImaxL = max3(ILm,ILp,IL);

      float dLR =  max3( 0, IL - ImaxR, IminR - IL);
      float dRL =  max3( 0, IR - ImaxL, IminL - IR);

      float BT = __min(dLR, dRL);
      return fabs(BT);
}


// birchfield and tomasi absolute differences
inline float computeC_BTAD( Point p, Point q, const struct Img &u, const struct Img &v) {
	if( !check_inside_image(p,u) ) return INFINITY;
	if( !check_inside_image(q,v) ) return INFINITY;
	float val = 0;
	for(int t=0;t<u.nch;t++)  {
      val += BTAD(p,q,t,u,v);
   }
	return val;
}
// birchfield and tomasi squared differences
inline float computeC_BTSD( Point p, Point q, const struct Img &u, const struct Img &v) {
   if( !check_inside_image(p,u) ) return INFINITY;
   if( !check_inside_image(q,v) ) return INFINITY;
   float val = 0;
   for(int t=0;t<u.nch;t++)  {
      float x = BTAD(p,q,t,u,v);
      val += x*x;
   }
   return val;
}

// Clipped NCC
// NOTE: window size = 3x3
inline float computeC_clippedNCC( Point p, Point q, const struct Img &u, const struct Img &v)
{
   int hwindow = CENSUS_NCC_WIN()/2;
	float r = 0;
   float NCC = 0;
	for (int t = 0; t < u.nch; t++)
	{
      float mu1  = 0; float mu2  = 0;
      float s1   = 0; float s2   = 0;
      float prod = 0;
      int n = 0;
		for (int i = -hwindow; i <= hwindow; i++)
		for (int j = -hwindow; j <= hwindow; j++)
		{
			float v1 = valnan(u, p + Point(i, j), t);
			float v2 = valnan(v, q + Point(i, j), t);
         if (isnan(v1) || isnan(v2)) return INFINITY;
         mu1+=v1;    mu2+=v2;
         s1 +=v1*v1; s2 +=v2*v2; prod+=v1*v2;
         n++;
      }
      mu1/=n; mu2/=n;
      s1 /=n;  s2/=n; prod/=n;

      NCC += (prod - mu1*mu2) / sqrt( __max(0.0000001,(s1 - mu1*mu1)*(s2 - mu2*mu2)) );
	}
   float clippedNCC = u.nch - __max(0,__min(NCC,u.nch));
	return clippedNCC*64;
}



//// global table of all the cost functions
struct distance_functions{
   cost_t f;
   const char *name;
} global_table_of_distance_functions[] = {
         #define REGISTER_FUNCTIONN(x,xn) {x, xn}
         REGISTER_FUNCTIONN(computeC_AD,"ad"),
         REGISTER_FUNCTIONN(computeC_SD,"sd"),
         REGISTER_FUNCTIONN(computeC_census_on_preprocessed_images,"census"),
         REGISTER_FUNCTIONN(computeC_clippedNCC,"ncc"),
         REGISTER_FUNCTIONN(computeC_BTAD,"btad"),
         REGISTER_FUNCTIONN(computeC_BTSD,"btsd"),
         REGISTER_FUNCTIONN(computeC_AD_sub,"ad_sub"),
         #undef REGISTER_FUNCTIONN
         {NULL, ""},
};
int get_distance_index(const char *name) {
   int r=0; // default cost function is computeC_AD (first in table distance_functions)
   for(int i=0; global_table_of_distance_functions[i].f; i++)
      if (strcmp (name,global_table_of_distance_functions[i].name)==0) 
         r=i;
   return r;
}


//// global table of the prefilter names
const char* global_table_of_prefilters[] = {
                                             "none",
                                             "census", 
                                             "sobelx",
                                             "gblur",
                                              NULL,
                                           };
int get_prefilter_index(const char *name) {
   int r=0;
   for(int i=0; global_table_of_prefilters[i]; i++)
      if (strcmp (name,global_table_of_prefilters[i])==0) 
         r=i;
   return r;
}




//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
#ifdef DVEC_ALLOCATION_HACK

#include "dvec2.cc"

struct costvolume_t {
   int npix;
   int ndata;
   struct Dvec *vectors;
   float       *alldata;


   inline Dvec  operator[](int i) const  { 
      return this->vectors[i];
   }
   inline Dvec& operator[](int i)        { 
      return this->vectors[i];
   }

   costvolume_t() {
      int npix=0;
      int ndata=0;
      this->vectors = NULL;
      this->alldata = NULL;
   }

   costvolume_t(const struct costvolume_t &src)
   {
      npix =src.npix;
      ndata=src.ndata;
      vectors = (struct Dvec*) malloc(sizeof(struct Dvec)*npix);
      alldata = (float*)       calloc(ndata, sizeof(float));
      float *baseptr = alldata;
      memcpy(vectors, src.vectors, sizeof(struct Dvec)*npix);
      memcpy(alldata, src.alldata, sizeof(float)*ndata);

      int size = 0;
      for (int i=0; i<npix; i++)  {
         vectors[i].data = baseptr + size;
         size += (int)((int) vectors[i].max - (int) vectors[i].min + 1);
      }
   }

   ~costvolume_t(void)
   {
         if(this->vectors!=NULL) free(this->vectors);
         if(this->alldata!=NULL) free(this->alldata);
   }

};

struct costvolume_t allocate_costvolume (struct Img min, struct Img max) 
{
   int npix=min.nx*min.ny;
   int size=0;
   std::vector< int > pos(npix);
   for (int i=0; i<npix; i++)  {
      pos[i] = size;
      size += (int)((int) max[i] - (int) min[i] + 1);
   }

   struct costvolume_t cv;
   cv.npix =npix;
   cv.ndata=size;
   cv.vectors = (struct Dvec*) malloc(sizeof(struct Dvec)*npix);
   cv.alldata = (float*)       calloc(size, sizeof(float));
   float *baseptr = cv.alldata;

#pragma omp parallel for
   for (int i=0;i< npix;i++) {
      cv.vectors[i].init(min[i], max[i], baseptr + pos[i]);
   }

   return cv;
}



//////////////////////////////////////////////
//////////////////////////////////////////////
#else // DVEC_ALLOCATION_HACK
//////////////////////////////////////////////
////////    USING VECTORS     ////////////////

#include "dvec.cc"

struct costvolume_t {
   std::vector< Dvec > vectors;
   inline Dvec  operator[](int i) const  { return this->vectors[i];}
   inline Dvec& operator[](int i)        { return this->vectors[i];}
};


struct costvolume_t allocate_costvolume (struct Img min, struct Img max) 
{
   struct costvolume_t cv;
   cv.vectors = std::vector< Dvec >(min.npix);
   for (int i=0;i< min.npix;i++) {
      cv[i].init(min[i], max[i]);
   }

   return cv;
}


#endif // DVEC_ALLOCATION_HACK
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

struct costvolume_t allocate_and_fill_sgm_costvolume (struct Img &in_u, // source (reference) image                      
                                                      struct Img &in_v, // destination (match) image                  
                                                      struct Img &dminI,// per pixel max&min disparity
                                                      struct Img &dmaxI,
                                                      char* prefilter,        // none, sobel, census(WxW)
                                                      char* distance,         // census, l1, l2, ncc(WxW), btl1, btl2
                                                      float truncDist)        // truncated differences
{
   int nx = in_u.nx; 
   int ny = in_u.ny;
   int nch= in_u.nch;

   struct Img u(in_u);
   struct Img v(in_v);

   // 0. pick the prefilter and cost functions
   int distance_index  = get_distance_index(distance);
   int prefilter_index = get_prefilter_index(prefilter);
   cost_t cost = global_table_of_distance_functions[distance_index].f;

   // 1. parameter consistency check
   if (distance_index == get_distance_index("census") || prefilter_index == get_prefilter_index("census")) {
       if (TSGM_DEBUG()) printf("costvolume: changing both distance and prefilter to CENSUS\n");
       distance_index  = get_distance_index("census");
       prefilter_index = get_prefilter_index("census");
   }
   if (TSGM_DEBUG()) printf("costvolume: selecting distance  %s\n", global_table_of_distance_functions[distance_index].name);
   if (TSGM_DEBUG()) printf("costvolume: selecting prefilter %s\n", global_table_of_prefilters[prefilter_index]);
   if (TSGM_DEBUG()) printf("costvolume: truncate distances at %f\n", truncDist);

   // 2. apply prefilters if needed
   if (prefilter_index == get_prefilter_index("census")) {
      int winradius = CENSUS_NCC_WIN() / 2;
      if (TSGM_DEBUG()) printf("costvolume: applying census with window of size %d\n", winradius*2+1);
      u = census_transform(in_u, winradius);
      v = census_transform(in_v, winradius);
   }
   if (prefilter_index == get_prefilter_index("sobelx")) {
      if (TSGM_DEBUG()) printf("costvolume: applying sobel filter\n" );
      float sobel_x[] = {-1,0,1, -2,0,2, -1,0,1};
      u = apply_filter(in_u, sobel_x, 3, 3, 1);
      v = apply_filter(in_v, sobel_x, 3, 3, 1);
   }
   if (prefilter_index == get_prefilter_index("gblur")) {
      if (TSGM_DEBUG()) printf("costvolume: applying gblur(s=1) filter\n" );
      u = gblur_truncated(in_u, 1.0);
      v = gblur_truncated(in_v, 1.0);
   }

   // 3. allocate the cost volume 
   struct costvolume_t CC = allocate_costvolume(dminI, dmaxI);

   // 4. apply it 
   #pragma omp parallel for
   for(int jj=0; jj<ny; jj++) for(int ii=0; ii<nx; ii++)
   {
      int pidx = (ii + jj*nx);
      int allinvalid = 1;

      for(int o=CC[pidx].min;o<=CC[pidx].max;o++) 
      {
         Point p(ii,jj);      // current point on left image
         Point q = p + Point(o,0); // other point on right image
         // 4.1 compute the cost 
         float e = truncDist * u.nch;
         if (check_inside_image(q, v)) 
            e = cost(p, q, u, v);
         // 4.2 truncate the cost (if needed)
         e = __min(e, truncDist * u.nch);
         // 4.3 store it in the costvolume
         CC[pidx].set_nolock(o, e); // pragma omp critic is inside set
         if(isfinite(e)) allinvalid=0;
      }
      // SAFETY MEASURE: If there are no valid hypotheses for this pixel 
      // (ie all hypotheses fall outside the target image or are invalid in some way)
      // then the cost must be set to 0, for all the available hypotheses
      // Leaving inf would be propagated and invalidate the entire solution 
      if (allinvalid) {
         for(int o=CC[pidx].min;o<=CC[pidx].max;o++) 
         {
            Point p(ii,jj);      // current point on left image
            Point q = p + Point(o,0); // other point on right image
            CC[pidx].set_nolock(o, 0); // pragma omp critic is inside set
         }
      }
   }
   return CC;
}


#endif //COSTVOLUME_H_
