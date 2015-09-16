/* Copyright 2014, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include <numeric>
#include <algorithm>
#include <vector>
#include <cstring>
#include "assert.h"

#include "smartparameter.h"

//// a structure to wrap images 
#include "img_interp.h"
#include "point.h"

#include "img_tools.h"

// // not used here but generally useful 
// typedef std::vector<float> FloatVector;


SMART_PARAMETER(TSGM_DEBUG,0)

/********************** COSTVOLUME *****************************/

#include "mgm_costvolume.h"
struct costvolume_t allocate_and_fill_sgm_costvolume (struct Img &in_u, // source (reference) image                      
                                                      struct Img &in_v, // destination (match) image                  
                                                      struct Img &dminI,// per pixel max&min disparity
                                                      struct Img &dmaxI,
                                                      char* prefilter,        // none, sobel, census(WxW)
                                                      char* distance,         // census, l1, l2, ncc(WxW), btl1, btl2
                                                      float truncDist);        // truncated differences

/********************** MGM *****************************/

#include "mgm_core.cc"
struct costvolume_t mgm(struct costvolume_t CC, const struct Img &in_w, 
                        const struct Img &dminI, const struct Img &dmaxI, 
                        struct Img *out, struct Img *outcost, 
                        const float P1, const float P2, const int NDIR, const int MGM, 
                        const int USE_FELZENSZWALB_POTENTIALS, // USE SGM(0) or FELZENSZWALB(1) POTENTIALS
                        int SGM_FIX_OVERCOUNT);                // fix the overcounting in SGM following (Drory etal. 2014)

#include "mgm_weights.h"
struct Img compute_mgm_weights(struct Img &u, float aP, float aThresh);

/********** SOLUTION REFINEMENT AND ENERGY COMPUTATION ***********/

#include "mgm_refine.h"
void subpixel_refinement_sgm(struct costvolume_t &S,       // modifies out and outcost
                             std::vector<float > &out,  
                             std::vector<float > &outcost, 
                             char *refinement); //none, vfit, parabola, cubic, parabolaOCV

#include "mgm_print_energy.h"
void print_solution_energy(const struct Img &in_u, std::vector<float > &disp, 
                           struct costvolume_t &CC, float P1, float P2);

/********************** OTHERSTUFF  *****************************/


void leftright_test(struct Img &dx, struct Img &Rdx, float threshold=1)
{
   int nc = dx.ncol;
   int nr = dx.nrow;
   int Rnc = Rdx.ncol;
   int Rnr = Rdx.nrow;

   for(int y=0;y<nr;y++)
      for(int x=0;x<nc;x++) {
         int i=x+y*nc;
         int Lx,Rx;
         Lx = round(x+ dx[i]);
         if( (Lx)<Rnc && (Lx)>=0 ){
            int Lidx = Lx + y*Rnc;
            float Rx = Lx + Rdx[Lidx];
            if ( fabs(Rx-x) > threshold) {
               dx[i]  = NAN;
            }
         }else {
            dx[i]  = NAN;
         }
      }

}

void leftright_test_bleyer(struct Img &dx, struct Img &Rdx)
   // warps the pixels of the right image to the left, if no pixel in the
   // left image receives a contribution then it is marked as occluded
{
   int nc = dx.ncol;
   int nr = dx.nrow;
   int Rnc = Rdx.ncol;
   int Rnr = Rdx.nrow;

   struct Img occL(nc,nr);
   for(int i=0;i<nr*nc;i++) occL[i]=0;

   for(int y=0;y<Rnr;y++)
      for(int x=0;x<Rnc;x++) {
         int i=x+y*Rnc;
         int Lx = round(x+ Rdx[i]);
         if( (Lx)<nc && (Lx)>=0 ){
            occL[Lx + y*nc] = 255;
         }
      }

   for(int i=0;i<nr*nc;i++) 
      if(occL[i]==0) dx[i] = NAN;

}


std::pair<float, float> update_dmin_dmax(struct Img outoff, struct Img *dminI, struct Img *dmaxI, int slack=3, int radius=2) {
   struct Img dminI2(*dminI);
   struct Img dmaxI2(*dmaxI);
   int nx = outoff.nx;
   int ny = outoff.ny;

   // global (finite) min and max
   std::pair<float,float>gminmax = image_minmax(outoff);
   float gmin = gminmax.first; float gmax = gminmax.second;

   if (slack<0) slack = -slack;
   int r=radius;

   for (int j=0;j<ny;j++)
   for (int i=0;i<nx;i++)
   {
      float dmin = INFINITY; float dmax = -INFINITY;
      for (int dj=-r;dj<=r;dj++)
      for (int di=-r;di<=r;di++)
      {
         float v = valneumann(outoff, i+di, j+dj);
         if (isfinite(v)) {
            dmin = fmin( dmin, v - slack );
            dmax = fmax( dmax, v + slack );
         } else {
            dmin = fmin( dmin, gmin - slack );
            dmax = fmax( dmax, gmax + slack );
         }
      }
      if (isfinite(dmin)) { 
         dminI2[i+j*nx] = dmin; dmaxI2[i+j*nx] = dmax; 
      }

   }

   *dminI = dminI2;
   *dmaxI = dmaxI2;
   return std::pair<float, float> (gmin, gmax);
}


// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
static char *pick_option(int *c, char ***v, char *o, char *d)
{
   int argc = *c;
   char **argv = *v;
   int id = d ? 1 : 0;
   for (int i = 0; i < argc - id; i++)
      if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, o)) {
	 char *r = argv[i + id] + 1 - id;
	 *c -= id + 1;
	 for (int j = i; j < argc - id; j++)
	    (*v)[j] = (*v)[j + id + 1];
	 return r;
      }
   return d;
}



/*MGM*/


SMART_PARAMETER(TSGM,4);
SMART_PARAMETER(TSGM_FIX_OVERCOUNT,1);
SMART_PARAMETER(TSGM_2LMIN,0);
SMART_PARAMETER(USE_TRUNCATED_LINEAR_POTENTIALS,0);

SMART_PARAMETER(WITH_MGM2,0);

SMART_PARAMETER(TSGM_ITER,1)
SMART_PARAMETER(TESTLRRL,1)
SMART_PARAMETER(MEDIAN,0)



int main(int argc, char* argv[]) 
{
	/* patameter parsing - parameters*/
	if(argc<4)
	{
		fprintf (stderr, "too few parameters\n");
		fprintf (stderr, "   usage: %s  [-r dmin -R dmax] [-m dminImg -M dmaxImg] [-O NDIR: 2, (4), 8, 16] u v out [cost [backflow]]\n",argv[0]);
		fprintf (stderr, "        [-P1 (8) -P2 (32)]: sgm regularization parameters P1 and P2\n");
		fprintf (stderr, "        [-p PREFILT(none)]: prefilter = {none|census|sobelx|gblur} (census is WxW)\n");
		fprintf (stderr, "        [-t      DIST(ad)]: distance = {census|ad|sd|ncc|btad|btsd}  (ncc is WxW, bt is Birchfield&Tomasi)\n");
		fprintf (stderr, "        [-truncDist (inf)]: truncate distances at nch*truncDist  (default INFINITY)\n");
		fprintf (stderr, "        [-s  SUBPIX(none)]: subpixel refinement = {none|vfit|parabola|cubic}\n");
		fprintf (stderr, "        [-aP1         (1)]: multiplier factors of P1 and P2 when\n");
		fprintf (stderr, "        [-aP2         (1)]:    \\sum |I1 - I2|^2 < nch*aThresh^2\n");
		fprintf (stderr, "        [-aThresh     (5)]: Threshold for the multiplier factor (default 5)\n");
		fprintf (stderr, "        [-l   FILE (none)]: write disparity without LR test (default none)\n");
		fprintf (stderr, "        ENV: CENSUS_NCC_WIN=3   : size of the window for census and NCC\n");
		fprintf (stderr, "        ENV: TESTLRRL=1   : lrrl\n");
		fprintf (stderr, "        ENV: MEDIAN=0     : radius of the median filter postprocess\n");
		fprintf (stderr, "        ENV: TSGM=4       : regularity level\n");
		fprintf (stderr, "        ENV: TSGM_ITER=1  : iterations\n");
		fprintf (stderr, "        ENV: TSGM_FIX_OVERCOUNT=1   : fix overcounting of the data term in the energy\n");
		fprintf (stderr, "        ENV: TSGM_DEBUG=0 : prints debug informtion\n");
		fprintf (stderr, "        ENV: TSGM_2LMIN=0 : use the improved TSGM cost only for TSGM=2. Overrides TSGM value\n");
		fprintf (stderr, "        ENV: USE_TRUNCATED_LINEAR_POTENTIALS=0 : use the Felzenszwalb-Huttenlocher\n");
		fprintf (stderr, "                          : truncated linear potential (when=1). P1 and P2 change meaning\n");
		fprintf (stderr, "                          : The potential they describe becomes:  V(p,q) = min(P2,  P1*|p-q|)\n");
		return 1;
	}
	
	
	//read the parameters
	int i = 1;
   char *in_min_disp_file = pick_option(&argc, &argv, (char*) "m", (char*) "");
   char *in_max_disp_file = pick_option(&argc, &argv, (char*) "M", (char*) "");
   int dmin = atoi(pick_option(&argc, &argv, (char*) "r", (char*) "-30"));
   int dmax = atoi(pick_option(&argc, &argv, (char*) "R", (char*) "30"));
   int NDIR  = atoi(pick_option(&argc, &argv, (char*) "O", (char*) "4"));
   float P1  = atof(pick_option(&argc, &argv, (char*) "P1", (char*) "8"));
   float P2  = atof(pick_option(&argc, &argv, (char*) "P2", (char*) "32"));
   float aP1 = atof(pick_option(&argc, &argv, (char*) "aP1", (char*) "1"));
   float aP2 = atof(pick_option(&argc, &argv, (char*) "aP2", (char*) "1"));
   float aThresh = atof(pick_option(&argc, &argv, (char*) "aThresh", (char*) "5"));

   char* distance  = pick_option(&argc, &argv, (char*) "t", (char*) "ad");   //{census|ad|sd|ncc|btad|btsd}
   char* prefilter = pick_option(&argc, &argv, (char*) "p", (char*) "none"); //{none|census|sobelx}
   char* refine    = pick_option(&argc, &argv, (char*) "s", (char*) "none"); //{none|vfit|parabola|cubic}
   float truncDist = atof(pick_option(&argc, &argv, (char*) "truncDist",  (char*) "inf"));
   char *nolr_disp_file = pick_option(&argc, &argv, (char*) "l", (char*) "");

	char* f_u     = (argc>i) ? argv[i] : NULL;      i++;
	char* f_v     = (argc>i) ? argv[i] : NULL;      i++;
	char* f_out   = (argc>i) ? argv[i] : NULL;      i++;
	char* f_cost  = (argc>i) ? argv[i] : NULL;      i++;
	char* f_back  = (argc>i) ? argv[i] : NULL;      i++;
	
	printf("%d %d\n", dmin, dmax);
	
	
	// read input
	struct Img u = iio_read_vector_split(f_u);
	struct Img v = iio_read_vector_split(f_v);

   remove_nonfinite_values_Img(u, 0);
   remove_nonfinite_values_Img(v, 0);

   struct Img dminI(u.nx, u.ny);
   struct Img dmaxI(u.nx, u.ny);
   for(int i=0;i<u.npix;i++) {dminI[i]=dmin; dmaxI[i]=dmax;}

   if(strcmp (in_min_disp_file,"")!=0 ){
   	dminI = iio_read_vector_split(in_min_disp_file);
	   dmaxI = iio_read_vector_split(in_max_disp_file);
      // sanity check for nans
      remove_nonfinite_values_Img(dminI, dmin);
      remove_nonfinite_values_Img(dmaxI, dmax);

      // more hacks to prevent produce due to bad inputs (min>=max)
      for (int i=0;i<u.npix;i++) {
         if (dmaxI[i] < dminI[i] + 1) dmaxI[i] = ceil(dminI[i] + 1);
      }
   }
	

	P1 = P1*u.nch; //8
	P2 = P2*u.nch; //32
	
	// call
	struct Img outoff  = Img(u.nx, u.ny);
	struct Img outcost = Img(u.nx, u.ny);

   // variables for LR
	struct Img outoffR  = Img(v.nx, v.ny);
	struct Img outcostR = Img(v.nx, v.ny);
   struct Img dminRI(v.nx, v.ny);
   struct Img dmaxRI(v.nx, v.ny);
   for(int i = 0; i < v.npix; i++) {dminRI[i] = -dmax; dmaxRI[i] = -dmin;}



   struct Img u_w = compute_mgm_weights(u, aP2, aThresh); // missing aP1 !! TODO
   struct Img v_w = compute_mgm_weights(v, aP2, aThresh);


   struct costvolume_t CC = allocate_and_fill_sgm_costvolume (u, v, dminI, dmaxI, prefilter, distance, truncDist);
   for(int i = 0; i < TSGM_ITER(); i++) {
      struct costvolume_t S = WITH_MGM2() ?  
                                 mgm2(CC, u_w, dminI, dmaxI, &outoff, &outcost, P1, P2, 
                                    NDIR, TSGM(), USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) :
                                 mgm(CC, u_w, dminI, dmaxI, &outoff, &outcost, P1, P2, 
                                    NDIR, TSGM(), USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) ;
      print_solution_energy(u, outoff.data, CC, P1, P2);
      // call subpixel refinement  (modifies out and outcost)
      subpixel_refinement_sgm(S, outoff.data, outcost.data, refine);
      std::pair<float,float>gminmax = update_dmin_dmax(outoff, &dminI, &dmaxI);
      remove_nonfinite_values_Img(dminI, gminmax.first);
      remove_nonfinite_values_Img(dmaxI, gminmax.second);
//       char name[200]; sprintf(name, "/tmp/%02d.tif", i); // DEBUG
//	      iio_write_vector_split(name, outoff); // DEBUG
//         // dump disp range
//       struct Img rr = Img(dmaxI);
//       for(int i=0;i<rr.npix;i++) rr[i] -= dminI[i];
//	      iio_write_vector_split(name, rr); // DEBUG
   }
   if(MEDIAN()) outoff = median_filter(outoff,MEDIAN());


	// save the disparity without LR
	if( ! strcmp (nolr_disp_file, "") == 0 ) 
      iio_write_vector_split(nolr_disp_file, outoff);

   
   if(TESTLRRL()) {
      struct costvolume_t CC = allocate_and_fill_sgm_costvolume (v, u, dminRI, dmaxRI, prefilter, distance, truncDist);
      for(int i = 0; i < TSGM_ITER(); i++) {
         struct costvolume_t S = WITH_MGM2() ?  
                                    mgm2(CC, v_w, dminRI, dmaxRI, &outoffR, &outcostR, P1, P2, 
                                       NDIR, TSGM(), USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) :
                                    mgm(CC, v_w, dminRI, dmaxRI, &outoffR, &outcostR, P1, P2, 
                                       NDIR, TSGM(), USE_TRUNCATED_LINEAR_POTENTIALS(), TSGM_FIX_OVERCOUNT()) ;
         print_solution_energy(v, outoffR.data, CC, P1, P2);
         // call subpixel refinement  (modifies out and outcost)
         subpixel_refinement_sgm(S, outoffR.data, outcostR.data, refine);
         std::pair<float,float>gminmax = update_dmin_dmax(outoffR, &dminRI, &dmaxRI);
         remove_nonfinite_values_Img(dminRI, gminmax.first);
         remove_nonfinite_values_Img(dmaxRI, gminmax.second);
      }
      if(MEDIAN()) outoffR = median_filter(outoffR,MEDIAN());
      Img tmpL(outoff);
      Img tmpR(outoffR);
      leftright_test(outoffR, tmpL); // R-L
      leftright_test(outoff, tmpR);  // L-R
   }


	
	// save the disparity
	struct Img out = Img(u.nx, u.ny);
	for(int i=0;i<u.nx*u.ny;i++) out.data[i]=outoff[i];
	
	// generate the backprojected image
	struct Img syn = Img(u.nx, u.ny, u.nch);
	for(int x=0;x<u.nx;x++)
		for(int y=0;y<u.ny;y++){
			Point p(x,y);
			Point q = Point(outoff[x+u.nx*y],0);
			for(int c=0;c<u.nch;c++)
				if( check_inside_image(p+q, v) ) 
					syn.data[x+y*u.nx + c*u.npix] = v.data[x+q.x+(y+q.y)*v.nx + c*v.npix];
				else 
					syn.data[x+y*u.nx+ c*u.npix] = u.data[x+y*u.nx+ c*u.npix];
		}
	
	iio_write_vector_split(f_out, out);
	if(f_cost) iio_write_vector_split(f_cost, outcost);
	if(f_back) iio_write_vector_split(f_back, syn);
	
	return 0;
}
