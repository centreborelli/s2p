/* You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
/* Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include <numeric>
#include <algorithm>
#include <vector>
#include <cstring>
#include "assert.h"

//// structures to wrap images and points
#include "img.h"
#include "img_tools.h"
#include "point.h"

// from img_tools.h
//inline int check_inside_image(const Point p, const struct Img &u) {
//   int nx = u.nx;
//   int ny = u.ny;
//   float x = p.x;
//   float y = p.y;
//   if(x>=0 && y>=0 && x<nx && y<ny) return 1;
//   else return 0;
//}


/********************** COSTVOLUME *****************************/
#include "mgm_costvolume.h"

//struct costvolume_t allocate_costvolume (struct Img min, struct Img max) {
//   struct costvolume_t cv;
//   cv.vectors = std::vector< Dvec >(min.npix);
//   for (int i=0;i< min.npix;i++) 
//      cv[i].init(min[i], max[i]);
//   return cv;
//}
struct costvolume_t allocate_costvolume (struct Img min, struct Img max);


/********************** MGM *****************************/

#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))
// fast alternatives to: __min(a,__min(b,c)) 
// fastestest ?
#define fmin3_(x, y, z) \
   (((x) < (y)) ? (((z) < (x)) ? (z) : (x)) : (((z) < (y)) ? (z) : (y)))
// fast and easy to understand
static inline float fmin3(float a, float b, float c)
{
   float m = a;
   if (m > b) m = b;
   if (m > c) m = c;
   return m;
}


// intervening points p,q,r 
// faster variant for the case 2
// THIS IS THE SIMPLEST MGM WEIGHT UPDATE FUNCTION
inline void update_cost2(Dvec &Lp, Dvec &CCp, Dvec &Lq, Dvec &Lr, const float P1, const float P2) {

            float min1L_all = Lq.get_minvalue();
            float min2L_all = Lr.get_minvalue();

            for(int o=Lp.min;o<=Lp.max;o++) {
               float C    = CCp[o];  // the matching cost for p <-> p+d

               float vL0  = Lq[o];	// the neighbor has the same label
               float vLP1 = __min( Lq[o-1],  Lq[o+1]) + P1;		// the neighbour has a similar (+-1) label
               float vLP2 = min1L_all + P2;  // the minimum label of the neighbour

               float v2L0  = Lr[o];	// the neighbor has the same label
               float v2LP1 = __min( Lr[o-1],  Lr[o+1]) + P1;		// the neighbour has a similar (+-1) label
               float v2LP2 = min2L_all + P2;  // the minimum label of the neighbour

               float edge_potentials = 0;
               edge_potentials += (fmin3(vL0 , vLP1 , vLP2 ) - min1L_all) / 2;
               edge_potentials += (fmin3(v2L0, v2LP1, v2LP2) - min2L_all) / 2;


               Lp.set_nolock(o,  C + edge_potentials);
            }

}


// intervening points p,q,r 
// THIS FUNCTION CONSIDERS 4 NEIGHBORS AND WEIGHTED EDGES
inline void update_costW(Dvec &Lp, Dvec &CCp, Dvec &Lq, Dvec &Lr, Dvec &Ls, Dvec &Lt, const float P1, const float P2,
      const float DeltaI1, const float DeltaI2, const float DeltaI3, const float DeltaI4, const int howmany) {

            float minL_all  = INFINITY;
            float min2L_all = INFINITY;
            float min3L_all = INFINITY;
            float min4L_all = INFINITY;

            minL_all = Lq.get_minvalue();
            if (howmany >= 2) min2L_all = Lr.get_minvalue();
            if (howmany >= 3) min3L_all = Ls.get_minvalue();
            if (howmany >= 4) min4L_all = Lt.get_minvalue();

            for(int o=Lp.min;o<=Lp.max;o++) {
               float C    = CCp[o];  // the matching cost for p <-> p+d
               //float C    = computeC( p,p+Point(o,0), u,v);  // the matching cost for p <-> p+d // SLOWER
               float edge_potentials = 0;
               
               float vL0  = Lq[o];	// the neighbor has the same label
               float vLP1 = __min( Lq[o-1],  Lq[o+1]) + P1*DeltaI1;		// the neighbour has a similar (+-1) label
               float vLP2 = minL_all + P2*DeltaI1;  // the minimum label of the neighbour
               edge_potentials += fmin3(vL0 , vLP1 , vLP2 ) - minL_all;

               if (howmany >= 2) {
                  float v2L0  = Lr[o];	// the neighbor has the same label
                  float v2LP1 = __min( Lr[o-1],  Lr[o+1]) + P1*DeltaI2;		// the neighbour has a similar (+-1) label
                  float v2LP2 = min2L_all + P2*DeltaI2;  // the minimum label of the neighbour

                  edge_potentials += fmin3(v2L0, v2LP1, v2LP2) - min2L_all;
               }

               if (howmany >= 3) {
                  float v3L0  = Ls[o];	// the neighbor has the same label
                  float v3LP1 = __min( Ls[o-1],  Ls[o+1]) + P1*DeltaI3;		// the neighbour has a similar (+-1) label
                  float v3LP2 = min3L_all + P2*DeltaI3;  // the minimum label of the neighbour

                  edge_potentials += fmin3(v3L0, v3LP1, v3LP2) - min3L_all;
               }

               if (howmany >= 4) {
                  float v4L0  = Lt[o];	// the neighbor has the same label
                  float v4LP1 = __min( Lt[o-1],  Lt[o+1]) + P1*DeltaI4; // the neighbour has a similar (+-1) label
                  float v4LP2 = min4L_all + P2*DeltaI4;  // the minimum label of the neighbour
                  edge_potentials += fmin3(v4L0, v4LP1, v4LP2) - min4L_all;
               }

               Lp.set_nolock(o,  C + edge_potentials / howmany);
            }

}



// compute in place the min convolution of vector M[] (of lenght mm)
// with the distance function with slope P1 and truncated at P2 (may be INFINITY)
// minMall is the minimum of the input M[] which is needed for the truncated distance
// minMall may be a value lower than the values stored in M[]
static void minConvTruncatedLinear(float M[], const int mm, const float minMall, const float P1, const float P2) {
   // forward pass
   for(int o=1; o<mm; o++) 
      M[o] = __min(M[o-1] + P1, M[o]);
   // backward pass
   for(int o=mm-2; o>=0; o--)  
      M[o] = __min(M[o+1] + P1, M[o]);
   // truncated distance
   if (P2 < INFINITY)
      for(int o=0; o<mm; o++) 
         M[o] = __min(M[o], minMall + P2);
}


static void FixBounrady_for_minConvTruncatedLinear(const float I[], const int imin, const int imax, float M[], const int mmin, const int mmax, const float P1) {
   // handle boundary cases (left)
   if (imin < mmin) {
      float T = I[0];
      for(int o=imin+1;o<=mmin;o++) {
         float Inext = o<=imax ? I[o-imin]: INFINITY;
         T = __min(T + P1, Inext);
      }
      M[0] = __min(M[0], T);
   }

   // handle boundary cases (right)
   if (imax > mmax) {
      float T = I[imax-imin];
      for(int o=imax-1;o>=mmax;o--) {
         float Inext = o>=imin ? I[o-imin]: INFINITY;
         T = __min(T + P1, Inext);
      }
      M[mmax-mmin] = __min(M[mmax-mmin], T);
   }
}




// intervening points p,q,r 
// faster variant for the case 2
// Adaptation of the Felzenszwalb-Huttenlocher message passing for the truncated linear model
// see: "Efficient Belief Propagation for Early Vision"
// P1 and P2 ARE USED WITH A DIFFERENT MEANING
// HERE THE COST IS:   V(p,q) = min(P2,  P1*|p-q|)
inline void update_cost2_trunclinear(Dvec &Lp, Dvec &CCp, Dvec &Lq, Dvec &Lr, const float P1, const float P2) {

            float min1L_all = Lq.get_minvalue();
            float min2L_all = Lr.get_minvalue();

            float M1[Lp.max-Lp.min+1];
            float M2[Lp.max-Lp.min+1];
            int mm = Lp.min;

            // initialize copying the values
            for(int o=Lp.min;o<=Lp.max;o++)  M1[o-mm] = Lq[o];
            FixBounrady_for_minConvTruncatedLinear(&(Lq.data[0]), Lq.min, Lq.max, M1, Lp.min, Lp.max, P1);
            minConvTruncatedLinear(M1, Lp.max-Lp.min+1, min1L_all, P1, P2);

            for(int o=Lp.min;o<=Lp.max;o++)  M2[o-mm] = Lr[o];
            FixBounrady_for_minConvTruncatedLinear(&(Lr.data[0]), Lr.min, Lr.max, M2, Lp.min, Lp.max, P1);
            minConvTruncatedLinear(M2, Lp.max-Lp.min+1, min2L_all, P1, P2);

            for(int o=Lp.min;o<=Lp.max;o++) {
               Lp.set_nolock(o, CCp[o] + (M1[o-mm] -min1L_all + M2[o-mm] - min2L_all)/2);
            }

}



// intervening points p,q,r,s,t
// Adaptation of the Felzenszwalb-Huttenlocher message passing for the truncated linear model
// see: "Efficient Belief Propagation for Early Vision"
// P1 and P2 ARE USED WITH A DIFFERENT MEANING
// HERE THE COST IS:   V(p,q) = min(P2,  P1*|p-q|)
// THIS FUNCTION CONSIDERS 4 NEIGHBORS AND WEIGHTED EDGES
inline void update_costW_trunclinear(Dvec &Lp, Dvec &CCp, Dvec &Lq, Dvec &Lr, Dvec &Ls, Dvec &Lt, const float P1, const float P2, 
      const float DeltaI1, const float DeltaI2, const float DeltaI3, const float DeltaI4, const int howmany) {

            float min1L_all = INFINITY;
            float min2L_all = INFINITY;
            float min3L_all = INFINITY;
            float min4L_all = INFINITY;

            min1L_all = Lq.get_minvalue();
            if (howmany >= 2) min2L_all = Lr.get_minvalue();
            if (howmany >= 3) min3L_all = Ls.get_minvalue();
            if (howmany >= 4) min4L_all = Lt.get_minvalue();

            int mm = Lp.min;
            int NN = Lp.max-Lp.min+1;
            float M1[NN];
            float M2[NN];
            float M3[NN];
            float M4[NN];

            // initialize copying the values
            for(int o=Lp.min;o<=Lp.max;o++) M1[o-mm] = Lq[o];
            minConvTruncatedLinear(M1, NN, min1L_all, P1*DeltaI1, P2*DeltaI1);


            if (howmany >= 2) {
               // initialize copying the values
               for(int o=Lp.min;o<=Lp.max;o++) M2[o-mm] = Lr[o];
               minConvTruncatedLinear(M2, NN, min2L_all, P1*DeltaI2, P2*DeltaI2);
            }

            if (howmany >= 3) {
               // initialize copying the values
               for(int o=Lp.min;o<=Lp.max;o++) M3[o-mm] = Ls[o];
               minConvTruncatedLinear(M3, NN, min3L_all, P1*DeltaI3, P2*DeltaI3);
            }

            if (howmany >= 4) {
               // initialize copying the values
               for(int o=Lp.min;o<=Lp.max;o++) M4[o-mm] = Lt[o];
               minConvTruncatedLinear(M4, NN, min4L_all, P1*DeltaI4, P2*DeltaI4);
            }

            // compute the cost
            for(int o=Lp.min;o<=Lp.max;o++) {
               float edge_potentials              = M1[o-mm] - min1L_all;
               if (howmany >= 2) edge_potentials += M2[o-mm] - min2L_all;
               if (howmany >= 3) edge_potentials += M3[o-mm] - min3L_all;
               if (howmany >= 4) edge_potentials += M4[o-mm] - min4L_all;
               Lp.set_nolock(o, CCp[o] + edge_potentials/howmany);
            }

}


inline void update_cost2Lmin(Dvec &Lp, Dvec &CCp, Dvec &Lq, Dvec &Lr, float P1, float P2) {
   // The profiles of the 1D cost functions are of this form
   //
   //                    P2 --------------------------
   //                   /
   //      \           /
   //       P1       P1
   //         \     /
   //          \   /
   //   -2  -1 --0---1---2---3----
   //
   // When combined in 2D these functions lead to the following cases (with P1<P2)
   // produces an odd shape shape in 2D
   //                             |           
   //                             |          
   //                            (2)P2     P2+P1    ...      P2+P2   
   //                             |           |
   //                             |           |               .
   //                             |           |               . 
   //                             |           |
   //                P1+P1 ---   (1)P1 ---  P1+P1   ...      P1+P2
   //                  |          |           |               |
   //                  |          |           |               |
   //               (-1)P1 ===   (0)   ===  (1)P1 ======... (2)P2 .......
   //                  |          |           |       
   //                  |          |           |       
   //                P1+P1 ---    P1   ---  P1+P1  -----
   //                             |
   //                             |
   //  
   //  We would like to try something more isotropic like
   //                             |           
   //                             |          
   //                            (2)P2     P2+P1    ...      P2+P2   
   //                             |           |
   //                             |           |               .
   //                             |           |               . 
   //                             |           |
   //                P1+P1 ---  (1).3P1 ---  P1+P1   ...      P1+P2
   //                  |          |           |               |
   //                  |          |           |               |
   //               (-1).3P1 === (0)   ===  (1).3P1======... (2)P2 .......
   //                  |          |           |       
   //                  |          |           |       
   //                P1+P1 ---   .3P1   ---  P1+P1  -----
   //                             |
   //                             |
   //  
   //  But the net effect is rather inperceptible

            float min1L_all = Lq.get_minvalue();
            float min2L_all = Lr.get_minvalue();

            // MUST TEST 9 CONFIGURATIONS!
            for(int o=Lp.min;o<=Lp.max;o++) {
               float C    = CCp[o];  // the matching cost for p <-> p+d
               //float C    = computeC( p,p+Point(o,0), u,v);  // the matching cost for p <-> p+d // SLOWER
               
               float edge_potentials = 0;
               
               float vL0  = Lq[o];	// the neighbor has the same label
               float vLP1 = __min( Lq[o-1],  Lq[o+1]) + P1;		// the neighbour has a similar (+-1) label
               float vLP2 = min1L_all + P2;  // the minimum label of the neighbour

               float v2L0  = Lr[o];	// the neighbor has the same label
               float v2LP1 = __min( Lr[o-1],  Lr[o+1]) + P1;		// the neighbour has a similar (+-1) label
               float v2LP2 = min2L_all + P2;  // the minimum label of the neighbour

               edge_potentials = fmin3( 
                  fmin3( 
                  vL0  + v2LP1   -0.7*P1 ,    // +-1,0
                  vLP1 + v2L0    -0.7*P1  ,    // 0,+-1
                  vLP1 + v2LP1  //-2*P1+2*sqrt(2.0)*P1      // +-1,+-1
                  ),
                  fmin3( 
                  vL0  + v2L0    ,   // 0,0
                  vLP1 + v2LP2   ,//-P1 +P2 ,    // +-1,others
                  vLP2 + v2LP1   //-P1 +P2       // others,+-1
                  ),
                  fmin3( 
                  vLP2 + v2LP2      ,    // others,others
                  vL0  + v2LP2  , //+P2  ,   // others,0
                  vLP2 + v2L0   //+P2      // 0,others
                  )
               ) / 2 - (min1L_all + min2L_all)/2;

               Lp.set_nolock(o,  C + edge_potentials);

            }

}




struct Pass_setup {
   int row_major;
   Point dir1;
   Point dir2;
   Point dir3;
   Point dir4;
   int inc_x;  // 1: ascending x  0: descending x
   int inc_y;  // 1: ascending y  0: descending y
   // 0: use dir1 if y odd dir2 otherwise
   Pass_setup(Point d1, Point d2, Point d3, Point d4, int ix, int iy, int rm) {
      dir1=d1; dir2=d2; dir3=d3; dir4=d4;
      inc_x=ix;inc_y=iy;
      row_major=rm;
   }
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



// mgm returns the "aggregated" cost volume, out, and outcost without any other refinement
struct costvolume_t mgm(struct costvolume_t CC, const struct Img &in_w, 
                        const struct Img &dminI, const struct Img &dmaxI, 
                        struct Img *out, struct Img *outcost, 
                        const float P1, const float P2, const int NDIR, const int MGM, 
                        const int USE_FELZENSZWALB_POTENTIALS = 0, // USE SGM(0) or FELZENSZWALB(1) POTENTIALS
                        int SGM_FIX_OVERCOUNT = 1)                 // fix the overcounting in SGM following (Drory etal. 2014)
{

   int nx = dminI.nx;
   int ny = dminI.ny;

   // check the content of in_w is it all 1?
   int USE_IMAGE_DEPENDENT_WEIGHTS=0;
   for (int i=0; i<in_w.ncol*in_w.nrow*in_w.nch; i++) 
      if (in_w[i] != 1.0) USE_IMAGE_DEPENDENT_WEIGHTS = 1;
   if (USE_IMAGE_DEPENDENT_WEIGHTS) printf(" USING IMAGE DEPENDENT WEIGHTS\n");

   // run SGM              // ALLOCATED AND INITIALIZED TO 0 (THIS IS THE costvolume THAT IS RETURNED!)
   struct costvolume_t S = allocate_costvolume(dminI, dmaxI);

   std::vector<Pass_setup > direct;
   //                         PASSES
   //
   // O: first pixel in the scan
   // c: current pixel in the scan
   // - or |: processed pixels (according to the scan order)
   // 1,2,3,4: considered neighbours in the corresponding order 
   //
   //      (I)              (II)           (III)             (IV)
   //
   // O-----------                      O | | | |         | | | | | O 
   // ------------           c--1--     | | | | 4         1 3 | | | | 
   // ------------     ---4--2--3--     | | | | |         | | | | | | 
   // ------------     ------------     | | | | 2 o       o 2 | | | | 
   // -3--2--4----     ------------     | | | | | |         | | | | | 
   // -1--c            ------------     | | | | 3 1         4 | | | | 
   //                  -----------O     | | | | | |         | | | | | 
   // 
   //
   //                         NEIGHBORS
   //
   //            (-1,-1)       (0,-1)        (1,-1) 
   //                            |
   //                            |
   //            (-1,0)   ---    o    ---    (1,0)
   //                            | 
   //                            |
   //            (-1,1)        (0,1)         (1,1) 
   //
   // with MGM == 1  only the neighbor #1 is considered
   // with MGM == 2  #1 and #2 are considered
   // with MGM == 4  #1 to #4 are considered
   
   
   // horizontal and vertical
   direct.push_back( Pass_setup(Point(-1,0)  , Point(0,-1)  , Point(-1,-1) , Point(1,-1)  ,1,1,1) ); // (I)
   direct.push_back( Pass_setup(Point(1,0)   , Point(0,1)   , Point(1,1)   , Point(-1,1)  ,0,0,1) ); // (II)
   direct.push_back( Pass_setup(Point(0,1)   , Point(-1,0)  , Point(-1,1)  , Point(-1,-1) ,1,0,0) ); // (III)
   direct.push_back( Pass_setup(Point(0,-1)  , Point(1,0)   , Point(1,-1)  , Point(1,1)   ,0,1,0) ); // (IV)
   // diagonals 45º                                          
   direct.push_back( Pass_setup(Point(-1,-1) , Point(1,-1)  , Point(0,-1)  , Point(1,0)   ,0,1,1) );
   direct.push_back( Pass_setup(Point(1,-1)  , Point(1,1)   , Point(1,0)   , Point(0,1)   ,0,0,0) );
   direct.push_back( Pass_setup(Point(1,1)   , Point(-1,1)  , Point(0,1)   , Point(-1,0)  ,1,0,1) );
   direct.push_back( Pass_setup(Point(-1,1)  , Point(-1,-1) , Point(-1,0)  , Point(0,-1)  ,1,1,0) );
   // 22.5º
   // ...  

   // translate pass directions to edges encoded in channels of the image w
   // THIS IS TIED TO THE INFORMATION IN THE VECTOR direct
   // THE SAME VECTORS COULD BE OBTAINED PROGRAMATICALLY FROM direct 
   //   str::vector < std::pair<Point, int> > dir_to_idx;
   //   for(int i=0;i<7;i++)
   //      dir_to_idx.push_back( std::pair<Point, int> (direct[i].dir1 ,i ));
   int pass_to_channel_1[] = {0,1,2,3,4,5,6,7};
   int pass_to_channel_2[] = {3,2,0,1,5,6,7,4};
   int pass_to_channel_3[] = {4,6,7,5,3,1,2,0};
   int pass_to_channel_4[] = {5,7,4,6,1,2,0,3};

   // local Lr (could be implemented with a couple of line buffers)
   struct costvolume_t Lr(CC);

   for(int pass=0;pass<NDIR;pass++)
   {
      printf("%d", pass); fflush(stdout);
      Pass_setup dir = direct[pass];

      // reset the values of Lr for this passage
      #pragma omp parallel for
      for(int pidx=0; pidx<nx*ny; pidx++) 
         for(int o=Lr[pidx].min;o<=Lr[pidx].max;o++) 
            Lr[pidx].set_nolock(o, CC[pidx][o] ); // no omp critic is inside set_nolock

      int maxii = nx, maxjj = ny;
      if( !dir.row_major ) { maxii = ny; maxjj = nx; }
      

      // scan in the horizontal direction left to right
      for(int ii=0; ii<maxii+2*maxjj; ii++) {
      #pragma omp parallel for schedule(static,1)
      for(int jj=0; jj<maxjj; jj++) 
      {
         // ensure diagonal scan (slope 2)
         int x=ii -2*jj, y=jj;
         if(x < 0 || x >= maxii) continue;

         int maxnx = maxii, maxny = maxjj;
         // swap the indices if we are in column major 
         #define SWAPi(a,b) {int swap=a;a=b;b=swap;}
         if (!dir.row_major) {
            SWAPi(x, y);
            SWAPi(maxnx, maxny);
         }

         // reverse the direction
         if(dir.inc_x==0) x = (maxnx-1)-x;
         if(dir.inc_y==0) y = (maxny-1)-y;

         Point p(x,y);				   // current point
         Point pr  = p + dir.dir1;  // dir1 neighbor
         Point pr2 = p + dir.dir2;  // dir2 neighbor
         Point pr3 = p + dir.dir3;  // dir3 neighbor
         Point pr4 = p + dir.dir4;  // dir4 neighbor

         // base index of the neighbor
         int pidx   = (p.x +p.y *nx);
         int pridx  = (pr.x+pr.y*nx);
         int pr2idx = (pr2.x+pr2.y*nx);
         int pr3idx = (pr3.x+pr3.y*nx);
         int pr4idx = (pr4.x+pr4.y*nx);

         if (!check_inside_image(pr ,dminI)) continue;
         if (!check_inside_image(pr2,dminI)) continue;
         if (!check_inside_image(pr3,dminI)) continue;
         if (!check_inside_image(pr4,dminI)) continue;

         int TSGM_2LMIN = 0;
         if(TSGM_2LMIN>0) { // THIS IS A LEGACY FEATURE
            //   update_cost2L2(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
            update_cost2Lmin(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
         } 
         else if(USE_IMAGE_DEPENDENT_WEIGHTS) {      // IMAGE DEPENDENT WEIGHTS

            #define val(u, p, ch)  u.data[(p.x) + (u.nx)*(p.y) + (ch)*(u.npix)]
            float DeltaI1 = val(in_w, p, pass_to_channel_1[pass]);
            float DeltaI2 = val(in_w, p, pass_to_channel_2[pass]);
            float DeltaI3 = val(in_w, p, pass_to_channel_3[pass]);
            float DeltaI4 = val(in_w, p, pass_to_channel_4[pass]);
            #undef val
            if(USE_FELZENSZWALB_POTENTIALS>0) 
               update_costW_trunclinear(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2, 
                     DeltaI1, DeltaI2, DeltaI3, DeltaI4, MGM);
            else
               update_costW(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2,
                  DeltaI1, DeltaI2, DeltaI3, DeltaI4, MGM);
         }
         else {                                 // WITHOUT IMAGE DEPENDENT WEIGHTS
            if(USE_FELZENSZWALB_POTENTIALS>0) {
               if(MGM==2)
                  update_cost2_trunclinear(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
               else 
                  update_costW_trunclinear(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2, 
                        1.0, 1.0, 1.0, 1.0, MGM);
            }
            else if(MGM==2) 
               update_cost2(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
            else
               update_costW(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2, 
                     1.0, 1.0, 1.0, 1.0, MGM);
         }

      }
   }

      // accumulate S for the current orientation
      #pragma omp parallel for
      for(int i=0; i<nx*ny; i++) {
         for(int o=Lr[i].min;o<=Lr[i].max;o++) {
            S[i].increment_nolock(o, Lr[i][o]); // pragma omp critic is inside set
         }
      }
   }


   // WTA 
   #pragma omp parallel for
   for(int i=0;i<nx*ny;i++) {
      float minP;
      float minL=INFINITY;
      for(int o=S[i].min;o<=S[i].max;o++) {
         // overcounting correction (Drory etal. 2014) 
         if (SGM_FIX_OVERCOUNT==1)
            S[i].set_nolock(o, S[i][o] - (NDIR -1) * CC[i][o]);

         if(isfinite(S[i][o]))
         if(minL > S[i][o]) {
            minL = S[i][o];
            minP = o;
         }
      }	  
      (*out)[i] = minP;
      (*outcost)[i] = minL;
   }

   // return the aggregated costvolume
   return S;
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////






// mgm returns the "aggregated" cost volume, out, and outcost without any other refinement
struct costvolume_t mgm2(struct costvolume_t CC, const struct Img &in_w, 
                        const struct Img &dminI, const struct Img &dmaxI, 
                        struct Img *out, struct Img *outcost, 
                        const float P1, const float P2, const int NDIR, const int MGM, 
                        const int USE_FELZENSZWALB_POTENTIALS = 0, // USE SGM(0) or FELZENSZWALB(1) POTENTIALS
                        int SGM_FIX_OVERCOUNT = 1)                 // fix the overcounting in SGM following (Drory etal. 2014)
{

   int nx = dminI.nx;
   int ny = dminI.ny;

   // check the content of in_w is it all 1?
   int USE_IMAGE_DEPENDENT_WEIGHTS=0;
   for (int i=0; i<in_w.ncol*in_w.nrow*in_w.nch; i++) 
      if (in_w[i] != 1.0) USE_IMAGE_DEPENDENT_WEIGHTS = 1;
   if (USE_IMAGE_DEPENDENT_WEIGHTS) printf(" USING IMAGE DEPENDENT WEIGHTS\n");

   // run SGM              // ALLOCATED AND INITIALIZED TO 0 (THIS IS THE costvolume THAT IS RETURNED!)
   struct costvolume_t S = allocate_costvolume(dminI, dmaxI);

   std::vector<Pass_setup > direct;
   //                         PASSES
   //
   // O: first pixel in the scan
   // c: current pixel in the scan
   // - or |: processed pixels (according to the scan order)
   // 1,2,3,4: considered neighbours in the corresponding order 
   //
   //      (I)              (II)           (III)             (IV)
   //
   // O-----------                      O | | | |         | | | | | O 
   // ------------           c--1--     | | | | 4         1 3 | | | | 
   // ------------     ---4--2--3--     | | | | |         | | | | | | 
   // ------------     ------------     | | | | 2 o       o 2 | | | | 
   // -3--2--4----     ------------     | | | | | |         | | | | | 
   // -1--c            ------------     | | | | 3 1         4 | | | | 
   //                  -----------O     | | | | | |         | | | | | 
   // 
   //
   //                         NEIGHBORS
   //
   //            (-1,-1)       (0,-1)        (1,-1) 
   //                            |
   //                            |
   //            (-1,0)   ---    o    ---    (1,0)
   //                            | 
   //                            |
   //            (-1,1)        (0,1)         (1,1) 
   //
   // with MGM == 1  only the neighbor #1 is considered
   // with MGM == 2  #1 and #2 are considered
   // with MGM == 4  #1 to #4 are considered
   
   
   // horizontal and vertical
   direct.push_back( Pass_setup(Point(-1,0)  , Point(0,-1)  , Point(-1,-1) , Point(1,-1)  ,1,1,1) ); // (I)
   direct.push_back( Pass_setup(Point(1,0)   , Point(0,1)   , Point(1,1)   , Point(-1,1)  ,0,0,1) ); // (II)
   direct.push_back( Pass_setup(Point(0,1)   , Point(-1,0)  , Point(-1,1)  , Point(-1,-1) ,1,0,0) ); // (III)
   direct.push_back( Pass_setup(Point(0,-1)  , Point(1,0)   , Point(1,-1)  , Point(1,1)   ,0,1,0) ); // (IV)
   // diagonals 45º                                          
   direct.push_back( Pass_setup(Point(-1,-1) , Point(1,-1)  , Point(0,-1)  , Point(1,0)   ,0,1,1) );
   direct.push_back( Pass_setup(Point(1,-1)  , Point(1,1)   , Point(1,0)   , Point(0,1)   ,0,0,0) );
   direct.push_back( Pass_setup(Point(1,1)   , Point(-1,1)  , Point(0,1)   , Point(-1,0)  ,1,0,1) );
   direct.push_back( Pass_setup(Point(-1,1)  , Point(-1,-1) , Point(-1,0)  , Point(0,-1)  ,1,1,0) );
   // 22.5º
   // ...  

   // translate pass directions to edges encoded in channels of the image w
   // THIS IS TIED TO THE INFORMATION IN THE VECTOR direct
   // THE SAME VECTORS COULD BE OBTAINED PROGRAMATICALLY FROM direct 
   //   str::vector < std::pair<Point, int> > dir_to_idx;
   //   for(int i=0;i<7;i++)
   //      dir_to_idx.push_back( std::pair<Point, int> (direct[i].dir1 ,i ));
   int pass_to_channel_1[] = {0,1,2,3,4,5,6,7};
   int pass_to_channel_2[] = {3,2,0,1,5,6,7,4};
   int pass_to_channel_3[] = {4,6,7,5,3,1,2,0};
   int pass_to_channel_4[] = {5,7,4,6,1,2,0,3};

   #pragma omp parallel for
   for(int pass=0;pass<NDIR;pass++)
   {
      printf("%d", pass); fflush(stdout);
      Pass_setup dir = direct[pass];

      // local Lr (could be implemented with a couple of line buffers)
      struct costvolume_t Lr(CC);

      int maxii = nx, maxjj = ny;
      if( !dir.row_major ) { maxii = ny; maxjj = nx; }
      

      // scan in the horizontal direction left to right
      for(int jj=0; jj<maxjj; jj++) {
      for(int ii=0; ii<maxii; ii++)
      {
         int x=ii, y=jj;

         int maxnx = maxii, maxny = maxjj;
         // swap the indices if we are in column major 
         #define SWAPi(a,b) {int swap=a;a=b;b=swap;}
         if (!dir.row_major) {
            SWAPi(x, y);
            SWAPi(maxnx, maxny);
         }

         // reverse the direction
         if(dir.inc_x==0) x = (maxnx-1)-x;
         if(dir.inc_y==0) y = (maxny-1)-y;

         Point p(x,y);				   // current point
         Point pr  = p + dir.dir1;  // dir1 neighbor
         Point pr2 = p + dir.dir2;  // dir2 neighbor
         Point pr3 = p + dir.dir3;  // dir3 neighbor
         Point pr4 = p + dir.dir4;  // dir4 neighbor

         // base index of the neighbor
         int pidx   = (p.x +p.y *nx);
         int pridx  = (pr.x+pr.y*nx);
         int pr2idx = (pr2.x+pr2.y*nx);
         int pr3idx = (pr3.x+pr3.y*nx);
         int pr4idx = (pr4.x+pr4.y*nx);

         if (!check_inside_image(pr ,dminI)) continue;
         if (!check_inside_image(pr2,dminI)) continue;
         if (!check_inside_image(pr3,dminI)) continue;
         if (!check_inside_image(pr4,dminI)) continue;

         int TSGM_2LMIN = 0;
         if(TSGM_2LMIN>0) { // THIS IS A LEGACY FEATURE (DISABLED)
            //   update_cost2L2(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
            update_cost2Lmin(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
         } 
         else if(USE_IMAGE_DEPENDENT_WEIGHTS) {      // IMAGE DEPENDENT WEIGHTS

            #define val(u, p, ch)  u.data[(p.x) + (u.nx)*(p.y) + (ch)*(u.npix)]
            float DeltaI1 = val(in_w, p, pass_to_channel_1[pass]);
            float DeltaI2 = val(in_w, p, pass_to_channel_2[pass]);
            float DeltaI3 = val(in_w, p, pass_to_channel_3[pass]);
            float DeltaI4 = val(in_w, p, pass_to_channel_4[pass]);
            #undef val
            if(USE_FELZENSZWALB_POTENTIALS>0) 
               update_costW_trunclinear(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2, 
                     DeltaI1, DeltaI2, DeltaI3, DeltaI4, MGM);
            else
               update_costW(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2,
                  DeltaI1, DeltaI2, DeltaI3, DeltaI4, MGM);
         }
         else {                                 // WITHOUT IMAGE DEPENDENT WEIGHTS
            if(USE_FELZENSZWALB_POTENTIALS>0) {
               if(MGM==2)
                  update_cost2_trunclinear(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
               else 
                  update_costW_trunclinear(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2, 
                        1.0, 1.0, 1.0, 1.0, MGM);
            }
            else if(MGM==2) 
               update_cost2(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], P1, P2);
            else
               update_costW(Lr[pidx], CC[pidx], Lr[pridx], Lr[pr2idx], Lr[pr3idx], Lr[pr4idx], P1, P2, 
                     1.0, 1.0, 1.0, 1.0, MGM);
         }

      }
   }

      // accumulate S for the current orientation
      #pragma omp critical
      {
      for(int i=0; i<nx*ny; i++) {
         for(int o=Lr[i].min;o<=Lr[i].max;o++) {
            S[i].increment_nolock(o, Lr[i][o]); // pragma omp critic is inside set
         }
      }
      }
   }


   // WTA 
   #pragma omp parallel for
   for(int i=0;i<nx*ny;i++) {
      float minP;
      float minL=INFINITY;
      for(int o=S[i].min;o<=S[i].max;o++) {
         // overcounting correction (Drory etal. 2014) 
         if (SGM_FIX_OVERCOUNT==1)
            S[i].set_nolock(o, S[i][o] - (NDIR -1) * CC[i][o]);

         if(isfinite(S[i][o]))
         if(minL > S[i][o]) {
            minL = S[i][o];
            minP = o;
         }
      }	  
      (*out)[i] = minP;
      (*outcost)[i] = minL;
   }

   // return the aggregated costvolume
   return S;
}



