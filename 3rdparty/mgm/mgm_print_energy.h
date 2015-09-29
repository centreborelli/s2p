#ifndef MGM_PRINT_ENERGY_H_
#define MGM_PRINT_ENERGY_H_
/******************** PRINT ENERGY ********************/

#include "img.h"
#include "img_tools.h"
#include "point.h"

// EVALUATE THE ENERGY OF THE CURRENT SOLUTION
// type==0 : truncated L1
// type==1 : L1
// type==2 : L2
float evaluate_energy_4connected(const struct Img &u, std::vector<float > &outdisp, struct costvolume_t &CC, struct Img *E_out, float P1, float P2, int type) {
   int nx=u.nx;
   int ny=u.ny;

   // EVALUATE THE ENERGY OF THE CURRENT SOLUTION
   float ENERGY = 0;
   float ENERGYL2 = 0;
   float ENERGYtrunc = 0;

   struct Img E(nx,ny);
   struct Img EL2(nx,ny);
   struct Img Etrunc(nx,ny);

   for(int y=0; y<ny; y++)
   for(int x=0; x<nx; x++) 
   {
      Point p(x,y);				   // current point
      int pidx   = (p.x +p.y *nx);

      E[pidx] = 0;
      EL2[pidx] = 0;
      Etrunc[pidx] = 0;

      float G     =0; // contribution for the current point 
      float GL2   =0; // contribution for the current point 
      float Gtrunc=0; // contribution for the current point 

      // DATA TERM
      int o = outdisp[pidx];
      G      += CC[pidx][o];
      GL2      += CC[pidx][o];
      Gtrunc += CC[pidx][o];

      // EDGE POTENTIALS
      Point directions[] = {Point(-1,0) , Point(0,1), 
                            Point(1,0)  , Point(0,-1),
                            Point(-1,0) };

      int N = 4; 
      for(int t=0;t<N;t++) {
         Point pr = p + directions[t]; 
         Point pq = p + directions[t+1]; // needed for L2
         if (!check_inside_image(pr ,u)) continue;
         if (!check_inside_image(pq ,u)) continue;
         int pridx = (pr.x+pr.y*nx);
         int pqidx = (pq.x+pq.y*nx);
         int oor   = outdisp[pridx];
         int ooq   = outdisp[pqidx];

         G   += fabs(oor - o)/N;

         GL2 += sqrt( (oor-o)*(oor-o) + (ooq-o)*(ooq-o))/N;

         if(fabs(oor - o) <= 1) Gtrunc += P1/N;
         else Gtrunc += P2/N;
      }

      ENERGY += G;
      E[pidx] = G;

      ENERGYL2 += GL2;
      EL2[pidx] = GL2;

      ENERGYtrunc += Gtrunc;
      Etrunc[pidx] = G;

   }
   if(type==1) {
      *E_out=E;
      return ENERGY;
   }
   if(type==2) {
      *E_out=EL2;
      return ENERGYL2;
   }
   *E_out=Etrunc;
   return ENERGYtrunc;
   

}



void print_solution_energy(const struct Img &in_u, std::vector<float > &disp, struct costvolume_t &CC, float P1, float P2) {
   if (TSGM_DEBUG()) {
      // DEBUG INFORMATION
      struct Img E;
      printf(" ENERGY L1trunc: %.9e\t", evaluate_energy_4connected(in_u,disp,CC,&E,P1,P2,0));
      iio_write_vector_split((char*)"/tmp/ENERGY_L1trunc.tif", E);
      printf("L1: %.9e\t", evaluate_energy_4connected(in_u,disp,CC,&E,P1,P2,1));
      printf("L2: %.9e\n", evaluate_energy_4connected(in_u,disp,CC,&E,P1,P2,2));
   }
   else {
      printf("\n");
   }
}


#endif //MGM_PRINT_ENERGY_H_
