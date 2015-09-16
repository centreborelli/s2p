#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))
// fast alternatives to: __min(a,__min(b,c)) 
// fastestest ?
#define fmin3_(x, y, z) \
   (((x) < (y)) ? (((z) < (x)) ? (z) : (x)) : (((z) < (y)) ? (z) : (y)))
// fast and easy to understand
//static inline float fmin3(float a, float b, float c)
//{
//   float m = a;
//   if (m > b) m = b;
//   if (m > c) m = c;
//   return m;
//}

inline float fastexp(float x) {
   int result = static_cast<int>(12102203 * x) + 1065353216;
   result *= result > 0;
   std::memcpy(&x, &result, sizeof(result));
   return x;
}

inline float deltaImage(const struct Img &u, const Point p, const Point q)
{
   float d = 0;
   for (int c = 0; c < u.nch; c++){
      float diff = val(u, p, c) - val(u, q, c);
//      d += diff > 0 ? diff : -diff;
      d += diff * diff;
//      d = __max(d, fabs(val(u, p, c) - val(u, q, c))) ;
   }
	return d/u.nch;
}

inline float ws( float DeltaI, float aP3, float Thresh) {

   if (fabs(DeltaI) < Thresh*Thresh) return aP3;
   else return 1;

//// implements weights from "simple but effective tree structure for dynamic programming-based stereo matching" Blayer, Gelautz
//   float T=30; 
//   float P3=4;
//   if (fabs(DeltaI) < T) return P3;
//   else return 1;
//
//// implements image adaptive weights from formula 1 of "Efficient High-Resolution Stereo Matching using Local Plane Sweeps"
//   float sigmaI=128;
//   float alpha=10;
//   return (8 + alpha * fastexp( - fabs(DeltaI) / sigmaI))/32.0;
//   return 8.0/32.0 + 1.0/__max(fabs(DeltaI),0.0001);
//
//   return 1.;
}




// For a pixel p each image of the stack correspond the weights to 
// int neighbouring pixels: W, E, S, N, (NW, NE, SE, SW)
struct Img compute_mgm_weights(struct Img &u, float aP, float aThresh) 
{
   int nx = u.nx;
   int ny = u.ny;
   // load the edge weights (or initialize them to 1: TODO)
   struct Img w(nx, ny, 8);
   Point scans[] = {Point(-1,0), Point(1,0), Point(0,1), Point(0,-1), Point(-1,-1), Point(1,-1), Point(1,1), Point(-1,1)};

   for (int o = 0; o < 8 ; o++) 
   for (int j = 0; j < ny; j++)
   for (int i = 0; i < nx; i++)
      {
         float wvalue = 1.0;
         Point p(i,j);				   // current point
         Point pr  = p + scans[o];   // neighbor
         if (check_inside_image(pr ,u))  {
            float Delta = deltaImage(u,p,pr);
            wvalue = ws(Delta, aP, aThresh);
         }
         w[i + j*nx + o*nx*ny] = wvalue;
      }
   return w;
}
