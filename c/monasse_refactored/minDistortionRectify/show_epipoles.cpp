#include <iostream>

#include "libNumerics/numerics.h"
#include "libNumerics/matrix.h"
#include "libMatch/match.h"
#include "libNumerics/homography.h"
#include "libIO/draw.h"
#include "libIO/io_png.h"
#include "libLWImage/LWImage.h"

#include "minDistortionRectify/utils.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>


static const color CYAN={0,255,255};
static const color RED={255,0,0};
static const int CROSS_HALF_LENGTH = 3; ///< Half-length of cross for SIFT
static const int LINE_RESOLUTION=500.;
static const int MAX_NUMBER_PLOTS=10.;

double max3abs(double a, double b, double c){
  return std::max(std::abs(a),std::max(std::abs(b),std::abs(c)));
}

/// Draw a point at the image
static void drawPoint(unsigned char* data, int w, int h,
                     double dx, double dy, struct color c){

  if(dx<0||dx>=w||dy<0||dy>=h) return;
  int x = dx;
  int y = dy;
  data[x + y*w + 0*w*h] = c.r;
  data[x + y*w + 1*w*h] = c.g;
  data[x + y*w + 2*w*h] = c.b;
}

///  Draw a generic line in image. This line is represent by its epipolar
/// representation \f$l = [l_u,l_v,l_w]\f$.
///  \a data must have 3 non-interlaced channels.
void draw_line(unsigned char* data, int w, int h,
               double lu, double lv, double lw, struct color c)
{

  // Max value in the lines coordinate = 1
  double p =1./max3abs(lu,lv,lw);
  lu*=p;
  lv*=p;
  lw*=p;

  // Ceils to zero small values
  if(std::abs(lu)<1e-5) lu=0;
  if(std::abs(lv)<1e-5) lv=0;
  if(std::abs(lw)<1e-5) lw=0;

  if(lv != 0){
    // line [ y = f(x) = (-x*lu -lw)/lv ]
    double dx = (double)w/LINE_RESOLUTION;
    for(double x=0;x<w; x+= dx){
      drawPoint(data,w,h,lrint(x),lrint((-x*lu-lw)/lv),c );
    }
  }
  else 
    if (lu != 0){
      // line [x = f(y) = (-y*lv -lw )/lu]
      double dy = (double)h/LINE_RESOLUTION;
      for(double y=0;y<h; y+= dy){
        drawPoint(data,w,h,lrint((-y*lv-lw)/lu) ,lrint(y), c);
      }

    }
}

/// Usage: show_epipoles in.png out.png match.txt left|right H.txt F.txt
/// Take as background image @in.png, superimpose lines and small crosses 
/// at left or right points in @match.txt to which the
/// homography given as a matrix in file @H.txt (Matlab text format) is applied.
/// The ouput image file (PNG format) is @out.png. The argument before @H.txt
/// must be the string "left" or "right". The fundamental Matrix @F must also
/// be provided so that the epipoles may be computed. This matrix is such that
/// \f$ m'^TFm=0\f$ where \f$m,m'\f$ are correspondence points from left and
/// right image respectively.
int main(int argc, const char *argv[])
{

  if(argc < 7 && argc > 8) {
    std::cerr << "Usage: " << argv[0]
        << " in.png out.png match.txt left|right H.txt F.txt nmax(default=10)" <<std::endl;
    return 1;
  }

  // Read image
  size_t nx, ny;
  LWImage<unsigned char> in(0,0,0);
  in.data = read_png_u8_rgb(argv[1], &nx, &ny);
  in.w = static_cast<int>(nx); in.h = static_cast<int>(ny); in.comps=3;
  if(! in.data) {
    std::cerr << "Error reading image " << argv[1] << std::endl;
    return 1;
  }

  // Read correspondences
  std::vector<Match> match;
  if(! loadMatch(argv[3],match)) {
    std::cerr << "Failed reading " << argv[3] << std::endl;
    return 1;
  }

  // Read left/right image
  bool bLeft = (strcmp(argv[4],"left")==0);
  if(!bLeft && strcmp(argv[4], "right")!=0) {
    std::cerr << "Error: arg 4 must be 'left' or 'right'" << std::endl;
    return 1;
  }

  // Read homography
  libNumerics::Homography H;
  std::ifstream f(argv[5]);
  if((f >> H.mat()).fail()) {
    std::cerr << "Error reading homography file " << argv[2] << std::endl;
    return 1;
  }

  //Reading F matrix
  std::ifstream ffile(argv[6]);
  libNumerics::matrix<double> F(3,3);
  ffile >> F;

  int nplots = MAX_NUMBER_PLOTS;
  if (argc==8){
    nplots = atoi(argv[7]);
    if(nplots < 0){
      std::cerr << "Error converting max number of plots."<<std::endl;
      return 1;
    }
  }

  //Computing epipole
  if (!bLeft) 
    F=F.t();
  libNumerics::vector<double> epipole = rect::computeEpipole(F);
  epipole = H.mat()*epipole;
  double a = max3abs(epipole(0),epipole(1),epipole(2));
  epipole *= 1./a;
  libNumerics::matrix<double> ce =  rect::getCrossProdMatrix(epipole);

  libNumerics::vector<double> x(3);
  libNumerics::vector<double> l(3);
  // Drawing SIFT and EPipolar lines
  std::vector<Match>::const_iterator it=match.begin();
  for(; it != match.end(); ++it) {
    if(--nplots < 0) 
      break;
    x= libNumerics::vector<double>(it->x1, it->y1,1);
    if(! bLeft) {
      x= libNumerics::vector<double>(it->x2, it->y2,1);
    }
    x=H.mat()*x;
    x*=1./x(2);
    l=ce*x;


    int ix = static_cast<int>(std::floor(x(0)+0.5));
    int iy = static_cast<int>(std::floor(x(1)+0.5));
    if(0<=ix && ix<in.w && 0<=iy && iy<in.h){
      draw_cross(in.data, in.w, in.h, ix, iy, CROSS_HALF_LENGTH, RED);
      draw_line(in.data, in.w, in.h, l(0),l(1),l(2), CYAN);
    }
  }

  // Write image
  if(write_png_u8(argv[2], in.data, in.w, in.h, in.comps) != 0) {
    std::cerr << "Error writing file " << argv[3] << std::endl;
    return 1;
  }
  free(in.data);
  return 0;
}
