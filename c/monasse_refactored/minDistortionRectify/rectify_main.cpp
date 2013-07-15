#include <unistd.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdio>

#include <minDistortionRectify/rectify_mindistortion.hpp>
#include <minDistortionRectify/utils.h>

using namespace rect;

void trace(const StereoModel& model, int phase_number, std::string phase_name,
        StereoHomography sh,StereoHomography result, bool fix) {
    if (fix)
        model.computeFitOutput(result);
    std::cout<<phase_name << ":"<<std::endl;
    std::cout<<"LH"<<phase_number<<" = "<< sh.getHLeft().mat() <<std::endl;;
    std::cout<<"RH"<<phase_number<<" = "<< sh.getHRight().mat() <<std::endl;;
    std::cout<<"LS"<<phase_number<<" = " << result.getHLeft().mat()<<std::endl;;
    std::cout<<"RS"<<phase_number<<" = " << result.getHRight().mat() <<std::endl;;
    std::cout<<std::endl;
}

/// Computes the rectification of the 2 images and prints on standard output
/// each Matrix corresponding to each Homography, as well as a Homography computed
/// from the previous homography so we can visualize each step of the process.
/// if parameter fix==true, it will scale fix the dimensions of the output
/// images.
void compRectifVerbose(
    const Scene& scene_left,
    const Scene& scene_right,
    const matrix<flnum> & F,
    Homography& Hl,
    Homography& Hr,
    bool fix = true
    ){

  StereoHomography cur;
  StereoHomography H;
  StereoModel model(StereoRelation(F), scene_left, scene_right);
  trace(model, 0,"Original State:", H, cur, false);

  // Projective Transform Phase
  H = model.computeProjective(cur) ;
  trace(model,1,"Projective Transform Phase",H , cur, fix);

  // Similarity Transform Phase
  H = model.computeSimilarity(cur);
  trace(model, 2,"Similarity Transform Phase",H , cur, fix);

  // Shearing Transform Phase
  H = model.computeShearing(cur);
  trace(model, 3,"Shearing Transform Phase", H , cur, fix);

  // Fixing Output Format Phase:
  H = model.computeFitOutput(cur);
  trace(model, 4,"Fixing Output Format Phase", H , cur, fix);

  Hl = cur.getHLeft();
  Hr = cur.getHRight();

  std::pair<double,double> dist = model.computeDistortion(cur);
  std::cerr << "Left Image Projective distortion = " << dist.first  << std::endl;
  std::cerr << "Right Image Projective distortion = "<< dist.second << std::endl;
  std::cerr << "Total Distortion = " << dist.first + dist.second   << std::endl;

}


/// Compute rectifying homographies.
/// Rectification invariant degrees of freedom are computed so as to keep image
/// centers fixed.
/// (Wrapper for the StereoModel class)
void compRectif(
    const Scene& scene_left,
    const Scene& scene_right,
    const matrix<flnum> & F,
    Homography& Hl,
    Homography& Hr
    )
{
  StereoModel model(StereoRelation(F), scene_left, scene_right);
  StereoHomography sH = model.rectify();
  Hl = sH.getHLeft();
  Hr = sH.getHRight();
  std::pair<double,double> dist = model.computeDistortion(sH);
  std::cerr << "Left Image Projective distortion = " << dist.first  << std::endl;
  std::cerr << "Rigth Image Projective distortion = "<< dist.second << std::endl;
  std::cerr << "Total Distortion = " << dist.first + dist.second   << std::endl;
}

/// Compute and print min and max disparity
void printDisparity(const std::vector<Match>& match,
                    const libNumerics::Homography& Hl,
                    const libNumerics::Homography& Hr)
{
  std::vector<Match>::const_iterator it=match.begin();
  double min=DBL_MAX, max=-DBL_MAX;
  for(; it != match.end(); ++it) {
    double xl=it->x1, yl=it->y1;
    Hl(xl,yl);
    double xr=it->x2, yr=it->y2;
    Hr(xr,yr);
    xr -= xl;
    if(xr < min)
      min = xr;
    if(xr > max)
      max = xr;
  }
  std::cout << "Disparity: "
      << (int)floor(min) << " " << (int)ceil(max) << std::endl;
}

/// Usage: recify_mindistortion F.txt w h [match.txt] Hl Hr
/// Take as input a file with the Fundamental matrix F and the image dimensions,
/// @w and @h, and output the homographies to apply to left (@Hl) and right
/// (@Hr) images in the form of \f$3\times 3\f$ matrices, stored in Matlab
/// format. F is such that if m belongs to image on the left, and m' on the
/// image on the right, we have m'^T*F*m = 0 .
/// If the matching points list is also given (optional), the distortion is
/// computed in function of these points, and not on all the image points.
/// Outputs on the standard output the homographies corresponding to each phase
/// LH and RH corresponding to the homographies to be applied to the left and
/// right image respectively.
/// Also outputs the homographies corresponding to the cumulative phases. LS and
/// RS analogously. These last are filtered to fit the output format well, by
/// conservating the sum of image areas.
int main(int argc, char** argv)
{

  if(argc != 6 && argc != 7 ) {
    std::cerr << "Usage: " << argv[0] << " F.txt w h [match.txt] Hl Hr" <<std::endl;
    return 1;
  }

  //Reading F matrix
  std::ifstream ffile(argv[1]);
  matrix<flnum> F(3,3);
  ffile >> F;

  int w=0,h=0;
  if(! (std::istringstream(argv[2]) >> w).eof()) w=0;
  if(! (std::istringstream(argv[3]) >> h).eof()) h=0;
  if(w <=0 || h <= 0) {
    std::cerr << "Wrong dimensions of image" << std::endl;
    return 1;
  }

  bool isVerbose = true;

  Scene* scene_left = NULL;
  Scene* scene_right = NULL;
  // Using matches
  int usesMatches=0;
  std::vector<Match> match;
  if (argc == 7){
    usesMatches=1;
    if(! loadMatch(argv[4],match)) {
      std::cerr << "Failed reading " << argv[7] << std::endl;
      return 1;
    }
    typedef  InterestPointsScene::point_t pdd;
    std::vector< pdd > points_left;
    std::vector< pdd > points_right;
    for (unsigned int i = 0; i < match.size(); i++) {
      points_left.push_back(pdd (match[i].x1, match[i].y1) );
      points_right.push_back(pdd (match[i].x2, match[i].y2) );
    } // i
    scene_left  = new InterestPointsScene(w,h,points_left);
    scene_right = new InterestPointsScene(w,h,points_right);
  }else{
    scene_left  = new Scene(w,h);
    scene_right = new Scene(w,h);
  }

  libNumerics::Homography Hl, Hr;
  if (isVerbose){
    compRectifVerbose(*scene_left, *scene_left, F, Hl, Hr);
  }else{
    compRectif(*scene_left, *scene_right, F, Hl, Hr);
  }

  if (!match.empty())
    printDisparity(match, Hl, Hr);

  std::ofstream f1(argv[4+usesMatches]), f2(argv[5+usesMatches]);

  if((f1 << Hl.mat() << std::endl).fail()) {
    std::cerr << "Error writing file " << argv[4+usesMatches] << std::endl;
    return 1;
  }

  if((f2 << Hr.mat() << std::endl).fail()) {
    std::cerr << "Error writing file " << argv[5+usesMatches] << std::endl;
    return 1;
  }

  delete scene_right;
  delete scene_left;

  return 0;
}

