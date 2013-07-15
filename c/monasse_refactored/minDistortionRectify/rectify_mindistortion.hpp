#ifndef RECTIFY_MINDISTORTION
#define RECTIFY_MINDISTORTION
#include <utility>
#include "libNumerics/matrix.h"
#include "libNumerics/numerics.h"
#include "libNumerics/homography.h"
#include "libMatch/match.h"


using namespace libNumerics;

namespace rect{

/// This class describes the image dimension and it's interest points
/// for the distortion computation. The default interest points are 
/// all the points in the image.
class Scene{
 public:
  /// A Scene is defined simply by its dimensions.
  Scene(int w, int h);

  /// Computes P x P^T 
  virtual matrix<flnum> getPPT() const;

  /// Computes pc x pc^T
  virtual matrix<flnum> getpcpcT() const;

  /// Transform from 0,1 coordinate to 0..w, 0..h coordinate.
  vector<flnum> getPoint(double x, double y) const;

  /// Destructor
  virtual ~Scene();

 protected:
  int w;
  int h;
};

/// Subclass for which we explicetly declare the interest points for
/// the distortion calculation
class InterestPointsScene: public Scene{
 public:
  typedef std::pair<flnum,flnum> point_t;
  /// A scene is defined by its dimension and its interest points (for which
  ///  we want to minimize the distortion)
  InterestPointsScene(int w, int h, const std::vector<point_t>& match);

  /// Computes P x P^T 
  virtual matrix<flnum> getPPT() const;

  /// Computes pc x pc^T
  virtual matrix<flnum> getpcpcT() const;

  /// Destructor
  virtual ~InterestPointsScene();

 private:
  std::vector<point_t> points;
  matrix<flnum> PPT;
  matrix<flnum> pcpcT;
};


/// Model for the Stereo Relationship
class StereoRelation{
 public:

  /// The Stereo relationship is determined by the Fundamental Matrix 
  StereoRelation(const matrix<flnum> & F);

  /// Returns the epipole of the left scene (i.e., such that F*e=0)
  vector<flnum> getEpipoleLeft() const;

  /// Returns the epipole of the left scene (i.e., such that e'^T*F=0)
  vector<flnum> getEpipoleRight() const;

  /// Returns the matrix corresponding to the cross product with the left epipole of F
  const matrix<flnum>& getEpipoleLeftCrossProduct() const;

  /// Returns the Fundamental matrix
  const matrix<flnum> & getF();

 private:
  const matrix<flnum> F;
  const vector<flnum> e_left;
  const vector<flnum> e_right;
  const matrix<flnum> e_left_cross_product;
};

/// Stores two homographies: One for the left image, and another for the rigth
/// image. Homographies are set on the contructor only. 
class StereoHomography{
 public:

  /// Default Constructor. Both Homographies are initalized as the identity
  /// homography.
  StereoHomography();

  /// Explicit contructor where we set each homography.
  StereoHomography( Homography H_l ,Homography H_r);

  /// Performs the composition operation on each Homography
  StereoHomography operator*(const StereoHomography& rhs) const;

  /// Returns H_l, the Homography to be applied on the left image
  const Homography& getHRight() const;
  
  /// Returns H_r, the Homography to be applied on the right image
  const Homography& getHLeft() const ;

 private:
  Homography H_l;
  Homography H_r;
};

/// This is a wrapper for the full caracterization of our problem.
/// It consists on two scenes and a stereo relationship. 
/// It contains the implementation of the rectification algorithm step.
class StereoModel{
 public:

  /// Constructor. The rectification problem is defined by 2 images and a stereo
  /// relationship.
  StereoModel(StereoRelation stereo_lelation, const Scene& scene_left, const Scene& scene_right);

  /// Rectifies the model. Will perform each of the corresponding phases
  /// and return the StereoHomography corresponding to the final rectification.
  StereoHomography rectify();

  /// Computes and returns the Projective transforms of the Homography Decomposition.
  /// Applies to input Stereo Homographies.
  StereoHomography computeProjective(StereoHomography& sH);

  /// Computes and returns the Similarity transforms of the Homography decomposition.
  /// Applies to input Stereo Homographies.
  StereoHomography computeSimilarity(StereoHomography& sH);

  /// Computes and returns the Shearing transforms of the Homography decomposition.
  /// Applies to input Stereo Homographies.
  StereoHomography computeShearing(StereoHomography& sH);

  /// Computes an algorithm to fix the output of the Stereo Homography
  /// such that the sum of areas are conserved.
  /// Applies to input Stereo Homographies.
  StereoHomography computeFitOutput(StereoHomography& sH) const;

  /// Returns the distortion of each side of a StereoHomography according to the criteria of the model.
  std::pair<double,double> computeDistortion(const StereoHomography& sH) const;

 private:

  /// Initializes the from parameters of the model
  void initData();

  /// Computes the Shearing Homography component for a single image
  void findShearingHomography(Homography prev, Scene scene, Homography & H);

  StereoRelation rel;
  const Scene& scene_right;
  const Scene& scene_left;

  //For use in the projective component calculation and distortion 
  //computation
  matrix<flnum> A_l;
  matrix<flnum> B_l;
  matrix<flnum> A_r;
  matrix<flnum> B_r;

};


}

#endif /* end of include guard: RECTIFY_MINDISTORTION */

