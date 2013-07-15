#define _USE_MATH_DEFINES // For Windows
#include "libNumerics/numerics.h"
#include "libNumerics/matrix.h"
#include "libNumerics/homography.h"
#include "libMatch/match.h"

#include "minDistortionRectify/rectify_mindistortion.hpp"
#include "minDistortionRectify/utils.h"
#include <cmath>
#include <cfloat>

using namespace libNumerics;

namespace rect {

Scene::Scene(int _w, int _h) :
    w(_w),
    h(_h) {
    }

Scene::~Scene(){
}

matrix<flnum> Scene::getPPT() const{
  matrix<flnum> ret(3, 3);
  ret = 0;
  ret(0, 0) = w * w - 1;
  ret(1, 1) = h * h - 1;
  ret *= (w * h / 12.);
  return ret;
}

matrix<flnum> Scene::getpcpcT() const{
  matrix<flnum> ret(3, 3);
  int w1 = w - 1;
  int h1 = h - 1;
  ret = 0;
  ret(0, 0) = w1*w1;
  ret(1, 0) = ret(0, 1) = w1*h1;
  ret(2, 0) = ret(0, 2) = 2 * w1;
  ret(1, 1) = h1*h1;
  ret(2, 1) = ret(1, 2) = 2 * h1;
  ret(2, 2) = 4;
  ret *= 0.25;
  return ret;
}

vector<flnum> Scene::getPoint(double x, double y) const {
  return vector<flnum > (w*x, h*y, 1.);
}

InterestPointsScene::InterestPointsScene(int w, int h, const std::vector< point_t >& _points):
    Scene(w,h),
    points(_points),
    PPT(3,3),
    pcpcT(3,3)
{

  int n = points.size();
  matrix<flnum> P(3,n);
  matrix<flnum> pc(3,1);
  pc = 0.;
  for (int i = 0; i < n; i++) {
    pc(0) += points[i].first;
    pc(1) += points[i].second;
    pc(2)++;
  } // i
  pc/=n;
  pcpcT = pc*pc.t();
  P=0;
  for (int i = 0; i < n; i++) {
    P(0,i) = points[i].first - pc(0);
    P(1,i) = points[i].second  - pc(1);
    P(2,i) = 0.;
  } // i
  PPT = P*P.t();
}

matrix<flnum> InterestPointsScene::getPPT() const{
  return PPT;
}

matrix<flnum> InterestPointsScene::getpcpcT() const{
  return pcpcT;
}

InterestPointsScene::~InterestPointsScene(){
}



StereoRelation::StereoRelation(const matrix<flnum> & _F) :
    F(_F),
    e_left(computeEpipole(F)),
    e_right(computeEpipole(F.t())),
    e_left_cross_product(getCrossProdMatrix(e_left)) {
      vector<flnum> Fe = F*e_left;
      vector<flnum> CrossE = e_left_cross_product*e_left;

    }

const matrix<flnum> & StereoRelation::getF() {
  return F;
}

vector<flnum> StereoRelation::getEpipoleLeft() const {
  return e_left;
}

vector<flnum> StereoRelation::getEpipoleRight() const {
  return e_right;
}

const matrix<flnum>& StereoRelation::getEpipoleLeftCrossProduct() const {
  return e_left_cross_product;
}


StereoModel::StereoModel(StereoRelation _rel,const Scene& _scene_right, const Scene& _scene_left) :
    rel(_rel),
    scene_right(_scene_right),
    scene_left(_scene_left),
    A_l(2, 2),
    B_l(2, 2),
    A_r(2, 2),
    B_r(2, 2) {
    }


StereoHomography::StereoHomography() {
  H_l.setId();
  H_r.setId();
}

StereoHomography::StereoHomography(Homography _H_l, Homography _H_r) :
    H_l(_H_l),
    H_r(_H_r) {
    }

StereoHomography StereoHomography::operator*(const StereoHomography& rhs) const {
  return StereoHomography(H_l * rhs.getHLeft(), H_r * rhs.getHRight());
}

const Homography& StereoHomography::getHLeft() const {
  return H_l;
}

const Homography& StereoHomography::getHRight() const {
  return H_r;
}


StereoHomography StereoModel::computeProjective(StereoHomography& sH){

  initData();

  vector<flnum> z(2);

  minimize_sum_polyratio(A_l, B_l, A_r, B_r, z, false);

  matrix<flnum> ce = rel.getEpipoleLeftCrossProduct();
  matrix<flnum> F = rel.getF();

  vector<flnum> z3(z(0), z(1), 0);

  vector<flnum> w_l = ce*z3;
  vector<flnum> w_r = F*z3;
  w_l *= 1. / w_l(2);
  w_r *= 1. / w_r(2);

  Homography H_r;
  H_r.setId();
  H_r.mat()(2, 0) = w_r(0);
  H_r.mat()(2, 1) = w_r(1);

  Homography H_l;
  H_l.setId();
  H_l.mat()(2, 0) = w_l(0);
  H_l.mat()(2, 1) = w_l(1);

  sH = StereoHomography(H_l, H_r);
  return sH;
}

void StereoModel::initData() {

  matrix<flnum> F = rel.getF();
  matrix<flnum> ce = rel.getEpipoleLeftCrossProduct();

  A_l = (ce.t() * scene_left.getPPT() * ce).copy(0, 1, 0, 1);
  B_l = (ce.t() * scene_left.getpcpcT() * ce).copy(0, 1, 0, 1);

  A_r = (F.t() * scene_right.getPPT() * F).copy(0, 1, 0, 1);
  B_r = (F.t() * scene_right.getpcpcT() * F).copy(0, 1, 0, 1);
}

std::pair<double, double> StereoModel::computeDistortion(const StereoHomography& sH) const {

  vector<flnum> w_l = sH.getHLeft().mat().t().col(2);
  vector<flnum> w_r = sH.getHRight().mat().t().col(2);

  std::pair<double, double> ret;

  double p1 = dot(w_l, scene_left.getPPT() * w_l);
  double q1 = dot(w_l, scene_left.getpcpcT() * w_l);
  ret.first = p1 / q1;

  double p2 = dot(w_r, scene_right.getPPT() * w_r);
  double q2 = dot(w_r, scene_right.getpcpcT() * w_r);
  ret.second = p2 / q2;

  return ret;
}


StereoHomography StereoModel::computeSimilarity(StereoHomography& sH){

  const Homography Hp_l = sH.getHLeft();
  const Homography Hp_r = sH.getHRight();
  const matrix<flnum> F = rel.getF();

  vector<flnum> w_l(Hp_l.mat()(2, 0), Hp_l.mat()(2, 1), 1);
  vector<flnum> w_r(Hp_r.mat()(2, 0), Hp_r.mat()(2, 1), 1);

  Homography H_l, H_r;
  H_l.mat() = 0;
  H_l.mat()(0, 0) = F(2, 1) - w_l(1) * F(2, 2);
  H_l.mat()(0, 1) = -F(2, 0) + w_l(0) * F(2, 2);
  H_l.mat()(1, 0) = F(2, 0) - w_l(0) * F(2, 2);
  H_l.mat()(1, 1) = F(2, 1) - w_l(1) * F(2, 2);
  H_l.mat()(1, 2) = F(2, 2);
  H_l.mat()(2, 2) = 1;

  H_r.mat() = 0;
  H_r.mat()(0, 0) = -F(1, 2) + w_r(1) * F(2, 2);
  H_r.mat()(0, 1) = F(0, 2) - w_r(0) * F(2, 2);
  H_r.mat()(1, 0) = -F(0, 2) + w_r(0) * F(2, 2);
  H_r.mat()(1, 1) = -F(1, 2) + w_r(1) * F(2, 2);
  H_r.mat()(2, 2) = 1;

  //Finds the last y-translation such that the least y-value is 0.
  double least_y_value = DBL_MAX;
  std::vector< vector<flnum> > corners_r;
  corners_r.push_back(scene_right.getPoint(0, 0));
  corners_r.push_back(scene_right.getPoint(0, 1));
  corners_r.push_back(scene_right.getPoint(1, 0));
  corners_r.push_back(scene_right.getPoint(1, 1));

  std::vector< vector<flnum> > corners_l;
  corners_l.push_back(scene_left.getPoint(0, 0));
  corners_l.push_back(scene_left.getPoint(0, 1));
  corners_l.push_back(scene_left.getPoint(1, 0));
  corners_l.push_back(scene_left.getPoint(1, 1));


  for (int i = 0; i < 4; i++) {
    vector<flnum> v = (H_r * Hp_r).mat() * corners_r[i];
    v(1) /= v(2);
    least_y_value = (v(1) < least_y_value) ? v(1) : least_y_value;
  } // i

  for (int i = 0; i < 4; i++) {
    vector<flnum> v = (H_l * Hp_l).mat() * corners_l[i];
    v(1) /= v(2);
    least_y_value = (v(1) < least_y_value) ? v(1) : least_y_value;
  } // i
  double vcp = -least_y_value;
  H_l.mat()(1, 2) = F(2, 2) + vcp;
  H_r.mat()(1, 2) = vcp;
  StereoHomography Hsimi = StereoHomography(H_l, H_r);
  sH=Hsimi*sH;
  return Hsimi;
}

StereoHomography StereoModel::computeShearing(StereoHomography& sH){
  Homography H_l, H_r;
  findShearingHomography(sH.getHLeft(), scene_left, H_l);
  findShearingHomography(sH.getHRight(),scene_right, H_r);
  StereoHomography sH_shearing = StereoHomography(H_l, H_r);
  sH=sH_shearing*sH;
  return sH_shearing;
}

void StereoModel::findShearingHomography(Homography prev, Scene scene, Homography & H) {

  vector<flnum> a = scene.getPoint(.5, 0);
  vector<flnum> b = scene.getPoint(1., .5);
  vector<flnum> c = scene.getPoint(.5, 1.);
  vector<flnum> d = scene.getPoint(0., .5);

  vector<flnum> ah = prev.mat() * a;
  ah *= 1. / ah(2);
  vector<flnum> bh = prev.mat() * b;
  bh *= 1. / bh(2);
  vector<flnum> ch = prev.mat() * c;
  ch *= 1. / ch(2);
  vector<flnum> dh = prev.mat() * d;
  dh *= 1. / dh(2);

  vector<flnum> x = bh - dh;
  vector<flnum> y = ch - ah;

  double w = b(0);
  double h = c(1);
  double as = (h * h * x(1) * x(1) + w * w * y(1) * y(1)) / (h * w * (x(1) * y(0) - x(0) * y(1)));
  double bs = (h * h * x(0) * x(1) + w * w * y(0) * y(1)) / (h * w * (x(0) * y(1) - x(1) * y(0)));

  //We prefer positive a
  if (as < 0) {
    as = -as;
    bs = -bs;
  }

  H.setId();
  H.mat()(0, 0) = as;
  H.mat()(0, 1) = bs;
}

StereoHomography  StereoModel::rectify(){
  StereoHomography ret;
  computeProjective(ret);
  computeSimilarity(ret);
  computeShearing(ret);
  computeFitOutput(ret);
  return ret;
}

StereoHomography StereoModel::computeFitOutput(StereoHomography& sH) const{

  vector<flnum> a_l = scene_left.getPoint(0, 0);
  vector<flnum> b_l = scene_left.getPoint(0., 1);
  vector<flnum> c_l = scene_left.getPoint(1, 0);
  vector<flnum> d_l = scene_left.getPoint(1., 1);

  vector<flnum> a_r = scene_right.getPoint(0, 0);
  vector<flnum> b_r = scene_right.getPoint(0., 1);
  vector<flnum> c_r = scene_right.getPoint(1, 0);
  vector<flnum> d_r = scene_right.getPoint(1., 1);

  vector<flnum> pa_l = sH.getHLeft().mat() * a_l;
  pa_l *= 1. / pa_l(2);
  vector<flnum> pb_l = sH.getHLeft().mat() * b_l;
  pb_l *= 1. / pb_l(2);
  vector<flnum> pc_l = sH.getHLeft().mat() * c_l;
  pc_l *= 1. / pc_l(2);
  vector<flnum> pd_l = sH.getHLeft().mat() * d_l;
  pd_l *= 1. / pd_l(2);

  vector<flnum> pa_r = sH.getHRight().mat() * a_r;
  pa_r *= 1. / pa_r(2);
  vector<flnum> pb_r = sH.getHRight().mat() * b_r;
  pb_r *= 1. / pb_r(2);
  vector<flnum> pc_r = sH.getHRight().mat() * c_r;
  pc_r *= 1. / pc_r(2);
  vector<flnum> pd_r = sH.getHRight().mat() * d_r;
  pd_r *= 1. / pd_r(2);

  //Area of the 2 pictures before the transformation
  double S0 = (getAreaP(b_l - a_l, c_l - a_l) + getAreaP(d_l - b_l, d_l - c_l) +
               getAreaP(b_r - a_r, c_r - a_r) + getAreaP(d_r - b_r, d_r - c_r));

  //Area of the 2 pictures after the transformation
  double S1 = (getAreaP(pb_l - pa_l, pc_l - pa_l) + getAreaP(pd_l - pb_l, pd_l - pc_l) +
               getAreaP(pb_r - pa_r, pc_r - pa_r) + getAreaP(pd_r - pb_r, pd_r - pc_r));

  double r = sqrt(S0 / S1);
  Homography zoom;
  zoom.setZoom(r, r);

  pa_l = zoom.mat() * pa_l;
  pb_l = zoom.mat() * pb_l;
  pc_l = zoom.mat() * pc_l;
  pd_l = zoom.mat() * pd_l;

  pa_r = zoom.mat() * pa_r;
  pb_r = zoom.mat() * pb_r;
  pc_r = zoom.mat() * pc_r;
  pd_r = zoom.mat() * pd_r;

  //lower x-coordinate 0 on each image
  //lower v-coordinate 0 on both  image (must translate simultaneously)
  double least_v = DBL_MAX;
  double least_u_r = DBL_MAX;
  double least_u_l = DBL_MAX;

  least_v = (pa_l(1) < least_v) ? pa_l(1) : least_v;
  least_v = (pb_l(1) < least_v) ? pb_l(1) : least_v;
  least_v = (pc_l(1) < least_v) ? pc_l(1) : least_v;
  least_v = (pd_l(1) < least_v) ? pd_l(1) : least_v;

  least_v = (pa_r(1) < least_v) ? pa_r(1) : least_v;
  least_v = (pb_r(1) < least_v) ? pb_r(1) : least_v;
  least_v = (pc_r(1) < least_v) ? pc_r(1) : least_v;
  least_v = (pd_r(1) < least_v) ? pd_r(1) : least_v;

  least_u_l = (pa_l(0) < least_u_l) ? pa_l(0) : least_u_l;
  least_u_l = (pb_l(0) < least_u_l) ? pb_l(0) : least_u_l;
  least_u_l = (pc_l(0) < least_u_l) ? pc_l(0) : least_u_l;
  least_u_l = (pd_l(0) < least_u_l) ? pd_l(0) : least_u_l;

  least_u_r = (pa_r(0) < least_u_r) ? pa_r(0) : least_u_r;
  least_u_r = (pb_r(0) < least_u_r) ? pb_r(0) : least_u_r;
  least_u_r = (pc_r(0) < least_u_r) ? pc_r(0) : least_u_r;
  least_u_r = (pd_r(0) < least_u_r) ? pd_r(0) : least_u_r;

  Homography Translation_r;
  Homography Translation_l;
  Translation_r.setTrans(-least_u_r, -least_v);
  Translation_l.setTrans(-least_u_l, -least_v);
  Homography H_r = Translation_r*zoom;
  Homography H_l = Translation_l*zoom;
  StereoHomography sH_fix = StereoHomography(H_l, H_r);
  sH=sH_fix*sH;
  return sH_fix;
}


}

