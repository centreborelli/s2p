#ifndef KEYPOINT_H_INCLUDED
#define KEYPOINT_H_INCLUDED


//! Global includes
#include <cstdlib>


//! Local includes
#include "ScaleSpace.h"
#include "../Utilities/Parameters.h"


/**
 * @brief keypoint structure, related to a keypoint.
 *
 * stores SIFT output (keypoint position, scale, orientation, feature vector...)
 * stores intermediary results (orientation histogram, summarizes...)
 **/
class KeyPoint {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    KeyPoint();


    /**
     * @brief Copy constructor.
     **/
    KeyPoint(
      const KeyPoint& i_keyPoint);


    /**
     * @brief Overload constructor.
     **/
    KeyPoint(
      const size_t i_nbOri,
      const size_t i_nbHist,
      const size_t i_nbBins);


    /**
     * @brief Operator overload.
     **/
    KeyPoint& operator=(
      const KeyPoint& i_keyPoint);


    /**
     * @brief Default destructor.
     **/
    ~KeyPoint();


    /**
     * @brief Setters.
     **/
    void setX       (const float i_x) {m_x        = i_x;}
    void setY       (const float i_y) {m_y        = i_y;}
    void setSigma   (const float i_s) {m_sigma    = i_s;}
    void setTheta   (const float i_t) {m_theta    = i_t;}
    void setO       (const int   i_o) {m_o        = i_o;}
    void setS       (const int   i_s) {m_s        = i_s;}
    void setI       (const int   i_i) {m_i        = i_i;}
    void setJ       (const int   i_j) {m_j        = i_j;}
    void setVal     (const float i_v) {m_val      = i_v;}
    void setEdgeResp(const float i_e) {m_edgeResp = i_e;}

    /**
     * @brief Getters.
     **/
    float  getX       () const {return m_x       ;}
    float  getY       () const {return m_y       ;}
    float  getSigma   () const {return m_sigma   ;}
    float  getTheta   () const {return m_theta   ;}
    int    getO       () const {return m_o       ;}
    int    getS       () const {return m_s       ;}
    int    getI       () const {return m_i       ;}
    int    getJ       () const {return m_j       ;}
    float  getVal     () const {return m_val     ;}
    size_t getNbOri   () const {return m_nbOri   ;}
    size_t getNbHist  () const {return m_nbHist  ;}
    size_t getNbBins  () const {return m_nbBins  ;}
    float  getEdgeResp() const {return m_edgeResp;}
    float  getDescr   (const size_t p_n) const {return m_descriptors[p_n];}
    float  getOriHist (const size_t p_n) const {return m_oriHist    [p_n];}
    float* getPtrHist () {return m_oriHist;}
    float* getPtrDescr() {return m_descriptors;}


    /**
     * @brief Extract principal orientations from gradient orientation histogram,
     *        and update the histogram.
     *
     * @param o_principalOrientations: pointer to tab storing principal orientations
     *
     * @return the number of principal orientation (>=1).
     */
    int extractPrincipalOrientations(
      const int p_nbBins,
      const float p_threshold,
      float* o_principalOrientations);


    /**
     * @brief Extract feature vector and update the current keypoint.
     **/
    void extractFeatureVector(
      const ScaleSpace* i_sx,
      const ScaleSpace* i_sy,
      const Parameters& p_params);


    /**
     * @brief Threshold and quantize the descriptors.
     *
     * @param p_n: number of value to deal with.
     **/
    void thresholdAndQuantizeFeatureVector(
      const size_t p_n);


  private:


    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


  //! Data members
  private:
    float m_x; // coordinates
    float m_y;
    float m_sigma; // level of blur (it includes the assumed image blur)
    float m_theta; // orientation

    int m_o; // discrete coordinates
    int m_s;
    int m_i;
    int m_j;

    float m_val; // normalized operator value (independant of the scalespace sampling)
    float m_edgeResp; // edge response
    size_t m_nbHist;     // number of histograms in each direction
    size_t m_nbOri;      // number of bins per histogram
    float* m_descriptors;

    size_t m_nbBins; // number of bins in the orientation histogram
    float* m_oriHist; // gradient orientation histogram
};
#else
class KeyPoint;



#endif // KEYPOINT_H_INCLUDED
