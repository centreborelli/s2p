#ifndef LIBSIFT_H_INCLUDED
#define LIBSIFT_H_INCLUDED


//! Global includes
#include <list>


//! Local includes
#include "KeyPoint.h"
#include "ScaleSpace.h"
#include "../LibImages/LibImages.h"
#include "../Utilities/Parameters.h"
#include "../LibSSE/LibSSE.h"
#include "../Utilities/Time.h"


/**
 * @brief Class that implements the Sift method.
 **/
class Sift {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    Sift();

    /**
     * @brief Surcharged constructor.
     **/
    Sift(
      const Parameters& p_params);


    /**
     * @brief Copy constructor.
     **/
    Sift(
      const Sift& i_sift);


    /**
     * @brief Operator overload.
     **/
    Sift& operator=(
      const Sift& i_sift);


    /**
     * @brief Default destructor.
     **/
    ~Sift();


    /**
     * @brief Compute the list of keypoints from the input image.
     **/
     void computeKeyPoints(
      const Image& i_im);


    /**
     * @brief Accumulate orientation histogram.
     *
     * @param io_key: current KeyPoint to update.
     **/
    void accumulateOrientationHistogram(
      KeyPoint* io_key);


    /**
     * @brief Get the keypoints.
     **/
    void getKeyPoints(
      std::vector<KeyPoint>& o_keyPoints) const;


    /**
     * @brief Write the keypoints into a file.
     **/
    void writeKeyPoints(
      std::string const& i_fileName) const;


  private:


    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


    /**
     * @brief Initialization
     **/
    void init(
      const Image& i_im);


    /**
     * @brief Determine the number of octaves
     **/
    size_t getNbOctaves() const;


    /**
     * @brief Adapt the threshold to the scalespace discretization.
     **/
    float convertThreshold() const;


    /**
     * @brief Compute the Gaussian scalespace for SIFT.
     **/
    void computeScaleSpace();


    /**
     * @brief Compute the difference of Gaussians.
     **/
    void computeDoG();


    /**
     * @brief Find 3D discrete extrema.
     **/
    void find3dDiscreteExtrema();


    /**
     * @brief Discard keypoints with low response.
     **/
    void discardKeyPointsWithLowResponse(
      const float p_threshold);


    /**
     * @brief Interpolate keypoints positions.
     **/
    void interpolateKeyPointsPosition();


    /**
     * @brief Compute edge response.
     **/
    void computeEdgeResponse();


    /**
     * @brief Discard keypoints on edge.
     **/
    void discardKeyPointsOnEdge();


    /**
     * @brief Discard keypoints near the border.
     **/
    void discardKeyPointsNearBorder();


    /**
     * @brief Compute the scalespace gradient.
     **/
    void computeScaleSpaceGradient();


    /**
     * @brief Attribute the orientation for all the keypoints.
     **/
    void attributeKeyPointsOrientations();


    /**
     * @brief Attribute the descriptors for all the keypoints.
     **/
    void attributeKeyPointsDescriptors();


  //! Data members
  public:
    //! List of keypoints,
    std::list<KeyPoint*>* m_keyPoints;

  private:

    //! Set of parameters
    Parameters* m_params;

    //! Image (gray version)
    Image* m_im;
    size_t m_width;
    size_t m_height;

    //! Scale spaces
    ScaleSpace* m_s;
    ScaleSpace* m_d;
    ScaleSpace* m_sx;
    ScaleSpace* m_sy;

    //! Miscellaneous
    Time* m_time;
};
#else
class Sift;

#endif // LIBSIFT_H_INCLUDED
