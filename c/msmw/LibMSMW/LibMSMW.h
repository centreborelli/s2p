#ifndef LIBMSMW_H_INCLUDED
#define LIBMSMW_H_INCLUDED


//! Global includes


//! Local includes
#include "../LibImages/LibImages.h"
#include "../Utilities/Parameters.h"


/**
 * @brief Class to perform the MSMW stereo algorithm.
 **/

class MSMW {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    MSMW();


    /**
     * @brief Copy constructor.
     **/
    MSMW(
      const MSMW& i_msmw);


    /**
     * @brief Operator overload.
     **/
    MSMW& operator=(
      const MSMW& i_msmw);


    /**
     * @brief Default destructor.
     **/
    ~MSMW();


    /**
     * @brief Initialize the images and the parameters.
     **/
    void init(
      const Image& i_im1,
      const Image& i_im2,
      const Parameters& p_params);


    /**
     * @brief Run the MSMW estimator.
     **/
    void run(
      Image& o_imDisp1,
      Image& o_imDisp2,
      Image& o_imDist1,
      Image& o_imDist2,
      Image& o_imMask1,
      Image& o_imMask2);


  private:
    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


    /**
     * @brief Recursive call to the MSMW multiscale algorithm.
     **/
    void runMultiscaleChainRecursive();


    /**
     * @brief This function implements multi-window for a current scale.
     **/
    void runChainMultiWindow();


    /**
     * @brief Run the chain for a given window.
     **/
    void runPixelChain();


    /**
     * @brief Pixelian correlation over one direction.
     *
     * @param o_imSelf: will contains the distance for self similarity;
     * @param p_isFirst: if true, considers the left images (1), otherwise
     *                   considers right images (2).
     **/
    void runPixelChain1D(
      Image& o_imSelf,
      const bool p_isFirst);


    /**
     * @brief Check the pixelian reciprocity on both disparity images.
     **/
    void checkPixelianReciprocity();


    /**
     * @brief Update both mask to remove isolated points.
     **/
    void removeIsolatedPoints();


    /**
     * @brief Perform grain filter. Update both mask at once.
     **/
    void performGrainFilter();


    /**
     * @brief Perform MinDist: consider the rejections of the current mask.
     **/
    void checkMinDist();


    /**
     * @brief Check strobe and self similarity effect.
     *
     * @param i_imSelf: self distance image;
     * @param p_isFirst: if true, will work on left images (1), otherwise on
     *                   right images (2).
     **/
    void checkSelf(
      const Image& i_imSelf,
      const bool p_isFirst);


    /**
     * @brief Update Dmin and Dmax for both images.
     **/
    void updateWindowBoundary(
      const Image& p_window);


  //! Data members
  private:

    //! Images
    Image* m_im1;
    Image* m_im2;
    Image* m_imDmin1;
    Image* m_imDmax1;
    Image* m_imDmin2;
    Image* m_imDmax2;
    Image* m_imDisp1;
    Image* m_imDisp2;
    Image* m_imDist1;
    Image* m_imDist2;
    Image* m_imMask1;
    Image* m_imMask2;

    //! Convenience images: precomputed translation
    Image* m_imTrans1;
    Image* m_imTrans2;
    Image* m_imStrobe11;
    Image* m_imStrobe12;
    Image* m_imStrobe21;
    Image* m_imStrobe22;

    //! Windows
    Image* m_windows;

    //! Size
    size_t m_channels;
    size_t m_height;
    size_t m_width;

    //! Parameters
    Parameters* m_params;
    size_t m_currentScale;
    size_t m_nbWindows;
};

#else
class MSMW;

#endif // LIBMSMW_H_INCLUDED
