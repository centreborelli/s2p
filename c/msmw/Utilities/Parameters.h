#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED


//! Global includes
#include <iostream>
#include <string.h>


//! Local includes
#include "../LibImages/LibImages.h"


/**
 * @brief Class for the parameters of the main function.
 **/
class Parameters {

  //! Methods
  public:

    /**
     * @brief Default constructor.
     **/
    Parameters();

    /**
     * @brief Copy constructor.
     **/
    Parameters(
      const Parameters& i_params);

    /**
     * @brief Operator overload.
     **/
    Parameters& operator=(
      const Parameters& i_params);

    /**
     * @brief Default destructor.
     **/
    ~ Parameters();

    /**
     * @brief Getters.
     **/
    char*  inpLeft           () const {return m_inpLeft            ;}
    char*  inpRight          () const {return m_inpRight           ;}
    char*  outDispL          () const {return m_outDispL           ;}
    char*  outDispR          () const {return m_outDispR           ;}
    char*  outMaskL          () const {return m_outMaskL           ;}
    char*  outMaskR          () const {return m_outMaskR           ;}
    float  minDisp           () const {return m_minDisp            ;}
    float  maxDisp           () const {return m_maxDisp            ;}
    size_t orientation       () const {return m_orientation        ;}
    size_t nbScales          () const {return m_nbScales           ;}
    size_t x                 () const {return m_x                  ;}
    size_t y                 () const {return m_y                  ;}
    bool   verbose           () const {return m_verbose            ;}
    Image* window            () const {return m_window             ;}
    size_t grainArea         () const {return m_grainArea          ;}
    float valueRecip         () const {return m_valueRecip         ;}
    float valueRemoveIsolated() const {return m_valueRemoveIsolated;}
    float valueMinDist       () const {return m_valueMinDist       ;}
    float dmin1              () const {return m_dmin1              ;}
    float dmax1              () const {return m_dmax1              ;}
    float dmin2              () const {return m_dmin2              ;}
    float dmax2              () const {return m_dmax2              ;}
    size_t dist              () const {return m_dist               ;}

    /**
     * @brief Read the input arguments. Detect the "-h" option, and print the
     *        informations. Otherwise, get the input arguments, if valid.
     **/
    int checkArgs(
      const int i_argc,
      char** i_argv);

    /**
     * @brief Update the parameters according to the current scale.
     **/
    void update(
      const size_t p_currentScale);

  //! Miscellaneous functions
  private:

    /**
     * @brief Contain the copyright.
     **/
    void printCopyright() const;

    /**
     * @brief Print the synopsis of this program.
     **/
    void printSynopsis() const;

    /**
     * @brief Print the description of this program.
     **/
    void printDescription() const;

    /**
     * @brief Print the input parameters list.
     **/
    void printInput() const;

    /**
     * @brief Print the output parameters list.
     **/
    void printOutput() const;

    /**
     * @brief Print the optional parameters list.
     **/
    void printOptional() const;

    /**
     * @brief Print the signature.
     **/
    void printSignature() const;

    /**
     * @brief Print a line of a specific m_size maximum number of characters,
     *        with a pad start.
     **/
    void printLine(
      const std::string &i_sentence,
      const std::string &i_pad) const;

    /**
     * @brief Print a line of a specific m_size maximum number of characters,
     *        and add a word in the last position.
     **/
    void printWord(
      const std::string &i_line,
      const std::string &i_word,
      const std::string &i_pad = "") const;

    /**
     * @brief Release the memory.
     **/
    void releaseMemory();

  //! Data members
  private:

    //! Misellaneaous
    size_t m_size;          // Maximum number of characters to print by line

    //! Parameters set by users
    char*  m_inpLeft;    // Path of the input left image        -il [%s]
    char*  m_inpRight;   // Path of the input right image       -ir [%s]
    char*  m_outDispL;   // Path of the output left disparity   -dl [%s]
    char*  m_outDispR;   // Path of the output right disparity  -dr [%s]
    char*  m_outMaskL;   // Path of the output left mask        -kl [%s]
    char*  m_outMaskR;   // Path of the output right mask       -kr [%s]
    float  m_minDisp;    // minimum displacement                -m  [%f]
    float  m_maxDisp;    // maximum displacement                -M  [%f]
    size_t m_orientation;// Indicates nb of orientations to use -W  [%d]
    size_t m_nbScales;   // Number of scales                    -n  [%d]
    size_t m_x;          // Window x size                       -x  [%d]
    size_t m_y;          // Window y size                       -y  [%d]
    size_t m_dist;       // Distance to use                     -p  [%d]
    bool   m_verbose;    // Activate the verbose mode           -v

    //! Hard-coded parameters
    Image* m_window;
    size_t m_grainArea;
    float m_valueRecip;
    float m_valueRemoveIsolated;
    float m_valueMinDist;
    float m_dmin1;
    float m_dmax1;
    float m_dmin2;
    float m_dmax2;
};
#else
class Parameters;

#endif // PARAMETERS_H_INCLUDED
