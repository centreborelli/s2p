#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED


//! Global includes
#include <iostream>
#include <string.h>


//! Local includes


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
     * @brief Surcharged constructor.
     **/
    Parameters(
            const bool adjustSize,
            const size_t oWidth,
            const size_t oHeight,
            const bool verbose) {
        m_adjustSize = adjustSize;
        m_oWidth = oWidth;
        m_oHeight = oHeight;
        m_verbose = verbose;
    }

    /**
     * @brief Default destructor.
     **/
    ~ Parameters();

    /**
     * @brief Copy constructor.
     **/
    Parameters(
      const Parameters& i_params);

    /**
     * @brief Operator overload.
     **/
    Parameters& operator=(const Parameters& i_params) {
      new (this) Parameters(i_params);
      return *this;
    }

    /**
     * @brief Getters.
     **/
    char*  inpName   () const {return m_inpName   ;}
    char*  inpHomo   () const {return m_inpHomo   ;}
    char*  outName   () const {return m_outName   ;}
    bool   adjustSize() const {return m_adjustSize;}
    size_t oWidth    () const {return m_oWidth    ;}
    size_t oHeight   () const {return m_oHeight   ;}
    bool   verbose   () const {return m_verbose   ;}

    /**
     * @brief Setters.
     **/
    void setAdjustSize(const bool p_adjustSize) {m_adjustSize = p_adjustSize;}

    /**
     * @brief Read the input arguments. Detect the "-h" option, and print the
     *        informations. Otherwise, get the input arguments, if valid.
     **/
    int checkArgs(
      const int i_argc,
      char** i_argv);


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

  //! Data members
  private:

    //! Misellaneaous
    size_t m_size;          // Maximum number of characters to print by line

    //! Parameters
    char*  m_inpName;    // Path of the input image            -i [%s]
    char*  m_inpHomo;    // Path to the homography             -t [%s]
    char*  m_outName;    // Path of the ouput image            -o [%s]
    bool   m_adjustSize; // Activate the automatic size        -a
    size_t m_oWidth;     // Width of the output image          -c [%d]
    size_t m_oHeight;    // Height of the output image         -l [%d]
    bool   m_verbose;    // Activate the verbose mode          -v
};
#else
class Parameters;

#endif // PARAMETERS_H_INCLUDED
