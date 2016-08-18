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
    char*  input    () const {return m_input    ;}
    char*  outPath  () const {return m_outPath  ;}
    bool   verbose  () const {return m_verbose  ;}
    size_t nbOct    () const {return m_nbOct    ;}
    size_t nbSpo    () const {return m_nbSpo    ;}
    size_t nbHist   () const {return m_nbHist   ;}
    size_t nbBins   () const {return m_nbBins   ;}
    size_t nbOri    () const {return m_nbOri    ;}
    size_t maxIter  () const {return m_maxIter  ;}
    float  sigmaMin () const {return m_sigmaMin ;}
    float  deltaMin () const {return m_deltaMin ;}
    float  sigmaIn  () const {return m_sigmaIn  ;}
    float  dog      () const {return m_dog      ;}
    float  edge     () const {return m_edge     ;}
    float  lambdaOri() const {return m_lambdaOri;}
    float  lambdaDes() const {return m_lambdaDes;}
    float  t        () const {return m_t        ;}

    /**
     * @brief Setters.
     **/
    void setDefaultValues();
    void set_thresh_dog(float t) {m_dog = t;}
    void set_noct(int n) {m_nbOct = (size_t) n;}
    void set_nspo(int n) {m_nbSpo = (size_t) n;}

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

    /**
     * @brief Release the memory.
     **/
    void releaseMemory();

  //! Data members
  private:

    //! Misellaneaous
    size_t m_size;       // Maximum number of characters to print by line

    //! Parameters set by users
    char*  m_input;      // Path of the input image             -i [%s]
    char*  m_outPath;    // Path of the output folder           -o [%s]
    bool   m_verbose;    // Activate the verbose mode           -v

    //! Sift specific parameters
    size_t m_nbOct;      // Number of octave                    -no [%d]
    size_t m_nbSpo;      //                                     -ns [%d]
    size_t m_nbHist;     // Number of histogram                 -nh [%d]
    size_t m_nbBins;     // Number of bins                      -nb [%d]
    size_t m_nbOri;      // Number of orientation               -nr [%d]
    size_t m_maxIter;    // Maximum number of iterations        -m  [%d]
    float m_sigmaMin;    //                                     -sm [%f]
    float m_deltaMin;    //                                     -dm [%f]
    float m_sigmaIn;     //                                     -si [%f]
    float m_dog;         //                                     -d  [%f]
    float m_edge;        //                                     -e  [%f]
    float m_lambdaOri;   //                                     -lo [%f]
    float m_lambdaDes;   //                                     -ld [%f]
    float m_t;           //                                     -t  [%f]
};
#else
class Parameters;

#endif // PARAMETERS_H_INCLUDED
