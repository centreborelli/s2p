/**
 * @file Parameters.cpp
 *
 * @brief Parameters class: will initialize and contains the parameters used
 *        for the main function.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


//! Global includes
#include <stdlib.h>
#include <sstream>
#include <vector>


//! Local includes
#include "Parameters.h"


using namespace std;


//! Default constructor
Parameters::Parameters(
  ) :
  //! Miscellaneous
  m_size     (43 + 2 * 15),  // Number of characters per line to print

  //! Parameters set by users
  m_input    (NULL),         // Path of the input image             -i [%s]
  m_outPath  (NULL),         // Path of the output folder           -o [%s]
  m_verbose  (false),        // Verbose option                      -v

  //! Sift specific parameters
  m_nbOct    (8),            // Number of octave                    -no [%d]
  m_nbSpo    (3),            //                                     -ns [%d]
  m_nbHist   (4),            // Number of histogram                 -nh [%d]
  m_nbBins   (36),           // Number of bins                      -nb [%d]
  m_nbOri    (8),            // Number of orientation               -nr [%d]
  m_maxIter  (5),            // Maximum number of iterations        -m  [%d]
  m_sigmaMin (0.8f),         //                                     -sm [%f]
  m_deltaMin (0.5f),         //                                     -dm [%f]
  m_sigmaIn  (0.5f),         //                                     -si [%f]
  m_dog      (0.013333333f), //                                     -d  [%f]
  m_edge     (10),           //                                     -e  [%f]
  m_lambdaOri(1.5f),         //                                     -lo [%f]
  m_lambdaDes(6.f),          //                                     -ld [%f]
  m_t        (0.8f) {       //                                      -t  [%f]

}


//! Default destructor
Parameters::~Parameters(
  ) {

  //! Release the memory
  releaseMemory();
}


//! Operator overload.
Parameters& Parameters::operator=(
  const Parameters& i_params) {

  if (&i_params == this) {
    return *this;
  }

  releaseMemory();

  new (this) Parameters(i_params);
  return *this;
}


//! Copy constructor.
Parameters::Parameters(
  const Parameters& i_params) :
  //! Miscalleneous
  m_size(i_params.m_size),

  //! Parameters set by users
  m_input  (i_params.input  ()),
  m_outPath(i_params.outPath()),
  m_verbose(i_params.verbose()),

  //! Sift specific parameters
  m_nbOct    (i_params.m_nbOct),
  m_nbSpo    (i_params.m_nbSpo),
  m_nbHist   (i_params.m_nbHist),
  m_nbBins   (i_params.m_nbBins),
  m_nbOri    (i_params.m_nbOri),
  m_maxIter  (i_params.m_maxIter),
  m_sigmaMin (i_params.m_sigmaMin),
  m_deltaMin (i_params.m_deltaMin),
  m_sigmaIn  (i_params.m_sigmaIn),
  m_dog      (i_params.m_dog),
  m_edge     (i_params.m_edge),
  m_lambdaOri(i_params.m_lambdaOri),
  m_lambdaDes(i_params.m_lambdaDes),
  m_t        (i_params.m_t) {

}


//! Read the input arguments.
int Parameters::checkArgs(
  const int i_argc,
  char** i_argv) {

  //! Release memory
  releaseMemory();

  //! If no arguments entered
  if (i_argc <= 1) {
    cout << "No arguments detected. Type -h or -help for some help." << endl;
    return EXIT_FAILURE;
  }

  //! Detect the parameter entered
  for (int n = 1; n < i_argc; n++) {

    //! Get the current argument
    const string sarg = i_argv[n];

    //! If help is required
    if (sarg.find("-h") == 0 || sarg.find("-help") == 0) {

      //! Copyright
      this->printCopyright();

      //! Synopsis
      this->printSynopsis();

      //! Description
      this->printDescription();

      //! Input parameters
      this->printInput();

      //! Output parameters
      this->printOutput();

      //! Other options
      this->printOptional();

      //! Signature
      this->printSignature();

      exit(EXIT_FAILURE);
    }

    //! Input path
    if (sarg == "-i") {
      if (n + 1 < i_argc) {
        m_input = i_argv[n + 1];
      }
    }

    //! Output path
    if (sarg == "-o") {
      if (n + 1 < i_argc) {
        m_outPath = i_argv[n + 1];
      }
    }

    //! Number of octaves
    if (sarg == "-no") {
      if (n + 1 < i_argc) {
        m_nbOct = atoi(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-ns") {
      if (n + 1 < i_argc) {
        m_nbSpo = atoi(i_argv[n + 1]);
      }
    }

    //! Number of histograms
    if (sarg == "-nh") {
      if (n + 1 < i_argc) {
        m_nbHist = atoi(i_argv[n + 1]);
      }
    }

    //! Number of bins
    if (sarg == "-nb") {
      if (n + 1 < i_argc) {
        m_nbBins = atoi(i_argv[n + 1]);
      }
    }

    //! Number of orientations
    if (sarg == "-no") {
      if (n + 1 < i_argc) {
        m_nbOri = atoi(i_argv[n + 1]);
      }
    }

    //! Maximum number of iterations
    if (sarg == "-m") {
      if (n + 1 < i_argc) {
        m_maxIter = atoi(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-sm") {
      if (n + 1 < i_argc) {
        m_sigmaMin = atof(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-dm") {
      if (n + 1 < i_argc) {
        m_deltaMin = atof(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-si") {
      if (n + 1 < i_argc) {
        m_sigmaIn = atof(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-d") {
      if (n + 1 < i_argc) {
        m_dog = atof(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-e") {
      if (n + 1 < i_argc) {
        m_edge = atof(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-lo") {
      if (n + 1 < i_argc) {
        m_lambdaOri = atof(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-ld") {
      if (n + 1 < i_argc) {
        m_lambdaDes = atof(i_argv[n + 1]);
      }
    }

    //!
    if (sarg == "-t") {
      if (n + 1 < i_argc) {
        m_t = atof(i_argv[n + 1]);
      }
    }

    //! Verbose option
    if (sarg == "-v") {
      m_verbose = true;
    }
  }

  //! Check that the mandatory arguments have been given
  if (m_input == NULL) {
    cout << "The input image path must be precised. Use \n";
    cout << "    -i path/image.ext" << endl;
    exit(EXIT_FAILURE);
  }
  if (m_outPath == NULL) {
    cout << "The output path must be precised. Use\n";
    cout << "    -o path/" << endl;
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}


//! Contain the copyright.
void Parameters::printCopyright(
  ) const {

  const string s0 = "* -------------------------------- *";
  const string s1 = "| Copyright (c) 2015-2016 by IPOL. |";
  const string s2 = "|       All Rights Reserved.       |";
  const string st = std::string((m_size - s0.size()) / 2, ' ');

  cout << endl;
  cout << st << s0 << st << endl;
  cout << st << s1 << st << endl;
  cout << st << s2 << st << endl;
  cout << st << s0 << st << endl;
  cout << endl;
}


//! Print the synopsis of this program.
void Parameters::printSynopsis(
  ) const {

  //! Header
  const string header = "SYNOPSIS:";

  //! Name of the program
  const string name = "./sift";

  //! Parameters
  string sentence = name;
  sentence += " [-i  input image path]";
  sentence += " [-o  output path]";
  sentence += " [-no nb of octaves]";
  sentence += " [-ns nb of spo]";
  sentence += " [-nh nb of histograms]";
  sentence += " [-nb nb of bins]";
  sentence += " [-nr nb of orientations]";
  sentence += " [-m  max iterations]";
  sentence += " [-sm sigma min]";
  sentence += " [-dm delta min]";
  sentence += " [-si sigma in]";
  sentence += " [-d  DoG]";
  sentence += " [-e  edge]";
  sentence += " [-lo lambda orientation]";
  sentence += " [-ld lambda descriptors]";
  sentence += " [-t  T]";
  sentence += " [-v verbose]";

  //! Print the synopsis
  cout << header << endl << endl;

  //! Print the content
  this->printLine(sentence, "  ");

  //! Finish the synopsis
  cout << endl << endl;
}


//! Print the description of this program.
void Parameters::printDescription(
  ) const {

  //! Header
  const string header = "DESCRIPTION:";

  //! Content
  const string s1 = "This algorithm implements the SIFT method.";

  //! Print the description
  cout << header << endl << endl;

  //! Print the content
  this->printLine(s1, "  ");

  //! Finish the description
  cout << endl << endl;
}


//! Print the input parameters list.
void Parameters::printInput(
  ) const {

  //! Header
  const string header = "INPUT OPTIONS:";

  //! List of options
  string sentence;
  sentence += "-i [%s] : path to the input image\n";

  //! Print the header
  cout << header << endl << endl;

  //! Print the options
  this->printLine(sentence, "    ");

  //! Finish the list
  cout << endl << endl;
}


//! Print the output parameters list.
void Parameters::printOutput(
  ) const {

  //! Header
  const string header = "OUTPUT OPTIONS (results):";

  //! List of options
  string sentence;
  sentence += "-o [%s] : path to the output path\n";

  //! Print the header
  cout << header << endl << endl;
}


//! Print the optional parameters list.
void Parameters::printOptional(
  ) const {

  //! For convenience
  const string s4 = std::string(4, ' ');
  const string s7 = std::string(7, ' ');

  //! Print header
  this->printWord("OTHER OPTIONS:", "Default");
  cout << endl;

  //! Number of octaves
  this->printWord("-no (optional)", "8", s4);
  this->printLine("Number of octaves.", s7);
  cout << endl;

  //!
  this->printWord("-ns (optional)", "3", s4);
  this->printLine(".", s7);
  cout << endl;

  //! Number of histograms
  this->printWord("-nh (optional)", "4", s4);
  this->printLine("Number of histograms.", s7);
  cout << endl;

  //! Number of bins
  this->printWord("-nb (optional)", "36", s4);
  this->printLine("Number of bins.", s7);
  cout << endl;

  //! Number of orientations
  this->printWord("-nr (optional)", "8", s4);
  this->printLine("Number of orientations.", s7);
  cout << endl;

  //! Maximum number of iterations
  this->printWord("-m  (optional)", "5", s4);
  this->printLine("Maximum number of iterations.", s7);
  cout << endl;

  //! Sigma min
  this->printWord("-sm (optional)", "0.8", s4);
  this->printLine("Sigma min.", s7);
  cout << endl;

  //! Delta min
  this->printWord("-dm (optional)", "0.5", s4);
  this->printLine("Delta min.", s7);
  cout << endl;

  //! Sigma in
  this->printWord("-si (optional)", "0.5", s4);
  this->printLine("Sigma in.", s7);
  cout << endl;

  //! DoG
  this->printWord("-d  (optional)", "0.04/3.0", s4);
  this->printLine("DoG coefficient.", s7);
  cout << endl;

  //! Edge
  this->printWord("-e  (optional)", "10.0", s4);
  this->printLine("Edge.", s7);
  cout << endl;

  //! Lambda ori
  this->printWord("-lo (optional)", "1.5", s4);
  this->printLine("Lambda orientation.", s7);
  cout << endl;

  //! Lambda des
  this->printWord("-ld (optional)", "6.0", s4);
  this->printLine("Lambda descriptors.", s7);
  cout << endl;

  //! T
  this->printWord("-t  (optional)", "0.8", s4);
  this->printLine("T.", s7);
  cout << endl;

  //! Verbose
  this->printWord("-v  (optional)", "False", s4);
  this->printLine("Activate the verbose mode.", s7);
  cout << endl;

  //! Help
  this->printWord("-h  (optional)", "False", s4);
  this->printLine("Print the help.", s7);
  cout << endl << endl;
}


//! Print the signature.
void Parameters::printSignature(
  ) const {

  const string pad = std::string((m_size - 19) / 2, ' ');

  cout << pad << " ___ _          _ " << endl;
  cout << pad << "|_ _| |_   __ _| |" << endl;
  cout << pad << " | || '_ \\/ _` | |" << endl;
  cout << pad << "|___|_.__/\\__,_|_|" << endl;
  cout << endl;
}


//! Print a line of a specific m_size maximum number of characters
void Parameters::printLine(
  const std::string &i_sentence,
  const std::string &i_pad) const {

  //! Initializations
  istringstream iss(i_sentence);

  //! Get all the words contained in this sentence
  vector<string> line;
  do {
    string word;
    iss >> word;
    line.push_back(word);
  } while (iss);

  //! Print the line
  size_t nb = i_pad.size();
  cout << i_pad;
  for (size_t n = 0; n < line.size(); n++) {

    //! Get the number of characters contained in the current word
    const size_t nc = line[n].size();

    //! Print the current word
    if (nb + nc <= m_size) {
      cout << line[n];
      nb += nc;
    }
    else {
      nb = nc + i_pad.size();
      cout << endl;
      cout << i_pad << line[n];
    }

    //! Print the space character
    if (nb < m_size) {
      cout << " ";
      nb++;
    }
    else {
      cout << endl << i_pad;
      nb = i_pad.size();
    }
  }

  //! Finish the line
  cout << endl;
}


//! Print a line of a specific m_size maximum number of characters,
//! and add a word in the last position.
void Parameters::printWord(
  const std::string &i_line,
  const std::string &i_word,
  const std::string &i_pad) const {

  //! Check the number of characters
  const size_t nbL = i_line.size();
  const size_t nbW = i_word.size();
  const size_t nbP = i_pad.size();

  if (nbL + nbW + nbP < m_size) {
    const string s = std::string(m_size - nbL - nbW - nbP, ' ');
    cout << i_pad << i_line << s << i_word << endl;
  }
  else {
    const string s = std::string(m_size - nbP - nbW, ' ');
    cout << i_pad << i_line << endl;
    cout << i_pad << s << i_word << endl;
  }
}


//! Set the default values for the parameters
void Parameters::setDefaultValues() {

  m_nbOct     = 8;
  m_nbSpo     = 3;
  m_nbHist    = 4;
  m_nbBins    = 36;
  m_nbOri     = 8;
  m_maxIter   = 5;
  m_sigmaMin  = 0.8f;
  m_deltaMin  = 0.5f;
  m_sigmaIn   = 0.5f;
  m_dog       = 0.013333333f;
  m_edge      = 10.f;
  m_lambdaOri = 1.5f;
  m_lambdaDes = 6.f;
  m_t         = 0.8f;
}


//! Release the memory
void Parameters::releaseMemory() {
}


