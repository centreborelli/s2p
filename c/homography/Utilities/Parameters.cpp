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
  m_size         (43 + 2 * 15),

  //! Parameters
  m_inpName      (NULL),       // Path of the input image            -i [%s]
  m_inpHomo      (NULL),       // Path of the input homography       -t [%s]
  m_outName      (NULL),       // Path of the ouput image            -o [%s]
  m_adjustSize   (false),      // Activate the automatic size        -a
  m_oWidth       (0),          // Wanted width of te output image    -c [%d]
  m_oHeight      (0),          // Wanted height of the output image  -l [%d]
  m_verbose      (false)       // Activate verbose mode              -v
  {

}


//! Default destructor
Parameters::~Parameters(
  ) {
}


//! Copy constructor.
Parameters::Parameters(
  const Parameters& i_params) :
  //! Miscalleneous
  m_size(i_params.m_size),

  //! Parameters
  m_inpName      (i_params.inpName()),
  m_inpHomo      (i_params.inpHomo()),
  m_outName      (i_params.outName()),
  m_adjustSize   (i_params.adjustSize()),
  m_oWidth       (i_params.oWidth()),
  m_oHeight      (i_params.oHeight()),
  m_verbose      (i_params.verbose()) {

}


//! Read the input arguments.
int Parameters::checkArgs(
  const int i_argc,
  char** i_argv) {

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

      return EXIT_FAILURE;
    }

    //! Input path
    if (sarg.find("-i") == 0) {
      if (n + 1 < i_argc) {
        m_inpName = i_argv[n + 1];
      }
    }

    if (sarg.find("-t") == 0) {
      if (n + 1 < i_argc) {
        m_inpHomo = i_argv[n + 1];
      }
    }

    //! Output path
    if (sarg.find("-o") == 0) {
      if (n + 1 < i_argc) {
        m_outName = i_argv[n + 1];
      }
    }

    //! Ajust size
    if (sarg.find("-a") == 0) {
      m_adjustSize = true;
    }

    //! Output width
    if (sarg.find("-c") == 0) {
      if (n + 1 < i_argc) {
        m_oWidth = atoi(i_argv[n + 1]);
      }
    }

    //! Output height
    if (sarg.find("-l") == 0) {
      if (n + 1 < i_argc) {
        m_oHeight = atoi(i_argv[n + 1]);
      }
    }

    //! Verbose option
    if (sarg.find("-v") == 0) {
      m_verbose = true;
    }
  }

  //! Check that the mandatory arguments have been given
  if (m_inpName == NULL) {
    cout << "The input image path must be precised. Use \n";
    cout << "    -i path/image.ext" << endl;
    return EXIT_FAILURE;
  }
  if (m_inpHomo == NULL) {
    cout << "The input homography path must be precised. Use \n";
    cout << "    -t path/Homography.txt" << endl;
    return EXIT_FAILURE;
  }
  if (m_outName == NULL) {
    cout << "The output image path must be precised. Use \n";
    cout << "    -o path/image.ext" << endl;
    return EXIT_FAILURE;
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
  const string name = "./Homography";

  //! Parameters
  string sentence = name;
  sentence += " [-i input image path]";
  sentence += " [-t input homography path]";
  sentence += " [-o output image path]";
  sentence += " [-a adjust size]";
  sentence += " [-c output width]";
  sentence += " [-l output height]";
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
  const string s1 = "This algorithm implements a landscape evolution model.";

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
  sentence = "-t [%s] : path to the input homography";
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
  sentence += "-o [%s] : path to the output image";

  //! Print the header
  cout << header << endl << endl;

  //! Print the options
  this->printLine(sentence, "    ");

  //! Finish the list
  cout << endl << endl;
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

  //! Adjust size
  this->printWord("-a (optional)", "False", s4);
  this->printLine("Activate the automatic size adjustment.", s7);

  //! Output width
  this->printWord("-c (optional)", "0", s4);
  this->printLine("If > 0, will set the width to the output image.", s7);

  //! Output height
  this->printWord("-l (optional)", "0", s4);
  this->printLine("If > 0, will set the height to the output image.", s7);

  //! Verbose
  this->printWord("-v (optional)", "False", s4);
  this->printLine("Activate the verbose mode.", s7);
  cout << endl;

  //! Help
  this->printWord("-h (optional)", "False", s4);
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




