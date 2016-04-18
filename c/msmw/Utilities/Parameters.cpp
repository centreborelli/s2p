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

  //! Parameters set by users
  m_inpLeft    (NULL),   // Path of the input left image        -il [%s]
  m_inpRight   (NULL),   // Path of the input right image       -ir [%s]
  m_outDispL   (NULL),   // Path of the output left disparity   -dl [%s]
  m_outDispR   (NULL),   // Path of the output right disparity  -dr [%s]
  m_outMaskL   (NULL),   // Path of the output left mask        -kl [%s]
  m_outMaskR   (NULL),   // Path of the output right mask       -kr [%s]
  m_minDisp    (0),      // minimum displacement                -m  [%f]
  m_maxDisp    (0),      // maximum displacement                -M  [%f]
  m_orientation(5),      // Indicates nb of orientations to use -W  [%d]
  m_nbScales   (4),      // Number of scales                    -n  [%d]
  m_x          (9),      // Window x size                       -x  [%d]
  m_y          (9),      // Window y size                       -y  [%d]
  m_dist       (0),      // Distance to use                     -p  [%d]
  m_verbose    (false),  // Activate verbose mode               -v

  //! Hard-coded parameters
  m_window(NULL),
  m_grainArea(25),
  m_valueRecip(1.f),
  m_valueRemoveIsolated(0.25f),
  m_valueMinDist(1.f),
  m_dmin1(m_minDisp),
  m_dmax1(m_maxDisp),
  m_dmin2(-m_maxDisp),
  m_dmax2(-m_minDisp) {

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
  m_inpLeft    (i_params.inpLeft()),
  m_inpRight   (i_params.inpRight()),
  m_outDispL   (i_params.outDispL()),
  m_outDispR   (i_params.outDispR()),
  m_outMaskL   (i_params.outMaskL()),
  m_outMaskR   (i_params.outMaskR()),
  m_minDisp    (i_params.minDisp()),
  m_maxDisp    (i_params.maxDisp()),
  m_orientation(i_params.orientation()),
  m_nbScales   (i_params.nbScales()),
  m_x          (i_params.x()),
  m_y          (i_params.y()),
  m_dist       (i_params.dist()),
  m_verbose    (i_params.verbose()),

  //! Hard-coded parameters
  m_window(new Image(*i_params.m_window)),
  m_grainArea(i_params.m_grainArea),
  m_valueRecip(i_params.m_valueRecip),
  m_valueRemoveIsolated(i_params.m_valueRemoveIsolated),
  m_valueMinDist(i_params.m_valueMinDist),
  m_dmin1(i_params.m_dmin1),
  m_dmax1(i_params.m_dmax1),
  m_dmin2(i_params.m_dmin2),
  m_dmax2(i_params.m_dmax2) {

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

      return EXIT_FAILURE;
    }

    //! Input path
    if (sarg.find("-il") == 0) {
      if (n + 1 < i_argc) {
        m_inpLeft = i_argv[n + 1];
      }
    }

    if (sarg.find("-ir") == 0) {
      if (n + 1 < i_argc) {
        m_inpRight = i_argv[n + 1];
      }
    }

    //! Output path
    if (sarg.find("-dl") == 0) {
      if (n + 1 < i_argc) {
        m_outDispL = i_argv[n + 1];
      }
    }

    if (sarg.find("-dr") == 0) {
      if (n + 1 < i_argc) {
        m_outDispR = i_argv[n + 1];
      }
    }

    if (sarg.find("-kl") == 0) {
      if (n + 1 < i_argc) {
        m_outMaskL = i_argv[n + 1];
      }
    }

    if (sarg.find("-kr") == 0) {
      if (n + 1 < i_argc) {
        m_outMaskR = i_argv[n + 1];
      }
    }

    //! Minimum disparity range
    if (sarg.find("-m") == 0) {
      if (n + 1 < i_argc) {
        m_minDisp = atof(i_argv[n + 1]);
      }
    }

    //! Maximum disparity range
    if (sarg.find("-M") == 0) {
      if (n + 1 < i_argc) {
        m_maxDisp = atof(i_argv[n + 1]);
      }
    }

    //! Orientation
    if (sarg.find("-W") == 0) {
      if (n + 1 < i_argc) {
        m_orientation = atoi(i_argv[n + 1]);
      }
    }

    //! Number of scales
    if (sarg.find("-n") == 0) {
      if (n + 1 < i_argc) {
        m_nbScales = atoi(i_argv[n + 1]);
      }
    }

    //! Window X
    if (sarg.find("-x") == 0) {
      if (n + 1 < i_argc) {
        m_x = atoi(i_argv[n + 1]);
      }
    }

    //! Window Y
    if (sarg.find("-y") == 0) {
      if (n + 1 < i_argc) {
        m_y = atoi(i_argv[n + 1]);
      }
    }

    //! Patches distance to use
    if (sarg.find("-p") == 0) {
      if (n + 1 < i_argc) {
        m_dist = atoi(i_argv[n + 1]);
      }
    }

    //! Verbose option
    if (sarg.find("-v") == 0) {
      m_verbose = true;
    }
  }

  //! Check that the mandatory arguments have been given
  if (m_inpLeft == NULL) {
    cout << "The input left image path must be precised. Use \n";
    cout << "    -il path/imageLeft.ext" << endl;
    return EXIT_FAILURE;
  }
  if (m_inpRight == NULL) {
    cout << "The input right image path must be precised. Use \n";
    cout << "    -ir path/imageRight.ext" << endl;
    return EXIT_FAILURE;
  }
  if (m_outDispL == NULL) {
    cout << "The output left disparity path must be precised. Use \n";
    cout << "    -dl path/disparityLeft.txt" << endl;
    return EXIT_FAILURE;
  }
  if (m_outMaskL == NULL) {
    cout << "The output left mask path must be precised. Use \n";
    cout << "    -kl path/maskLeft.txt" << endl;
    return EXIT_FAILURE;
  }

  //! Automatic parameters
  m_dmin1 =  m_minDisp;
  m_dmax1 =  m_maxDisp;
  m_dmin2 = -m_maxDisp;
  m_dmax2 = -m_minDisp;

  if (m_window != NULL) {
    m_window->init(m_x, m_y);
  }
  else {
    m_window = new Image(m_x, m_y);
  }
  *m_window = 1.f;
  m_window->normalizeL1();

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
  const string name = "./msmw";

  //! Parameters
  string sentence = name;
  sentence += " [-il input left image path]";
  sentence += " [-ir input right image path]";
  sentence += " [-dl output left disparity path]";
  sentence += " [-dr output right disparity path]";
  sentence += " [-kl output left mask path]";
  sentence += " [-kr output right mask path]";
  sentence += " [-m  minimum disparity]";
  sentence += " [-M  maximum disparity]";
  sentence += " [-W  orientation]";
  sentence += " [-n  number of scales]";
  sentence += " [-x  window X]";
  sentence += " [-y  window Y]";
  sentence += " [-p  distance]";
  sentence += " [-v  verbose]";

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
  const string s1 = "This algorithm implements the MSMW method.";

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
  sentence += "-il [%s] : path to the input left image\n";

  //! Print the header
  cout << header << endl << endl;

  //! Print the options
  this->printLine(sentence, "    ");
  sentence = "-ir [%s] : path to the input right image\n";
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
  sentence += "-dl [%s] : path to the output left disparity\n";

  //! Print the header
  cout << header << endl << endl;

  //! Print the options
  this->printLine(sentence, "    ");
  sentence = "-kl [%s] : path to the output left mask\n";
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

  //! Output right disparity
  this->printWord("-dr (optional)", "NULL", s4);
  this->printLine("Path to the output right disparity.", s7);

  //! Output right mask
  this->printWord("-kr (optional)", "NULL", s4);
  this->printLine("Path to the output right mask.", s7);

  //! Minimum disparity
  this->printWord("-m  (optional)", "0", s4);
  this->printLine("Minimum value of the disparity range.", s7);

  //! Maximum disparity
  this->printWord("-M  (optional)", "0", s4);
  this->printLine("Maximum value of the disparity range.", s7);

  //! Orientation
  this->printWord("-W  (optional)", "5", s4);
  this->printLine("Number of orientations to use.", s7);

  //! Number of scales
  this->printWord("-n  (optional)", "4", s4);
  this->printLine("Number of scales to use.", s7);

  //! Window X size
  this->printWord("-x  (optional)", "9", s4);
  this->printLine("Window X size.", s7);

  //! Window Y size
  this->printWord("-y  (optional)", "9", s4);
  this->printLine("Window Y size.", s7);

  //! Patch distance
  this->printWord("-p  (optional)", "0", s4);
  this->printLine("If 0, then the mean of the patches will be removed before\
    applying the L2 distance, if 1 the classic L2 distance will be used.", s7);

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


//! Update the parameters according to the current scale.
void Parameters::update(
  const size_t p_currentScale) {

  //! Scale
  const float scale = 1 << (p_currentScale - 1);//powf(2.0, (float) p_currentScale - 1);

  //! Update the range
  m_dmin1 /= scale;
  m_dmax1 /= scale;
  m_dmin2 /= scale;
  m_dmax2 /= scale;

  //! Reduce the size of the grain at lower scales
  m_grainArea /= scale;
}


//! Release the memory
void Parameters::releaseMemory() {

  if (m_window != NULL) {
    delete m_window;
    m_window = NULL;
  }
}


