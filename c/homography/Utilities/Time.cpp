/**
 * @file Time.cpp
 *
 * @brief Class to handle the elapsed time. For now, only works on linux system.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


#include <iostream>

#include "Time.h"


using namespace std;


//! Default constructor
Time::Time() :

  m_time((struct timespec) {0, 0}) {

  //! Initialize the time
#ifdef __linux__
  clock_gettime(CLOCK_MONOTONIC, &m_time);
#endif
}



//! Compute an amount of time elapsed.
void Time::getTime(
  const char* p_name,
  const size_t p_nbChar) {
#ifdef __linux__
  //! Check the current time
  struct timespec finish;
  clock_gettime(CLOCK_MONOTONIC, &finish);

  //! Compute the elapsed time
  double elapsed = (finish.tv_sec - m_time.tv_sec) * 1000;
  elapsed += (finish.tv_nsec - m_time.tv_nsec) / 1000000.0;

  //! Determine the space to add at the end
  string sentence = p_name;
  sentence += string(max(0, int(p_nbChar) - int(sentence.size())), ' ');

  //! Print the result
  cout << sentence << ": (ms) = " << elapsed << endl;

  //! Start a new timer
  clock_gettime(CLOCK_MONOTONIC, &m_time);
#endif
}
