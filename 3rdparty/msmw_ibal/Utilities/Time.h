#ifndef TIME_H_INCLUDED
#define TIME_H_INCLUDED


//! Global includes
#include <ctime>


//! Local includes


/**
 * @brief Small class to handle the elapsed time. For now, only works on linux
 *        system (problem of portability of CLOCK_MONOTONIC).
 **/

class Time {

  //! Methods
  public:
    /**
     * @brief Default constructor.
     **/
    Time();


    /**
     * @brief Copy constructor.
     **/
    Time(
      const Time& i_time);


    /**
     * @brief Operator overload.
     **/
    Time& operator=(
      const Time& i_time);


    //! Default destructor
    ~Time();


    /**
     * @brief Compute an amount of time elapsed.
     *
     * @param io_start : when the clock started; Will be
     *        reinitialized just after;
     * @param p_name : name of the section to print;
     * @param p_nbChar: number of characters to print.
     *
     * @return none.
     **/
    void getTime(
      const char* p_name,
      const size_t p_nbChar = 35);


  private:


    /**
     * @brief Release the memory.
     **/
    void releaseMemory();


  //! Data members
  private:

    struct timespec* m_time;
};
#else
class Time;

#endif // TIME_H_INCLUDED
