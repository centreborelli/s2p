#ifndef TIME_H_INCLUDED
#define TIME_H_INCLUDED


#include <ctime>


/**
 * @brief Small class to handle the elapsed time. For now, only works on linux
 *        system (problem of portability of CLOCK_MONOTONIC).
 **/

class Time {

  public:
    /**
     * @brief Default constructor.
     **/
    Time();


    /**
     * @brief Compute an amount of time elapsed.
     *
     * @param p_name : name of the section to print;
     * @param p_nbChar: number of characters to print.
     *
     * @return none.
     **/
    void get_time(
            const char *p_name,
            const size_t p_nbChar = 35);


  private:

    struct timespec m_time;
};
#endif // TIME_H_INCLUDED
