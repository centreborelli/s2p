/**
 * @file Time.cpp
 *
 * @brief Class to handle the elapsed time.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com> (original linux version)
 * @author Carlo de Franchis <carlodef@gmail.com> (osx compatible version)
 **/


#include <iostream>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif


#include "Time.h"


using namespace std;

//! Portable way to get current time
static void portable_gettime(struct timespec *ts)
{
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, ts);
#endif
}


//! Default constructor
Time::Time() :

        m_time((struct timespec) {0, 0}) {

    //! Initialize the time
    portable_gettime(&m_time);
}




//! Compute an amount of time elapsed.
void Time::get_time(
        const char *p_name,
        const size_t p_nbChar) {
    //! Check the current time
    struct timespec finish;
    portable_gettime(&finish);

    //! Compute the elapsed time
    double elapsed = (finish.tv_sec - m_time.tv_sec) * 1000;
    elapsed += (finish.tv_nsec - m_time.tv_nsec) / 1000000.0;

    //! Determine the space to add at the end
    string sentence = p_name;
    sentence += string(max(0, int(p_nbChar) - int(sentence.size())), ' ');

    //! Print the result
    cout << sentence << ": (ms) = " << elapsed << endl;

    //! Start a new timer
    portable_gettime(&m_time);
}
