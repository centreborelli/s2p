#include <stdio.h>
#include <time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

// portable way to get current time
void portable_gettime(struct timespec *ts)
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


// print an amount of time elapsed and start a new timer
void print_elapsed_time(struct timespec *ts, char *title, int title_width)
{
    // get the current time
    struct timespec finish;
    portable_gettime(&finish);

    // compute the elapsed time
    double elapsed = (finish.tv_sec - ts->tv_sec) * 1000;
    elapsed += (finish.tv_nsec - ts->tv_nsec) / 1000000.0;

    // print the result
    printf("%-*s %7.3f (ms)\n", title_width, title, elapsed); // see K&R p154

    // start a new timer
    portable_gettime(ts);
}
