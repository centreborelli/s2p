/*
 * Copyright (c) 2011-2014, Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

/**
 * @file timing.c
 * @brief timing and profiling tools
 *
 * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 */

#define _POSIX_C_SOURCE 200112L /* ask for POSIX 2001 */

#ifndef USE_TIMING
#define USE_TIMING
#endif
#include "timing.h"

/*
 * OS DETECTION
 */

#if (defined(_WIN32) || defined(__WIN32__) \
     || defined(__TOS_WIN__) || defined(__WINDOWS__))
/* from http://sourceforge.net/p/predef/wiki/OperatingSystems/ */
#define TIMING_WINDOWS
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#elif (defined(__unix__) || defined(__unix))
/* from http://sourceforge.net/p/predef/wiki/Standards/ */
#include <unistd.h>
#if (defined(_POSIX_VERSION) && (_POSIX_VERSION >= 200112L))
#define TIMING_POSIX
#include <sys/time.h>
#endif                          /* POSIX test */
#endif                          /* Windows/Unix test */

/*
 * ARCHITECTURE DETECTION
 */

#if (defined(__amd64__) || defined(__amd64) || defined(_M_X64))
/* from http://sourceforge.net/p/predef/wiki/Architectures/ */
#define TIMING_AMD64

#elif (defined(__i386__) || defined(__i386) || defined(_M_IX86) \
       || defined(__X86__) || defined(_X86_) || defined(__I86__))
/* from http://predef.sourceforge.net/prearch.html#sec6 */
#define TIMING_I386
#endif


/****************************************************************
 * WALL CLOCK TIMER
 ****************************************************************/

/** number of wall clock counters */
#ifndef TIMING_WALLCLOCK_NB
#define TIMING_WALLCLOCK_NB 16
#endif

/** wall clock counter array, initialized to 0 */
unsigned long _timing_wallclock_counter[TIMING_WALLCLOCK_NB];

/**
 * @brief portable wallclock call
 */
#if defined(TIMING_POSIX)
/* use POSIX gettomeofday() with microsecond precision */
unsigned long _timing_wallclock() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (unsigned long)(t.tv_usec + t.tv_sec * 1000000);
}
#elif defined(TIMING_WINDOWS)
/* use Windows GetSystemTime() with millisecond precision */
unsigned long _timing_wallclock() {
    static SYSTEMTIME t;
    GetSystemTime(&t);
#define _UL unsigned long       /* temporary, for shorter lines */
    return (_UL)(1000 * ((_UL) t.wMilliseconds
                         + 1000 * ((_UL) t.wSecond
                                   + 60 * ((_UL) t.wMinute
                                           + 60 * ((_UL) t.wHour
                                                   + 24 * (_UL) t.wDay)))));
#undef _UL
}
#else
/* fall back to libc time() with second precision */
#warning wallclock uses a low-precision libc implementation
unsigned long _timing_wallclock() {
    time_t rawtime;
    struct tm *t;
    (void) time(&rawtime);
    t = localtime(&rawtime);
#define _UL unsigned long       /* temporary, for shorter lines */
    return (_UL)(1000000 * ((_UL) t->tm_sec
                            + 60 * ((_UL) t->tm_min +
                                    +60 * ((_UL) t->tm_hour +
                                           +24 * (_UL) t->tm_mday))));
#undef _UL
}
#endif

/****************************************************************
 * CPU CLOCK TIMER
 ****************************************************************/

/** number of CPU clock counters */
#ifndef TIMING_CPUCLOCK_NB
#define TIMING_CPUCLOCK_NB 16
#endif

/** CPU clock counter array, initialized to 0 (K&R2, p.86) */
unsigned long _timing_cpuclock_counter[TIMING_CPUCLOCK_NB];

/**
 * @brief portable cpuclock call
 */
#if defined(TIMING_POSIX)
/* use POSIX clock_gettime (reduced to microsecond precision) */
unsigned long _timing_cpuclock() {
    struct timespec tmp;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tmp);
    return (unsigned long)(tmp.tv_sec * 1000000 + tmp.tv_nsec / 1000);
}
#else
unsigned long _timing_cpuclock()
/* fall back to libc clock() (1/100s to 1/10000s precision) */
{
    return (unsigned long)(clock() * 1000000 / CLOCKS_PER_SEC);
}
#endif

/****************************************************************
 * CPU CYCLES COUNTER
 ****************************************************************/

/**
 * Counter functions use the work of Daniel J. Bernstein:
 *   http://ebats.cr.yp.to/cpucycles.html
 *   http://ebats.cr.yp.to/cpucycles-20060326.tar.gz
 */

#if (defined(__STDC__) && defined(__STDC_VERSION__) \
     && (__STDC_VERSION__ >= 199409L))
#define _LL long long
#else
#define _LL long
#endif

/** number of cycle counters */
#ifndef TIMING_CYCLE_NB
#define TIMING_CYCLE_NB 16
#endif

/** cycle counter array, initialized to 0 */
_LL _timing_cycle_counter[TIMING_CYCLE_NB];

/**
 * @brief portable cycles counter
 */
_LL _timing_cpucycles(void);
#if defined(TIMING_AMD64)
/** CPU cycles counter for amd64 */
_LL _timing_cpucycles() {
    unsigned _LL result;
    __asm__ volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax":"=a"
                     (result)::"%rdx");
    return result;
}

#elif defined(TIMING_I386)
/** CPU cycles counter for x86 */
_LL _timing_cpucycles() {
    _LL result;
    __asm__ volatile(".byte 15;.byte 49":"=A"(result));
    return result;
}

#else
/** dummy CPU cycles counter */
_LL _timing_cpucycles() {
    return 0;
}

#endif                          /* architecture detection */
