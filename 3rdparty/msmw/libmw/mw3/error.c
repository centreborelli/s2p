/**
 * @file error.c
 *
 * @version 1.14
 * @author Jacques Froment & Sylvain Parrino (1995-2005)
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>

#include "definitions.h"

#include "error.h"

/* TODO: move to config.h or definitions.h */
int mwerrcnt = 0;

void mwdebug(char *fmt, ...)
{
    if (mwdbg)
    {
        va_list marker;

        va_start(marker, fmt);
        fprintf(stderr, "<dbg> ");
        vfprintf(stderr, fmt, marker);
        va_end(marker);
    }
}

void mwerror(int code, int exit_code, char *fmt, ...)
{
    va_list marker;

    va_start(marker, fmt);

    switch (code)
    {
    case WARNING:
        fprintf(stderr, "megawave warning (%s) : ", mwname);
        vfprintf(stderr, fmt, marker);
        break;
    case ERROR:
        fprintf(stderr, "megawave error (%s) : ", mwname);
        vfprintf(stderr, fmt, marker);
        mwerrcnt++;
        break;
    case FATAL:
        fprintf(stderr, "megawave fatal (%s) : ", mwname);
        vfprintf(stderr, fmt, marker);
        fprintf(stderr, "Exit.\n");
        exit(exit_code);
        break;
    case INTERNAL:
        fprintf(stderr, "megawave internal (%s) : ", mwname);
        vfprintf(stderr, fmt, marker);
        fprintf(stderr, "Exit.\n");
        exit(exit_code);
        break;
    case USAGE:
        fprintf(stderr, "Bad parameter; use the '--help' option for details");
        exit(exit_code);
        break;
    default:
        mwerror(FATAL, 1, "Bad usage of mwerror : code %d is unknown\n",
                code);
        break;
    }
    va_end(marker);
}
