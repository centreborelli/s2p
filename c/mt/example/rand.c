/*
 * This file shows how to use mt.c. It is released in the public
 * domain and as such comes with no copyright requirement.
 *
 * compile with: cc rand.c mt.c
 */

#include <stdlib.h>
#include <stdio.h>

/* include the mt prototypes */
#include "mt.h"

int main(int argc, char *const *argv)
{
    int i;

    if (1 < argc)
        /* if a parameter is given, it is used for the initialisation */
        mt_init((unsigned long int) atol(argv[1]));
    else
        /*
         * without parameter, use the automatic initialization
         */
        mt_init_auto();

    /*
     * now we can get 1024 pseudo-random value
     */
    for (i = 0; i < 1024; i++)
        printf("%18.16f\n", mt_drand53());

    return EXIT_SUCCESS;
}
