/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fsignal.c

  Vers. 1.3
  Author : Jacques Froment
  Basic memory routines for the fsignal internal type

  Main changes :
  v1.3 (JF): added include <string> (Linux 2.6.12 & gcc 4.0.2)

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~  This file is part of the MegaWave2 system library ~~~~~~~~~~~~~~~
  MegaWave2 is a "soft-publication" for the scientific community. It has
  been developed for research purposes and it comes without any warranty.
  The last version is available at http://www.cmla.ens-cachan.fr/Cmla/Megawave
  CMLA, Ecole Normale Superieure de Cachan, 61 av. du President Wilson,
  94235 Cachan cedex, France. Email: megawave@cmla.ens-cachan.fr
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <stdlib.h>
#include <string.h>

#include "definitions.h"
#include "error.h"

#include "fsignal.h"

/* creates a new fsignal structure */

Fsignal mw_new_fsignal()
{
    Fsignal signal;

    if (!(signal = (Fsignal) (malloc(sizeof(struct fsignal)))))
    {
        mwerror(ERROR, 0, "[mw_new_fsignal] Not enough memory\n");
        return (NULL);
    }

    signal->size = 0;
    signal->allocsize = 0;
    signal->firstp = signal->lastp = 0;
    signal->param = 0.0;

    signal->scale = 1.0;
    signal->shift = 0.0;
    signal->gain = 1.0;
    signal->sgrate = 1.0;
    signal->bpsample = 8 * sizeof(float);

    strcpy(signal->cmt, "?");
    strcpy(signal->name, "?");

    signal->values = NULL;

    return (signal);
}

/* allocates the values array for N samples */

Fsignal mw_alloc_fsignal(Fsignal signal, int N)
{
    int mem_size;

    if (signal == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_alloc_fsignal] cannot alloc array : "
                "fsignal structure is NULL\n");
        return (NULL);
    }

    mem_size = N * sizeof(float);
    if (mem_size <= 0)
    {
        mwerror(ERROR, 0,
                "[mw_alloc_fsignal] Attempts to alloc a fsignal "
                "with null size\n");
        return NULL;
    }

    if (signal->values != NULL)
    {
        mwerror(ERROR, 0,
                "[mw_alloc_fsignal] Attempts to alloc a fsignal "
                "which is already allocated\n");
        return (NULL);
    }

    signal->values = (float *) malloc(mem_size);
    if (signal->values == NULL)
    {
        signal->size = 0;
        signal->allocsize = 0;
        mwerror(ERROR, 0, "[mw_alloc_fsignal] Not enough memory\n");
        return (NULL);
    }

    signal->size = N;
    signal->allocsize = mem_size;
    return (signal);
}

/* desallocate the array in the fsignal structure and the structure itself */

void mw_delete_fsignal(Fsignal signal)
{
    if (signal == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_delete_fsignal] cannot delete : "
                "fsignal structure is NULL\n");
        return;
    }

    if (signal->values != NULL)
        free(signal->values);
    signal->values = NULL;
    free(signal);
    signal = NULL;
}

/* Change the size of the allocated array */
/* May define the struct if not defined */
/* So you have to call it with signal = mw_change_fsignal(signal,...) */

Fsignal mw_change_fsignal(Fsignal signal, int N)
{
    int mem_size;

    if (signal == NULL)
        signal = mw_new_fsignal();
    if (signal == NULL)
        return (NULL);

    mem_size = N * sizeof(float);
    if (mem_size > signal->allocsize)
    {
        if (signal->values != NULL)
        {
            free(signal->values);
            signal->values = NULL;
        }
        if (mw_alloc_fsignal(signal, N) == NULL)
        {
            mw_delete_fsignal(signal);
            return (NULL);
        }
    }
    else
        signal->size = N;

    return (signal);
}

/* Clear the array of a fsignal with the value v */

void mw_clear_fsignal(Fsignal signal, float v)
{
    register float *ptr;
    register int l;

    if ((!signal) || (!signal->values))
    {
        mwerror(ERROR, 0,
                "[mw_clear_fsignal] NULL signal struct or NULL values array\n");
        return;
    }

    for (l = 1, ptr = signal->values; l <= signal->size; l++, ptr++)
        *ptr = v;
}

/* Copy the values array of a fsignal into another fsignal */

void mw_copy_fsignal_values(Fsignal in, Fsignal out)
{
    if ((!in) || (!out) || (!in->values) || (!out->values)
        || (in->size != out->size))
    {
        mwerror(ERROR, 0,
                "[mw_copy_fsignal_values] NULL input or output "
                "signal or signals of different sizes !\n");
        return;
    }

    memcpy(out->values, in->values, sizeof(float) * in->size);
}

/* Copy the header of a fsignal into another fsignal */

void mw_copy_fsignal_header(Fsignal in, Fsignal out)
{
    if ((!in) || (!out))
    {
        mwerror(ERROR, 0,
                "[mw_copy_fsignal_header] NULL input or output signal !\n");
        return;
    }
    out->firstp = in->firstp;
    out->lastp = in->lastp;
    out->param = in->param;
    out->scale = in->scale;
    out->shift = in->shift;
    out->gain = in->gain;
    out->sgrate = in->sgrate;
    out->bpsample = in->bpsample;
    strcpy(out->cmt, in->cmt);
}

/* Copy a fsignal into another fsignal */

void mw_copy_fsignal(Fsignal in, Fsignal out)
{
    if ((!in) || (!out) || (!in->values) || (!out->values)
        || (in->size != out->size))
    {
        mwerror(ERROR, 0,
                "[mw_copy_fsignal] NULL input or output signal "
                "or signals of different sizes !\n");
        return;
    }

    memcpy(out->values, in->values, sizeof(float) * in->size);
    mw_copy_fsignal_header(in, out);
}
