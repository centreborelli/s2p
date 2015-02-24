/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  list.c

  Vers. 1.5
  Author : Lionel Moisan
  Basic routines for the Flist(s) / Dlist(s) internal types

  Main changes :
  v??? (JF): not enough memory errors handling and alloc_size->max_size.
  v??? (LM): data fields are not ignored any more in copy,...
  v1.5 (JF): added include <string> (Linux 2.6.12 & gcc 4.0.2)
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

#include "list.h"

/* Memory growth rate for reallocations */
#define MW_LIST_ENLARGE_FACTOR 1.5

/******************************** FLIST ********************************/

/* Create a new flist structure */

Flist mw_new_flist(void)
{
    Flist l;

    if (!(l = (Flist) (malloc(sizeof(struct flist)))))
    {
        mwerror(ERROR, 0, "[mw_new_flist] Not enough memory\n");
        return (NULL);
    }
    l->size = l->max_size = l->dim = l->data_size = 0;
    l->values = NULL;
    l->data = NULL;

    return (l);
}

/* (Re)allocate the values[] array in a flist structure */

Flist mw_realloc_flist(Flist l, int n)
{
    float *ptr;

    if (!l)
    {
        mwerror(ERROR, 0, "[mw_realloc_flist] flist structure is NULL\n");
        return (NULL);
    }

    if (n < l->size)
        mwerror(WARNING, 0, "[mw_realloc_flist] data may have been lost\n");
    ptr = (float *) realloc(l->values, n * l->dim * sizeof(float));
    if (!ptr && n != 0)
    {
        mwerror(ERROR, 0, "[mw_realloc_flist] Not enough memory\n");
        if (l->values)
            free(l->values);
        free(l);
        return (NULL);
    }
    l->values = ptr;
    l->max_size = n;
    return (l);
}

/* Enlarge the values[] array in a flist structure */

Flist mw_enlarge_flist(Flist l)
{
    int n;

    if (!l)
    {
        mwerror(ERROR, 0, "[mw_enlarge_flist] flist structure is NULL\n");
        return (NULL);
    }

    n = (int) (1 + (float) l->max_size * MW_LIST_ENLARGE_FACTOR);
    return (mw_realloc_flist(l, n));
}

/* Configure and (re)allocate a flist structure */

Flist mw_change_flist(Flist l, int max_size, int size, int dimension)
{
    if (!l)
    {
        l = mw_new_flist();
        if (!l)
            return (NULL);
    }

    l->dim = dimension;
    l->size = size;
    if (size < 0 || size > max_size)
    {
        mwerror(ERROR, 0, "[mw_change_flist] Inconsistant size\n");
        return (l);
    }

    /* (re)allocate values[] */
    if (max_size > l->max_size)
        return (mw_realloc_flist(l, max_size));

    return (l);
}

/* Delete a flist structure */

void mw_delete_flist(Flist l)
{
    if (!l)
    {
        mwerror(ERROR, 0, "[mw_delete_flist] flist structure is NULL\n");
        return;
    }
    if (l->values)
        free(l->values);
    if (l->data_size > 0 && l->data)
        free(l->data);
    free(l);
}

/* Clear the values[] array in a flist structure */

void mw_clear_flist(Flist l, float v)
{
    int n;

    if (!l)
    {
        mwerror(ERROR, 0, "[mw_clear_flist] flist structure is NULL\n");
        return;
    }

    for (n = l->size * l->dim; n--;)
        l->values[n] = v;
}

/* Duplicate a flist structure */

Flist mw_copy_flist(Flist in, Flist out)
{
    int n;

    if (!in)
    {
        mwerror(ERROR, 0, "[mw_copy_flist] NULL input\n");
        return (NULL);
    }
    out = mw_change_flist(out, in->size, in->size, in->dim);
    if (out)
    {
        for (n = in->size * in->dim; n--;)
            out->values[n] = in->values[n];
        if (in->data_size && in->data)
        {
            if (out->data_size < in->data_size)
            {
                out->data = realloc(out->data, in->data_size);
                if (!out->data)
                {
                    mwerror(ERROR, 0, "[mw_copy_flist] Not enough memory\n");
                    return (out);
                }
            }
            memcpy(out->data, in->data, in->data_size);
            out->data_size = in->data_size;
        }
    }
    return (out);
}

/******************************** FLISTS ********************************/

/* Create a new flists structure */

Flists mw_new_flists(void)
{
    Flists ls;

    if (!(ls = (Flists) (malloc(sizeof(struct flists)))))
    {
        mwerror(ERROR, 0, "[mw_new_flists] Not enough memory\n");
        return (NULL);
    }
    ls->size = ls->max_size = ls->data_size = 0;
    ls->list = NULL;
    ls->data = NULL;

    return (ls);
}

/* (Re)allocate the list[] array in a flists structure */

Flists mw_realloc_flists(Flists ls, int n)
{
    Flist *ptr;
    int i;

    if (!ls)
    {
        mwerror(ERROR, 0, "[mw_realloc_flists] flists structure is NULL\n");
        return (NULL);
    }

    if (n < ls->size)
        mwerror(WARNING, 0, "[mw_realloc_flists] data may have been lost\n");
    ptr = (Flist *) realloc(ls->list, n * sizeof(Flist));
    if (!ptr && n != 0)
    {
        mwerror(ERROR, 0, "[mw_realloc_flists] Not enough memory\n");
        free(ls);
        return (NULL);
    }

    /* clear new allocated Flists */
    ls->list = ptr;
    for (i = ls->max_size; i < n; i++)
        ls->list[i] = NULL;
    ls->max_size = n;
    return (ls);
}

/* Enlarge the list[] array in a flists structure */

Flists mw_enlarge_flists(Flists ls)
{
    int n;

    if (!ls)
    {
        mwerror(ERROR, 0, "[mw_enlarge_flists] flists structure is NULL\n");
        return (NULL);
    }

    n = (int) (1 + (float) ls->max_size * MW_LIST_ENLARGE_FACTOR);
    return (mw_realloc_flists(ls, n));
}

/* Configure and (re)allocate a flists structure */

Flists mw_change_flists(Flists ls, int max_size, int size)
{
    if (!ls)
    {
        ls = mw_new_flists();
        if (!ls)
            return (NULL);
    }

    ls->size = size;
    if (size < 0 || size > max_size)
    {
        mwerror(ERROR, 0, "[mw_change_flists] Inconsistant size\n");
        return (ls);
    }

    /* (re)allocate list[] */
    if (max_size > ls->max_size)
        return (mw_realloc_flists(ls, max_size));

    return (ls);
}

/* Delete a flists structure (delete all associated Flist) */

void mw_delete_flists(Flists ls)
{
    int i;

    if (!ls)
    {
        mwerror(ERROR, 0, "[mw_delete_flists] flists structure is NULL\n");
        return;
    }
    for (i = ls->size; i--;)
        mw_delete_flist(ls->list[i]);
    if (ls->list)
        free(ls->list);
    if (ls->data_size > 0 && ls->data)
        free(ls->data);
    free(ls);
}

/* Duplicate a flists structure */

Flists mw_copy_flists(Flists in, Flists out)
{
    int n;

    if (!in)
    {
        mwerror(ERROR, 0, "[mw_copy_flists] NULL input\n");
        return (NULL);
    }
    out = mw_change_flists(out, in->size, in->size);
    if (out)
    {
        for (n = in->size; n--;)
            out->list[n] = mw_copy_flist(in->list[n], out->list[n]);
        if (in->data_size && in->data)
        {
            if (out->data_size < in->data_size)
            {
                out->data = realloc(out->data, in->data_size);
                if (!out->data)
                {
                    mwerror(ERROR, 0, "[mw_copy_flists] Not enough memory\n");
                    return (out);
                }
            }
            memcpy(out->data, in->data, in->data_size);
            out->data_size = in->data_size;
        }
    }
    return (out);
}

/******************************** DLIST ********************************/

/* Create a new dlist structure */

Dlist mw_new_dlist(void)
{
    Dlist l;

    if (!(l = (Dlist) (malloc(sizeof(struct dlist)))))
    {
        mwerror(ERROR, 0, "[mw_new_dlist] Not enough memory\n");
        return (NULL);
    }
    l->size = l->max_size = l->dim = l->data_size = 0;
    l->values = NULL;
    l->data = NULL;

    return (l);
}

/* (Re)allocate the values[] array in a dlist structure */

Dlist mw_realloc_dlist(Dlist l, int n)
{
    double *ptr;

    if (!l)
    {
        mwerror(ERROR, 0, "[mw_realloc_dlist] dlist structure is NULL\n");
        return (NULL);
    }

    if (n < l->size)
        mwerror(WARNING, 0, "[mw_realloc_dlist] data may have been lost\n");

    ptr = (double *) realloc(l->values, n * l->dim * sizeof(double));
    if (!ptr && n != 0)
    {
        mwerror(ERROR, 0, "[mw_realloc_dlist] Not enough memory\n");
        if (l->values)
            free(l->values);
        free(l);
        return (NULL);
    }
    l->values = ptr;
    l->max_size = n;
    return (l);
}

/* Enlarge the values[] array in a dlist structure */

Dlist mw_enlarge_dlist(Dlist l)
{
    int n;

    if (!l)
    {
        mwerror(ERROR, 0, "[mw_enlarge_dlist] dlist structure is NULL\n");
        return (NULL);
    }

    n = (int) (1 + (double) l->max_size * MW_LIST_ENLARGE_FACTOR);
    return (mw_realloc_dlist(l, n));
}

/* Configure and (re)allocate a dlist structure */

Dlist mw_change_dlist(Dlist l, int max_size, int size, int dimension)
{
    if (!l)
    {
        l = mw_new_dlist();
        if (!l)
            return (NULL);
    }

    l->dim = dimension;
    l->size = size;
    if (size < 0 || size > max_size)
    {
        mwerror(ERROR, 0, "[mw_change_dlist] Inconsistant size\n");
        return (l);
    }

    /* (re)allocate values[] */
    if (max_size > l->max_size)
        return (mw_realloc_dlist(l, max_size));

    return (l);
}

/* Delete a dlist structure */

void mw_delete_dlist(Dlist l)
{
    if (!l)
    {
        mwerror(ERROR, 0, "[mw_delete_dlist] dlist structure is NULL\n");
        return;
    }
    if (l->values)
        free(l->values);
    if (l->data_size > 0 && l->data)
        free(l->data);
    free(l);
}

/* Clear the values[] array in a dlist structure */

void mw_clear_dlist(Dlist l, double v)
{
    int n;

    if (!l)
    {
        mwerror(ERROR, 0, "[mw_clear_dlist] dlist structure is NULL\n");
        return;
    }

    for (n = l->size * l->dim; n--;)
        l->values[n] = v;
}

/* Duplicate a dlist structure */

Dlist mw_copy_dlist(Dlist in, Dlist out)
{
    int n;

    if (!in)
    {
        mwerror(ERROR, 0, "[mw_copy_dlist] NULL input\n");
        return (NULL);
    }
    out = mw_change_dlist(out, in->size, in->size, in->dim);
    if (out)
    {
        for (n = in->size * in->dim; n--;)
            out->values[n] = in->values[n];
        if (in->data_size && in->data)
        {
            if (out->data_size < in->data_size)
            {
                out->data = realloc(out->data, in->data_size);
                if (!out->data)
                {
                    mwerror(ERROR, 0, "[mw_copy_dlist] Not enough memory\n");
                    return (out);
                }
            }
            memcpy(out->data, in->data, in->data_size);
            out->data_size = in->data_size;
        }
    }
    return (out);
}

/******************************** DLISTS ********************************/

/* Create a new dlists structure */

Dlists mw_new_dlists(void)
{
    Dlists ls;

    if (!(ls = (Dlists) (malloc(sizeof(struct dlists)))))
    {
        mwerror(ERROR, 0, "[mw_new_dlists] Not enough memory\n");
        return (NULL);
    }
    ls->size = ls->max_size = ls->data_size = 0;
    ls->list = NULL;
    ls->data = NULL;

    return (ls);
}

/* (Re)allocate the list[] array in a dlists structure */

Dlists mw_realloc_dlists(Dlists ls, int n)
{
    Dlist *ptr;
    int i;

    if (!ls)
    {
        mwerror(ERROR, 0, "[mw_realloc_dlists] dlists structure is NULL\n");
        return (NULL);
    }

    if (n < ls->size)
        mwerror(WARNING, 0, "[mw_realloc_dlists] data may have been lost\n");

    ptr = (Dlist *) realloc(ls->list, n * sizeof(Dlist));
    if (!ptr && n != 0)
    {
        mwerror(ERROR, 0, "[mw_realloc_dlists] Not enough memory\n");
        free(ls);
        return (NULL);
    }

    /* clear new allocated Dlists */
    ls->list = ptr;
    for (i = ls->max_size; i < n; i++)
        ls->list[i] = NULL;
    ls->max_size = n;
    return (ls);
}

/* Enlarge the list[] array in a dlists structure */

Dlists mw_enlarge_dlists(Dlists ls)
{
    int n;

    if (!ls)
    {
        mwerror(ERROR, 0, "[mw_enlarge_dlists] dlists structure is NULL\n");
        return NULL;
    }

    n = (int) (1 + (double) ls->max_size * MW_LIST_ENLARGE_FACTOR);
    return (mw_realloc_dlists(ls, n));
}

/* Configure and (re)allocate a dlists structure */

Dlists mw_change_dlists(Dlists ls, int max_size, int size)
{
    if (!ls)
    {
        ls = mw_new_dlists();
        if (!ls)
            return (NULL);
    }

    ls->size = size;
    if (size < 0 || size > max_size)
    {
        mwerror(ERROR, 0, "[mw_change_dlists] Inconsistant size\n");
        return (ls);
    }

    /* (re)allocate list[] */
    if (max_size > ls->max_size)
        return (mw_realloc_dlists(ls, max_size));

    return (ls);
}

/* Delete a dlists structure (delete all associated Dlist) */

void mw_delete_dlists(Dlists ls)
{
    int i;

    if (!ls)
    {
        mwerror(ERROR, 0, "[mw_delete_dlists] dlists structure is NULL\n");
        return;
    }
    for (i = ls->size; i--;)
        mw_delete_dlist(ls->list[i]);
    if (ls->list)
        free(ls->list);
    if (ls->data_size > 0 && ls->data)
        free(ls->data);
    free(ls);
}

/* Duplicate a dlists structure */

Dlists mw_copy_dlists(Dlists in, Dlists out)
{
    int n;

    if (!in)
    {
        mwerror(ERROR, 0, "[mw_copy_dlists] NULL input\n");
        return (NULL);
    }
    out = mw_change_dlists(out, in->size, in->size);
    if (out)
    {
        for (n = in->size; n--;)
            out->list[n] = mw_copy_dlist(in->list[n], out->list[n]);
        if (in->data_size && in->data)
        {
            if (out->data_size < in->data_size)
            {
                out->data = realloc(out->data, in->data_size);
                if (!out->data)
                {
                    mwerror(ERROR, 0, "[mw_copy_dlists] Not enough memory\n");
                    return (out);
                }
            }
            memcpy(out->data, in->data, in->data_size);
            out->data_size = in->data_size;
        }
    }
    return (out);
}
