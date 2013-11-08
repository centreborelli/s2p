/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  shape.c

  Vers. 1.2
  Authors : Pascal Monasse, Frederic Guichard, Jacques Froment.
  Basic memory routines for the Point_plane, Shape and Shapes structures

  Main changes :
  v1.2 (JF): added include <string> (Linux 2.6.12 & gcc 4.0.2)
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

#include "shape.h"

/*--- Point_plane ---*/

/* creates a new point_plane structure */

Point_plane mw_new_point_plane(void)
{
    Point_plane point;

    if (!(point = (Point_plane) (malloc(sizeof(struct point_plane)))))
    {
        mwerror(ERROR, 0, "[mw_new_point_plane] Not enough memory\n");
        return (NULL);
    }
    point->x = point->y = -1;
    return (point);
}

/* Define the struct if it's not defined */

Point_plane mw_change_point_plane(Point_plane point)
{
    if (point == NULL)
        point = mw_new_point_plane();
    return (point);
}

/*--- Shape ---*/

/* Creates a new shape structure */

Shape mw_new_shape(void)
{
    Shape shape;

    if (!(shape = (Shape) (malloc(sizeof(struct shape)))))
    {
        mwerror(ERROR, 0, "[mw_new_shape] Not enough memory\n");
        return (NULL);
    }

    shape->inferior_type = 0;
    shape->value = 0.0;
    shape->open = 0;
    shape->area = 0;
    shape->removed = 0;
    shape->pixels = NULL;
    shape->boundary = NULL;
    shape->parent = shape->next_sibling = shape->child = NULL;
    shape->data = NULL;
    shape->data_size = 0;
    return (shape);
}

/* Define the struct if it's not defined */

Shape mw_change_shape(Shape sh)
{
    if (sh == NULL)
        sh = mw_new_shape();
    return (sh);
}

/* deallocate the shape structure */

void mw_delete_shape(Shape shape)
{
    if (shape == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_delete_shape] cannot delete : shape structure is NULL\n");
        return;
    }
    if (shape->pixels)
        free(shape->pixels);
    shape->pixels = NULL;
    shape->area = 0;
    if (shape->boundary)
        mw_delete_flist(shape->boundary);
    shape->boundary = NULL;
    if (shape->data)
        free(shape->data);
    free(shape);
    shape = NULL;
}

/* Get in the subtree of root sh a shape that is not removed, NULL if
 * all shapes are removed */

Shape mw_get_not_removed_shape(Shape sh)
{
    Shape NotSup = NULL;
    if ((sh == NULL) || (!sh->removed))
        return (sh);
    for (sh = sh->child; sh != NULL; sh = sh->next_sibling)
        if ((NotSup = mw_get_not_removed_shape(sh)) != NULL)
            break;
    return (NotSup);
}

/* Get the true parent, that is the greatest non removed ancestor */

Shape mw_get_parent_shape(Shape sh)
{
    if (sh == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_get_parent_shape] input shape structure is NULL\n");
        return (NULL);
    }
    if (sh->parent == NULL)     /* It is the root of the shape */
        return (NULL);

    do
        if ((sh = sh->parent) == NULL)
        {
            mwerror(ERROR, 0,
                    "[mw_get_parent_shape] the root of the shapes "
                    "is removed\n");
            return (NULL);
        }
    while (sh->removed);
    return (sh);
}

/* Get the first child, taking into account that some shapes are removed */

Shape mw_get_first_child_shape(Shape sh)
{
    Shape NotSup = NULL;
    if (sh == NULL)
    {
        mwerror(ERROR, 0, "[mw_get_first_child_shape] input shape is NULL\n");
        return (NULL);
    }
    for (sh = sh->child; sh != NULL; sh = sh->next_sibling)
        if ((NotSup = mw_get_not_removed_shape(sh)) != NULL)
            break;
    return (NotSup);
}

/* Get the next sibling, taking into account that some shapes are removed */

Shape mw_get_next_sibling_shape(Shape sh)
{
    Shape sh1 = NULL, sh2 = NULL;
    if (sh == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_get_next_sibling_shape] input shape is NULL\n");
        return (NULL);
    }
    /* First look at the siblings in the original tree */
    for (sh1 = sh->next_sibling; sh1 != NULL; sh1 = sh1->next_sibling)
        if ((sh2 = mw_get_not_removed_shape(sh1)) != NULL)
            return (sh2);
    if (sh->parent == NULL || !sh->parent->removed)
        return (NULL);
    /* The parent in the original tree is also the parent in the true
     * tree, nothing more to do */
    /* If not found, find the node in the original tree just before
     * the true parent */
    do
    {
        sh = sh->parent;
        /* Look at the siblings of this node */
        for (sh1 = sh->next_sibling; sh1 != NULL; sh1 = sh1->next_sibling)
            if ((sh2 = mw_get_not_removed_shape(sh1)) != NULL)
                return (sh2);
    }
    while (sh->parent->removed);
    return (NULL);
}

/* Return the smallest shape containing the given pixel */

Shape mw_get_smallest_shape(Shapes shs, int iX, int iY)
{
    Shape sh;

    if (shs == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_get_smallest_shape] input shapes structure is NULL\n");
        return (NULL);
    }
    sh = shs->smallest_shape[iY * shs->ncol + iX];
    if (sh == NULL)
    {
        mwerror(ERROR, 0, "[mw_get_smallest_shape] smallest shape is NULL\n");
        return (NULL);
    }
    if (sh->removed)
        sh = mw_get_parent_shape(sh);
    return (sh);
}

/*-- Shapes --*/

/* creates a new shapes structure */

Shapes mw_new_shapes(void)
{
    Shapes shapes;

    if (!(shapes = (Shapes) (malloc(sizeof(struct shapes)))))
    {
        mwerror(ERROR, 0, "[mw_new_shapes] Not enough memory\n");
        return (NULL);
    }
    shapes->the_shapes = NULL;
    shapes->smallest_shape = NULL;
    shapes->data = NULL;
    shapes->data_size = 0;
    shapes->nb_shapes = 0;
    shapes->nrow = shapes->ncol = 0;
    shapes->interpolation = 0;
    strcpy(shapes->cmt, "?");
    strcpy(shapes->name, "?");
    return (shapes);
}

/* Allocate the root of the shapes */

Shapes mw_alloc_shapes(Shapes shs, int nrow, int ncol, float value)
{
    int size, i;
    Shape root;

    if (shs == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_alloc_shapes] cannot alloc root : "
                "shapes structure is NULL\n");
        return (NULL);
    }

    size = nrow * ncol;
    if (size <= 0)
    {
        mwerror(ERROR, 0,
                "[mw_alloc_shapes] Attempts to alloc shapes with null size\n");
        return (NULL);
    }
    if ((shs->the_shapes != NULL) || (shs->smallest_shape != NULL))
    {
        mwerror(ERROR, 0,
                "[mw_alloc_shapes] Attempts to alloc root "
                "which is already allocated\n");
        return (NULL);
    }

    root = shs->the_shapes =
        (Shape) malloc((size + 1) * sizeof(struct shape));
    if (root == NULL)
    {
        mwerror(ERROR, 0, "[mw_alloc_shapes] Not enough memory\n");
        return (NULL);
    }
    root->inferior_type = 1;
    root->value = value;
    root->open = 1;
    root->area = size;
    root->removed = 0;
    root->pixels = NULL;
    root->boundary = NULL;
    root->parent = root->next_sibling = root->child = NULL;
    root->data = NULL;
    root->data_size = 0;

    shs->nb_shapes = 1;
    shs->nrow = nrow;
    shs->ncol = ncol;

    shs->smallest_shape = (Shape *) malloc(size * sizeof(Shape));
    if (shs->smallest_shape == NULL)
    {
        mwerror(ERROR, 0, "[mw_alloc_shapes] Not enough memory\n");
        free(root);
        return (NULL);
    }
    for (i = size - 1; i >= 0; i--)
        shs->smallest_shape[i] = root;

    return (shs);
}

/* Allocate the struct and the root if they are not defined */

Shapes mw_change_shapes(Shapes shs, int nrow, int ncol, float value)
{
    int i;

    if (shs == NULL)
        shs = mw_new_shapes();
    if (shs == NULL)
        return (NULL);

    /* Deallocate the shapes but the structure itself */
    for (i = shs->nb_shapes - 1; i > 0; i--)
        if (shs->the_shapes[i].boundary != NULL)
            mw_delete_flist(shs->the_shapes[i].boundary);

    if ((shs->the_shapes != NULL) && (shs->nb_shapes > 0))
        free(shs->the_shapes[0].pixels);
    if (shs->the_shapes != NULL)
        free(shs->the_shapes);
    if (shs->smallest_shape != NULL)
        free(shs->smallest_shape);

    shs = mw_alloc_shapes(shs, nrow, ncol, value);
    return (shs);
}

/* deallocate the shapes structure */

void mw_delete_shapes(Shapes shs)
{
    int i;

    if (shs == NULL)
    {
        mwerror(ERROR, 0,
                "[mw_delete_shapes] cannot delete : "
                "shapes structure is NULL\n");
        return;
    }

    for (i = shs->nb_shapes - 1; i > 0; i--)
        if (shs->the_shapes[i].boundary != NULL)
            mw_delete_flist(shs->the_shapes[i].boundary);

    if ((shs->the_shapes != NULL) && (shs->nb_shapes > 0))
        free(shs->the_shapes[0].pixels);
    if (shs->the_shapes != NULL)
        free(shs->the_shapes);
    if (shs->smallest_shape != NULL)
        free(shs->smallest_shape);
    free(shs);
}
