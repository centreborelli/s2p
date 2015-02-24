/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
 name = {flst_pixels};
 version = {"1.0"};
 author = {"Pascal Monasse"};
 function = {"Fill the field PIXELS of a tree of shapes.
             Module meant to be called from another one,
             not from the command line"};
 usage = {
    tree->pTree  "The tree of the interpolated image"
};
*/

#include <stdlib.h>
#include "mw3.h"
//#include "mw3-modules.h"
#include "mw3-modules-specific.h"

/* For each shape, find its number of proper pixels */
static void compute_proper_pixels(Shapes pTree, int *tabNbOfProperPixels)
{
    Shape pShape;
    int *pNbOfProperPixels;
    int i;

    /* 1) Initialize by the area */
    pShape = pTree->the_shapes + pTree->nb_shapes - 1;
    pNbOfProperPixels = tabNbOfProperPixels + pTree->nb_shapes - 1;
    for (i = pTree->nb_shapes - 1; i >= 0; i--)
        *pNbOfProperPixels-- = (pShape--)->area;
    /* 2) For each shape, substract its area to its parent */
    pShape = pTree->the_shapes + pTree->nb_shapes - 1;
    for (i = pTree->nb_shapes - 1; i > 0; i--, pShape--)
        tabNbOfProperPixels[pShape->parent - pTree->the_shapes] -=
            pShape->area;
}

/* Allocate the array of pixels of each shape. Thanks to the tree structure,
we allocate only memory for the pixels of the root, and other arrays are
just pointers */
static void allocate_pixels(Shapes pTree, int *tabNbOfProperPixels)
{
    Point_plane tabPixelsOfRoot;
    Shape pShape, *tabpShapesOfStack;
    int i, iIndex, iSizeOfStack;

    /* 1) Memory allocation */
    tabPixelsOfRoot = pTree->the_shapes[0].pixels = (Point_plane)
        malloc(pTree->nrow * pTree->ncol * sizeof(struct point_plane));
    if (tabPixelsOfRoot == NULL)
        mwerror(FATAL, 1, "allocate_pixels --> Allocation of pixels\n");
    tabpShapesOfStack = (Shape *) malloc(pTree->nb_shapes * sizeof(Shape));
    if (tabpShapesOfStack == NULL)
        mwerror(FATAL, 1, "allocate_pixels --> Allocation of stack\n");

    /* 2) Enumeration of the tree in preorder, using the stack */
    pShape = &pTree->the_shapes[0];
    iSizeOfStack = 0;
    i = 0;
    while (1)
        if (pShape != NULL)
        {
            /* Write pixels of pShape */
            pShape->pixels = &tabPixelsOfRoot[i];
            iIndex = pShape - pTree->the_shapes;
            i += tabNbOfProperPixels[iIndex];
            tabpShapesOfStack[iSizeOfStack++] = pShape; /* Push */
            pShape = pShape->child;
        }
        else
        {
            if (iSizeOfStack == 0)
                break;
            pShape = tabpShapesOfStack[--iSizeOfStack]->next_sibling;
            /* Pop */
        }
    free(tabpShapesOfStack);
}

/* Associate to each shape its array of pixels. Fills the field PIXELS of
the tree structure. From the command line, this function has no interest,
since that field is not saved to the file. It is meant to be called from
another module, when this field is needed */
void flst_pixels(Shapes pTree)
{
    Shape *ppShape;
    int *tabNbOfProperPixels;   /* For each shape,
                                 * its number of proper pixels */
    Point_plane pCurrentPoint;
    int i, j, iIndex;

    /* 1) Compute nb of proper pixels in each shape */
    if ((tabNbOfProperPixels =
         (int *) malloc(pTree->nb_shapes * sizeof(int))) == NULL)
        mwerror(FATAL, 1, "Allocation of array error");
    compute_proper_pixels(pTree, tabNbOfProperPixels);

    /* 2) Allocate array for the root and make others sub-arrays */
    allocate_pixels(pTree, tabNbOfProperPixels);

    /* 3) Fill the array */
    ppShape = pTree->smallest_shape + pTree->ncol * pTree->nrow - 1;
    for (i = pTree->nrow - 1; i >= 0; i--)
        for (j = pTree->ncol - 1; j >= 0; j--, ppShape--)
        {
            iIndex = (*ppShape) - pTree->the_shapes;
            pCurrentPoint =
                &(*ppShape)->pixels[--tabNbOfProperPixels[iIndex]];
            pCurrentPoint->x = j;
            pCurrentPoint->y = i;
        }
    free(tabNbOfProperPixels);
}
