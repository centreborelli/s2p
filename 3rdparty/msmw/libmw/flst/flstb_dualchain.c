/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
 name = {flstb_dualchain};
 version = {"1.0"};
 author = {"Pascal Monasse"};
 function = {"Find chain of dual pixels adjacent to a shape
             in bilinear interpolated image"};
 usage = {
    tree -> pTree "The tree of shapes",
    shape -> pShape "The shape whose boundary is to be computed",
    boundary <- pBoundary "The chain of adjacent dual pixels (output flist)",
    notused -> ctabtabSaddleValues "Saddle values"
};
*/

#include <stdlib.h>
#include <float.h>
#include "mw3.h"
//#include "mw3-modules.h"         /* for flst_reconstruct(), fsaddles(),
//                                 * flst_pixels() */
#include "mw3-modules-specific.h"

#ifndef POINT_T
#define POINT_T
typedef struct {
    float x, y;
} point_t;
#endif

#define NOT_A_SADDLE FLT_MAX    /* Value must be coherent with fsaddles.c */
#define is_a_saddle(x) ((x) != NOT_A_SADDLE)

#define EAST  0
#define NORTH 1
#define WEST  2
#define SOUTH 3
typedef struct {
    short int x, y;             /* Top left corner coordinates */
    unsigned char cDirection;   /* Direction along which we enter
                                 * the dual pixel */
} DualPixel;

/* Is the center of pixel (x+.5,y+.5) in the shape? */
static char point_in_shape(short int x, short int y, Shape pShape,
                           Shapes pTree)
{
    Shape pShapePoint = pTree->smallest_shape[y * pTree->ncol + x];
    return (pShapePoint->area <= pShape->area &&
            pShape->pixels <= pShapePoint->pixels &&
            pShapePoint->pixels < pShape->pixels + pShape->area);
}

/* If there is a saddle point in the dual pixel of top left corner (x,y) and
this saddle point belongs to pShape, return 1 */
static char saddle_inside(short int x, short int y, Shape pShape,
                          float **tabtabSaddleValues)
{
    float v = tabtabSaddleValues[y - 1][x - 1];

    if (!is_a_saddle(v))
        return 0;
    if (pShape->inferior_type)
        return (char) (pShape->value >= v);
    return (char) (pShape->value <= v);
}

#define TURN_LEFT(dir) \
dir = (dir==NORTH ? WEST :\
      (dir==WEST ? SOUTH :\
      (dir==SOUTH ? EAST : NORTH)))
#define TURN_RIGHT(dir) \
dir = (dir==NORTH ? EAST :\
      (dir==EAST ? SOUTH :\
      (dir==SOUTH ? WEST : NORTH)))

/* Find the dual pixel following pDualPixel as we follow the shape boundary */
static void find_next_dual_pixel(DualPixel * pDualPixel, Shape pShape,
                                 Shapes pTree, float **tabtabSaddleValues)
{
    char bLeftIn, bRightIn;

    switch (pDualPixel->cDirection)
    {
    case NORTH:
        bLeftIn =
            point_in_shape(pDualPixel->x - 1, pDualPixel->y - 1, pShape,
                           pTree);
        bRightIn =
            point_in_shape(pDualPixel->x, pDualPixel->y - 1, pShape, pTree);
        if (bLeftIn
            || saddle_inside(pDualPixel->x, pDualPixel->y, pShape,
                             tabtabSaddleValues))
            if (bRightIn)
            {
                ++pDualPixel->x;
                TURN_RIGHT(pDualPixel->cDirection);
            }
            else
                --pDualPixel->y;
        else
        {
            --pDualPixel->x;
            TURN_LEFT(pDualPixel->cDirection);
        }
        break;
    case WEST:
        bLeftIn =
            point_in_shape(pDualPixel->x - 1, pDualPixel->y, pShape, pTree);
        bRightIn =
            point_in_shape(pDualPixel->x - 1, pDualPixel->y - 1, pShape,
                           pTree);
        if (bLeftIn
            || saddle_inside(pDualPixel->x, pDualPixel->y, pShape,
                             tabtabSaddleValues))
            if (bRightIn)
            {
                --pDualPixel->y;
                TURN_RIGHT(pDualPixel->cDirection);
            }
            else
                --pDualPixel->x;
        else
        {
            ++pDualPixel->y;
            TURN_LEFT(pDualPixel->cDirection);
        }
        break;
    case SOUTH:
        bLeftIn = point_in_shape(pDualPixel->x, pDualPixel->y, pShape, pTree);
        bRightIn =
            point_in_shape(pDualPixel->x - 1, pDualPixel->y, pShape, pTree);
        if (bLeftIn
            || saddle_inside(pDualPixel->x, pDualPixel->y, pShape,
                             tabtabSaddleValues))
            if (bRightIn)
            {
                --pDualPixel->x;
                TURN_RIGHT(pDualPixel->cDirection);
            }
            else
                ++pDualPixel->y;
        else
        {
            ++pDualPixel->x;
            TURN_LEFT(pDualPixel->cDirection);
        }
        break;
    case EAST:
        bLeftIn =
            point_in_shape(pDualPixel->x, pDualPixel->y - 1, pShape, pTree);
        bRightIn =
            point_in_shape(pDualPixel->x, pDualPixel->y, pShape, pTree);
        if (bLeftIn
            || saddle_inside(pDualPixel->x, pDualPixel->y, pShape,
                             tabtabSaddleValues))
            if (bRightIn)
            {
                ++pDualPixel->y;
                TURN_RIGHT(pDualPixel->cDirection);
            }
            else
                ++pDualPixel->x;
        else
        {
            --pDualPixel->y;
            TURN_LEFT(pDualPixel->cDirection);
        }
        break;
    }
}

/* Find the boundary of the shape, which is closed */
static int find_closed_boundary(Shapes pTree, Shape pShape,
                                float **tabtabSaddleValues, Flist pBoundary)
{
    short int x0, y0;
    DualPixel dualPixel;
    point_t *pPoint = (point_t *) pBoundary->values;

    /* 1) Find initial dual pixel, with NORTH direction */
    x0 = pShape->pixels[0].x;
    y0 = pShape->pixels[0].y;
    do
        ++x0;
    while (point_in_shape(x0, y0, pShape, pTree));

    /* 2) Follow the boundary */
    dualPixel.cDirection = NORTH;
    dualPixel.x = x0;
    dualPixel.y = y0;
    do
    {
        pPoint[pBoundary->size].x = (float) dualPixel.x;
        pPoint[pBoundary->size++].y = (float) dualPixel.y;
        find_next_dual_pixel(&dualPixel, pShape, pTree, tabtabSaddleValues);
    }
    while (dualPixel.x != x0 || dualPixel.y != y0 ||
           dualPixel.cDirection != NORTH);
    return 0;
}

/* Find an initial point (to follow the boundary) at the border of the image */
static void initial_point_border(DualPixel * pDualPixel, Shape pShape,
                                 Shapes pTree)
{
    short int iWidth = (short int) pTree->ncol, iHeight =
        (short int) pTree->nrow;
    short int x, y;

    /* Right border */
    pDualPixel->cDirection = WEST;
    x = iWidth - 1;
    y = 0;
    if (point_in_shape(x, y++, pShape, pTree))
        while (y < iHeight && point_in_shape(x, y++, pShape, pTree));
    while (y < iHeight && !point_in_shape(x, y, pShape, pTree))
        ++y;
    if (y < iHeight)
    {
        pDualPixel->x = iWidth;
        pDualPixel->y = y;
        return;
    }
    /* Top border */
    pDualPixel->cDirection = SOUTH;
    x = 0;
    y = 0;
    if (point_in_shape(x++, y, pShape, pTree))
        while (x < iWidth && point_in_shape(x++, y, pShape, pTree));
    while (x < iWidth && !point_in_shape(x, y, pShape, pTree))
        ++x;
    if (x < iWidth)
    {
        pDualPixel->x = x;
        pDualPixel->y = 0;
        return;
    }
    /* Left border */
    pDualPixel->cDirection = EAST;
    x = 0;
    y = iHeight - 1;
    if (point_in_shape(x, y--, pShape, pTree))
        while (y >= 0 && point_in_shape(x, y--, pShape, pTree));
    while (y >= 0 && !point_in_shape(x, y, pShape, pTree))
        --y;
    if (y >= 0)
    {
        pDualPixel->x = 0;
        pDualPixel->y = y + 1;
        return;
    }
    /* Bottom border */
    pDualPixel->cDirection = NORTH;
    x = iWidth - 1;
    y = iHeight - 1;
    if (point_in_shape(x--, y, pShape, pTree))
        while (x >= 0 && point_in_shape(x--, y, pShape, pTree));
    while (x >= 0 && !point_in_shape(x, y, pShape, pTree))
        --x;
    if (x >= 0)
    {
        pDualPixel->x = x + 1;
        pDualPixel->y = iHeight;
        return;
    }
    mwerror(FATAL, 1, "initial_point_border --> Coherency?");
}

/* Find an open boundary */
static int find_open_boundary(Shapes pTree, Shape pShape,
                              float **tabtabSaddleValues, Flist pBoundary)
{
    DualPixel dualPixel;
    short int iWidth = (short int) pTree->ncol, iHeight =
        (short int) pTree->nrow;
    point_t *pPoint = (point_t *) pBoundary->values;

    initial_point_border(&dualPixel, pShape, pTree);
    do
    {
        pPoint[pBoundary->size].x = (float) dualPixel.x;
        pPoint[pBoundary->size++].y = (float) dualPixel.y;
        find_next_dual_pixel(&dualPixel, pShape, pTree, tabtabSaddleValues);
    }
    while (0 < dualPixel.x && dualPixel.x < iWidth &&
           0 < dualPixel.y && dualPixel.y < iHeight);
    pPoint[pBoundary->size].x = (float) dualPixel.x;    /* We store the exit */
    pPoint[pBoundary->size++].y = (float) dualPixel.y;
    return 0;
}

/* Find boundary of the shape */
void flstb_dualchain(Shapes pTree, Shape pShape, Flist pBoundary,
                     char *ctabtabSaddleValues)
{
    struct fimage imageSaddles, *pImage;
    char bBuildSaddleValues = 0;
    float **tabtabSaddleValues = (float **) ctabtabSaddleValues;

    if (pTree->interpolation != 1)
        mwerror(USAGE, 1, "Please apply to a *bilinear* tree");
    if (tabtabSaddleValues == NULL)
    {
        bBuildSaddleValues = (char) 1;
        imageSaddles.nrow = imageSaddles.ncol = 0;
        imageSaddles.gray = NULL;
        imageSaddles.allocsize = 0;
        if ((pImage = mw_new_fimage()) == NULL)
            mwerror(FATAL, 1, "Void image allocation error!");
        flst_reconstruct(pTree, pImage);
        fsaddles(pImage, &imageSaddles);
        mw_delete_fimage(pImage);
        if ((tabtabSaddleValues =
             mw_newtab_gray_fimage(&imageSaddles)) == NULL)
            mwerror(FATAL, 1, "Lines of saddle values: allocation error");
    }
    if (pTree->the_shapes[0].pixels == NULL)
        flst_pixels(pTree);

    if (mw_change_flist(pBoundary, pShape->area * 4, 0, 2) != pBoundary)
        mwerror(FATAL, 1, "List allocation error");

    if (pShape->open)
        find_open_boundary(pTree, pShape, tabtabSaddleValues, pBoundary);
    else
        find_closed_boundary(pTree, pShape, tabtabSaddleValues, pBoundary);
    mw_realloc_flist(pBoundary, pBoundary->size);

    if (bBuildSaddleValues)
    {
        free(tabtabSaddleValues[0]);
        free(tabtabSaddleValues);
    }
}
