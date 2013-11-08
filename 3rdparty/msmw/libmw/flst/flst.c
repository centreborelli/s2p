/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
 name = {flst};
 version = {"2.0"};
 author = {"Pascal Monasse"};
 function = {"Fast Level Sets Transform of an image"};
 usage = {
    'a':min_area->pMinArea  "argument of the grain filter",
    image->pImageInput      "Input fimage",
    tree<-pTree             "The tree of shapes"
};
*/

#include <stdlib.h>
#include <assert.h>
#include "mw3.h"
//#include "mw3-modules.h"
#include "mw3-modules-specific.h"

/* Optimization parameters. Probably should depend on image size, but these
values seem good for most images. */
#define INIT_MAX_AREA 10
#define STEP_MAX_AREA 8

/* A 'local configuration' of the pixel is the logical 'or' of these values,
stored in one byte. Corresponding bit is set if the neighbor is in region */
#define EAST         1
#define NORTH_EAST   2
#define NORTH        4
#define NORTH_WEST   8
#define WEST        16
#define SOUTH_WEST  32
#define SOUTH       64
#define SOUTH_EAST 128

/* Gives for each configuration of the neighborhood around the pixel the number
of new cc's of the complement created (sometimes negative, since
cc's can be deleted), when the pixel is appended to the region.
Configurations are stored on one byte.
tabPatterns[0]: set in 4-connectedness and complement in 8-connectedness.
tabPatterns[1]: set in 8-connectedness and complement in 4-connectedness. */
static char tabPatterns[2][256] =
    { {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
       0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
       0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
       0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
       0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 1,
       1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1,
       0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 1,
       0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
       0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
       0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
       0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
       0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
       0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0,
       1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0,
       0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0,
       0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1,
 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0,
 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 0, -1, 0, -1,
 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0,
 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 0, -1, 0, -1,
 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0,
 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0,
 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 2, 2, 1, 2, 1,
 1, 1, 2, 1, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0,
 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0,
 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 0, -1, 0, -1,
 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 0,
 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 0, -1, 0, -1}
};

/* The structure representing the neighborhood of a region. It is a priority
queue organized as a heap (a binary tree, where each node has a key larger,
resp. smaller, than its two children), stored in an array.
Index 1 is the root, its children are at index 2 and 3.
Therefore, the parent of node at k is at k/2 and its children at 2*k and 2*k+1.
Index 0 is only a temporary buffer to swap two nodes. */
enum TypeOfTree { AMBIGUOUS, MIN, MAX, INVALID };
typedef struct {
    struct point_plane point;   /* Neighbor pixel */
    float value;                /* Its gray level */
} Neighbor;
typedef struct {
    Neighbor *tabPoints;        /* The array of neighbors,
                                 * organized as a binary tree */
    int iNbPoints;              /* The size of the previous arrays */
    enum TypeOfTree type;       /* max- or min- oriented heap? */
    float otherBound;           /* Min gray level if max-oriented,
                                 * max if min-oriented */
} Neighborhood;

/* Structure to find connections between shapes. This is used when a monotone
section is extracted. The goal is to find the parent of its largest shape. The
gray level of the parent is known, so as a point of the parent, since we use
in fact an image of connections */
typedef struct {
    Shape shape;                /* Largest shape of the monotone section */
    float level;                /* Connection level */
} Connection;

/* Well, here are global variables. Ugly, but avoids to weigh the code,
since they are used almost everywhere */
static int iWidth, iHeight;
static int iMinArea, iMaxArea, iAreaImage, iHalfAreaImage, iPerimeterImage;
static int iExploration;        /* Used to avoid reinitializing images */
static Point_plane tabPointsInShape;
static int **tabtabVisitedNeighbors;    /* Exterior boundary */
static Shapes pGlobalTree;
static int iAtBorder;           /* #points meeting at border of image */

#define ARRAY_2D_ALLOC(array, iWidth, iHeight, type) \
(array) = (type**) malloc((iHeight)*sizeof(type*)); \
if((array) == NULL) mwerror(FATAL, 1, "ARRAY_2D_ALLOC --> memory"); \
(array)[0] = (type*) calloc((iWidth)*(iHeight), sizeof(type)); \
if((array)[0] == NULL) mwerror(FATAL, 1, "ARRAY_2D_ALLOC --> memory"); \
for(i = 1; i < (iHeight); i++) \
  (array)[i] = (array)[i-1]+(iWidth);

/* ------------------------------------------------------------------------
   --------- Functions to manage the neighborhood structure ---------------
   ------------------------------------------------------------------------ */

/* Reinitialise the neighborhood, so that it will be used for a new region */
static void reinit_neighborhood(Neighborhood * pNeighborhood,
                                enum TypeOfTree type)
{
    pNeighborhood->iNbPoints = 0;
    pNeighborhood->type = type;
}

/* To allocate the structure representing the neighborhood of a region */
static void init_neighborhood(Neighborhood * pNeighborhood, int init_iMaxArea)
{
    init_iMaxArea = 4 * (init_iMaxArea + 1);
    if (init_iMaxArea > iWidth * iHeight)
        init_iMaxArea = iWidth * iHeight;

    pNeighborhood->tabPoints =
        (Neighbor *) malloc((init_iMaxArea + 1) * sizeof(Neighbor));
    if (pNeighborhood->tabPoints == NULL)
        mwerror(FATAL, 1,
                "init_neighborhood --> neighbors allocation error\n");
    reinit_neighborhood(pNeighborhood, AMBIGUOUS);
}

/* Free the structure representing the neighborhood of a region */
static void free_neighborhood(Neighborhood * pNeighborhood)
{
    free(pNeighborhood->tabPoints);
}

#define ORDER_MAX(k,l) (tabPoints[k].value > tabPoints[l].value)
#define ORDER_MIN(k,l) (tabPoints[k].value < tabPoints[l].value)
#define SWAP(k,l) tabPoints[0] = tabPoints[k]; \
                  tabPoints[k] = tabPoints[l]; \
                  tabPoints[l] = tabPoints[0];
/* Put the last neighbor at a position so that we fix the heap */
static void fix_up(Neighborhood * pNeighborhood)
{
    Neighbor *tabPoints = pNeighborhood->tabPoints;
    int k = pNeighborhood->iNbPoints, l;

    if (pNeighborhood->type == MAX)
        while (k > 1 && ORDER_MAX(k, l = k >> 1))
        {
            SWAP(k, l);
            k = l;
        }
    else
        while (k > 1 && ORDER_MIN(k, l = k >> 1))
        {
            SWAP(k, l);
            k = l;
        }
}

#define ORDER_MAX2(k,l) (tabPoints[k].value >= tabPoints[l].value)
#define ORDER_MIN2(k,l) (tabPoints[k].value <= tabPoints[l].value)
/* Put the first neighbor at a position so that we fix the heap */
static void fix_down(Neighborhood * pNeighborhood)
{
    Neighbor *tabPoints = pNeighborhood->tabPoints;
    int N = pNeighborhood->iNbPoints, k = 1, l;

    if (pNeighborhood->type == MAX)
        while ((l = k << 1) <= N)
        {
            if (l < N && ORDER_MAX(l + 1, l))
                ++l;
            if (ORDER_MAX2(k, l))
                break;
            SWAP(k, l);
            k = l;
        }
    else
        while ((l = k << 1) <= N)
        {
            if (l < N && ORDER_MIN(l + 1, l))
                ++l;
            if (ORDER_MIN2(k, l))
                break;
            SWAP(k, l);
            k = l;
        }
}

/* Add the pixel (x,y), of gray-level VALUE, to the neighbor pixels */
static void add_neighbor(Neighborhood * pNeighborhood, short int x,
                         short int y, float value)
{
    Neighbor *pNewNeighbor;

    /* 0) Tag the pixel as 'visited neighbor' */
    tabtabVisitedNeighbors[y][x] = iExploration;
    /* 1) Update the fields TYPE and OTHERBOUND of PNEIGHBORHOOD */
    if (pNeighborhood->iNbPoints == 0)
        pNeighborhood->otherBound = value;
    else
        switch (pNeighborhood->type)
        {
        case MIN:
            if (value > pNeighborhood->otherBound)
                pNeighborhood->otherBound = value;
            else if (value < pNeighborhood->tabPoints[1].value)
                pNeighborhood->type = INVALID;
            break;
        case MAX:
            if (value < pNeighborhood->otherBound)
                pNeighborhood->otherBound = value;
            else if (value > pNeighborhood->tabPoints[1].value)
                pNeighborhood->type = INVALID;
            break;
        case AMBIGUOUS:
            if (value != pNeighborhood->tabPoints[1].value)
            {
                pNeighborhood->type =
                    (value < pNeighborhood->tabPoints[1].value) ? MAX : MIN;
                pNeighborhood->otherBound = value;
            }
            break;
        case INVALID:
            return;
        }
    if (pNeighborhood->type == INVALID)
        return;
    /* 2) Add the point in the heap and update it */
    pNewNeighbor = &pNeighborhood->tabPoints[++pNeighborhood->iNbPoints];
    pNewNeighbor->point.x = x;  /* Initialise the new neighbor point */
    pNewNeighbor->point.y = y;
    pNewNeighbor->value = value;
    fix_up(pNeighborhood);      /* Update the heap of neighbors */
}

/* Remove the neighbor at the top of the heap, that is the root of the tree. */
static void remove_neighbor(Neighborhood * pNeighborhood)
{
    Neighbor *pTop = &pNeighborhood->tabPoints[1];
    float value = pTop->value;

    if (pNeighborhood->type == INVALID)
        return;
    *pTop = pNeighborhood->tabPoints[pNeighborhood->iNbPoints--];
    if (pNeighborhood->iNbPoints == 0)
        return;
    fix_down(pNeighborhood);
    if (value != pTop->value && pTop->value == pNeighborhood->otherBound)
        pNeighborhood->type = AMBIGUOUS;
}

/* ------------------------------------------------------------------------
   --------- Allocations of structures used in the algorithm --------------
   ------------------------------------------------------------------------ */

/* Allocate image of the tags for visited pixels and the visited neighbors.
Do not be afraid about the parameters: pointers to 2-D arrays. */
static void init_image_of_visited_pixels(int ***ptabtabVisitedPixels)
{
    int i;

    ARRAY_2D_ALLOC(*ptabtabVisitedPixels, iWidth, iHeight, int);
    ARRAY_2D_ALLOC(tabtabVisitedNeighbors, iWidth, iHeight, int);
}

static void free_image_of_visited_pixels(int **tabtabVisitedPixels)
{
    free(tabtabVisitedPixels[0]);       /* Actually a 1-D array */
    free(tabtabVisitedPixels);

    free(tabtabVisitedNeighbors[0]);
    free(tabtabVisitedNeighbors);
}

/* Initialize the output image */
static void init_output_image(float *tabPixelsIn,
                              float ***ptabtabPixelsOutput)
{
    int i;

    *ptabtabPixelsOutput = (float **) malloc(iHeight * sizeof(float *));
    if (*ptabtabPixelsOutput == 0)
        mwerror(FATAL, 1, "init_output_image --> allocation error\n");
    for (i = 0; i < iHeight; i++)
        (*ptabtabPixelsOutput)[i] = tabPixelsIn + i * iWidth;
}

static void free_output_image(float **tabtabPixelsOutput)
{
    free(tabtabPixelsOutput);
}

static void init_region(int init_region_iMaxArea)
{
    tabPointsInShape =
        (Point_plane) malloc(init_region_iMaxArea *
                             sizeof(struct point_plane));
    if (tabPointsInShape == NULL)
        mwerror(FATAL, 1,
                "init_region --> impossible to allocate the array\n");
}

static void free_region(void)
{
    free(tabPointsInShape);
}

/* ------------------------------------------------------------------------
   --------- The core extraction algorithm --------------------------------
   ------------------------------------------------------------------------ */

/* Is pixel (x, y) a local minimum? */
static char is_local_min(float **ou, short int x, short int y,
                         char b8Connected)
            /* A 2-D array of the image */
{
    float v;
    char n = 0;

    v = ou[y][x];
    return (x == iWidth - 1 || (ou[y][x + 1] > v && ++n) || ou[y][x + 1] == v)
        && (x == 0 || (ou[y][x - 1] > v && ++n) || ou[y][x - 1] == v)
        && (y == iHeight - 1 || (ou[y + 1][x] > v && ++n)
            || ou[y + 1][x] == v) && (y == 0 || (ou[y - 1][x] > v && ++n)
                                      || ou[y - 1][x] == v)
        && (b8Connected == 0
            ||
            ((x == iWidth - 1 || y == 0 || (ou[y - 1][x + 1] > v && ++n)
              || ou[y - 1][x + 1] == v) && (x == iWidth - 1
                                            || y == iHeight - 1
                                            || (ou[y + 1][x + 1] > v && ++n)
                                            || ou[y + 1][x + 1] == v)
             && (x == 0 || y == iHeight - 1 || (ou[y + 1][x - 1] > v && ++n)
                 || ou[y + 1][x - 1] == v) && (x == 0 || y == 0
                                               || (ou[y - 1][x - 1] > v
                                                   && ++n)
                                               || ou[y - 1][x - 1] == v)))
        && n != 0;
}

/* Is pixel (x,y) a local maximum? */
static char is_local_max(float **ou, short int x, short int y,
                         char b8Connected)
            /* A 2-D array of the image */
{
    float v;
    char n = 0;

    v = ou[y][x];
    return (x == iWidth - 1 || (ou[y][x + 1] < v && ++n) || ou[y][x + 1] == v)
        && (x == 0 || (ou[y][x - 1] < v && ++n) || ou[y][x - 1] == v)
        && (y == iHeight - 1 || (ou[y + 1][x] < v && ++n)
            || ou[y + 1][x] == v) && (y == 0 || (ou[y - 1][x] < v && ++n)
                                      || ou[y - 1][x] == v)
        && (b8Connected == 0
            ||
            ((x == iWidth - 1 || y == 0 || (ou[y - 1][x + 1] < v && ++n)
              || ou[y - 1][x + 1] == v) && (x == iWidth - 1
                                            || y == iHeight - 1
                                            || (ou[y + 1][x + 1] < v && ++n)
                                            || ou[y + 1][x + 1] == v)
             && (x == 0 || y == iHeight - 1 || (ou[y + 1][x - 1] < v && ++n)
                 || ou[y + 1][x - 1] == v) && (x == 0 || y == 0
                                               || (ou[y - 1][x - 1] < v
                                                   && ++n)
                                               || ou[y - 1][x - 1] == v)))
        && n != 0;
}

/* Set pixels and saddle points in `tabPoints' at level newGray */
static void levelize(float **tabtabPixelsOutput, Point_plane tabPoints,
                     int iNbPoints, float newGray)
{
    int i;
    for (i = iNbPoints - 1; i >= 0; i--)
        tabtabPixelsOutput[tabPoints[i].y][tabPoints[i].x] = newGray;
}

/* Return, coded in one byte, the local configuration around the pixel (x,y) */
static unsigned char configuration(int **tabtabVisitedPixels, short int x,
                                   short int y)
{
    short int iMaxX = iWidth - 1, iMaxY = iHeight - 1;
    unsigned char cPattern = 0;

    if (x != 0)
    {
        if (tabtabVisitedPixels[y][x - 1] == iExploration)
            cPattern = WEST;
        if ((y == 0 && iAtBorder) ||
            (y > 0 && tabtabVisitedPixels[y - 1][x - 1] == iExploration))
            cPattern |= NORTH_WEST;
        if ((y == iMaxY && iAtBorder) ||
            (y < iMaxY && tabtabVisitedPixels[y + 1][x - 1] == iExploration))
            cPattern |= SOUTH_WEST;
    }
    else if (iAtBorder)
        cPattern = SOUTH_WEST | WEST | NORTH_WEST;

    if (x != iMaxX)
    {
        if (tabtabVisitedPixels[y][x + 1] == iExploration)
            cPattern |= EAST;
        if ((y == 0 && iAtBorder) ||
            (y > 0 && tabtabVisitedPixels[y - 1][x + 1] == iExploration))
            cPattern |= NORTH_EAST;
        if ((y == iMaxY && iAtBorder) ||
            (y < iMaxY && tabtabVisitedPixels[y + 1][x + 1] == iExploration))
            cPattern |= SOUTH_EAST;
    }
    else if (iAtBorder)
        cPattern |= SOUTH_EAST | EAST | NORTH_EAST;

    if ((y == 0 && iAtBorder) ||
        (y > 0 && tabtabVisitedPixels[y - 1][x] == iExploration))
        cPattern |= NORTH;

    if ((y == iMaxY && iAtBorder) ||
        (y < iMaxY && tabtabVisitedPixels[y + 1][x] == iExploration))
        cPattern |= SOUTH;

    return cPattern;
}

/* Insert a new shape and its siblings in the tree, with parent pParent */
static void insert_children(Shape pParent, Shape pNewChildToInsert)
{
    Shape pSibling = pNewChildToInsert;
    while (pSibling->next_sibling != NULL)
    {
        pSibling->parent = pParent;
        pSibling = pSibling->next_sibling;
    }
    pSibling->parent = pParent;
    pSibling->next_sibling = pParent->child;
    pParent->child = pNewChildToInsert;
}

static Shape new_shape(int iCurrentArea, float currentGrayLevel,
                       char bOfInferiorType, Shape pChild)
              /* Supposed to have no sibling. Can be NULL */
{
    Shape pNewShape = &pGlobalTree->the_shapes[pGlobalTree->nb_shapes++];

    pNewShape->inferior_type = bOfInferiorType;
    pNewShape->value = currentGrayLevel;
    pNewShape->open = (char) (iAtBorder != 0);
    pNewShape->area = iCurrentArea;
    pNewShape->removed = 0;
    pNewShape->pixels = NULL;
    pNewShape->data = NULL;
    pNewShape->boundary = (Flist) NULL;
    pNewShape->data_size = 0;
    /* Make links */
    pNewShape->parent = NULL;
    pNewShape->next_sibling = NULL;
    pNewShape->child = pChild;
    if (pChild != NULL)
        pChild->parent = pNewShape;
    return pNewShape;
}

/* Knowing that the last extracted shape contains the points, update,
for each one, the smallest shape containing it */
static void update_smallest_shapes(Point_plane tabPoints, int iNbPoints)
{
    int i, iIndex;
    Shape pNewShape, pRoot = &pGlobalTree->the_shapes[0];

    pNewShape = &pGlobalTree->the_shapes[pGlobalTree->nb_shapes - 1];
    for (i = iNbPoints - 1; i >= 0; i--)
    {
        iIndex = tabPoints[i].y * iWidth + tabPoints[i].x;
        if (pGlobalTree->smallest_shape[iIndex] == pRoot)
            pGlobalTree->smallest_shape[iIndex] = pNewShape;
    }
}

/* Find children of the last constructed monotone section, which is composed
of the interval between pSmallestShape and the last extracted shape. That is,
find shapes in other monotone sections whose parent is inside this interval */
static void connect(Point_plane tabPoints, int iNbPoints,
                    Connection * tabConnections, Shape pSmallestShape)
{
    int i, iIndex;
    Shape pShape, pParent;
    float level;

    for (i = iNbPoints - 1; i >= 0; i--)
    {
        iIndex = tabPoints[i].y * iWidth + tabPoints[i].x;
        pShape = tabConnections[iIndex].shape;
        if (pShape != NULL)
        {
            level = tabConnections[iIndex].level;
            pParent = pSmallestShape;
            while (pParent->value != level)
            {
                assert(pParent !=
                       &pGlobalTree->the_shapes[pGlobalTree->nb_shapes - 1]);
                pParent = pParent->parent;
            }
            insert_children(pParent, pShape);
            tabConnections[iIndex].shape = NULL;
        }
    }
}

/* Make a new connection structure at the given point */
static void new_connection(Point_plane pPoint, float level,
                           Connection * tabConnections)
{
    int iIndex;
    Shape pSibling, pShape =
        &pGlobalTree->the_shapes[pGlobalTree->nb_shapes - 1];

    iIndex = pPoint->y * iWidth + pPoint->x;
    if (tabConnections[iIndex].shape == NULL)
    {
        tabConnections[iIndex].shape = pShape;
        tabConnections[iIndex].level = level;
    }
    else
    {
        assert(tabConnections[iIndex].level == level);
        pSibling = tabConnections[iIndex].shape;
        while (pSibling->next_sibling != NULL)
            pSibling = pSibling->next_sibling;
        pSibling->next_sibling = pShape;
    }
}

/* Is the neighbor pixel already stored for this exploration? */
#define NEIGHBOR_NOT_STORED(x,y) (tabtabVisitedNeighbors[y][x] < iExploration)

/* Store the 4-neighbors of pixel (x,y) */
static void store_4neighbors(float **ou, short int x, short int y,
                             Neighborhood * pNeighborhood)
{
    if (x > 0 && NEIGHBOR_NOT_STORED(x - 1, y))
        add_neighbor(pNeighborhood, x - 1, y, ou[y][x - 1]);
    if (x < iWidth - 1 && NEIGHBOR_NOT_STORED(x + 1, y))
        add_neighbor(pNeighborhood, x + 1, y, ou[y][x + 1]);
    if (y > 0 && NEIGHBOR_NOT_STORED(x, y - 1))
        add_neighbor(pNeighborhood, x, y - 1, ou[y - 1][x]);
    if (y < iHeight - 1 && NEIGHBOR_NOT_STORED(x, y + 1))
        add_neighbor(pNeighborhood, x, y + 1, ou[y + 1][x]);
}

/* Store the diagonal neighbors of pixel (x,y) */
static void store_8neighbors(float **ou, short int x, short int y,
                             Neighborhood * pNeighborhood)
{
    if (x > 0)
    {
        if (y > 0 && NEIGHBOR_NOT_STORED(x - 1, y - 1))
            add_neighbor(pNeighborhood, x - 1, y - 1, ou[y - 1][x - 1]);
        if (y < iHeight - 1 && NEIGHBOR_NOT_STORED(x - 1, y + 1))
            add_neighbor(pNeighborhood, x - 1, y + 1, ou[y + 1][x - 1]);
    }
    if (++x < iWidth)
    {
        if (y > 0 && NEIGHBOR_NOT_STORED(x, y - 1))
            add_neighbor(pNeighborhood, x, y - 1, ou[y - 1][x]);
        if (y < iHeight - 1 && NEIGHBOR_NOT_STORED(x, y + 1))
            add_neighbor(pNeighborhood, x, y + 1, ou[y + 1][x]);
    }
}

/* Add the points in the neighborhood of gray level currentGrayLevel to the
region tabIsoPointsInShape and return 1 if a new shape is (maybe) detected.
This "maybe" is linked to `pIgnoreHoles', indicating if we can count the
holes. New points are added to `tabIsoPointsInShape'
from position `pCurrentArea'.
This value is changed at exit in case of success. */
static char add_iso_level(Point_plane tabIsoPointsInShape, int *pCurrentArea,
                          float currentGrayLevel,
                          Neighborhood * pNeighborhood, float **ou,
                          int **tabtabVisitedPixels, char *p8Connected,
                          char *pIgnoreHoles)
{
    short int x, y;
    Neighbor *pNeighbor;
    int iCurrentArea, iNbHoles;
    unsigned char cPattern;

    iNbHoles = 0;
    iCurrentArea = *pCurrentArea;
    pNeighbor = &pNeighborhood->tabPoints[1];
    do
    {                           /* 1) Neighbor is added to the region */
        x = pNeighbor->point.x;
        y = pNeighbor->point.y;
        tabIsoPointsInShape[iCurrentArea].x = x;
        tabIsoPointsInShape[iCurrentArea++].y = y;
        if (!*pIgnoreHoles)
        {
            cPattern = configuration(tabtabVisitedPixels, x, y);
            iNbHoles += tabPatterns[(int) *p8Connected][cPattern];
        }
        if (x == 0 || x == iWidth - 1 || y == 0 || y == iHeight - 1)
            ++iAtBorder;
        tabtabVisitedPixels[y][x] = iExploration;
        /* 2) Store new neighbors */
        store_4neighbors(ou, x, y, pNeighborhood);
        if (pNeighborhood->type == MAX)
        {
            if (!*p8Connected)
                *pIgnoreHoles = *p8Connected = (char) 1;
            store_8neighbors(ou, x, y, pNeighborhood);
        }
        remove_neighbor(pNeighborhood);
    }
    while (iCurrentArea <= iMaxArea &&
           pNeighbor->value == currentGrayLevel &&
           pNeighborhood->type != INVALID);

    if (iCurrentArea <= iMaxArea &&
        iAtBorder != iPerimeterImage &&
        (!iAtBorder || iCurrentArea <= iHalfAreaImage) &&
        pNeighborhood->type != INVALID && (*pIgnoreHoles || iNbHoles == 0))
    {
        *pCurrentArea = iCurrentArea;
        return (char) 1;
    }
    return 0;
}

/* Extract the terminal branch containing the point (x,y) */
static void find_terminal_branch(float **ou, int **tabtabVisitedPixels,
                                 short int x, short int y, char b8Connected,
                                 Neighborhood * pNeighborhood,
                                 Connection * tabConnections)
{
    float level;
    int iArea, iLastShapeArea;
    char bUnknownHoles;
    Shape pSmallestShape, pLastShape;

  restart_branch:
    ++iExploration;
    iArea = iAtBorder = 0;
    bUnknownHoles = 0;
    pSmallestShape = pLastShape = NULL;
    level = ou[y][x];
    reinit_neighborhood(pNeighborhood, b8Connected ? MAX : MIN);
    add_neighbor(pNeighborhood, x, y, level);
    while (add_iso_level(tabPointsInShape, &iArea,
                         level, pNeighborhood, ou, tabtabVisitedPixels,
                         &b8Connected, &bUnknownHoles) != 0)
    {
        if (bUnknownHoles)
        {
            assert(iArea != 0);
            if (pLastShape != NULL)
            {
                iArea = pLastShape->area;
                connect(tabPointsInShape, iArea, tabConnections,
                        pSmallestShape);
                new_connection(&tabPointsInShape[iArea], level,
                               tabConnections);
            }
            levelize(ou, tabPointsInShape, iArea, level);
            goto restart_branch;
        }
        if (iMinArea <= iArea)
        {                       /* Store new shape? */
            iLastShapeArea = (pLastShape == NULL) ? 0 : pLastShape->area;
            pLastShape = new_shape(iArea, level, !b8Connected, pLastShape);
            if (pSmallestShape == NULL)
                pSmallestShape = pLastShape;
            update_smallest_shapes(tabPointsInShape + iLastShapeArea,
                                   iArea - iLastShapeArea);
        }
        if (iAtBorder && iArea == iHalfAreaImage)
            break;
        bUnknownHoles = (char) (b8Connected
                                && pNeighborhood->type == AMBIGUOUS);
        if (bUnknownHoles)
            b8Connected = 0;
        level = pNeighborhood->tabPoints[1].value;
    }
    if (pLastShape != NULL)
    {
        connect(tabPointsInShape, iArea, tabConnections, pSmallestShape);
        if (iAtBorder && iArea == iHalfAreaImage)
            insert_children(pGlobalTree->the_shapes, pLastShape);
        else if (iArea != 0)
            new_connection(&tabPointsInShape[iArea], level, tabConnections);
    }
    levelize(ou, tabPointsInShape, iArea, level);
}

/* Scan the image, calling a procedure to extract terminal branch at each
   (not yet visited) local extremum */
static void scan(float **tabtabPixelsOutput, int **tabtabVisitedPixels,
                 Neighborhood * pNeighborhood, Connection * tabConnections)
{
    short int i, j;
    char b8Connected = 0;
    int iExplorationInit;

    iExplorationInit = iExploration;
    for (i = 0; i < iHeight; i++)
        for (j = 0; j < iWidth; j++)
            if (tabtabVisitedPixels[i][j] <= iExplorationInit &&
                (is_local_min(tabtabPixelsOutput, j, i, (char) 0) ||
                 (is_local_max(tabtabPixelsOutput, j, i, (char) 1)
                  && (b8Connected = 1) == 1)))
            {
                find_terminal_branch(tabtabPixelsOutput, tabtabVisitedPixels,
                                     j, i, b8Connected, pNeighborhood,
                                     tabConnections);
                b8Connected = 0;
            }
}

/* ------------------------------------------------------------------------
   --------- The main function --------------------------------------------
   ------------------------------------------------------------------------ */
/* The "Fast Level Set Transform" gives the tree of interiors of level lines
(named 'shapes') representing the image.
Only shapes of area >= *pMinArea are put in the tree. pMinArea==NULL means 1.
Output: *pTree is filled */
void flst(int *pMinArea, Fimage pImageInput, Shapes pTree)
{
    float **tabtabPixelsOutput; /* Array accessing pixels of output image */
    Neighborhood neighborhood;  /* The neighborhood of the current region */
    int **tabtabVisitedPixels;  /* Image of last visit for each pixel */
    Connection *tabConnections;
    int i;

    pGlobalTree = pTree;
    iWidth = pImageInput->ncol;
    iHeight = pImageInput->nrow;

    iAreaImage = iWidth * iHeight;
    iHalfAreaImage = (iAreaImage + 1) / 2;
    if (iWidth == 1)
        iPerimeterImage = iHeight;
    else if (iHeight == 1)
        iPerimeterImage = iWidth;
    else
        iPerimeterImage = (iWidth + iHeight - 2) * 2;
    iMinArea = (pMinArea != NULL) ? *pMinArea : 1;
    if (iMinArea > iAreaImage)
        mwerror(USAGE, 1, "min area > image");

    pGlobalTree =
        mw_change_shapes(pTree, iHeight, iWidth, pImageInput->gray[0]);
    pGlobalTree->interpolation = 0;
    tabConnections = (Connection *) malloc(iAreaImage * sizeof(Connection));
    if (tabConnections == NULL)
        mwerror(FATAL, 1, "Image of largest shapes: allocation error\n");
    for (i = iAreaImage - 1; i >= 0; i--)
        tabConnections[i].shape = NULL;

    init_output_image(pImageInput->gray, &tabtabPixelsOutput);
    init_image_of_visited_pixels(&tabtabVisitedPixels);
    init_neighborhood(&neighborhood, iAreaImage);
    init_region(iAreaImage);

    iExploration = 0;
    /* Loop with increasing iMaxArea: speed optimization. Result correct
     * with only one call to `scan' and iMaxArea = iAreaImage */
    iMaxArea = 0;
    do
    {
        if (iMaxArea == 0)
            iMaxArea = INIT_MAX_AREA;
        else
            iMaxArea <<= STEP_MAX_AREA;
        if (iMaxArea == 0 || iMaxArea >= iHalfAreaImage)
            /* iMaxArea==0: overflow */
            iMaxArea = iAreaImage - 1;
        scan(tabtabPixelsOutput, tabtabVisitedPixels, &neighborhood,
             tabConnections);
    }
    while (iMaxArea + 1 < iAreaImage);

    /* Make connections with root */
    pTree->the_shapes[0].value = tabtabPixelsOutput[0][0];
    for (i = iAreaImage - 1; i >= 0; i--)
        if (tabConnections[i].shape != NULL)
        {
            assert(tabConnections[i].level == pTree->the_shapes[0].value);
            insert_children(pTree->the_shapes, tabConnections[i].shape);
        }

    free(tabConnections);
    free_image_of_visited_pixels(tabtabVisitedPixels);
    free_region();
    free_neighborhood(&neighborhood);
    free_output_image(tabtabPixelsOutput);
}
