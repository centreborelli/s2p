/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
 name = {flst_reconstruct};
 version = {"1.0"};
 author = {"Pascal Monasse, Frederic Guichard"};
 function = {"Reconstruct an image from its Fast Level Sets Transform"};
 usage = {
    tree->pTree               "Input tree of shapes",
    image<-pFloatImageOutput  "Output fimage"
};
*/

#include "mw3.h"
//#include "mw3-modules.h"
#include "mw3-modules-specific.h"

void flst_reconstruct(Shapes pTree, Fimage pFloatImageOutput)
{
    int i;
    float *pOutputPixel;
    Shape *ppShapeOfPixel;
    Shape pShape;

    if (!mw_change_fimage(pFloatImageOutput, pTree->nrow, pTree->ncol))
        mwerror(FATAL, 1,
                "flst_reconstruct --> impossible "
                "to allocate the output image\n");

    pOutputPixel = pFloatImageOutput->gray;
    ppShapeOfPixel = pTree->smallest_shape;
    for (i = pTree->nrow * pTree->ncol - 1; i >= 0; i--)
    {
        pShape = *ppShapeOfPixel++;
        while (pShape->removed)
            pShape = pShape->parent;
        *pOutputPixel++ = pShape->value;
    }
}
