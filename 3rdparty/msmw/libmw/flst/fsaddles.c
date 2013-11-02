/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
 name = {fsaddles};
 version = {"1.0"};
 author = {"Pascal Monasse"};
 function = {"Position and value of saddle points in bilinear interpolation"};
 usage = {
    image->pImage           "Input fimage",
    saddles<-pSaddlesImage  "Output fimage of saddle points",
    nb<-fsaddles            "Number of saddle points"
};
*/

#include <stdlib.h>
#include <float.h>
#include "mw3.h"
//#include "mw3-modules.h"
#include "mw3-modules-specific.h"

/* The value set at a pixel that is not a saddle point. The most appropriate
value would be NAN ('Not A Number'); unfortunately, this constant is not
systematically defined in 'math.h'. */
#define NOT_A_SADDLE FLT_MAX

/* Compute the saddle values and return the number of saddle points */
int fsaddles(Fimage pImage, Fimage pSaddlesImage)
{
    int i, j, iNumber;
    double minDiag1, maxDiag1, minDiag2, maxDiag2, temp;
    float **tabtabPixels, **tabtabSaddleValues;

    if (mw_change_fimage(pSaddlesImage, pImage->nrow - 1, pImage->ncol - 1) ==
        NULL)
        mwerror(FATAL, 1, "Image allocation error");
    tabtabPixels = mw_newtab_gray_fimage(pImage);
    tabtabSaddleValues = mw_newtab_gray_fimage(pSaddlesImage);
    if (tabtabPixels == NULL || tabtabSaddleValues == NULL)
        mwerror(FATAL, 1, "Array of lines: allocation error");

    iNumber = 0;
    for (i = pImage->nrow - 2; i >= 0; i--)
        for (j = pImage->ncol - 2; j >= 0; j--)
        {
            if ((minDiag1 = tabtabPixels[i][j]) > (maxDiag1 =
                                                   tabtabPixels[i + 1][j +
                                                                       1]))
            {
                temp = minDiag1;
                minDiag1 = maxDiag1;
                maxDiag1 = temp;
            }
            if ((minDiag2 = tabtabPixels[i + 1][j]) > (maxDiag2 =
                                                       tabtabPixels[i][j +
                                                                       1]))
            {
                temp = minDiag2;
                minDiag2 = maxDiag2;
                maxDiag2 = temp;
            }
            if (minDiag1 > maxDiag2 || maxDiag1 < minDiag2)
            {
                temp = minDiag1 + maxDiag1 - (minDiag2 + maxDiag2);
                tabtabSaddleValues[i][j] =
                    (minDiag1 * maxDiag1 - minDiag2 * maxDiag2) / temp;
                ++iNumber;
            }
            else
                tabtabSaddleValues[i][j] = NOT_A_SADDLE;
        }
    free(tabtabPixels);
    free(tabtabSaddleValues);

    return iNumber;
}
