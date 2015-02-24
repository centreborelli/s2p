/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand 
  name = {fgrain}; 
  version = {"1.2"}; 
  author = {"Pascal Monasse, Frederic Guichard, G. Facciolo"}; 
  function = {"Grain filter of an image"}; 
  usage = { 
    'a': [min_area=20]-> pMinArea   "Min area of grains we keep", 
    image_in -> pFloatImageInput    "Input fimage", 
    image_out <- pFloatImageOutput  "Output fimage" 
    }; 
*/ 
/*----------------------------------------------------------------------
 v1.2 (04/2007): simplified header (LM)
 v1.3 (2013): portable version (GF)
----------------------------------------------------------------------*/

#include "mw3.h" 
#include "stdio.h"
#include "stdlib.h"

extern void flst();
extern void flst_reconstruct();
 
/* This removes the shapes from the tree associated to pFloatImageInput 
that are too small (threshold *pMinArea). As a consequence all the remaining 
shapes of pFloatImageOutput are of area larger or equal than *pMinArea */ 

void mw_fgrain(int *pMinArea, Fimage pFloatImageInput, Fimage pFloatImageOutput) 
{ 
  int i; 
  Shapes pTree; 
 
  if(mw_change_fimage(pFloatImageOutput, pFloatImageInput->nrow, 
		      pFloatImageInput->ncol) == NULL) 
    mwerror(FATAL, 1, 
	    "fgrain --> Not enough memory to allocate the output image"); 
  if((pTree = mw_new_shapes()) == NULL) 
    mwerror(FATAL, 1, 
	    "fgrain --> Not enough memory to allocate the tree of shapes"); 

  /* Compute the Level Sets Transform of the input image */ 
  flst(NULL, pFloatImageInput, pTree, NULL, NULL); 
 
  /* Kill too small grains. 
     Bound i>0 because it is forbidden to delete the root, at index 0 */ 
  for(i = pTree->nb_shapes-1; i > 0; i--) 
    if(pTree->the_shapes[i].area < *pMinArea) 
      pTree->the_shapes[i].removed = (char)1; 
 
  /* Reconstruct in pFloatImageOutput the modified tree of shapes */ 
  flst_reconstruct(pTree, pFloatImageOutput); 
 
  mw_delete_shapes(pTree); 
} 

/*in and out must be allocated*/
void fgrain(int MinArea, float *in, int nx, int ny, float *out) {
   int i;
   Fimage mwin = mw_change_fimage(NULL,ny,nx);
   Fimage mwout = mw_new_fimage();
   for(i=0;i<nx*ny;i++) mwin->gray[i] = in[i]; 

   mw_fgrain( &MinArea, mwin, mwout);

   for(i=0;i<nx*ny;i++) out[i] = mwout->gray[i]; 
   mw_delete_fimage(mwin);
   mw_delete_fimage(mwout);
}
