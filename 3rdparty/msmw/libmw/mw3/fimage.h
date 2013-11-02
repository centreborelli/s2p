/*
 * fimage.h
 */

#ifndef _FIMAGE_H_
#define _FIMAGE_H_

/* src/fimage.c */
Fimage mw_new_fimage(void);
Fimage mw_alloc_fimage(Fimage image, int nrow, int ncol);
void mw_delete_fimage(Fimage image);
Fimage mw_change_fimage(Fimage image, int nrow, int ncol);
float mw_getdot_fimage(Fimage image, int x, int y);
void mw_plot_fimage(Fimage image, int x, int y, float v);
void mw_draw_fimage(Fimage image, int a0, int b0, int a1, int b1, float c);
void mw_clear_fimage(Fimage image, float v);
void mw_copy_fimage(Fimage in, Fimage out);
float **mw_newtab_gray_fimage(Fimage image);

#endif /* !_FIMAGE_H_ */
