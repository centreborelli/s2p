/*
 * list.h
 */

#ifndef _LIST_H_
#define _LIST_H_

/* src/list.c */
Flist mw_new_flist(void);
Flist mw_realloc_flist(Flist l, int n);
Flist mw_enlarge_flist(Flist l);
Flist mw_change_flist(Flist l, int max_size, int size, int dimension);
void mw_delete_flist(Flist l);
void mw_clear_flist(Flist l, float v);
Flist mw_copy_flist(Flist in, Flist out);
Flists mw_new_flists(void);
Flists mw_realloc_flists(Flists ls, int n);
Flists mw_enlarge_flists(Flists ls);
Flists mw_change_flists(Flists ls, int max_size, int size);
void mw_delete_flists(Flists ls);
Flists mw_copy_flists(Flists in, Flists out);
Dlist mw_new_dlist(void);
Dlist mw_realloc_dlist(Dlist l, int n);
Dlist mw_enlarge_dlist(Dlist l);
Dlist mw_change_dlist(Dlist l, int max_size, int size, int dimension);
void mw_delete_dlist(Dlist l);
void mw_clear_dlist(Dlist l, double v);
Dlist mw_copy_dlist(Dlist in, Dlist out);
Dlists mw_new_dlists(void);
Dlists mw_realloc_dlists(Dlists ls, int n);
Dlists mw_enlarge_dlists(Dlists ls);
Dlists mw_change_dlists(Dlists ls, int max_size, int size);
void mw_delete_dlists(Dlists ls);
Dlists mw_copy_dlists(Dlists in, Dlists out);

#endif /* !_LIST_H_ */
