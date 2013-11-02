/*
 * fsignal.h
 */

#ifndef _FSIGNAL_H_
#define _FSIGNAL_H_

/* src/fsignal.c */
Fsignal mw_new_fsignal(void);
Fsignal mw_alloc_fsignal(Fsignal signal, int N);
void mw_delete_fsignal(Fsignal signal);
Fsignal mw_change_fsignal(Fsignal signal, int N);
void mw_clear_fsignal(Fsignal signal, float v);
void mw_copy_fsignal_values(Fsignal in, Fsignal out);
void mw_copy_fsignal_header(Fsignal in, Fsignal out);
void mw_copy_fsignal(Fsignal in, Fsignal out);

#endif /* !_FSIGNAL_H_ */
