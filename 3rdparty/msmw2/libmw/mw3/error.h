/*
 * error.h
 */

#ifndef _ERROR_H_
#define _ERROR_H_

/* src/error.c */
void mwdebug(char *fmt, ...);
void mwerror(int code, int exit_code, char *fmt, ...);

#endif /* !_ERROR_H_ */
