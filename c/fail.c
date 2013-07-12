#ifndef _FAIL_C
#define _FAIL_C

#include <stdio.h>
#include <stdlib.h>

#ifdef __linux
#  include <sys/types.h>
#  include <unistd.h>
static const char *emptystring = "";
static const char *myname(void)
{
#  define n 0x29a
	//const int n = 0x29a;
	static char buf[n];
	pid_t p = getpid();
	snprintf(buf, n, "/proc/%d/cmdline", p);
	FILE *f = fopen(buf, "r");
	if (!f) return emptystring;
	int c, i = 0;
	while ((c = fgetc(f)) != EOF && i < n) {
#  undef n
		buf[i] = c ? c : ' ';
		i += 1;
	}
	if (i) buf[i-1] = '\0';
	fclose(f);
	return buf;
}
#else
static const char *myname(void) { return ""; }
#endif//__linux

#ifdef DOTRACE
#include <execinfo.h>
#endif

#ifndef BACKTRACE_SYMBOLS
#define BACKTRACE_SYMBOLS 50
#endif

static void print_trace(FILE *f)
{
	(void)f;
#ifdef DOTRACE
	void *array[BACKTRACE_SYMBOLS];
	size_t size, i;
	char **strings;

	size = backtrace (array, BACKTRACE_SYMBOLS);
	strings = backtrace_symbols (array, size);

	fprintf (f, "Obtained %zu stack frames.\n", size);

	for (i = 0; i < size; i++)
		fprintf (f, "%s\n", strings[i]);

	free (strings);
#endif
}

#include <stdarg.h>

//static void fail(const char *fmt, ...) __attribute__((noreturn,format(printf,1,2)));
static void fail(const char *fmt, ...) __attribute__((noreturn));
static void fail(const char *fmt, ...)

{
	va_list argp;
	fprintf(stderr, "\nFAIL(\"%s\"): ", myname());
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
	print_trace(stderr);
#ifdef NDEBUG
	exit(-1);
#else//NDEBUG
	exit(*(volatile int *)0x43);
#endif//NDEBUG
}

#endif//_FAIL_C
