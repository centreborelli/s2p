#ifndef _XFOPEN_C
#define _XFOPEN_C

#include <stdio.h>
#include <string.h>

#include "fail.c"

static FILE *xfopen(const char *s, const char *p)
{
	FILE *f;

	if (!s) fail("trying to open a file with NULL name");

	if (0 == strcmp("-", s))
	{
		if (0 == strcmp("w", p))
			return stdout;
		else if (0 == strcmp("r", p))
			return stdin;
		else
			fail("unknown fopen mode \"%s\"", p);
	}
	if (0 == strcmp("--", s) && 0 == strcmp("w", p)) return stderr;

	f = fopen(s, p);
	if (f == NULL)
		fail("can not open file \"%s\" in mode \"%s\", // (%s)",
				s, p);//, strerror(errno));
	return f;
}

static void xfclose(FILE *f)
{
	if (f != stdout && f != stdin && f != stderr) {
		int r = fclose(f);
		if (r) fail("fclose error \"%s\"", "");//strerror(errno));
	}
}

#endif//_XFOPEN_C
