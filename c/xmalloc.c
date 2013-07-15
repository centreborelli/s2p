#ifndef _XMALLOC_C
#define _XMALLOC_C

#include <stdlib.h>

#include "fail.c"

static void *xmalloc(size_t size)
{
#ifndef NDEBUG
	{
		double sm = size / (0x100000 * 1.0);
		if (sm > 1000)
			fprintf(stderr, "WARNING: large malloc"
					" %zu bytes (%gMB)\n", size, sm);
	}
#endif
	if (size == 0)
		fail("xmalloc: zero size");
	void *new = malloc(size);
	if (!new)
	{
		double sm = size / (0x100000 * 1.0);
		fail("xmalloc: out of memory when requesting "
			"%zu bytes (%gMB)",//:\"%s\"",
			size, sm);//, strerror(errno));
	}
	return new;
}

inline // to avoid unused warnings
static void *xrealloc(void *p, size_t s)
{
	void *r = realloc(p, s);
	if (!r) fail("realloc failed");
	return r;
}

inline // to avoid unused warnings
static void xfree(void *p)
{
	if (!p)
		fail("thou shalt not free a null pointer!");
	free(p);
}

#endif//_XMALLOC_C
