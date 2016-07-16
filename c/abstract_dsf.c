#include <assert.h>
#include "abstract_dsf.h"

// API
void adsf_assert_consistency(int *t, int n)
{
	assert(n > 0);
	assert(t);
	for (int i = 0; i < n; i++) {
		assert(t[i] >= 0);
		assert(t[i] < n);
	}
}

// API
void adsf_begin(int *t, int n)
{
	for (int i = 0; i < n; i++)
		t[i] = i;
}

// API
int adsf_find(int *t, int n, int a)
{
	assert(a >= 0 && a < n);
	if (a != t[a])
		t[a] = adsf_find(t, n, t[a]);
	return t[a];
}

static int adsf_make_link(int *t, int n, int a, int b)
{
	if (a < b) { // arbitrary choice
		t[b] = a;
		return a;
	} else {
		t[a] = b;
		return b;
	}
}

// API
int adsf_union(int *t, int n, int a, int b)
{
	assert(a >= 0 && a < n);
	assert(b >= 0 && b < n);
	a = adsf_find(t, n, a);
	b = adsf_find(t, n, b);
	if (a != b)
		b = adsf_make_link(t, n, a, b);
	return b;
}

// API
//int adsf_number_of_classes(int *t, int n)
//{
//	adsf_assert_consistency(t, n);
//	error("not yet implemented");
//	return -1;
//}
