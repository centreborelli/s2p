
#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

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


#include <stdarg.h>

//static void error(const char *fmt, ...) __attribute__((noreturn,format(printf,1,2)));
static void error(const char *fmt, ...) __attribute__((noreturn));
static void error(const char *fmt, ...)

{
	va_list argp;
	fprintf(stderr, "\nERROR(\"%s\"): ", myname());
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
#ifdef NDEBUG
	exit(-1);
#else//NDEBUG
	//print_trace(stderr);
	exit(*(int *)0x43);
#endif//NDEBUG
}


#include <stdlib.h>

static void *xmalloc(size_t size)
{
	if (size == 0)
		error("xmalloc: zero size");
	void *new = malloc(size);
	if (!new)
	{
		double sm = size / (0x100000 * 1.0);
		error("xmalloc: out of memory when requesting "
			"%zu bytes (%gMB)",//:\"%s\"",
			size, sm);//, strerror(errno));
	}
	return new;
}

static void *xrealloc(void *p, size_t s)
{
	void *r = realloc(p, s);
	if (!r) error("realloc failed");
	return r;
}

static void xfree(void *p)
{
	if (!p)
		error("thou shalt not free a null pointer!");
	free(p);
}

static void *matrix_build(int w, int h, size_t n)
{
	size_t p = sizeof(void *);
	char *r = xmalloc(h*p + w*h*n);
	for (int i = 0; i < h; i++)
		*(void **)(r + i*p) = r + h*p + i*w*n;
	return r;
}

static double random_uniform(void)
{
	return rand()/(RAND_MAX+1.0);
}

static double random_normal(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y1 = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	//double y2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	return y1;
}

#include <math.h>

static void hsv_to_rgb_doubles(double *out, double *in)
{
	//assert_hsv(in);
	double r, g, b, h, s, v; r=g=b=h=s=v=0;
	h = in[0]; s = in[1]; v = in[2];
	if (s == 0)
		r = g = b = v;
	else {
		int H = fmod(floor(h/60),6);
		double p, q, t, f = h/60 - H;
		p = v * (1 - s);
		q = v * (1 - f*s);
		t = v * (1 - (1 - f)*s);
		switch (H) {
			case 6:
			case 0: r = v; g = t; b = p; break;
			case 1: r = q; g = v; b = p; break;
			case 2: r = p; g = v; b = t; break;
			case 3: r = p; g = q; b = v; break;
			case 4: r = t; g = p; b = v; break;
			case -1:
			case 5: r = v; g = p; b = q; break;
			default:
				fprintf(stderr, "H=%d\n", H);
				assert(false);
		}
	}
	out[0] = r; out[1] = g; out[2] = b;
	//assert_rgb(out);
}


#if 0
static void assert_rgb(double t[3])
{
	for (int i = 0; i < 3; i++)
		assert(t[i] >= 0 && t[i] <= 1);
}
#endif

static void assert_hsv(double t[3])
{
	if (t[0] < 0 || t[0] >= 360) error("queca %g\n", t[0]);
	assert(t[0] >= 0 && t[0] < 360);
	if (!(t[1] >= 0 && t[1] <= 1))
		error("CACA S = %g\n", t[1]);
	assert(t[1] >= 0 && t[1] <= 1);
	assert(t[2] >= 0 && t[2] <= 1);
}

static void rgb_to_hsv_doubles(double *out, double *in)
{
	//assert_rgb(in);
	double r, g, b, h, s, v, M, m;
	r = in[0]; g = in[1]; b = in[2];

	//printf("rgb %g,%g,%g...\n", r, g, b);

	if (g >= r && g >= b) {
		M = g;
		m = fmin(r, b);
		h = M == m ? 0 : 60*(b-r)/(M-m)+120;
	}
	else if (b >= g && b >= r) {
		M = b;
		m = fmin(r, b);
		h = M == m ? 0 : 60*(r-g)/(M-m)+240;
	}
	else {
		assert(r >= g && r >= b);
		M = r;
		if (g >= b) {
			m = b;
			h = M == m ? 0 : 60*(g-b)/(M-m)+0;
		} else {
			m = g;
			h = M == m ? 0 : 60*(g-b)/(M-m)+360;
		}
	}

	s = M == 0 ? 0 : (M - m)/M;
	v = M;
	h = fmod(h, 360);

	//printf("\thsv %g,%g,%g\n", h, s, v);
	out[0] = h; out[1] = s; out[2] = v;
	//assert_hsv(out);
}

// draw a segment between two points
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (abs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px+j*slope), j+py, e);
		}
	}
}

// draw a segment between two points (somewhat anti-aliased)
void traverse_segment_aa(int px, int py, int qx, int qy,
		void (*f)(int,int,float,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, 1.0, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment_aa(qx, qy, px, py, f, e);
	else {
		if (abs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py); slope /= (qx - px);
			assert(px < qx);
			assert(fabs(slope) <= 1);
			for (int i = 0; i <= qx-px; i++) {
				float exact = py + i*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(i+px, whole, 1-part, e);
				f(i+px, owhole, part, e);
			}
		} else { // vertical
			float slope = (qx - px); slope /= (qy - py);
			assert(abs(qy - py) >= abs(qx - px));
			assert(py < qy);
			assert(fabs(slope) <= 1);
			for (int j = 0; j <= qy-py; j++) {
				float exact = px + j*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(whole, j+py, 1-part, e);
				f(owhole, j+py, part, e);
			}
		}
	}
}

// draw a segment between two points (somewhat anti-aliased)
void traverse_segment_aa2(float px, float py, float qx, float qy,
		void (*f)(int,int,float,void*), void *e)
{
	//if (px == qx && py == qy)
	//	f(px, py, 1.0, e);
	//else
	if (qx + qy < px + py) // bad quadrants
		traverse_segment_aa2(qx, qy, px, py, f, e);
	else {
		if (fabs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py); slope /= (qx - px);
			assert(px < qx);
			assert(fabs(slope) <= 1);
			for (int i = 0; i <= qx-px; i++) {
				float exact = py + i*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(i+px, whole, 1-part, e);
				f(i+px, owhole, part, e);
			}
		} else { // vertical
			float slope = (qx - px); slope /= (qy - py);
			assert(abs(qy - py) >= abs(qx - px));
			assert(py < qy);
			assert(fabs(slope) <= 1);
			for (int j = 0; j <= qy-py; j++) {
				float exact = px + j*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(whole, j+py, 1-part, e);
				f(owhole, j+py, part, e);
			}
		}
	}
}

int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}
