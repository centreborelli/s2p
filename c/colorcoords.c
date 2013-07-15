#ifndef _COLORCOORDS_C
#define _COLORCOORDS_C


#ifndef BAD_MIN
#define BAD_MIN(a,b) (b)<(a)?(b):(a)
#endif

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


static void rgb_to_hsv_doubles(double *out, double *in)
{
	//assert_rgb(in);
	double r, g, b, h, s, v, M, m;
	r = in[0]; g = in[1]; b = in[2];

	//printf("rgb %g,%g,%g...\n", r, g, b);

	if (g >= r && g >= b) {
		M = g;
		m = BAD_MIN(r, b);
		h = M == m ? 0 : 60*(b-r)/(M-m)+120;
	}
	else if (b >= g && b >= r) {
		M = b;
		m = BAD_MIN(r, b);
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
#endif//_COLORCOORDS_C
