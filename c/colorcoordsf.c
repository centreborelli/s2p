#ifndef COLORCOORDSF_C
#define COLORCOORDSF_C


#ifndef BAD_MIN
#define BAD_MIN(a,b) (b)<(a)?(b):(a)
#endif

static void hsv_to_rgb_floats(float *out, float *in)
{
	//assert_hsv(in);
	float r, g, b, h, s, v; r=g=b=h=s=v=0;
	h = in[0]; s = in[1]; v = in[2];
	if (s == 0)
		r = g = b = v;
	else {
		int H = fmod(floor(h/60),6);
		float p, q, t, f = h/60 - H;
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


static void rgb_to_hsv_floats(float *out, float *in)
{
	//assert_rgb(in);
	float r, g, b, h, s, v, M, m;
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

// Commission Intérnationale de l'Éclairage 1931
static void rgb_to_xyz_floats(float *xyz, float *rgb)
{
	float n = 0.17697;
	xyz[0] = (0.49    * rgb[0] + 0.31   * rgb[1] + 0.2     * rgb[2]) / n;
	xyz[1] = (0.17697 * rgb[0] + 0.8124 * rgb[1] + 0.01063 * rgb[2]) / n;
	xyz[2] = (0       * rgb[0] + 0.01   * rgb[1] + 0.99    * rgb[2]) / n;
}

// Commission Intérnationale de l'Éclairage 1931
static void xyz_to_rgb_floats(float *rgb, float *xyz)
{
	rgb[0] =  0.41847 * xyz[0]    - 0.15866 * xyz[1]   - 0.082835 * xyz[2];
	rgb[1] = -0.091169 * xyz[0]   + 0.25243 * xyz[1]   + 0.015708 * xyz[2];
	rgb[2] =  0.00092090 * xyz[0] - 0.0025498 * xyz[1] + 0.1786 * xyz[2];
}

#endif//COLORCOORDSF_C
