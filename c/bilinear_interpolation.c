static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static float getsamplec(float *fx, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (l < 0) l = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	if (l >= pd) l = pd-1;
	float (*x)[w][pd] = (void*)fx;
	return x[j][i][l];
	//return x[(i+j*w)*pd + l];
}

static void bilinear_interpolation_vec_at(float *result,
		float *x, int w, int h, int pd,
		float p, float q)
{
	int ip = p;
	int iq = q;
	for (int l = 0; l < pd; l++) {
		float a = getsamplec(x, w, h, pd, ip  , iq  , l);
		float b = getsamplec(x, w, h, pd, ip+1, iq  , l);
		float c = getsamplec(x, w, h, pd, ip  , iq+1, l);
		float d = getsamplec(x, w, h, pd, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
		result[l] = r;
	}
}

static float bilinear_interpolation_at(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getsamplec(x, w, h, 1, ip  , iq  , 0);
	float b = getsamplec(x, w, h, 1, ip+1, iq  , 0);
	float c = getsamplec(x, w, h, 1, ip  , iq+1, 0);
	float d = getsamplec(x, w, h, 1, ip+1, iq+1, 0);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}
