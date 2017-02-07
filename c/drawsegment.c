#ifndef _DRAWSEGMENT_C
#define _DRAWSEGMENT_C

#include <assert.h>
#include <math.h>

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

#endif//_DRAWSEGMENT_C
