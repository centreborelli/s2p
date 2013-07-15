#ifndef _OK_LIST_H
#define _OK_LIST_H

/* data structure and operations to maintain a finite set of points belonging
 * to a finite set of disjoint squares (or cubes).  Points and regions are
 * identified only by their inidices. */


struct ok_list {
	int number_of_regions;
	int number_of_points;	// (including the removed points!)
	int *r;		// gives the representative point of each region
	int *p;		// gives the region that contains each point
	int *plist;	// linked lists of points (by regions)
	int *buf;	// buffer to return lists of points

	// optional geometrical data useful when the regions are
	// the rectangles of a grid:
	float x0[3];
	float dx[3];
	int nx[3];
	// even more optional: a pointer to the coordinates of the points
//	float (*px)[3];
};

void ok_init(struct ok_list *, int nr, int np);
void ok_free(struct ok_list *);
void ok_add_point(struct ok_list *, int r, int p);
void ok_remove_point(struct ok_list *, int p); // the point p will never be returned
int ok_which_region(struct ok_list *, int p); // returns the index of the region
int ok_which_points(struct ok_list *, int r); // returns the number of points
					      // (and fills buf)

//#ifdef USE_IMAGE_STRUCTURES
//#include "image3d.h"
void ok_init_grid(struct ok_list *, int);
int ok_add_geo_point(struct ok_list *, float [3], int);
						// just like add_point, but
						// the region index is
						// implicit from the
						// coordinates
						//
int ok_regionindex_neigs(struct ok_list *, float [3]);
//#endif /* USE_IMAGE_STRUCTURES */
void ok_svg_layer(void *, struct ok_list *);
void ok_hack_assert_consistency(struct ok_list *);


/* data structure for non-necessarily disjoint regions
 *
 * NOTE: this stuff is typically not needed, because it can be replaced
 * by a low-brow implementation that  add  the points to a single region
 * and forces the user to search for each point in "neighboring" regions.
 *
 * Since this is typically sufficient, this structure is left implemented.
 *
 * See the structure "square_grid" for a comfortable way to operate with
 * square grids.
 *
 *
 */

//struct mok_list {
//	int number_of_regions;
//	int number_of_points;
//	int maximum_occupancy;
//	/* ...  */
//	int *rbuf;
//	int *pbuf;
//};
//void mok_init(struct mok_list *, int nr, int np, int mocc);
//void mok_free(struct mok_list *);
//void mok_add_point(struct mok_list *, int r, int p);
//int mok_which_regions(struct mok_list *, int p);
//int mok_which_points(struct mok_list *, int r);


#endif /* _OK_LIST_H */
