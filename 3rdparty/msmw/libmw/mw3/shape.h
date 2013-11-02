/*
 * shape.h
 */

#ifndef _SHAPE_H_
#define _SHAPE_H_

/* src/shape.c */
Point_plane mw_new_point_plane(void);
Point_plane mw_change_point_plane(Point_plane point);
Shape mw_new_shape(void);
Shape mw_change_shape(Shape sh);
void mw_delete_shape(Shape shape);
Shape mw_get_not_removed_shape(Shape sh);
Shape mw_get_parent_shape(Shape sh);
Shape mw_get_first_child_shape(Shape sh);
Shape mw_get_next_sibling_shape(Shape sh);
Shape mw_get_smallest_shape(Shapes shs, int iX, int iY);
Shapes mw_new_shapes(void);
Shapes mw_alloc_shapes(Shapes shs, int nrow, int ncol, float value);
Shapes mw_change_shapes(Shapes shs, int nrow, int ncol, float value);
void mw_delete_shapes(Shapes shs);

#endif /* !_SHAPE_H_ */
