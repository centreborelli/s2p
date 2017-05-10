#ifndef _TRIANGULATION_H
#define _TRIANGULATION_H

#include <stdbool.h>
#include <stdlib.h>
#include "vvector.h"
#include "coordconvert.h"
#include "rpc.h"
#include <string.h>

// 3x3 matrix object 
struct mat3x3 {
    double val[3][3];
};

// struct that defines the required 
// information to process a pair of tiles 
typedef struct {
    
    int ID;
    int sight_slave;
    bool process;
    
    // RPC
    struct rpc rpc_master;
    struct rpc rpc_slave;
    
    // homographies
    struct mat3x3 H_ref;
    struct mat3x3 H_sec;
    
    // pointing corrections
    struct mat3x3 pointing_correc;
    
    // take into account
    // pointing corrections
    // H_secB = H_sec * A^(-1)
    struct mat3x3 H_secB;
    
    // inverse slave homographies
    struct mat3x3 H_invSecB;
    
    // disparity maps
    int nx,ny,nch;
    float *dispx;
    float *dispy;
    
    // mask
    float *mask;
    
    double q0[3]; // a point inside ref image
    double q1[3]; // a point inside a given slave image
    
} type_pair;

// A simple list of pairs of tiles.
typedef struct
{
    type_pair *data;
    int tot_size;
    int real_size;
} list_of_pairs;

// define sight object
typedef struct 
{
    // sight ID
    int ID;

    // the sight passes through points "s" ans "p"
    // (in ECEF coord.)
    double s[3]; 
    double p[3]; 
    // its unit direction vector is "v"
    double v[3]; 
    
    // distance from the optimal point to this sight
    double err;
    
    // indices 0 to 2 : optimal point in  ECEF coord
    // indices 3 to 5 : closest point in this sight
    // to optimal point in  ECEF coord.
    // Finally, err-...[i+3] - err-...[i] gives
    // the ith component of the smallest vector starting
    // from the optimal point and ending to a point in this sight
    double err_vec_ptopt_to_sight[6];
    
} type_sight;

// distance between a 3D point P and a line (defined by its
// normalized direction vector V and passing through a 3D point S) 
double dist_line_point3D(double *V,double *S,double *P);

// compute the vector VEC between a 3D point P0 and a line (
// passing through two 3D points P and S)
double vec_pt_to_line3D(double *P,double *S,double *P0,double *VEC);

// find closest 3D point from from a set of 3D lines
void find_point_opt(type_sight *sights_list, int N,
		    double *point_opt, double *outerr);

// compute the height of a point given its location inside two images
// geometric solution
double rpc_height_geo(list_of_pairs *list_pairs, 
		      int local_nb_sights, 
		      type_sight *sights_list);

#endif // _TRIANGULATION_H
