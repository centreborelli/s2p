#include "triangulation.h"

// distance between a 3D point P and a line (defined by its
// normalized direction vector V and passing through a 3D point S) 
double dist_line_point3D(double *V,double *S,double *P)
{
	double cross_prod[3];
	double sm[3];
	double norm;
		
	sm[0] = P[0] - S[0];
	sm[1] = P[1] - S[1];
	sm[2] = P[2] - S[2];

	VEC_CROSS_PRODUCT(cross_prod,V,sm);
	VEC_LENGTH(norm,cross_prod);
	
	return norm;
}

// compute the vector VEC between a 3D point P0 and a line (
// passing through two 3D points P and S)
double vec_pt_to_line3D(double *P,double *S,double *P0,double *VEC)
{
	double diff10[3];
	double diff21[3];
	double dot_prod_diff10_diff21;
	double t,norm_diff21;
	double VT[3];
	double diffVT0[3];
	double norm_diffVT0;
	
	VEC_DIFF(diff10,P,P0);
	VEC_DIFF(diff21,S,P);
	VEC_LENGTH(norm_diff21,diff21);
	VEC_DOT_PRODUCT(dot_prod_diff10_diff21,diff10,diff21);
	t=-dot_prod_diff10_diff21 / ( norm_diff21*norm_diff21 );
	VEC_COPY(VT,P)
	VEC_ACCUM(VT,t,diff21);
	VEC_DIFF(diffVT0,VT,P0);
	
	for(int t=0;t<3;t++)
	{
		VEC[t] = P0[t];
		VEC[t+3] = VT[t];
	}
	
	VEC_LENGTH(norm_diffVT0,diffVT0);
	return norm_diffVT0;
}

// find closest 3D point from from a set of 3D lines
void find_point_opt(type_sight *sights_list, int N,
		    double *point_opt, double *outerr)
{
	double prod1[3][3];
	double prod2[3];
	double Id[3][3];
	double mat_sum[3][3]={{0,0,0},
					   {0,0,0},
					   {0,0,0}};
	double prod_sum[3]={0,0,0};
	
	for(int i=0; i<N; i++)
	  {
	    OUTER_PRODUCT_3X3(prod1,sights_list[i].v,sights_list[i].v);
	    IDENTIFY_MATRIX_3X3(Id);
	    ACCUM_SCALE_MATRIX_3X3(Id,-1,prod1); // Now Id actually contains a mat diff
	    ACCUM_SCALE_MATRIX_3X3(mat_sum,1,Id); // accumulates it into matrix 'mat_sum'
	    MAT_DOT_VEC_3X3(prod2,Id,sights_list[i].s); // Id still contains a mat diff
	    VEC_ACCUM(prod_sum,1,prod2); // accumulates the above result
	  }
	
	double mat_sum_inv[3][3];
	double det;
	INVERT_3X3(mat_sum_inv,det,mat_sum);
	
	// find the optimal point
	MAT_DOT_VEC_3X3(point_opt,mat_sum_inv,prod_sum);
}

//build the set of sights
void build_set_sights(type_sight *sights_list,list_of_pairs *list_pairs)
{
	int N_pairs = list_pairs->real_size;
	double alt1=3000.;
	double alt2=45.;

	bool first_sight_found = false;

	// pid (pair id), to loop over pairs
	// ind, insertion variable
	int pid=0, ind=0; 
	while(pid<N_pairs)
	{
		if (list_pairs->data[pid].process)
		{
			int ID;
			struct rpc *rpc_to_use;
			double X1,Y1,Z1,X2,Y2,Z2,q[3];
			if (!first_sight_found)
			{
				rpc_to_use = &list_pairs->data[pid].rpc_master;
				ID = 1;
				q[0] = list_pairs->data[pid].q0[0];
				q[1] = list_pairs->data[pid].q0[1];
				q[2] = list_pairs->data[pid].q0[2];
				first_sight_found = true;
				// look at this pair twice :
				// 1st : get the first sight
				// 2nd : get the slave sight
				// so : pid-- to go backwards
				pid--; 
			}
			else
			{
				rpc_to_use = &list_pairs->data[pid].rpc_slave;
				ID = list_pairs->data[pid].sight_slave;
				q[0] = list_pairs->data[pid].q1[0];
				q[1] = list_pairs->data[pid].q1[1];
				q[2] = list_pairs->data[pid].q1[2];
			}
			
			double point1[2],point2[2];
			eval_rpc(point2, rpc_to_use, q[0], q[1], alt2);
			eval_rpc(point1, rpc_to_use, q[0], q[1], alt1);

			geotedic_to_ECEF(point1[0],point1[1],alt1,&X1,&Y1,&Z1);
			geotedic_to_ECEF(point2[0],point2[1],alt2,&X2,&Y2,&Z2);
			
			sights_list[ind].p[0] = X1;
			sights_list[ind].p[1] = Y1;
			sights_list[ind].p[2] = Z1;
			sights_list[ind].s[0] = X2;
			sights_list[ind].s[1] = Y2;
			sights_list[ind].s[2] = Z2;
			sights_list[ind].v[0] = X2-X1;
			sights_list[ind].v[1] = Y2-Y1;
			sights_list[ind].v[2] = Z2-Z1;
			VEC_NORMALIZE(sights_list[ind].v);
			sights_list[ind].ID = ID;
			
			ind++;
		}
		pid++;
	}
}

// compute the height of a point given its location inside two images
// (geometric solution)
double rpc_height_geo(list_of_pairs *list_pairs,
		      int local_nb_sights, 
		      type_sight *sights_list)
{
	// build the set of sights
	build_set_sights(sights_list, list_pairs);

	double *outerr = (double *) malloc(local_nb_sights*sizeof(double));
	
	double point_opt[3];
	
	int i,j;
	double lgt,lat,h;
	// final estimation, without the outliers
	find_point_opt(sights_list, local_nb_sights,
		       point_opt, outerr);
	
	// Errors, defined as the distance
	// between each viewing line
	// and the optimal point
	for(int i=0; i<local_nb_sights; i++)
	  sights_list[i].err = vec_pt_to_line3D(sights_list[i].p,
						sights_list[i].s,
						point_opt,
						sights_list[i].err_vec_ptopt_to_sight);

	// compute altitude h
	h = get_altitude_from_ECEF(point_opt[0],
				   point_opt[1],
				   point_opt[2]);

	// clean mem
	free(outerr);
	
	return h;
}
