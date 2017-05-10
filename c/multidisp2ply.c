#include "pickopt.c"

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <libgen.h>

#include "vvector.h"
#include "../3rdparty/iio/iio.h"
#include "rpc.h"
#include "triangulation.h"
#include "coordconvert.h"
#include "read_matrix.c"
#include "parsenumbers.c"

void utm_alt_zone(double *out, double lat, double lon, int zone);
void utm_zone(int *zone, bool *northp, double lat, double lon);

unsigned char test_little_endian(void)
{
    int x = 1;
    return (*(char*) & (x) == 1);
}

void write_ply_header(FILE* f, bool ascii, int npoints, int zone, bool hem,
        bool colors, bool normals)
{
    if (!ascii)
        if (!test_little_endian())
              fail("BINARY PLY NOT SUPPORTED ON BIG ENDIAN SYSTEMS!\n");

    fprintf(f, "ply\n");
    if (ascii)
        fprintf(f, "format ascii 1.0\n");
    else
        fprintf(f, "format binary_little_endian 1.0\n");

    fprintf(f, "comment created by S2P\n");
    if (zone >= 0)
        fprintf(f, "comment projection: UTM %i%s\n", zone, (hem ? "N" : "S"));
    fprintf(f, "element vertex %d\n", npoints);
    fprintf(f, "property double x\n");
    fprintf(f, "property double y\n");
    fprintf(f, "property double z\n");
    if (normals) {
        fprintf(f, "property double nx\n");
        fprintf(f, "property double ny\n");
        fprintf(f, "property double nz\n");
    }
    if (colors) {
        fprintf(f, "property uchar red\n");
        fprintf(f, "property uchar green\n");
        fprintf(f, "property uchar blue\n");
    }
    fprintf(f, "end_header\n");
}

static void parse_utm_string(int *zone, bool *hem, const char *s)
{
  if (0 == strcmp(s, "no_utm_zone")) {
    *zone = -1;
    return;
  }
  char hem_string[FILENAME_MAX];
  if (2 == sscanf(s, "%02d%s", zone, hem_string)) {
    // hem_string must be equal to "N" or "S"
    if (hem_string[1] == '\0') {
      if (hem_string[0] == 'N' || hem_string[0] == 'S') {
	*hem = (hem_string[0] == 'N');
	return;
      }
    }
  }
  fprintf(stderr, "zone: %d\themisphere: %s\n", *zone, hem_string);
  fprintf(stderr, "incorrect value for --utm-zone."
	  " It must look like '27N'\n");
  *zone = -1;
}

static void alloc_full_outputs (float ***img_vec_tab,
				float ***rpj_vec_tab,
				int **nb_views,
				int ***img_selected_views,
				FILE **obj_file,
				int nb_sights,
				int width,
				int height,
				const char *tile_dir)
{
  // track the vector from an opt. point
  // to a given view
  *img_vec_tab = (float **) malloc(nb_sights*sizeof( float * ));
  for(int i=0;i<nb_sights;i++)
    (*img_vec_tab)[i] = (float *) malloc(3*width*height*sizeof( float ));

  // reproject the above vectors
  // into original geometry
  *rpj_vec_tab = (float **) malloc(nb_sights*sizeof( float * ));
  for(int i=0;i<nb_sights;i++)
    (*rpj_vec_tab)[i] = (float *) malloc(3*width*height*sizeof( float ));
  
  // we wish we could know
  // the number of views used for each pixel
  // as well as which views have been used for each pixel
  *nb_views = (int *) calloc(width*height, sizeof(int)); // view <=> sight  
  *img_selected_views = (int **) malloc(nb_sights*sizeof(int * ));
  for(int i=0;i<nb_sights;i++)
    (*img_selected_views)[i] = (int *) calloc(width*height, sizeof(int));  
  
  // let's vizualise sights !
  char buf[1000];
  sprintf(buf,"%s/sights.obj",tile_dir);
  *obj_file = fopen(buf,"w");
  fprintf(*obj_file, "o sights for tile %s\n", tile_dir);

  // init with NAN
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
      for (int i=0; i<nb_sights; i++)
	for (int t=0; t<3; t++)
	  {
	    (*img_vec_tab)[i][width*3*y+3*x+t] = NAN;
	    (*rpj_vec_tab)[i][width*3*y+3*x+t] = NAN;
	  }
}

static void update_full_outputs(float ***img_vec_tab,
				float ***rpj_vec_tab,
				int **nb_views,
				int ***img_selected_views,
				FILE **obj_file,
				int local_nb_sights,
				int width,
				int height,
				const char *tile_dir,
				float ***errMap_tab,
				type_sight *sights_list,
				int x,
				int y,
				list_of_pairs *list_pairs,
				unsigned int *vertex_index)
{
  int posH = x + width*y;
  for (int s = 0 ; s < local_nb_sights; s++)
    {
      // err by sight : dist from opt point to sight i
      // (index 0 of errMap_tab is dedicated to rms error)
      (*errMap_tab)[sights_list[s].ID][posH] = sights_list[s].err;
      
      // vec err by sight : vector from opt point to sight i
      // (sights numbered from 1 to N so index = ID-1)
      for(int t=0; t<3; t++)
	(*img_vec_tab)[sights_list[s].ID-1][width*3*y+3*x+t] =
	  sights_list[s].err_vec_ptopt_to_sight[t+3] - sights_list[s].err_vec_ptopt_to_sight[t];
      
      // same vectors as above, but reprojected into
      // ref image (so two components : drow and dcol)
      // (sights numbered from 1 to N so index = ID-1)
      double lgt1,lat1,alt1,pos1[2];
      double lgt2,lat2,alt2,pos2[2];
      ECEF_to_lgt_lat_alt(sights_list[s].err_vec_ptopt_to_sight[0],
			  sights_list[s].err_vec_ptopt_to_sight[1],
			  sights_list[s].err_vec_ptopt_to_sight[2],
			  &lgt1,&lat1,&alt1);
      ECEF_to_lgt_lat_alt(sights_list[s].err_vec_ptopt_to_sight[3],
			  sights_list[s].err_vec_ptopt_to_sight[4],
			  sights_list[s].err_vec_ptopt_to_sight[5],
			  &lgt2,&lat2,&alt2);
      
      eval_rpci(pos1, &list_pairs->data[0].rpc_master, lgt1, lat1, alt1);
      eval_rpci(pos2, &list_pairs->data[0].rpc_master, lgt2, lat2, alt2);
      
      for(int t=0; t<2; t++)
	(*rpj_vec_tab)[sights_list[s].ID-1][width*3*y+3*x+t] = pos2[t]-pos1[t];
      
      // nb of sights used for current pixel
      (*nb_views)[posH] = local_nb_sights;
      
      // tells which sight has been used for current pixel
      // (sights numbered from 1 to N so index = ID-1)
      (*img_selected_views)[sights_list[s].ID-1][posH]=1;
      
      // build obj
      double objpoint1[3],objpoint2[3];
      for(int t=0; t<3; t++)
	{
	  objpoint1[t] = sights_list[s].err_vec_ptopt_to_sight[t+3] + 0.5*sights_list[s].v[t];
	  objpoint2[t] = sights_list[s].err_vec_ptopt_to_sight[t+3] - 0.5*sights_list[s].v[t];
	}
      fprintf((*obj_file), "v %f %f %f\n", objpoint1[0],objpoint1[1],objpoint1[2]);
      fprintf((*obj_file), "v %f %f %f\n", objpoint2[0],objpoint2[1],objpoint2[2]);
      fprintf((*obj_file), "l %d %d\n", vertex_index,vertex_index+1);
      vertex_index += 2;
    }
}

static void free_full_outputs(float ***img_vec_tab,
			      float ***rpj_vec_tab,
			      int **nb_views,
			      int ***img_selected_views,
			      FILE **obj_file,
			      int nb_sights,
			      int width,
			      int height,
			      const char *tile_dir)
{
  char buf[1000];
  sprintf(buf, "%s/nb_sights.tif",tile_dir);
  iio_save_image_int(buf, *nb_views, width, height);
  free(*nb_views);
  
  for(int i=0;i<nb_sights;i++)
    {
      sprintf(buf,"%s/selected_sight_%d.tif",tile_dir,i+1);
      iio_save_image_int(buf, (*img_selected_views)[i], width, height);
      sprintf(buf,"%s/rpc_err_vec_sight_%d.tif",tile_dir,i+1);
      iio_save_image_float_vec(buf, (*img_vec_tab)[i], width, height, 3);
      sprintf(buf,"%s/rpc_err_rpjvec_sight_%d.tif",tile_dir,i+1);
      iio_save_image_float_vec(buf, (*rpj_vec_tab)[i], width, height, 3);
    }

  for(int i=0;i<nb_sights;i++)
    {
      free((*img_selected_views)[i]);
      free((*img_vec_tab)[i]);
      free((*rpj_vec_tab)[i]);
    }

  if (nb_sights > 1)
    {
      free(*img_selected_views);
      free(*img_vec_tab);
      free(*rpj_vec_tab);
    }
}



static void help(char *s)
{
  
  fprintf(stderr, "\t usage: %s out.ply N --disp[1..N] disp.tif"
	  "--rpc_ref rpc_ref.xml --rpc_sec[1..N] rpc_sec.xml [--color colors.png] "
	  "[--utm-zone ZONE] [--ascii]"
	  "[--lon-m l0] [--lon-M lf] [--lat-m l0] [--lat-M lf] "
	  "[--col-m x0] [--col-M xf] [--row-m y0] [--row-M yf]\n", s);
}

int main_disp_to_heights(int c, char *v[])
{
  if (c < 3)
    {
      help(v[0]);
      return 1;
    }

  // get the number of pairs
  int N_pair = atoi(v[2]);

  // utm zone and hemisphere: true for 'N' and false for 'S'
  int zone;
  bool hem;
  const char *utm_string = pick_option(&c, &v, "-utm-zone", "no_utm_zone");
  parse_utm_string(&zone, &hem, utm_string);

  // ascii flag
  bool ascii = pick_option(&c, &v, "-ascii", NULL);

  // longitude-latitude bounding box
  double lon_m = atof(pick_option(&c, &v, "-lon-m", "-inf"));
  double lon_M = atof(pick_option(&c, &v, "-lon-M", "inf"));
  double lat_m = atof(pick_option(&c, &v, "-lat-m", "-inf"));
  double lat_M = atof(pick_option(&c, &v, "-lat-M", "inf"));

  // x-y bounding box
  double col_m = atof(pick_option(&c, &v, "-col-m", "-inf"));
  double col_M = atof(pick_option(&c, &v, "-col-M", "inf"));
  double row_m = atof(pick_option(&c, &v, "-row-m", "-inf"));
  double row_M = atof(pick_option(&c, &v, "-row-M", "inf"));

  bool full_outputs = false;

  // read input rpc models
  char *rpc_ref_string = pick_option(&c, &v, "-rpc_ref", "");

  // initialize list_pairs
  int nb_pairs = N_pair;
  int nb_sights = N_pair+1;
  list_of_pairs list_pairs;
  list_pairs.real_size = nb_pairs;
  list_pairs.tot_size = nb_pairs;
  list_pairs.data = (type_pair *) malloc(list_pairs.real_size*sizeof(type_pair));

  int width = (col_M - col_m);
  int height = (row_M - row_m);

  for (int i = 0; i < N_pair ; i++)
    {
      char buf[100];
      list_pairs.data[i].ID = i+1;
      list_pairs.data[i].sight_slave = i+2;
	
      int *nx = &list_pairs.data[i].nx;
      int *ny = &list_pairs.data[i].ny;
      int *nch = &list_pairs.data[i].nch;

      sprintf(buf, "-disp%d", i+1);  
      const char *disp_string = pick_option(&c, &v, buf, "");
      list_pairs.data[i].dispx = iio_read_image_float_split(disp_string, nx, ny, nch);
      list_pairs.data[i].dispy = list_pairs.data[i].dispx + (*nx)*(*ny);
	
      sprintf(buf, "-rpc_sec%d", i+1);
      char *rpc_sec_string = pick_option(&c, &v, buf, "");
      read_rpc_file_xml(&list_pairs.data[i].rpc_master, rpc_ref_string);
      read_rpc_file_xml(&list_pairs.data[i].rpc_slave, rpc_sec_string);
    }

  // open color images if provided
  int wwc, hhc, pd;
  uint8_t *clr = NULL;
  const char *color_string = pick_option(&c, &v, "-color", "");
  clr = iio_read_image_uint8_vec(color_string, &wwc, &hhc, &pd);

  // need at least one pair to do something
  if (list_pairs.real_size == 0)
    return 1;
    
  
  char* plyfilename = v[1];
  const char *tile_dir = dirname(strdup(plyfilename));   

  float *heightMap = malloc(width*height*sizeof(float));
  double *ecef = malloc(3*width*height*sizeof(double));
    
  int size_of_fout_err_tab;
  if (full_outputs)
    size_of_fout_err_tab = nb_sights+1;
  else
    size_of_fout_err_tab = 1;
    
  float **errMap_tab = (float **) malloc((size_of_fout_err_tab)*sizeof(float *));
  for(int i=0; i<size_of_fout_err_tab; i++)
    errMap_tab[i] = (float *) calloc(width*height, sizeof(float));
    
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
      {
	heightMap[x + width*y] = NAN;
	for(int t=0; t<3; t++)
	  ecef[width*3*y+3*x+t] = NAN;
      }
	  
  /* ----------------------------------------
  /* full_output initialisation
  /* ---------------------------------------- */
  float** img_vec_tab, **rpj_vec_tab;
  int *nb_views, **img_selected_views;
  FILE *obj_file;
  unsigned int vertex_index = 1;

  if (full_outputs)
    {
      alloc_full_outputs(&img_vec_tab, &rpj_vec_tab,
			 &nb_views, &img_selected_views, &obj_file,
			 nb_sights, width, height, tile_dir);
    }

  int npoints = 0;
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
      {
	int local_nb_sights = N_pair+1;
	int posH = x + width*y;

	// for each pair (pid == pair id)
	for(int pid=0;pid<list_pairs.real_size;pid++)
	  {
	    list_pairs.data[pid].q0[0] = x+col_m;
	    list_pairs.data[pid].q0[1] = y+row_m;
	    list_pairs.data[pid].q0[2] = 1.;

	    double posx = list_pairs.data[pid].dispx[posH];
	    double posy = list_pairs.data[pid].dispy[posH];
	     
	    list_pairs.data[pid].q1[0] = posx;
	    list_pairs.data[pid].q1[1] = posy;
	    list_pairs.data[pid].q1[2] = 1.;

	    if (isnan(posx) || isnan(posy))
	      {
		list_pairs.data[pid].process = false;
		local_nb_sights --;
	      }
	    else
	      {
		list_pairs.data[pid].process = true;
	      }
	  }

	if (local_nb_sights >= 2)
	  {
	    // Altitude
	    double hg;

	    // list of sights
	    type_sight *sights_list = (type_sight *) malloc(local_nb_sights*sizeof(type_sight));

	    // Compute the related height & rpc_err
	    hg = rpc_height_geo(&list_pairs,
				local_nb_sights,
				sights_list);

	    // 1) errors
	    double sum = 0.0;
	    double nb_elt = 0.0;

	    for(int s=0; s<local_nb_sights; s++)
	      {
		sum += pow(sights_list[s].err, 2.0);
		nb_elt++;
	      }

	    if (full_outputs)
	      {
		update_full_outputs(&img_vec_tab, &rpj_vec_tab,
				    &nb_views, &img_selected_views, &obj_file,
				    local_nb_sights, width, height, tile_dir,
				    &errMap_tab, sights_list,
				    x, y, &list_pairs, &vertex_index);
	      }
		 
	    // r.m.s. error
	    sum = sqrt(sum/nb_elt);
		  
	    // Output the results in original geometry
		  
	    // * height
	    heightMap[posH] = hg;
		  
	    // * ecef
	    // - sights_list is a list of "type_sight" element (see c/rpc.h)
	    // - err_vec_ptopt_to_sight is a field made up of a 6 double tab
	    // - the three first elements correspond to the optimal point in  ECEF coord
	    // - each sight is related to the optimal point, so we only need the first one.
	    for(int t=0; t<3; t++)
	      ecef[width*3*y+3*x+t] = sights_list[0].err_vec_ptopt_to_sight[t];
	      
	    // * r.m.s. error (dedicated to index 0)
	    errMap_tab[0][posH] = sum;
	    npoints ++;

	    // if not defined, utm zone is that of the first point
	    if (zone < 0)
	      {
		// positions
		double X[3];
		ECEF_to_lgt_lat_alt(ecef[width*3*y+3*x], 
				    ecef[width*3*y+3*x+1], 
				    ecef[width*3*y+3*x+2],
				    &X[0],&X[1],&X[2]);

		 // check with lon/lat bounding box
		 if (X[0] < lon_m || X[0] > lon_M || X[1] < lat_m || X[1] > lat_M)
		   continue;

		 utm_zone(&zone, &hem, X[1], X[0]);
	       }

	     free(sights_list);
	   }
       }

   // print header for ply file
   FILE *ply_file = fopen(plyfilename, "w");
   if (ply_file)
     {
       write_ply_header(ply_file, ascii, npoints, zone, hem, (bool) clr, false);
     }
   else
     {
       fprintf(stderr, "Can't open this file: %s\n", v[1]);
       return 1;
     }
   
   for (int y = 0; y < height; y++)
     for (int x = 0; x < width; x++) 
       {
	 int pix = y*width + x;

	 bool ok = true;
	 for(int t = 0; t < 3; t++)
	   if (isnan(ecef[width*3*y+3*x+t]))
	     ok=false;

	 if (ok) 
	   {
	     // positions
	     double X[3];
	     ECEF_to_lgt_lat_alt(ecef[width*3*y+3*x], 
				 ecef[width*3*y+3*x+1], 
				 ecef[width*3*y+3*x+2],
				 &X[0],&X[1],&X[2]);

	     // check with lon/lat bounding box
	     if (X[0] < lon_m || X[0] > lon_M || X[1] < lat_m || X[1] > lat_M)
	       continue;

	     // convert (lon, lat, alt) to utm
	     utm_alt_zone(X, X[1], X[0], zone);
	       
	     // colorization: if greyscale, copy the grey level on each channel
	     uint8_t rgb[3];
	     if (clr) {
	       for (int k = 0; k < pd; k++) rgb[k] = clr[k + pd*pix];
	       for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
	     }
	       
	     // write to ply
	     if (ascii) {
	       fprintf(ply_file, "%a %a %a ", X[0], X[1], X[2]);
	       if (clr)
		 fprintf(ply_file, "%d %d %d", rgb[0], rgb[1], rgb[2]);
	       fprintf(ply_file, "\n");
	     } else {
	       double XX[3] = {X[0], X[1], X[2]};
	       fwrite(XX, sizeof(double), 3, ply_file);
	       if (clr) {
		 unsigned char C[3] = {rgb[0], rgb[1], rgb[2]};
		 fwrite(rgb, sizeof(unsigned char), 3, ply_file);
	       }		 
	     }
	   }  
       }

   // * save images
   char buf[1000];
   sprintf(buf,"%s/height_map.tif",tile_dir);
   iio_save_image_float(buf, heightMap, width, height);
   sprintf(buf,"%s/ecef_coord.tif",tile_dir);
   iio_save_image_double_vec(buf, ecef, width, height, 3);

   sprintf(buf,"%s/rpc_err_rms_allsights.tif",tile_dir);
   iio_save_image_float(buf, errMap_tab[0], width, height);
  
   for(int i=1;i<size_of_fout_err_tab;i++)
     {
       sprintf(buf,"%s/rpc_err_norm_sight_%d.tif",tile_dir,i);
       iio_save_image_float(buf, errMap_tab[i], width, height);
     }

   if (full_outputs)
     {
       free_full_outputs(&img_vec_tab, &rpj_vec_tab,
			 &nb_views, &img_selected_views, &obj_file,
			 nb_sights, width, height, tile_dir);
     }

   // clean mem
   free(heightMap);
   free(ecef);

   for(int i=0;i<list_pairs.real_size;i++)
     {
       free(list_pairs.data[i].dispx);
     }
   free(list_pairs.data);  

   for(int i=0;i<size_of_fout_err_tab;i++)
     free(errMap_tab[i]);
   if (size_of_fout_err_tab > 1)
     free(errMap_tab);

   return 0;
}

int main(int c, char *v[])
{
  return main_disp_to_heights(c, v);
}
