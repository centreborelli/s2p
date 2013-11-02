#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "libpm_io.h"




float*  _load_pm(const char *file, int &channels, int &width, int &height)
{
	
	
	
	int isize,j;
	char need_flipping = 0;
	pmpic *pm;
	FILE *fp;
	
	pm = (pmpic*)malloc(sizeof(pmpic));
	if (pm == NULL)
	{ printf("ERROR::_load_pm()::Not enough memory to load any PM image !\n"); return 0;}
	
	memset(pm,0,sizeof(pmpic));
	
	/*
	 pm->pm_id = PM_MAGICNO;
	 pm->pm_np = pm->pm_nband = 1;
	 pm->pm_nrow = pm->pm_ncol = 512;
	 pm->pm_form = PM_C;
	 pm->pm_cmtsize = 0;
	 pm->pm_image = NULL;
	 pm->pm_cmt = NULL;
	 */
	
	if (!(fp = fopen(file, "rb")))
	{
		fclose(fp);
		free(pm);
		return NULL;
	}
	
	
	if (fread(pm,PM_IOHDR_SIZE,1,fp) == 0)
	{
		free(pm);
		fclose(fp);
		return NULL;
	}
	
	
	if (pm->pm_id != PM_MAGICNO)
	{
		_mw_in_flip_b4(pm->pm_id);
		if (pm->pm_id == PM_MAGICNO)
			need_flipping = 1;
	}
	
	if (pm->pm_id != PM_MAGICNO)
	{
		free(pm);
		fclose(fp);
		return NULL;
	}
	
	
	if (need_flipping)
	{
		_mw_in_flip_b4(pm->pm_np);
		_mw_in_flip_b4(pm->pm_nrow);
		_mw_in_flip_b4(pm->pm_ncol);
		_mw_in_flip_b4(pm->pm_nband);
		_mw_in_flip_b4(pm->pm_form);
		_mw_in_flip_b4(pm->pm_cmtsize);
	}
	
	
	width = pm->pm_ncol;
	height = pm->pm_nrow;
	channels = pm->pm_np;
	
	float *data = new float[channels * width * height];
	
    isize = pm_isize(pm);
    if (pm->pm_form == PM_C)  // char image. convert to float
	{
		
		unsigned char *data_8b = new unsigned char [pm->pm_nrow * pm->pm_ncol * pm->pm_np];
		if (fread(data_8b,channels*width*height*sizeof(unsigned char),1,fp) == 0)
		{
			free(pm);
			fclose(fp);
			return NULL;
		}
		
		
		for( j = 0; j < width*height*channels; j--) data[j] = (float) data_8b[j];
		
		delete[] data_8b;
		
	}else if (pm->pm_form == PM_F)	// float image
	{
		
		
		if (fread(data,channels*width*height*sizeof(float),1,fp) == 0)
		{
			free(pm);
			fclose(fp);
			return NULL;
		}
		
		
	}else 
	{
		
		free(pm);
		fclose(fp);
		return NULL;
		
	}
	
	
	
	fclose(fp);
	free(pm);
	return data;
	
}








int _create_pm(const char *file, int channels, int width, int height, float *data)
{
	
	
	FILE *fp;
	pmpic   *pm;
	int isize;
	
	if (!(fp = fopen(file, "w")))
	{
		//printf("ERROR::_create_pm:: Cannot create the file \"%s\"\n",file);
		return 0;
	}
	
	
	pm = (pmpic*)malloc(sizeof(pmpic));
	if (pm == NULL) {  printf("ERROR::_create_pm::Not enough memory to create any PM image !\n"); return 0;}
	
	memset(pm,0,sizeof(pmpic));
	
	pm->pm_id = PM_MAGICNO;
	pm->pm_np = channels;
	pm->pm_nband = 1;
	pm->pm_nrow = height;
	pm->pm_ncol = width;
	pm->pm_form = PM_F;
	pm->pm_cmtsize = strlen("nocomments");
	pm->pm_image = NULL;
	pm->pm_cmt = NULL;
	
	
	if (fwrite(pm,PM_IOHDR_SIZE,1,fp) == 0)
	{
		free(pm);
		fclose(fp);
		return 0;
	}
	
	isize = pm_isize(pm);
	if (fwrite(data,channels*width*height*sizeof(float),1,fp) == 0) 
	{
				free(pm);
				fclose(fp);
				return 0;
			
	}
	
		
	if (pm->pm_cmtsize)
	{
		if (fwrite("nocomments",pm->pm_cmtsize,1,fp) == 0)
		{			
			fclose(fp);
			free(pm);
			return(0);
			
		}
	}
	
	
	free(pm);
	fclose(fp);
	return 1;

}


/*
int _create_pm(const char *file, int nchannels, int width, int height, float **tab, char *type)
{


     FILE *fp;
     pmpic   *pm;
     int isize,i,j;

     if (!(fp = fopen(file, "w")))
     {
	  printf("ERROR::_create_pm:: Cannot create the file \"%s\"\n",file);
	  return 0;
     }

     pm = (pmpic*)malloc(sizeof(pmpic));
     if (pm == NULL) {  printf("ERROR::_create_pm::Not enough memory to create any PM image !\n"); return 0;}

     memset(pm,0,sizeof(pmpic));

     pm->pm_id = PM_MAGICNO;
     pm->pm_np = nchannels;
     pm->pm_nband = 1;
     pm->pm_nrow = height;
     pm->pm_ncol = width;

     if (strcmp(type,"PM_C") == 0)    {pm->pm_form = PM_C;}
     else if (strcmp(type,"PM_S") == 0)     {pm->pm_form = PM_S;}
     else if (strcmp(type,"PM_F") == 0)     {pm->pm_form = PM_F;}
     else {printf("ERROR::_create_pm():: type must be char, short or float \n");}

     pm->pm_cmtsize = strlen("nocomments");
     pm->pm_image = NULL;
     pm->pm_cmt = NULL;


     if (fwrite(pm,PM_IOHDR_SIZE,1,fp) == 0)
     {
	  printf("ERROR::_create_pm:: writing pm header in file \"%s\" !\n",file);
	  free(pm);
	  fclose(fp);
	  return 0;
     }

     isize = pm_isize(pm);

     if (pm->pm_form == PM_C)
	{

		unsigned char *gray = new unsigned char [pm->pm_nrow * pm->pm_ncol];

		for(i=0; i < nchannels; i++)
		{
			for(j=pm->pm_nrow * pm->pm_ncol-1; j>=0; j--) gray[j] = (unsigned char) PMMAX(0.0, PMMIN(tab[i][j], 255.0)) ;

			 if (fwrite(gray,width*height*sizeof(unsigned char),1,fp) == 0) 
     			{
				printf("ERROR::_create_pm():: while writing file \"%s\" !\n",file);
				free(pm);
				fclose(fp);
				return 0;
     			}

		}

		delete[] gray;

	}
	else if (pm->pm_form == PM_S)
	{

		unsigned short *gray = new unsigned short [pm->pm_nrow * pm->pm_ncol];

		for(i=0; i < nchannels; i++)
		{
			for(j=pm->pm_nrow * pm->pm_ncol-1; j>=0; j--) gray[j] = (unsigned short) PMMAX(0.0, PMMIN(tab[i][j], 65535.0)) ;

			 if (fwrite(gray,width*height*sizeof(unsigned short),1,fp) == 0) 
     			{
				printf("ERROR::_create_pm():: while writing file \"%s\" !\n",file);
				free(pm);
				fclose(fp);
				return 0;
     			}

		}

		delete[] gray;


	} else
	{


		for(i=0; i < nchannels; i++)
		{
			 if (fwrite(tab[i],width*height*sizeof(float),1,fp) == 0) 
     			{
				printf("Error while writing file \"%s\" !\n",file);
				free(pm);
				fclose(fp);
				return 0;
     			}
		}

	}



	if (pm->pm_cmtsize)
	{
		if (fwrite("nocomments",pm->pm_cmtsize,1,fp) == 0)
		{			
			fclose(fp);
			free(pm);
			return(0);
			
		}
	}


     free(pm);
     fclose(fp);
     return 1;


}
*/