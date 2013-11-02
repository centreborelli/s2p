/*
 * pm_io.h
 */

#ifndef _PM_IO_H_
#define _PM_IO_H_


#define PMMAX(i,j) ( (i)<(j) ? (j):(i) )
#define PMMIN(i,j) ( (i)<(j) ? (i):(j) )


//static char temp_char;
#define flipl2(p) temp_char = *p ;  *p = *(p+3); *(p+3) = temp_char; \
                  temp_char = *(p+1); *(p+1) = *(p+2); *(p+2) = temp_char   
#define flipl(x) flipl2((char *)(&x))
#define flips(p) temp_char = *p; *p = *(p+1); *(p+1) = temp_char



/* Get flipped value */
/* invert 2 bytes data : 01 -> 10 */
#define _mw_get_flip_b2(b2) ( (( ((unsigned short) (b2)) & 0xff00)>>8) + \
                              ((((unsigned short) (b2)) & 0x00ff)<<8) )

/* invert 4 bytes data : 0123 -> 3210 */
#define _mw_get_flip_b4(b4) ( (( ((unsigned long) (b4)) & 0xff000000)>>24)+ \
                              ((((unsigned long) (b4)) & 0x00ff0000)>>8) + \
                              ((((unsigned long) (b4)) & 0x0000ff00)<<8) + \
			      ((((unsigned long) (b4)) & 0x000000ff)<<24) ) 

/* In-place flipping */
/* invert 2 bytes data : 01 -> 10 */
#define _mw_in_flip_b2(b2)  (b2) = _mw_get_flip_b2(b2)

/* invert 4 bytes data : 0123 -> 3210 */
#define _mw_in_flip_b4(b4)  (b4) = _mw_get_flip_b4(b4)


/* For float, need to perform the flip by means of an unsigned long buffer 
   fptr : input has to be a pointer to a float
   flip_float : buffer of type unsigned long * (to be defined in the calling func.)
   output : *fptr is the flipped float value.
*/
#define _mw_in_flip_float(fptr) flip_float=(unsigned long *)(fptr); \
                                  _mw_in_flip_b4(*flip_float); 


/* For double=(u1,u2), performs (flip(u2),flip(u1)) where u1,u2 are
   long (b4 = 32 bits). 
   b4 : input has to be a pointer to a unsigned long.
   flip_double : buffer of type unsigned long (to be defined in the calling func.)
   output : *((double *) b4) is the flipped double value.   
*/

#define _mw_in_flip_double(b4)   _mw_in_flip_b4(*(b4)); \
                                 _mw_in_flip_b4(*((b4)+1)); \
                                 flip_double=*((b4)+1);\
                                 *((b4)+1)=*(b4); \
                                 *(b4)=flip_double;




#define	PM_MAGICNO	0x56494557		/* Hex for VIEW */
#define PM_MAXNELM	1024
#define PM_MAXNBAND	1024
#define PM_NOSHIFT	0
#define PM_SHIFT	1

#define PM_RED		0xff
#define PM_GREEN	0xff00
#define PM_BLUE		0xff0000
#define PM_ALPHA	0xff000000
#define PM_BW		0xffffffff
#define PM_COLOR_DEFAULT  0xffffffff

#define	PM_CMAX		0xff
#define PM_SMAX		0x7fff
#define	PM_IMAX		0x7fffffff
#define PM_FMAX		1.7E38

#define PM_AVG		0x1
#define PM_MODE		0x2

#define PM_IOHDR_SIZE	(sizeof(pmpic)-(2*sizeof(char*)))

#define	pm_max(pm)	((pm)->pm_form == PM_C ? PM_CMAX :		  \
				(pm)->pm_form == PM_S ? PM_SMAX :	  \
					(pm)->pm_form == PM_I ? PM_IMAX : \
						 PM_FMAX)

/*
 * The bottom byte of the format specifies the number of bytes in the 
 * element.  The top bit specifies that the element is in floating point
 * format.
 */
#define	PM_A		0x8000
#define	PM_C		0x8001
#define	PM_S		0x8002
#define	PM_I		0x8004
#define PM_F		0xc004

/* Added to record the color model */
#define PM_F_RGB	0xc004
#define PM_F_YUV	0xc104
#define PM_F_HSI	0xc204
#define PM_F_HSV	0xc304

/*
 * Number of bytes in the format.
 */
#define pm_nbpfm(fm)	((fm)&0xff)

/*
 * Convert PM_X format to an array index.
 *	(PM_A,0) (PM_C,1) (PM_S,2) (PM_F,3) (PM_I,4)
 */
#define pm_fmtoindex(fm)	(pm_nbpfm(fm)-(((fm)>>14)&01))

/*
 * Select and call function from array of function pointers.
 *	eg pm_sel(pm->pm_form,pm_x_add)(arg1,arg2,arge)
 */
#define	pm_sel(fm,fn)	(*((fn)[pm_fmtoindex(fm)]))

/*
 * Find a pixel given pointer to an image, number of rows, number of
 * columns number of bytes per element, the desired plane, row and
 * column. 
 */
#define pm_index(p,nr,nc,nb,cp,cr,cc)					\
        ( ((unsigned_char*)p)+((cp)*(nr)*(nc)*(nb))+((cr)*(nc)*(nb))+((cc)*(nb)) )

/*
 * pm_iindex is provided for backward compatability
 */
#define pm_iindex(a,nc,cr,cc)   ((a)+((cr)*(nc))+(cc))

/*
 * Use the above macro, but dig the info out of the pm header.
 */
#define pm_pindex(pm,cp,cr,cc)						\
	( pm_index(	(pm)->pm_image,					\
			(pm)->pm_nrow,					\
			(pm)->pm_ncol,					\
			pm_nbpe(pm),					\
			cp,						\
			cr,						\
			cc)						\
	)

/*
 * Calculate the number of bytes per element.
 */
#define pm_nbpe(pm)	(pm_nbpfm((pm)->pm_form)*(pm)->pm_nband)

/*
 * Number of elements per plane
 */
#define pm_nelmpp(p)	((p)->pm_ncol * (p)->pm_nrow)

/*
 * Number of bands times number of elements per plane
 */
#define pm_nbtelmpp(p)  (pm_nelmpp(p) * (p)->pm_nband)

/*
 * Number of elements in the image
 */
#define pm_nelm(p)	(pm_nelmpp(p) * (p)->pm_np)

/*
 * Number of bands times number of elements in the image
 */
#define pm_btelm(p)     (pm_nelm(p) * (p)->pm_nband)

/*
 * Size of the plane in bytes.
 */
#define pm_psize(p)	(pm_nelmpp(p) * pm_nbpe(p))

/*
 * Size of the image in bytes.
 */
#define pm_isize(p)	((p)->pm_np * pm_psize(p))

typedef struct {
     int pm_id;		/* Magic number for pm format files.	*/
     int pm_np;		/* Number of planes. Normally 1.	*/
     int pm_nrow;	/* Number of rows. Typically 512.	*/
     int pm_ncol;	/* Number of columns. Typically 512.	*/
     int pm_nband;	/* Number of bands. Humans use only 1.	*/
     int pm_form;	/* Pixel format.			*/
     int pm_cmtsize;	/* Number comment bytes. Includes NULL. */
     unsigned char * pm_image;	/* The image itself.		*/
     char *pm_cmt;	/* Description of operations performed.	*/
} pmpic;

/*
 * Complex format.
 */
typedef struct {
	float	r;
	float	i;
} pm_comp;

#define	PM_EBASE	100
#define PM_EMALLOC	101
#define PM_EBADPARAM	102
#define PM_EBADPIC	103
#define PM_EBADFORM	104
#define PM_EBADMAGIC	105
#define PM_ENULLPIC	106	/* picture given was NULL		   */
#define PM_EBADPLANES   107	/* invalid # of planes chosen for format   */
#define PM_EBADBANDS 	108	/* invalid # of bands chosen for format	   */
#define PM_EBADSIZE 	109	/* # of rows/cols and x offsets, y offsets */
				/*	too big for ikonas		   */
#define PM_EBADCOLORS 	110	/* invalid number of colors chosen	   */
				/*	for format			   */
#define PM_EBADCOLORPLANE 111	/* invalid color plane entered		   */


#define PM_NERROR	12
#define PM_ERROR(e)	(((e) < PM_EBASE || (e) > (PM_EBASE + PM_NERROR)) ? \
				0 : (e) - PM_EBASE)



float*  _load_pm(const char *file, int &channels, int &width, int &height);

int _create_pm(const char *file, int channels, int width, int height, float *data);


#endif /* !_PM_IO_H_ */
