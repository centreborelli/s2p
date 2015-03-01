/**
 * @file mw3.h
 *
 * @author Nicolas Limare (2009)
 *
 * mw3 API header
 */

#ifndef _MW3_H_
#define _MW3_H_

#include <stdio.h>

/* definitions.h */
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846   /* pi */
#endif
#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923   /* pi/2 */
#endif
#ifndef M_PI_4
#define M_PI_4         0.78539816339744830962   /* pi/4 */
#endif

#define WARNING  0
#define ERROR    1
#define FATAL    2
#define USAGE    3
#define INTERNAL 4

#define MODEL_RGB 0
#define MODEL_YUV 1
#define MODEL_HSI 2
#define MODEL_HSV 3

#define _MW_DATA_ASCII_FILE_HEADER "MegaWave2 - DATA ASCII file -\n"

#define mw_not_an_argument 1.0e9
#define mw_not_a_magnitude -1.0

#define mw_mtype_size 20
#define mw_ftype_size 20
#define mw_cmtsize 255
#define mw_namesize 255

#define SIZE_OF_MW2_BIN_TYPE_ID 6

#define _mw_get_flip_b2(b2) (((((unsigned short) (b2)) & 0xff00)>>8) +  \
                             ((((unsigned short) (b2)) & 0x00ff)<<8) )
#define _mw_get_flip_b4(b4) (((((unsigned long) (b4)) & 0xff000000)>>24)+ \
                             ((((unsigned long) (b4)) & 0x00ff0000)>>8) + \
                             ((((unsigned long) (b4)) & 0x0000ff00)<<8) + \
                             ((((unsigned long) (b4)) & 0x000000ff)<<24) )

#define _mw_in_flip_b2(b2)  do                  \
     {                                          \
          (b2) = _mw_get_flip_b2(b2);           \
     }                                          \
     while (0)
#define _mw_in_flip_b4(b4)  do                  \
     {                                          \
          (b4) = _mw_get_flip_b4(b4);           \
     }                                          \
     while (0)

#define _mw_in_flip_float(b4) do                \
     {                                          \
          unsigned long * flip_float;           \
          flip_float = (unsigned long *)(&b4);  \
          _mw_in_flip_b4(* flip_float);         \
     }                                          \
     while (0)
#define _mw_in_flip_double(b4) do               \
     {                                          \
          unsigned long * flip_float, tmp;      \
          flip_float = (unsigned long *)(&b4);  \
          _mw_in_flip_b4(* flip_float);         \
          _mw_in_flip_b4(* (flip_float + 1));   \
          tmp = * (flip_float + 1);             \
          * (flip_float + 1)= * (flip_float);   \
          * (flip_float) = tmp;                 \
     }                                          \
     while (0)

#define mw_bzero(s, n) do      \
     {                         \
          memset((s), 0, (n)); \
     }                         \
     while (0)

typedef struct cimage {
    int nrow;              /**< number of rows (dy)    */
    int ncol;              /**< number of columns (dx) */
    int allocsize;         /**< size allocated (in bytes) for the gray plane */
    unsigned char *gray;   /**< the gray level plane (may be NULL)           */

    float scale;             /**< scale of the picture */
    char cmt[mw_cmtsize];    /**< comments */
    char name[mw_namesize];
                             /**< name of the image */

    /* defines the signifiant part of the picture */
    int firstcol;     /**< first col not affected by left side effect */
    int lastcol;      /**< last col not affected by right side effect */
    int firstrow;     /**< first row not aff. by upper side effect    */
    int lastrow;      /**< last row not aff. by lower side effect     */

    /* for use in movies */
    struct cimage *previous;   /**< pointer to the prev image */
    struct cimage *next;       /**< pointer to the next image */
} *Cimage;

typedef struct fimage {
    int nrow;       /**< number of rows (dy)    */
    int ncol;       /**< number of columns (dx) */
    int allocsize;
                    /**< size allocated (in bytes) for the gray plane */
    float *gray;    /**< the gray level plane (may be NULL)           */

    float scale;             /**< scale of the picture */
    char cmt[mw_cmtsize];    /**< comments */
    char name[mw_namesize];
                             /**< name of the image */

    /* defines the signifiant part of the picture */
    int firstcol;     /**< first col not affected by left side effect */
    int lastcol;      /**< last col not affected by right side effect */
    int firstrow;     /**< first row not aff. by upper side effect    */
    int lastrow;      /**< last row not aff. by lower side effect     */

    /* for use in movies */
    struct fimage *previous;   /**< pointer to the prev image */
    struct fimage *next;       /**< pointer to the next image */
} *Fimage;

typedef struct ccimage {
    int nrow;               /**< number of rows (dy)    */
    int ncol;               /**< number of columns (dx) */
    int allocsize;          /**< size allocated (in bytes) for the gray plane */

    unsigned char model;    /**< model of the color system           */
    unsigned char *red;     /**< the red   level plane (may be NULL) */
    unsigned char *green;   /**< the green level plane (may be NULL) */
    unsigned char *blue;    /**< the blue  level plane (may be NULL  */

    float scale;             /**< scale of the picture */
    char cmt[mw_cmtsize];    /**< comments             */
    char name[mw_namesize];
                             /**< name of the image    */

    /* defines the signifiant part of the picture */
    int firstcol;     /**< first col not affected by left side effect */
    int lastcol;      /**< last col not affected by right side effect */
    int firstrow;     /**< first row not aff. by upper side effect    */
    int lastrow;      /**< last row not aff. by lower side effect     */

    /* for use in movies */
    struct ccimage *previous;   /**< pointer to the prev image */
    struct ccimage *next;       /**< pointer to the next image */
} *Ccimage;

typedef struct cfimage {
    int nrow;       /**< number of rows (dy)    */
    int ncol;       /**< number of columns (dx) */
    int allocsize;
                    /**< size allocated (in bytes) for the gray plane */

    unsigned char model;  /**< model of the color system           */
    float *red;           /**< the red   level plane (may be NULL) */
    float *green;         /**< the green level plane (may be NULL) */
    float *blue;          /**< the blue  level plane (may be NULL) */

    float scale;             /**< scale of the picture */
    char cmt[mw_cmtsize];    /**< comments */
    char name[mw_namesize];
                             /**< name of the image */

    /* defines the signifiant part of the picture */
    int firstcol;     /**< first col not affected by left side effect */
    int lastcol;      /**< last col not affected by right side effect */
    int firstrow;     /**< first row not aff. by upper side effect    */
    int lastrow;      /**< last row not aff. by lower side effect     */

    /* for use in movies */
    struct cfimage *previous;   /**< pointer to the prev image */
    struct cfimage *next;       /**< pointer to the next image */
} *Cfimage;

typedef struct point_curve {
    int x, y;
               /**< coordinates of the point */

    /* for use in curve */
    struct point_curve *previous;   /**< pointer to the prev point */
    struct point_curve *next;       /**< pointer to the next point */
} *Point_curve;

typedef struct curve {
    Point_curve first;  /** pointer to the first point of the curve */

    /* For use in curves */
    struct curve *previous;   /** pointer to the prev curve */
    struct curve *next;       /** pointer to the next curve */
} *Curve;
typedef struct curves {
    char cmt[mw_cmtsize];    /** comments                   */
    char name[mw_namesize];
                             /** name of the set            */
    Curve first;             /** pointer to the first curve */
} *Curves;

typedef struct point_fcurve {
    float x, y;  /**< coordinates of the point */

    /* for use in curve */
    struct point_fcurve *previous;   /**< pointer to the prev point */
    struct point_fcurve *next;       /**< pointer to the next point */
} *Point_fcurve;

typedef struct fcurve {
    Point_fcurve first;  /** pointer to the first point of the fcurve */

    /* For use in curves */
    struct fcurve *previous;   /** pointer to the prev curve */
    struct fcurve *next;       /** pointer to the next curve */
} *Fcurve;
typedef struct fcurves {
    char cmt[mw_cmtsize];    /** comments                   */
    char name[mw_namesize];
                             /** name of the set            */
    Fcurve first;             /** pointer to the first curve */
} *Fcurves;

typedef struct point_dcurve {
    double x, y;  /**< coordinates of the point */

    /* for use in curve */
    struct point_dcurve *previous;   /**< pointer to the prev point */
    struct point_dcurve *next;       /**< pointer to the next point */
} *Point_dcurve;

typedef struct dcurve {
    Point_dcurve first;  /** pointer to the first point of the fcurve */

    /* For use in curves */
    struct dcurve *previous;   /** pointer to the prev curve */
    struct dcurve *next;       /** pointer to the next curve */
} *Dcurve;
typedef struct dcurves {
    char cmt[mw_cmtsize];    /** comments                   */
    char name[mw_namesize];
                             /** name of the set            */
    Dcurve first;             /** pointer to the first curve */
} *Dcurves;

typedef struct polygon {
    /* the number of elements is given by nb_channels */
    int nb_channels;    /**< number of channels */
    float *channel;     /**< tab to the channel */
    Point_curve first;  /**< pointer to the first point of the curve */

    /* for use in polygons only */
    struct polygon *previous;   /**< pointer to the prev poly. (may be NULL) */
    struct polygon *next;       /**< pointer to the next poly. (may be NULL) */
} *Polygon;
typedef struct polygons {
    char cmt[mw_cmtsize];    /**< comments                     */
    char name[mw_namesize];
                             /**< name of the set              */
    Polygon first;           /**< pointer to the first polygon */
} *Polygons;

typedef struct fpolygon {
    /* the number of elements is given by nb_channels */
    int nb_channels;     /**< number of channels */
    float *channel;      /**< tab to the channel */
    Point_fcurve first;  /**< pointer to the first point of the curve */

    /* for use in polygons only */
    struct fpolygon *previous;   /**< pointer to the prev poly. (may be NULL) */
    struct fpolygon *next;       /**< pointer to the next poly. (may be NULL) */
} *Fpolygon;
typedef struct fpolygons {
    char cmt[mw_cmtsize];    /**< comments                     */
    char name[mw_namesize];
                             /**< name of the set              */
    Fpolygon first;          /**< pointer to the first polygon */
} *Fpolygons;

typedef struct point_type {
    unsigned char type;  /**< Type of the point
                          *   0 : regular point
                          *   1 : point in the image's border
                          *   2 : T-junction
                          *   3 : Tau-junction
                          *   4 : X-junction
                          *   5 : Y-junction
                          */
    struct point_type *previous;   /**< pointer to the prev point */
    struct point_type *next;       /**< pointer to the next point */
} *Point_type;

typedef struct morpho_line {
    Point_curve first_point;  /**< pointer to the first point
                               *   of the morpho_line curve                 */
    Point_type first_type;    /**< pointer to the first Point_type          */
    float minvalue;           /**< min gray level value of this morpho line */
    float maxvalue;           /**< max gray level value of this morpho line */
    unsigned char open;       /**< 0 if the morpho line is closed
                               *   opened otherwise                 */
    float data;               /**< user-defined data field (saved)  */
    void *pdata;              /**< user-defined data field
                               *   pointer to something (not saved) */

    /* for use in mimage */
    struct morpho_sets *morphosets;   /**< pointer to the morpho sets */
    unsigned int num;                 /**< morpho line number         */
    struct morpho_line *previous;     /**< pointer to the prev m.l.   */
    struct morpho_line *next;         /**< pointer to the next m.l.   */
} *Morpho_line;

typedef struct fmorpho_line {
    Point_fcurve first_point;  /**< pointer to the first point
                                *   of the morpho_line curve                 */
    Point_type first_type;     /**< pointer to the first Point_type          */
    float minvalue;            /**< min gray level value of this morpho line */
    float maxvalue;            /**< max gray level value of this morpho line */
    unsigned char open;        /**< 0 if the morpho line is closed
                                *   opened otherwise                 */
    float data;                /**< user-defined data field (saved)  */
    void *pdata;               /**< user-defined data field
                                *   pointer to something (not saved) */

    /* for use in mimage */
    struct fmorpho_line *previous;     /**< pointer to the prev m.l.   */
    struct fmorpho_line *next;         /**< pointer to the next m.l.   */
} *Fmorpho_line;

typedef struct hsegment {
    int xstart;
                 /**< left  x-coordinate of the segment */
    int xend;    /**< right x-coordinate of the segment */
    int y;       /**< y-coordinate of the segment       */
    struct hsegment *previous;   /**< pointer to the prev segment */
    struct hsegment *next;       /**< pointer to the next segment */
} *Hsegment;

typedef struct morpho_set {
    unsigned int num;        /**< morpho set number                */
    Hsegment first_segment;  /**< pointer to the first segment     */
    Hsegment last_segment;   /**< pointer to the last segment      */
    float minvalue;          /**< min gray level value of this set */
    float maxvalue;          /**< max gray level value of this set */
    unsigned char stated;    /**< 1 if this m.s. has already been stated,
                              *   0 otherwise                              */
    int area;                /**< area of the set
                              *   (number of pixels belonging to this set) */
    struct morpho_sets *neighbor;   /**< pointer to a chain
                                     *   of neighbor morpho sets */
} *Morpho_set;
typedef struct morpho_sets {
    Morpho_set morphoset;           /**< pointer to the current morpho set */
    struct morpho_sets *previous;   /**< pointer to the prev morpho sets   */
    struct morpho_sets *next;       /**< pointer to the next morpho sets   */
    /* for use in mimage */
    struct morpho_line *morpholine;    /**< pointer to the morpho line */
} *Morpho_sets;

typedef struct mimage {
    char cmt[mw_cmtsize];    /**< comments                          */
    char name[mw_namesize];
                             /**< name of the set                   */
    int nrow;                /**< number of rows (dy)               */
    int ncol;                /**< number of columns (dx)            */
    float minvalue;          /**< min gray level value in the image */
    float maxvalue;          /**< max gray level value in the image */
    Morpho_line first_ml;    /**< pointer to the first morpho line  */
    Fmorpho_line first_fml;  /**< Pointer to the first morpho line  */
    Morpho_sets first_ms;    /**< Pointer to the first morpho sets  */
} *Mimage;

#include <float.h>
#define MORPHO_INFTY FLT_MAX

typedef struct color {
    unsigned char model;  /**< model of the colorimetric system   */
    float red;            /**< the red   value if model=MODEL_RGB */
    float green;          /**< the green value if model=MODEL_RGB */
    float blue;           /**< the blue  value if model=MODEL_RGB */
} Color;

typedef struct cmorpho_line {
    Point_curve first_point;  /**< pointer to the first point
                               *   of the cmorpho_line curve        */
    Point_type first_type;    /**< pointer to the first Point_type  */
    Color minvalue;           /**< min color of this cmorpho line   */
    Color maxvalue;           /**< max color of this cmorpho line   */
    unsigned char open;       /**< 0 if the cmorpho line is closed
                               *   opened otherwise                 */
    float data;               /**< user-defined data field (saved)  */
    void *pdata;              /**< user-defined data field
                               *   pointer to something (not saved) */

    /* for use in mimage */
    struct cmorpho_sets *cmorphosets;   /**< pointer to the cmorpho sets */
    unsigned int num;                   /**< cmorpho line number         */
    struct cmorpho_line *previous;      /**< pointer to the prev cm.l.   */
    struct cmorpho_line *next;          /**< pointer to the next cm.l.   */
} *Cmorpho_line;

typedef struct cfmorpho_line {
    Point_fcurve first_point;  /**< pointer to the first point
                                *   of the cmorpho_line curve        */
    Point_type first_type;     /**< pointer to the first Point_type  */
    Color minvalue;            /**< min color of this cmorpho line   */
    Color maxvalue;            /**< max color of this cmorpho line   */
    unsigned char open;        /**< 0 if the cmorpho line is closed
                                *   opened otherwise                 */
    float data;                /**< user-defined data field (saved)  */
    void *pdata;               /**< user-defined data field
                                *   pointer to something (not saved) */

    /* for use in mimage */
    struct cfmorpho_line *previous;    /**< pointer to the prev cm.l.  */
    struct cfmorpho_line *next;        /**< pointer to the next cm.l.  */
} *Cfmorpho_line;

typedef struct cmorpho_set {
    unsigned int num;        /**< cmorpho set number               */
    Hsegment first_segment;  /**< pointer to the first segment     */
    Hsegment last_segment;   /**< pointer to the last segment      */
    Color minvalue;          /**< min color of this set            */
    Color maxvalue;          /**< max color of this set            */
    unsigned char stated;    /**< 1 if this cm.s. has already been stated,
                              *   0 otherwise                               */
    int area;                /**< area of the set
                              *   (number of pixels belonging to this set)  */
    struct cmorpho_sets *neighbor;   /**< pointer to a chain
                                      *   of neighbor cmorpho sets */
} *Cmorpho_set;
typedef struct cmorpho_sets {
    Cmorpho_set cmorphoset;          /**< pointer to the current morpho set */
    struct cmorpho_sets *previous;   /**< pointer to the prev morpho sets   */
    struct cmorpho_sets *next;       /**< pointer to the next morpho sets   */
    /* for use in mimage */
    struct cmorpho_line *cmorpholine;    /**< pointer to the morpho line */
} *Cmorpho_sets;

/** morphological cimage */
typedef struct cmimage {
    char cmt[mw_cmtsize];     /**< comments                          */
    char name[mw_namesize];   /**< name of the set                   */
    int nrow;                 /**< number of rows (dy)               */
    int ncol;                 /**< number of columns (dx)            */
    Color minvalue;           /**< min color in the image            */
    Color maxvalue;           /**< max color in the image            */
    Cmorpho_line first_ml;    /**< pointer to the first cmorpho line */
    Cfmorpho_line first_fml;  /**< Pointer to the first cmorpho line */
    Cmorpho_sets first_ms;    /**< Pointer to the first cmorpho sets */
} *Cmimage;

typedef struct flist {
    int size;        /**< size (number of elements)                     */
    int max_size;    /**< currently allocated size (number of ELEMENTS) */
    int dim;         /**< dimension (number of components per element)  */

    float *values;   /**< values = size * dim array
                      *   nth element = values[n*dim+i], i=0..dim-1 */

    int data_size;
                     /**< size of data[] in bytes    */
    void *data;      /**< user defined field (saved) */
} *Flist;
typedef struct flists {
    char cmt[mw_cmtsize];    /**< comments */
    char name[mw_namesize];
                             /**< name     */

    int size;                /**< size (number of elements) */
    int max_size;            /**< currently alloc. size (number of ELEMENTS) */

    Flist *list;             /**< array of Flist             */

    int data_size;           /**< size of data[] in bytes    */
    void *data;              /**< user defined field (saved) */
} *Flists;

typedef struct dlist {
    int size;        /**< size (number of elements)                     */
    int max_size;    /**< currently allocated size (number of ELEMENTS) */
    int dim;         /**< dimension (number of components per element)  */

    double *values;   /**< values = size * dim array
                      *   nth element = values[n*dim+i], i=0..dim-1 */

    int data_size;
                     /**< size of data[] in bytes    */
    void *data;      /**< user defined field (saved) */
} *Dlist;
typedef struct dlists {
    char cmt[mw_cmtsize];    /**< comments */
    char name[mw_namesize];
                             /**< name     */

    int size;                /**< size (number of elements) */
    int max_size;            /**< currently alloc. size (number of ELEMENTS) */

    Dlist *list;             /**< array of Dlist             */

    int data_size;           /**< size of data[] in bytes    */
    void *data;              /**< user defined field (saved) */
} *Dlists;

typedef struct cmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the image                     */
    Cimage first;            /**< pointer to the first image            */
} *Cmovie;

typedef struct ccmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the movie                     */
    Ccimage first;           /**< pointer to the first image            */
} *Ccmovie;

typedef struct fmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the image                     */
    Fimage first;            /**< pointer to the first image            */
} *Fmovie;

typedef struct cfmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the image                     */
    Cfimage first;           /**< pointer to the first image            */
} *Cfmovie;

typedef struct fsignal {
    int size;        /**< number of samples                              */
    int allocsize;
                     /**< size allocated (in bytes) for the values plane */
    float *values;   /**< the samples                                    */

    float scale;     /**< scale of the signal                            */
    float shift;     /**< shifting of the signal with respect to zero    */
    float gain;      /**< gain of the signal                             */
    float sgrate;    /**< sampling rate of the signal                    */
    int bpsample;    /**< number of bits per sample for audio drivers    */

    char cmt[mw_cmtsize];    /**< comments          */
    char name[mw_namesize];
                             /**< name of the image */

    /* Defines the signifiant part of the signal : */
    int firstp;      /**< first point not aff. by left side effect */
    int lastp;       /**< last point not aff. by right side effect */
    float param;     /**< distance between two succesive uncorrelated points */
} *Fsignal;

#define mw_max_nlevel 20
#define mw_max_nvoice 50
#define mw_max_nfilter_1d 4
#define mw_max_norient 5
#define mw_max_nfilter 6

typedef struct wtrans1d {
    char cmt[mw_cmtsize];    /**< comments             */
    char name[mw_namesize];
                             /**< name of the wtrans1d */

    int type;
                /**< type of the wtrans1d performed */
    int edges;
                /**< type of the edges statments    */

    char filter_name[mw_namesize][mw_max_nfilter_1d];
     /**< define the Names (prefix) of the filters used
      * orthogonal case: 1 Fsignal "name.ri"
      * biorthogonal case: 2 Fsignal "name1.ri" & "name2.ri"
      * dyadic case: 4 Fsignal H1,G1,H2,G2
      */

    int size;
               /**< Size of the original signal which is in average[0][0] */

    Fsignal A[mw_max_nlevel + 1][mw_max_nvoice];
     /**< average or low-pass signals
      * A[l][v] is the low-pass signal at the octave #l and at the
      * voice #v, that is the signal at the scale 2^(l+v/nvoice) if
      * the wavelet is complex, A represents the modulus only
      */

    Fsignal AP[mw_max_nlevel + 1][mw_max_nvoice];
     /**< phase of the average
      * used in case of complex wavelet only
      */

    Fsignal D[mw_max_nlevel + 1][mw_max_nvoice];
     /**< detail or wavelet coefficients signals
      * D[l][v] is the difference signal at the octave #l and at
      * the voice #v. If the wavelet is complex, D represents the
      * modulus value only. Use DP to get the phase
      */

    Fsignal DP[mw_max_nlevel + 1][mw_max_nvoice];
     /**< phase of the Detail or wavelet coefficients signals
      * used in case of complex wavelet only.
      */

    int nlevel;
                  /**< number of levels (octave) for this decomposition   */
    int nvoice;
                  /**< number of voices per octave for this decomposition */

    int complex;
                  /**< 1 if the wavelet is complex that is, if P[][] is used */
    int nfilter;
                  /**< Number of filters used to compute the decomposition   */
} *Wtrans1d;

typedef struct wtrans2d {
    char cmt[mw_cmtsize];    /**< comments             */
    char name[mw_namesize];
                             /**< name of the wtrans2d */

    int type;
                /**< type of the wtrans2d performed */
    int edges;
                /**< type of the edges statments    */

    char filter_name[mw_namesize][mw_max_nfilter];
     /**< define the Names (prefix) of the filters used
      * orthogonal case: 1 Fsignal "name.ri"
      * biorthogonal case: 2 Fsignal "name1.ri" & "name2.ri"
      * if filters for edges, 1 or 2 Fimage with same prefix "name"
      * and extension "edge" instead of "ri"
      * dyadic case: 6 Fsignal H1,G1,K1,H2,G2,K2
      */

    int nrow, ncol;
                     /**< size of the original image which is in images[0][0] */

    Fimage images[mw_max_nlevel + 1][mw_max_norient + 1];
     /**< images[j][d] is the wavelet decomposition at the level #j
      * (i.e. at the scale 2^j in the dyadic case) and at the
      * direction d of the original image images[0][0]
      */

    int nlevel;
                  /**< number of levels for this decomposition             */
    int norient;
                  /**< number of orientations for this decomposition       */
    int nfilter;
                  /**< number of filters used to compute the decomposition */
} *Wtrans2d;

typedef struct wpack2d {
    char cmt[mw_cmtsize];    /**< comments            */
    char name[mw_namesize];
                             /**< name of the wpack2d */

    Fsignal signal1;     /**< impulse response of the filter 'h'         */
    Fsignal signal2;     /**< impulse response of the filter '\tilde h',
                          * for biorthogonal wavelet packet
                          */
    int level;           /**< decomposition level (calculated)           */
    Cimage tree;         /**< decomposition tree                         */
    Fimage *images;      /**< array for output images
                          * (containing the wavelet packet coefficients)
                          */
    int img_array_size;
                         /**< number of elements in * images             */

    int img_ncol;
                   /**< number of colums in the image
                    * before the decomposition
                    */
    int img_nrow;
                   /**< number of rows in the image
                    * before the decomposition
                    */
    struct wpack2d *previous;   /**< pointer to the prev wpack2d */
    struct wpack2d *next;       /**< pointer to the next wpack2d */
} *Wpack2d;

typedef struct vpoint_wmax {
    int x, y;                  /**< location of the point in the image */
    float mag[mw_max_nlevel];  /**< magnitudes at scale 2^(n+1)        */
    float arg[mw_max_nlevel];  /**< arguments at scale 2^(n+1)         */
    float argp;                /**< prediction, if any done            */
    struct vpoint_wmax *previous, *next;  /**< link to others vpoints */
} *Vpoint_wmax;

typedef struct vchain_wmax {
    int size;           /**< number of vpoint in the vchain */
    Vpoint_wmax first;  /**< beginning of the vchain        */
    struct vchain_wmax *previous, *next; /**< link to others vchains */
} *Vchain_wmax;

typedef struct vchains_wmax {
    char cmt[mw_cmtsize];
                          /**< comments        */
    char name[mw_namesize];
                          /**< name of the set */

    int size;    /**< number of vchain in the vchains */
    int ref_level;
                 /**< reference level
                  * that is, octave for which the location (x,y)
                  * corresponds really to maxima
                  * others are wavelet coefficients only
                  */
    int nlevel;
                  /**< total number of octaves in the wavelet decomp */
    int nrow, ncol;
                  /**< size of the picture where (x,y) belong         */
    struct vchain_wmax *first;
                             /**< first vchain_wmax */
} *Vchains_wmax;

typedef struct point_plane {
    short x, y;  /**< coordinates of the point */
} *Point_plane;

typedef struct shape {
    char inferior_type;
                         /**< indicates if it is extracted from a superior
                          * or inferior level set
                          */
    float value;   /**< limiting gray-level of the level set                 */
    char open;     /**< indicates if the shape meets the border of the image */
    int area;      /**< area of the shape = area of the cc of level set
                    * + areas of the holes
                    */
    char removed;
                   /**< indicates whether the shape exists or not        */
    Point_plane pixels;  /**< the array of pixels contained in the shape */
    Flist boundary;      /**< the boundary curve defining the shape      */

    struct shape *parent, *next_sibling, *child;
     /**< data to include it in a tree
      * it has a parent (the smallest containing shape),
      * children (the largest contained shapes,
      * whose first is pChild and the others are its siblings),
      * and siblings (the other children of its parent)
      */

    int data_size;
                    /**< size of data[] in bytes */
    void *data;     /**< user defined field (saved) */
} *Shape;

typedef struct shapes {
    char cmt[mw_cmtsize];    /**< comments                            */
    char name[mw_namesize];
                             /**< name of the set                     */
    int nrow;                /**< number of rows (dy) of the image    */
    int ncol;                /**< number of columns (dx) of the image */
    int interpolation;       /**< interpolation used for the level lines:
                              * - 0 nearest neighbor
                              * - 1 bilinear
                              */
    Shape the_shapes;  /**< array of the shapes.
                        * the root of the tree is at index 0
                        */

    int nb_shapes;
                    /**< the number of shapes
                     * (the size of the array the_shapes)
                     */
    Shape *smallest_shape;   /**< image giving for each pixel
                              * the smallest shape containing it
                              */
    int data_size;
                    /**< size of data[] in bytes    */
    void *data;     /**< user defined field (saved) */
} *Shapes;

typedef struct rawdata {
    int size;                   /* number of samples */
    unsigned char *data;        /* data field        */
} *Rawdata;

extern char *mw_native_ftypes[];
extern char *mw_ftypes_description[];
extern char *mw_type_conv_out[];
extern char *mw_type_conv_in[];

extern char *mwname;
extern char *mwgroup;

#define mw_orthogonal 1
#define mw_biorthogonal 2
#define mw_dyadic 3
#define mw_continuous 4

#define mw_edges_zeropad 1      /* image is Zero-padded (no statment) */
#define mw_edges_periodic 2     /* image is made periodic */
#define mw_edges_mirror 3       /* mirror effect */
#define mw_edges_wadapted 4     /* adapted wavelet on edges */

#define mw_max_nlevel 20

#define mw_sample 0             /* Low-pass images */

#define mw_horizontal 1
#define mw_vertical 2

#define mw_diagonal 3

#define mw_magnitude 3
#define mw_argument 4

#define mw_orthogonal 1
#define mw_biorthogonal 2
#define mw_dyadic 3

#define mw_edges_zeropad 1      /* image is Zero-padded (no statment) */
#define mw_edges_periodic 2     /* image is made periodic */
#define mw_edges_mirror 3       /* mirror effect */
#define mw_edges_wadapted 4     /* adapted wavelet on edges */

extern int mwdbg;

/* src/ccimage.c */
Ccimage mw_new_ccimage(void);
Ccimage mw_alloc_ccimage(Ccimage image, int nrow, int ncol);
void mw_delete_ccimage(Ccimage image);
Ccimage mw_change_ccimage(Ccimage image, int nrow, int ncol);
void mw_getdot_ccimage(Ccimage image, int x, int y, unsigned char *r,
                       unsigned char *g, unsigned char *b);
void mw_plot_ccimage(Ccimage image, int x, int y, unsigned char r,
                     unsigned char g, unsigned char b);
void mw_draw_ccimage(Ccimage image, int a0, int b0, int a1, int b1,
                     unsigned char r, unsigned char g, unsigned char b);
void mw_clear_ccimage(Ccimage image, unsigned char r, unsigned char g,
                      unsigned char b);
void mw_copy_ccimage(Ccimage in, Ccimage out);
unsigned char **mw_newtab_red_ccimage(Ccimage image);
unsigned char **mw_newtab_green_ccimage(Ccimage image);
unsigned char **mw_newtab_blue_ccimage(Ccimage image);

/* src/ccimage_io.c */
Ccimage _mw_ccimage_load_native(char *NomFic, char *Type);
short _mw_ccimage_create_native(char *NomFic, Ccimage image, char *Type);
Ccimage _mw_ccimage_load_image(char *NomFic, char *Type);
short _mw_ccimage_create_image(char *NomFic, Ccimage image, char *Type);

/* src/cfimage.c */
Cfimage mw_new_cfimage(void);
Cfimage mw_alloc_cfimage(Cfimage image, int nrow, int ncol);
void mw_delete_cfimage(Cfimage image);
Cfimage mw_change_cfimage(Cfimage image, int nrow, int ncol);
void mw_getdot_cfimage(Cfimage image, int x, int y, float *r, float *g,
                       float *b);
void mw_plot_cfimage(Cfimage image, int x, int y, float r, float g, float b);
void mw_draw_cfimage(Cfimage image, int a0, int b0, int a1, int b1, float r,
                     float g, float b);
void mw_clear_cfimage(Cfimage image, float r, float g, float b);
void mw_copy_cfimage(Cfimage in, Cfimage out);
float **mw_newtab_red_cfimage(Cfimage image);
float **mw_newtab_green_cfimage(Cfimage image);
float **mw_newtab_blue_cfimage(Cfimage image);

/* src/cfimage_io.c */
Cfimage _mw_cfimage_load_native(char *NomFic, char *Type);
short _mw_cfimage_create_native(char *NomFic, Cfimage image, char *Type);
Cfimage _mw_cfimage_load_image(char *NomFic, char *Type);
short _mw_cfimage_create_image(char *NomFic, Cfimage image, char *Type);

/* src/cimage.c */
Cimage mw_new_cimage(void);
Cimage mw_alloc_cimage(Cimage image, int nrow, int ncol);
void mw_delete_cimage(Cimage image);
Cimage mw_change_cimage(Cimage image, int nrow, int ncol);
unsigned char mw_getdot_cimage(Cimage image, int x, int y);
void mw_plot_cimage(Cimage image, int x, int y, unsigned char v);
void mw_draw_cimage(Cimage image, int a0, int b0, int a1, int b1,
                    unsigned char c);
void mw_clear_cimage(Cimage image, unsigned char v);
void mw_copy_cimage(Cimage in, Cimage out);
unsigned char **mw_newtab_gray_cimage(Cimage image);
unsigned char mw_isitbinary_cimage(Cimage image);

/* src/cimage_io.c */
Cimage _mw_cimage_load_megawave1(char *NomFic, char *Type);
short _mw_cimage_create_megawave1(char *NomFic, Cimage image, char *Type);
Cimage _mw_cimage_load_native(char *NomFic, char *Type);
short _mw_cimage_create_native(char *NomFic, Cimage image, char *Type);
Cimage _mw_cimage_load_image(char *NomFic, char *Type);
short _mw_cimage_create_image(char *NomFic, Cimage image, char *Type);

/* src/fimage.c */
Fimage mw_new_fimage(void);
Fimage mw_alloc_fimage(Fimage image, int nrow, int ncol);
void mw_delete_fimage(Fimage image);
Fimage mw_change_fimage(Fimage image, int nrow, int ncol);
float mw_getdot_fimage(Fimage image, int x, int y);
void mw_plot_fimage(Fimage image, int x, int y, float v);
void mw_draw_fimage(Fimage image, int a0, int b0, int a1, int b1, float c);
void mw_clear_fimage(Fimage image, float v);
void mw_copy_fimage(Fimage in, Fimage out);
float **mw_newtab_gray_fimage(Fimage image);

/* src/fimage_io.c */
Fimage _mw_fimage_load_megawave1(char *NomFic, char *Type);
short _mw_fimage_create_megawave1(char *NomFic, Fimage image, char *Type);
Fimage _mw_fimage_load_native(char *NomFic, char *Type);
short _mw_fimage_create_native(char *NomFic, Fimage image, char *Type);
Fimage _mw_fimage_load_image(char *NomFic, char *Type);
short _mw_fimage_create_image(char *NomFic, Fimage image, char *Type);

/* src/ccmovie.c */
Ccmovie mw_new_ccmovie(void);
void mw_delete_ccmovie(Ccmovie movie);
Ccmovie mw_change_ccmovie(Ccmovie movie);

/* src/ccmovie_io.c */
Ccmovie _mw_ccmovie_load_movie_old_format(char *NomFic, char *Type);
Ccmovie _mw_ccmovie_load_native(char *fname, char *Type);
Ccmovie _mw_ccmovie_load_movie(char *NomFic, char *Type);
short _mw_ccmovie_create_movie(char *NomFic, Ccmovie movie, char *Type);

/* src/cfmovie.c */
Cfmovie mw_new_cfmovie(void);
void mw_delete_cfmovie(Cfmovie movie);
Cfmovie mw_change_cfmovie(Cfmovie movie);

/* src/cfmovie_io.c */
Cfmovie _mw_cfmovie_load_movie_old_format(char *NomFic, char *Type);
Cfmovie _mw_cfmovie_load_native(char *fname, char *Type);
Cfmovie _mw_cfmovie_load_movie(char *NomFic, char *Type);
short _mw_cfmovie_create_movie(char *NomFic, Cfmovie movie, char *Type);

/* src/cmovie.c */
Cmovie mw_new_cmovie(void);
void mw_delete_cmovie(Cmovie movie);
Cmovie mw_change_cmovie(Cmovie movie);

/* src/cmovie_io.c */
Cmovie _mw_cmovie_load_movie_old_format(char *NomFic, char *Type);
Cmovie _mw_cmovie_load_native(char *fname, char *Type);
Cmovie _mw_cmovie_load_movie(char *NomFic, char *Type);
short _mw_cmovie_create_movie(char *NomFic, Cmovie movie, char *Type);

/* src/fmovie.c */
Fmovie mw_new_fmovie(void);
void mw_delete_fmovie(Fmovie movie);
Fmovie mw_change_fmovie(Fmovie movie);

/* src/fmovie_io.c */
Fmovie _mw_fmovie_load_movie_old_format(char *NomFic, char *Type);
Fmovie _mw_fmovie_load_movie(char *fname, char *Type);
short _mw_fmovie_create_movie(char *NomFic, Fmovie movie, char *Type);

/* src/cmimage.c */
Cmorpho_line mw_new_cmorpho_line(void);
Cmorpho_line mw_change_cmorpho_line(Cmorpho_line ll);
void mw_delete_cmorpho_line(Cmorpho_line cmorpho_line);
Cmorpho_line mw_copy_cmorpho_line(Cmorpho_line in, Cmorpho_line out);
unsigned int mw_length_cmorpho_line(Cmorpho_line cmorpho_line);
unsigned int mw_num_cmorpho_line(Cmorpho_line ml_first);
Cfmorpho_line mw_new_cfmorpho_line(void);
Cfmorpho_line mw_change_cfmorpho_line(Cfmorpho_line ll);
void mw_delete_cfmorpho_line(Cfmorpho_line cfmorpho_line);
Cfmorpho_line mw_copy_cfmorpho_line(Cfmorpho_line in, Cfmorpho_line out);
unsigned int mw_length_cfmorpho_line(Cfmorpho_line cfmorpho_line);
Cmorpho_set mw_new_cmorpho_set(void);
Cmorpho_set mw_change_cmorpho_set(Cmorpho_set is);
void mw_delete_cmorpho_set(Cmorpho_set cmorpho_set);
Cmorpho_set mw_copy_cmorpho_set(Cmorpho_set in, Cmorpho_set out);
unsigned int mw_length_cmorpho_set(Cmorpho_set cmorpho_set);
Cmorpho_sets mw_new_cmorpho_sets(void);
Cmorpho_sets mw_change_cmorpho_sets(Cmorpho_sets is);
void mw_delete_cmorpho_sets(Cmorpho_sets cmorpho_sets);
Cmorpho_sets mw_copy_cmorpho_sets(Cmorpho_sets in, Cmorpho_sets out);
unsigned int mw_length_cmorpho_sets(Cmorpho_sets cmorpho_sets);
unsigned int mw_num_cmorpho_sets(Cmorpho_sets mss_first);
void mw_cmorpho_sets_clear_stated(Cmorpho_sets mss_first);
Cmimage mw_new_cmimage(void);
Cmimage mw_change_cmimage(Cmimage mi);
void mw_delete_cmimage(Cmimage cmimage);
Cmimage mw_copy_cmimage(Cmimage in, Cmimage out);
unsigned int mw_length_ml_cmimage(Cmimage cmimage);
unsigned int mw_length_fml_cmimage(Cmimage cmimage);
unsigned int mw_length_ms_cmimage(Cmimage cmimage);

/* src/cmimage_io.c */
Cmorpho_line _mw_load_cml_mw2_cml(char *fname);
Cmorpho_line _mw_cmorpho_line_load_native(char *fname, char *type);
Cmorpho_line _mw_load_cmorpho_line(char *fname, char *type);
void _mw_write_cml_mw2_cml(FILE * fp, Cmorpho_line ll, unsigned int nml);
short _mw_create_cml_mw2_cml(char *fname, Cmorpho_line ll);
short _mw_cmorpho_line_create_native(char *fname, Cmorpho_line ll,
                                     char *Type);
short _mw_create_cmorpho_line(char *fname, Cmorpho_line ll, char *Type);
Cfmorpho_line _mw_load_cfml_mw2_cfml(char *fname);
Cfmorpho_line _mw_cfmorpho_line_load_native(char *fname, char *type);
Cfmorpho_line _mw_load_cfmorpho_line(char *fname, char *type);
void _mw_write_cfml_mw2_cfml(FILE * fp, Cfmorpho_line fll);
short _mw_create_cfml_mw2_cfml(char *fname, Cfmorpho_line fll);
short _mw_cfmorpho_line_create_native(char *fname, Cfmorpho_line ll,
                                      char *Type);
short _mw_create_cfmorpho_line(char *fname, Cfmorpho_line ll, char *Type);
Cmorpho_set _mw_read_cms_mw2_cms(char *fname, FILE * fp, int need_flipping);
Cmorpho_set _mw_load_cms_mw2_cms(char *fname);
Cmorpho_set _mw_cmorpho_set_load_native(char *fname, char *type);
Cmorpho_set _mw_load_cmorpho_set(char *fname, char *type);
void _mw_write_cms_mw2_cms(FILE * fp, Cmorpho_set is);
short _mw_create_cms_mw2_cms(char *fname, Cmorpho_set is);
short _mw_cmorpho_set_create_native(char *fname, Cmorpho_set ms, char *Type);
short _mw_create_cmorpho_set(char *fname, Cmorpho_set ms, char *Type);
Cmorpho_sets _mw_read_cmss_mw2_cmss(char *fname, FILE * fp,
                                    int need_flipping);
Cmorpho_sets _mw_load_cmss_mw2_cmss(char *fname);
Cmorpho_sets _mw_cmorpho_sets_load_native(char *fname, char *type);
Cmorpho_sets _mw_load_cmorpho_sets(char *fname, char *type);
void _mw_write_cmss_mw2_cmss(FILE * fp, Cmorpho_sets iss);
short _mw_create_cmss_mw2_cmss(char *fname, Cmorpho_sets iss);
short _mw_cmorpho_sets_create_native(char *fname, Cmorpho_sets mss,
                                     char *Type);
short _mw_create_cmorpho_sets(char *fname, Cmorpho_sets mss, char *Type);
Cmimage _mw_load_cmimage_mw2_cmimage(char *fname);
Cmimage _mw_cmimage_load_native(char *fname, char *type);
Cmimage _mw_load_cmimage(char *fname, char *type);
short _mw_create_cmimage_mw2_cmimage(char *fname, Cmimage cmimage);
short _mw_cmimage_create_native(char *fname, Cmimage cmimage, char *Type);
short _mw_create_cmimage(char *fname, Cmimage cmimage, char *Type);

/* src/mimage.c */
Point_type mw_new_point_type(void);
Point_type mw_change_point_type(Point_type point);
void mw_delete_point_type(Point_type point);
Point_type mw_copy_point_type(Point_type in, Point_type out);
Morpho_line mw_new_morpho_line(void);
Morpho_line mw_change_morpho_line(Morpho_line ll);
void mw_delete_morpho_line(Morpho_line morpho_line);
Morpho_line mw_copy_morpho_line(Morpho_line in, Morpho_line out);
unsigned int mw_length_morpho_line(Morpho_line morpho_line);
unsigned int mw_num_morpho_line(Morpho_line ml_first);
Fmorpho_line mw_new_fmorpho_line(void);
Fmorpho_line mw_change_fmorpho_line(Fmorpho_line ll);
void mw_delete_fmorpho_line(Fmorpho_line fmorpho_line);
Fmorpho_line mw_copy_fmorpho_line(Fmorpho_line in, Fmorpho_line out);
unsigned int mw_length_fmorpho_line(Fmorpho_line fmorpho_line);
Hsegment mw_new_hsegment(void);
Hsegment mw_change_hsegment(Hsegment segment);
void mw_delete_hsegment(Hsegment segment);
Morpho_set mw_new_morpho_set(void);
Morpho_set mw_change_morpho_set(Morpho_set is);
void mw_delete_morpho_set(Morpho_set morpho_set);
Morpho_set mw_copy_morpho_set(Morpho_set in, Morpho_set out);
unsigned int mw_length_morpho_set(Morpho_set morpho_set);
Morpho_sets mw_new_morpho_sets(void);
Morpho_sets mw_change_morpho_sets(Morpho_sets is);
void mw_delete_morpho_sets(Morpho_sets morpho_sets);
Morpho_sets mw_copy_morpho_sets(Morpho_sets in, Morpho_sets out);
unsigned int mw_length_morpho_sets(Morpho_sets morpho_sets);
unsigned int mw_num_morpho_sets(Morpho_sets mss_first);
void mw_morpho_sets_clear_stated(Morpho_sets mss_first);
Mimage mw_new_mimage(void);
Mimage mw_change_mimage(Mimage mi);
void mw_delete_mimage(Mimage mimage);
Mimage mw_copy_mimage(Mimage in, Mimage out);
unsigned int mw_length_ml_mimage(Mimage mimage);
unsigned int mw_length_fml_mimage(Mimage mimage);
unsigned int mw_length_ms_mimage(Mimage mimage);

/* src/mimage_io.c */
Morpho_line _mw_load_ml_mw2_ml(char *fname);
Morpho_line _mw_morpho_line_load_native(char *fname, char *type);
Morpho_line _mw_load_morpho_line(char *fname, char *type);
void _mw_write_ml_mw2_ml(FILE * fp, Morpho_line ll, unsigned int nml);
short _mw_create_ml_mw2_ml(char *fname, Morpho_line ll);
short _mw_morpho_line_create_native(char *fname, Morpho_line ll, char *Type);
short _mw_create_morpho_line(char *fname, Morpho_line ll, char *Type);
Fmorpho_line _mw_load_fml_mw2_fml(char *fname);
Fmorpho_line _mw_fmorpho_line_load_native(char *fname, char *type);
Fmorpho_line _mw_load_fmorpho_line(char *fname, char *type);
void _mw_write_fml_mw2_fml(FILE * fp, Fmorpho_line fll);
short _mw_create_fml_mw2_fml(char *fname, Fmorpho_line fll);
short _mw_fmorpho_line_create_native(char *fname, Fmorpho_line ll,
                                     char *Type);
short _mw_create_fmorpho_line(char *fname, Fmorpho_line ll, char *Type);
Morpho_set _mw_read_ms_mw2_ms(char *fname, FILE * fp, int need_flipping);
Morpho_set _mw_load_ms_mw2_ms(char *fname);
Morpho_set _mw_morpho_set_load_native(char *fname, char *type);
Morpho_set _mw_load_morpho_set(char *fname, char *type);
void _mw_write_ms_mw2_ms(FILE * fp, Morpho_set is);
short _mw_create_ms_mw2_ms(char *fname, Morpho_set is);
short _mw_morpho_set_create_native(char *fname, Morpho_set ms, char *Type);
short _mw_create_morpho_set(char *fname, Morpho_set ms, char *Type);
Morpho_sets _mw_read_mss_mw2_mss(char *fname, FILE * fp, int need_flipping);
Morpho_sets _mw_load_mss_mw2_mss(char *fname);
Morpho_sets _mw_morpho_sets_load_native(char *fname, char *type);
Morpho_sets _mw_load_morpho_sets(char *fname, char *type);
void _mw_write_mss_mw2_mss(FILE * fp, Morpho_sets iss);
short _mw_create_mss_mw2_mss(char *fname, Morpho_sets iss);
short _mw_morpho_sets_create_native(char *fname, Morpho_sets mss, char *Type);
short _mw_create_morpho_sets(char *fname, Morpho_sets mss, char *Type);
Mimage _mw_load_mimage_mw2_mimage(char *fname);
Mimage _mw_mimage_load_native(char *fname, char *type);
Mimage _mw_load_mimage(char *fname, char *type);
short _mw_create_mimage_mw2_mimage(char *fname, Mimage mimage);
short _mw_mimage_create_native(char *fname, Mimage mimage, char *Type);
short _mw_create_mimage(char *fname, Mimage mimage, char *Type);

/* src/curve.c */
Point_curve mw_new_point_curve(void);
Point_curve mw_change_point_curve(Point_curve point);
void mw_delete_point_curve(Point_curve point);
Point_curve mw_copy_point_curve(Point_curve in, Point_curve out);
Curve mw_new_curve(void);
Curve mw_change_curve(Curve cv);
void mw_delete_curve(Curve curve);
Curve mw_copy_curve(Curve in, Curve out);
unsigned int mw_length_curve(Curve curve);
Curves mw_new_curves(void);
Curves mw_change_curves(Curves cvs);
void mw_delete_curves(Curves curves);
unsigned int mw_length_curves(Curves curves);
unsigned int mw_npoints_curves(Curves curves);

/* src/curve_io.c */
Point_curve _mw_point_curve_load_native(char *fname, char *type);
Curve _mw_load_curve_mw2_curve(char *fname);
Curve _mw_curve_load_native(char *fname, char *type);
Curve _mw_load_curve(char *fname, char *type);
short _mw_create_curve_mw2_curve(char *fname, Curve cv);
short _mw_curve_create_native(char *fname, Curve cv, char *Type);
short _mw_create_curve(char *fname, Curve cv, char *Type);
Curves _mw_load_curves_mw2_curves(char *fname);
Curves _mw_curves_load_native(char *fname, char *type);
Curves _mw_load_curves(char *fname, char *type);
short _mw_create_curves_mw2_curves(char *fname, Curves cvs);
short _mw_curves_create_native(char *fname, Curves cvs, char *Type);
short _mw_create_curves(char *fname, Curves cvs, char *Type);

/* src/dcurve.c */
Point_dcurve mw_new_point_dcurve(void);
Point_dcurve mw_change_point_dcurve(Point_dcurve point);
void mw_delete_point_dcurve(Point_dcurve point);
Point_dcurve mw_copy_point_dcurve(Point_dcurve in, Point_dcurve out);
Dcurve mw_new_dcurve(void);
Dcurve mw_change_dcurve(Dcurve cv);
void mw_delete_dcurve(Dcurve dcurve);
Dcurve mw_copy_dcurve(Dcurve in, Dcurve out);
unsigned int mw_length_dcurve(Dcurve dcurve);
Dcurves mw_new_dcurves(void);
Dcurves mw_change_dcurves(Dcurves cvs);
void mw_delete_dcurves(Dcurves dcurves);
unsigned int mw_length_dcurves(Dcurves dcurves);
unsigned int mw_npoints_dcurves(Dcurves dcurves);

/* src/dcurve_io.c */
Point_dcurve _mw_point_dcurve_load_native(char *fname, char *type);
Dcurve _mw_load_dcurve_mw2_dcurve(char *fname);
Dcurve _mw_dcurve_load_native(char *fname, char *type);
Dcurve _mw_load_dcurve(char *fname, char *type);
short _mw_create_dcurve_mw2_dcurve(char *fname, Dcurve cv);
short _mw_dcurve_create_native(char *fname, Dcurve cv, char *Type);
short _mw_create_dcurve(char *fname, Dcurve cv, char *Type);
Dcurves _mw_load_dcurves_mw2_dcurves_1_00(char *fname);
Dcurves _mw_load_dcurves_mw2_dcurves(char *fname);
Dcurves _mw_dcurves_load_native(char *fname, char *type);
Dcurves _mw_load_dcurves(char *fname, char *type);
short _mw_create_dcurves_mw2_dcurves(char *fname, Dcurves cvs);
short _mw_dcurves_create_native(char *fname, Dcurves cvs, char *Type);
short _mw_create_dcurves(char *fname, Dcurves cvs, char *Type);

/* src/fcurve.c */
Point_fcurve mw_new_point_fcurve(void);
Point_fcurve mw_change_point_fcurve(Point_fcurve point);
void mw_delete_point_fcurve(Point_fcurve point);
Point_fcurve mw_copy_point_fcurve(Point_fcurve in, Point_fcurve out);
Fcurve mw_new_fcurve(void);
Fcurve mw_change_fcurve(Fcurve cv);
void mw_delete_fcurve(Fcurve fcurve);
Fcurve mw_copy_fcurve(Fcurve in, Fcurve out);
unsigned int mw_length_fcurve(Fcurve fcurve);
Fcurves mw_new_fcurves(void);
Fcurves mw_change_fcurves(Fcurves cvs);
void mw_delete_fcurves(Fcurves fcurves);
unsigned int mw_length_fcurves(Fcurves fcurves);
unsigned int mw_npoints_fcurves(Fcurves fcurves);

/* src/fcurve_io.c */
Point_fcurve _mw_point_fcurve_load_native(char *fname, char *type);
Fcurve _mw_load_fcurve_mw2_fcurve(char *fname);
Fcurve _mw_fcurve_load_native(char *fname, char *type);
Fcurve _mw_load_fcurve(char *fname, char *type);
short _mw_create_fcurve_mw2_fcurve(char *fname, Fcurve cv);
short _mw_fcurve_create_native(char *fname, Fcurve cv, char *Type);
short _mw_create_fcurve(char *fname, Fcurve cv, char *Type);
Fcurves _mw_load_fcurves_mw2_fcurves_1_00(char *fname);
Fcurves _mw_load_fcurves_mw2_fcurves(char *fname);
Fcurves _mw_fcurves_load_native(char *fname, char *type);
Fcurves _mw_load_fcurves(char *fname, char *type);
short _mw_create_fcurves_mw2_fcurves(char *fname, Fcurves cvs);
short _mw_fcurves_create_native(char *fname, Fcurves cvs, char *Type);
short _mw_create_fcurves(char *fname, Fcurves cvs, char *Type);

/* src/fpolygon.c */
Fpolygon mw_new_fpolygon(void);
Fpolygon mw_alloc_fpolygon(Fpolygon fpolygon, int nc);
Fpolygon mw_change_fpolygon(Fpolygon poly, int nc);
void mw_delete_fpolygon(Fpolygon fpolygon);
unsigned int mw_length_fpolygon(Fpolygon fpoly);
Fpolygons mw_new_fpolygons(void);
Fpolygons mw_change_fpolygons(Fpolygons poly);
void mw_delete_fpolygons(Fpolygons fpolygons);
unsigned int mw_length_fpolygons(Fpolygons fpolys);

/* src/fpolygon_io.c */
Fpolygon _mw_load_fpolygon_a_fpoly(char *fname);
Fpolygon _mw_fpolygon_load_native(char *fname, char *type);
Fpolygon _mw_load_fpolygon(char *fname, char *type);
short _mw_create_fpolygon_a_fpoly(char *fname, Fpolygon fpoly);
short _mw_fpolygon_create_native(char *fname, Fpolygon fpoly, char *Type);
short _mw_create_fpolygon(char *fname, Fpolygon fpoly, char *Type);
Fpolygons _mw_load_fpolygons_a_fpoly(char *fname);
Fpolygons _mw_fpolygons_load_native(char *fname, char *type);
Fpolygons _mw_load_fpolygons(char *fname, char *type);
short _mw_create_fpolygons_a_fpoly(char *fname, Fpolygons poly);
short _mw_fpolygons_create_native(char *fname, Fpolygons fpoly, char *Type);
short _mw_create_fpolygons(char *fname, Fpolygons fpoly, char *Type);

/* src/polygon.c */
Polygon mw_new_polygon(void);
Polygon mw_alloc_polygon(Polygon polygon, int nc);
Polygon mw_change_polygon(Polygon poly, int nc);
void mw_delete_polygon(Polygon polygon);
unsigned int mw_length_polygon(Polygon poly);
Polygons mw_new_polygons(void);
Polygons mw_change_polygons(Polygons poly);
void mw_delete_polygons(Polygons polygons);
unsigned int mw_length_polygons(Polygons polys);

/* src/polygon_io.c */
Polygon _mw_load_polygon_a_poly(char *fname);
Polygon _mw_polygon_load_native(char *fname, char *type);
Polygon _mw_load_polygon(char *fname, char *type);
short _mw_create_polygon_a_poly(char *fname, Polygon poly);
short _mw_polygon_create_native(char *fname, Polygon poly, char *Type);
short _mw_create_polygon(char *fname, Polygon poly, char *Type);
Polygons _mw_load_polygons_a_poly(char *fname);
Polygons _mw_polygons_load_native(char *fname, char *type);
Polygons _mw_load_polygons(char *fname, char *type);
short _mw_create_polygons_a_poly(char *fname, Polygons poly);
short _mw_polygons_create_native(char *fname, Polygons poly, char *Type);
short _mw_create_polygons(char *fname, Polygons poly, char *Type);

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

/* src/shape_io.c */
Shape _mw_load_mw2_shape(char *fname);
Shape _mw_load_shape(char *fname, char *Type);
void _mw_write_mw2_shape(FILE * fp, Shape sh, int iparent);
short _mw_create_mw2_shape(char *fname, Shape sh);
short _mw_create_shape(char *fname, Shape sh, char *Type);
Shapes _mw_load_mw2_shapes_1_00(char *fname);
Shapes _mw_load_mw2_shapes(char *fname);
Shapes _mw_load_shapes(char *fname, char *Type);
short _mw_create_mw2_shapes(char *fname, Shapes shs);
short _mw_create_shapes(char *fname, Shapes shs, char *Type);

/* src/fsignal.c */
Fsignal mw_new_fsignal(void);
Fsignal mw_alloc_fsignal(Fsignal signal, int N);
void mw_delete_fsignal(Fsignal signal);
Fsignal mw_change_fsignal(Fsignal signal, int N);
void mw_clear_fsignal(Fsignal signal, float v);
void mw_copy_fsignal_values(Fsignal in, Fsignal out);
void mw_copy_fsignal_header(Fsignal in, Fsignal out);
void mw_copy_fsignal(Fsignal in, Fsignal out);

/* src/fsignal_io.c */
Fsignal _mw_load_fsignal_ascii(char *fname, Fsignal signal);
short _mw_create_fsignal_ascii(char *fname, Fsignal signal);
Fsignal _mw_load_fsignal(char *fname, char *type, Fsignal signal);
short _mw_create_fsignal(char *fname, Fsignal signal, char *type);

/* src/list.c */
Flist mw_new_flist(void);
Flist mw_realloc_flist(Flist l, int n);
Flist mw_enlarge_flist(Flist l);
Flist mw_change_flist(Flist l, int max_size, int size, int dimension);
void mw_delete_flist(Flist l);
void mw_clear_flist(Flist l, float v);
Flist mw_copy_flist(Flist in, Flist out);
Flists mw_new_flists(void);
Flists mw_realloc_flists(Flists ls, int n);
Flists mw_enlarge_flists(Flists ls);
Flists mw_change_flists(Flists ls, int max_size, int size);
void mw_delete_flists(Flists ls);
Flists mw_copy_flists(Flists in, Flists out);
Dlist mw_new_dlist(void);
Dlist mw_realloc_dlist(Dlist l, int n);
Dlist mw_enlarge_dlist(Dlist l);
Dlist mw_change_dlist(Dlist l, int max_size, int size, int dimension);
void mw_delete_dlist(Dlist l);
void mw_clear_dlist(Dlist l, double v);
Dlist mw_copy_dlist(Dlist in, Dlist out);
Dlists mw_new_dlists(void);
Dlists mw_realloc_dlists(Dlists ls, int n);
Dlists mw_enlarge_dlists(Dlists ls);
Dlists mw_change_dlists(Dlists ls, int max_size, int size);
void mw_delete_dlists(Dlists ls);
Dlists mw_copy_dlists(Dlists in, Dlists out);

/* src/list_io.c */
Flist _mw_read_mw2_flist(char *fname, FILE * fp, int need_flipping);
Flist _mw_load_mw2_flist(char *fname);
Flist _mw_flist_load_native(char *fname, char *type);
Flist _mw_load_flist(char *fname, char *type);
int _mw_write_mw2_flist(FILE * fp, Flist lst);
short _mw_create_mw2_flist(char *fname, Flist lst);
short _mw_flist_create_native(char *fname, Flist lst, char *type);
short _mw_create_flist(char *fname, Flist lst, char *type);
Flists _mw_load_mw2_flists(char *fname);
Flists _mw_flists_load_native(char *fname, char *type);
Flists _mw_load_flists(char *fname, char *type);
short _mw_create_mw2_flists(char *fname, Flists lsts);
short _mw_flists_create_native(char *fname, Flists lsts, char *type);
short _mw_create_flists(char *fname, Flists lsts, char *type);
Dlist _mw_load_mw2_dlist(char *fname);
Dlist _mw_dlist_load_native(char *fname, char *type);
Dlist _mw_load_dlist(char *fname, char *type);
int _mw_write_mw2_dlist(FILE * fp, Dlist lst);
short _mw_create_mw2_dlist(char *fname, Dlist lst);
short _mw_dlist_create_native(char *fname, Dlist lst, char *type);
short _mw_create_dlist(char *fname, Dlist lst, char *type);
Dlists _mw_load_mw2_dlists(char *fname);
Dlists _mw_dlists_load_native(char *fname, char *type);
Dlists _mw_load_dlists(char *fname, char *type);
short _mw_create_mw2_dlists(char *fname, Dlists lsts);
short _mw_dlists_create_native(char *fname, Dlists lsts, char *type);
short _mw_create_dlists(char *fname, Dlists lsts, char *type);

/* src/wave_io.c */
Fsignal _mw_fsignal_load_wave_pcm(char *fname, Fsignal signal,
                                  int need_flipping);
short _mw_fsignal_create_wave_pcm(char *fname, Fsignal signal);

/* src/wmax2d.c */
Vpoint_wmax mw_new_vpoint_wmax(void);
Vpoint_wmax mw_change_vpoint_wmax(Vpoint_wmax vpoint);
void mw_delete_vpoint_wmax(Vpoint_wmax vpoint);
Vchain_wmax mw_new_vchain_wmax(void);
Vchain_wmax mw_change_vchain_wmax(Vchain_wmax vchain);
void mw_delete_vchain_wmax(Vchain_wmax vchain);
Vchains_wmax mw_new_vchains_wmax(void);
Vchains_wmax mw_change_vchains_wmax(Vchains_wmax vchains);
void mw_delete_vchains_wmax(Vchains_wmax vchains);
Vpoint_wmax mw_copy_vpoint_wmax(Vpoint_wmax vpoint1, Vpoint_wmax vpoint0);
Vchain_wmax mw_copy_vchain_wmax(Vchain_wmax vchain1, Vchain_wmax vchain0);
int mw_give_nlevel_vchain(Vchain_wmax vchain);

/* src/wmax2d_io.c */
Vchains_wmax _mw_load_vchains_wmax(char *fname);
short _mw_create_vchains_wmax(char *fname, Vchains_wmax vchains);
Vchain_wmax _mw_load_vchain_wmax(char *fname);
short _mw_create_vchain_wmax(char *fname, Vchain_wmax vchain);

/* src/wpack2d.c */
int mw_bandsize_wpack2d(int initial_size, int current_level);
Wpack2d mw_new_wpack2d(void);
int mw_checktree_wpack2d(Cimage tree);
Wpack2d mw_alloc_wpack2d(Wpack2d pack, Cimage tree, Fsignal signal1,
                         Fsignal signal2, int start_nrow, int start_ncol);
void mw_delete_wpack2d(Wpack2d pack);
Wpack2d mw_change_wpack2d(Wpack2d pack, Cimage tree, Fsignal signal1,
                          Fsignal signal2, int start_nrow, int start_ncol);
void mw_copy_wpack2d(Wpack2d in, Wpack2d out, int new_tree_size);
void mw_clear_wpack2d(Wpack2d pack, float v);
void mw_prune_wpack2d(Wpack2d in, Wpack2d out, Cimage tree);

/* src/wpack2d_io.c */
Wpack2d _mw_load_wpack2d_ascii(char *fname);
Wpack2d _mw_wpack2d_load_native(char *fname, char *type);
Wpack2d _mw_load_wpack2d(char *fname, char *type);
short _mw_create_wpack2d_ascii(char *fname, Wpack2d pack);
short _mw_wpack2d_create_native(char *fname, Wpack2d pack, char *Type);
short _mw_create_wpack2d(char *fname, Wpack2d pack, char *Type);

/* src/wtrans1d.c */
Wtrans1d mw_new_wtrans1d(void);
void *_mw_alloc_wtrans1d(Wtrans1d wtrans, int level, int voice, int size,
                         int complex, int sampling, int use_average);
void *mw_alloc_ortho_wtrans1d(Wtrans1d wtrans, int level, int size);
void *mw_alloc_biortho_wtrans1d(Wtrans1d wtrans, int level, int size);
void *mw_alloc_dyadic_wtrans1d(Wtrans1d wtrans, int level, int size);
void *mw_alloc_continuous_wtrans1d(Wtrans1d wtrans, int level, int voice,
                                   int size, int complex);
void mw_delete_wtrans1d(Wtrans1d wtrans);

/* src/wtrans1d_io.c */
Wtrans1d _mw_load_wtrans1d_header(char *fname);
short _mw_create_wtrans1d_header(char *fname, Wtrans1d wtrans);
void *_mw_wtrans1d_load_signal_wtrans(char *fname, char *type,
                                      Wtrans1d wtrans, Fsignal(*S)[50],
                                      char *Sname);
Wtrans1d _mw_wtrans1d_load_wtrans(char *fname, char *type);
void _mw_wtrans1d_create_signal_wtrans(char *fname, char *type,
                                       Wtrans1d wtrans, Fsignal(*S)[50],
                                       char *Sname);
short _mw_wtrans1d_create_wtrans(char *fname, Wtrans1d wtrans, char *type);

/* src/wtrans2d.c */
Wtrans2d mw_new_wtrans2d(void);
void *_mw_alloc_wtrans2d_norient(Wtrans2d wtrans, int level, int nrow,
                                 int ncol, int Norient, int sampling);
void *mw_alloc_ortho_wtrans2d(Wtrans2d wtrans, int level, int nrow, int ncol);
void *mw_alloc_biortho_wtrans2d(Wtrans2d wtrans, int level, int nrow,
                                int ncol);
void *mw_alloc_dyadic_wtrans2d(Wtrans2d wtrans, int level, int nrow,
                               int ncol);
void mw_delete_wtrans2d(Wtrans2d wtrans);

/* src/wtrans2d_io.c */
Wtrans2d _mw_load_wtrans2d_header(char *fname);
short _mw_create_wtrans2d_header(char *fname, Wtrans2d wtrans);
Wtrans2d _mw_wtrans2d_load_wtrans(char *fname, char *type);
short _mw_wtrans2d_create_wtrans(char *fname, Wtrans2d wtrans, char *type);

/* src/error.c */
void mwdebug(char *fmt, ...);
void mwerror(int code, int exit_code, char *fmt, ...);

/* src/rawdata.c */
Rawdata mw_new_rawdata(void);
Rawdata mw_alloc_rawdata(Rawdata rd, int newsize);
void mw_delete_rawdata(Rawdata rd);
Rawdata mw_change_rawdata(Rawdata rd, int newsize);
void mw_copy_rawdata(Rawdata in, Rawdata out);

/* src/rawdata_io.c */
Rawdata _mw_load_rawdata(char *fname);
short _mw_create_rawdata(char *fname, Rawdata rd);

/* src/bmp_io.c */
FILE *_mw_read_bmp_header(char *fname, unsigned int *nx, unsigned int *ny,
                          unsigned int *offset, unsigned int *size,
                          unsigned int *planes, unsigned int *bitcount,
                          unsigned int *compression);
Cimage _mw_cimage_load_bmp(char *file);
short _mw_cimage_create_bmp(char *file, Cimage image);
Ccimage _mw_ccimage_load_bmp(char *file);
short _mw_ccimage_create_bmp(char *file, Ccimage image);

/* src/epsf_io.c */
Cimage _mw_cimage_load_epsf(char *fname);
short _mw_cimage_create_epsf(char *fname, Cimage image);

/* src/gif_io.c */
Cimage _mw_cimage_load_gif(char *fname);
short _mw_cimage_create_gif(char *fname, Cimage image);

/* src/jpeg_io.c */
Cimage _mw_cimage_load_jpeg(char *fname);
Ccimage _mw_ccimage_load_jpeg(char *fname);
short _mw_cimage_create_jpeg(char *fname, Cimage image, char *Quality);
short _mw_ccimage_create_jpeg(char *fname, Ccimage image, char *Quality);

/* src/pgm_io.c */
int _mw_pgm_get_next_item(FILE * fp, char *comment);
Cimage _mw_cimage_load_pgma(char *file);
short _mw_cimage_create_pgma(char *file, Cimage image);
Cimage _mw_cimage_load_pgmr(char *file);
short _mw_cimage_create_pgmr(char *file, Cimage image);

/* src/pm_io.c */
Cimage _mw_cimage_load_pm(char *file);
short _mw_cimage_create_pm(char *file, Cimage image);
Fimage _mw_fimage_load_pm(char *file);
short _mw_fimage_create_pm(char *file, Fimage image);
Ccimage _mw_ccimage_load_pm(char *file);
short _mw_ccimage_create_pm(char *file, Ccimage image);
Cfimage _mw_cfimage_load_pm(char *file);
short _mw_cfimage_create_pm(char *file, Cfimage image);

/* src/ppm_io.c */
Ccimage _mw_ccimage_load_ppmr(char *file);
short _mw_ccimage_create_ppmr(char *file, Ccimage image);

/* src/ps_io.c */
Cimage _mw_cimage_load_ps(char *fname);
short _mw_cimage_create_ps(char *fname, Cimage image);

/* src/tiff_io.c */
Ccimage _mw_ccimage_load_tiff(char *fname);
Cimage _mw_cimage_load_tiff(char *fname);
short _mw_ccimage_create_tiff(char *fname, Ccimage image);
short _mw_cimage_create_tiff(char *fname, Cimage image);

/* src/ascii_file.c */
int _mw_fascii_search_string(FILE * fp, char *str);
void _mw_remove_first_spaces(char *s);
void _mw_basename(char *s, char *bname);
void _mw_dirbasename(char *s, char *dname, char *bname);
int _mw_fascii_get_field(FILE * fp, char *fname, char *field_name,
                         char *str_control, void *ptr);
int _mw_fascii_get_optional_field(FILE * fp, char *fname, char *field_name,
                                  char *str_control, void *ptr);
FILE *_mw_open_data_ascii_file(char *fname);
FILE *_mw_create_data_ascii_file(char *fname);

/* src/basic_conv.c */
void _mw_float_to_uchar(register float *ptr_float,
                        register unsigned char *ptr_uchar, int N,
                        char *data_name);
void _mw_uchar_to_float(register unsigned char *ptr_uchar,
                        register float *ptr_float, int N);
void _mw_1x24XV_to_3x8_ucharplanes(register unsigned char *ptr,
                                   register unsigned char *ptr1,
                                   register unsigned char *ptr2,
                                   register unsigned char *ptr3, int N);
void _mw_3x8_to_1x24XV_ucharplanes(register unsigned char *ptr1,
                                   register unsigned char *ptr2,
                                   register unsigned char *ptr3,
                                   register unsigned char *ptr, int N);
Cimage mw_fimage_to_cimage(Fimage image_fimage, Cimage image);
Fimage mw_cimage_to_fimage(Cimage image_cimage, Fimage image);
Ccimage mw_cfimage_to_ccimage(Cfimage image_cfimage, Ccimage image);
Cfimage mw_ccimage_to_cfimage(Ccimage image_ccimage, Cfimage image);
Fimage mw_cfimage_to_fimage(Cfimage image_cfimage, Fimage image);
Cfimage mw_fimage_to_cfimage(Fimage image_fimage, Cfimage image);
Fimage mw_ccimage_to_fimage(Ccimage image_ccimage, Fimage image);
Cimage mw_ccimage_to_cimage(Ccimage image_ccimage, Cimage image);
Ccimage mw_cimage_to_ccimage(Cimage image_cimage, Ccimage image);
Polygons mw_curves_to_polygons(Curves curves, Polygons polys);
Fpolygons mw_fcurves_to_fpolygons(Fcurves fcurves, Fpolygons fpolys);
Fcurves mw_curves_to_fcurves(Curves curves, Fcurves fcurves);
Curves mw_fcurves_to_curves(Fcurves fcurves, Curves curves);
Fcurves mw_dcurves_to_fcurves(Dcurves dcurves, Fcurves fcurves);
Dcurves mw_fcurves_to_dcurves(Fcurves fcurves, Dcurves dcurves);
Fpolygons mw_polygons_to_fpolygons(Polygons polys, Fpolygons fpolys);
Polygons mw_fpolygons_to_polygons(Fpolygons fpolys, Polygons polys);
Fcurves mw_fpolygons_to_fcurves(Fpolygons fpolys, Fcurves fcurves);
Curves mw_polygons_to_curves(Polygons polys, Curves curves);
Polygon mw_curve_to_polygon(Curve curve, Polygon poly);
Fpolygon mw_fcurve_to_fpolygon(Fcurve fcurve, Fpolygon fpoly);
Point_fcurve mw_point_curve_to_point_fcurve(Point_curve pcurve,
                                            Point_fcurve first_point);
Fcurve mw_curve_to_fcurve(Curve curve, Fcurve fcurve);
Point_curve mw_point_fcurve_to_point_curve(Point_fcurve pfcurve,
                                           Point_curve first_point);
Curve mw_fcurve_to_curve(Fcurve fcurve, Curve curve);
Point_fcurve mw_point_dcurve_to_point_fcurve(Point_dcurve pcurve,
                                             Point_fcurve first_point);
Fcurve mw_dcurve_to_fcurve(Dcurve dcurve, Fcurve fcurve);
Point_dcurve mw_point_fcurve_to_point_dcurve(Point_fcurve pfcurve,
                                             Point_dcurve first_point);
Dcurve mw_fcurve_to_dcurve(Fcurve fcurve, Dcurve dcurve);
Fpolygon mw_polygon_to_fpolygon(Polygon poly, Fpolygon fpoly);
Polygon mw_fpolygon_to_polygon(Fpolygon fpoly, Polygon poly);
Fcurve mw_fpolygon_to_fcurve(Fpolygon fpoly, Fcurve fcurve);
Curve mw_polygon_to_curve(Polygon poly, Curve curve);
Curve mw_curves_to_curve(Curves curves, Curve curve);
Curves mw_curve_to_curves(Curve curve, Curves curves);
Fcurve mw_fcurves_to_fcurve(Fcurves fcurves, Fcurve fcurve);
Fcurves mw_fcurve_to_fcurves(Fcurve fcurve, Fcurves fcurves);
Morpho_line mw_curve_to_morpho_line(Curve curve, Morpho_line ll);
Curve mw_morpho_line_to_curve(Morpho_line ll, Curve cv);
Fmorpho_line mw_fcurve_to_fmorpho_line(Fcurve fcurve, Fmorpho_line fll);
Fcurve mw_fmorpho_line_to_fcurve(Fmorpho_line fll, Fcurve cv);
Fmorpho_line mw_morpho_line_to_fmorpho_line(Morpho_line ll, Fmorpho_line fll);
Morpho_line mw_fmorpho_line_to_morpho_line(Fmorpho_line fll, Morpho_line ll);
Mimage mw_morpho_line_to_mimage(Morpho_line ll, Mimage mimage);
Morpho_line mw_mimage_to_morpho_line(Mimage mimage, Morpho_line ll);
Mimage mw_fmorpho_line_to_mimage(Fmorpho_line fll, Mimage mimage);
Fmorpho_line mw_mimage_to_fmorpho_line(Mimage mimage, Fmorpho_line fll);
Curves mw_mimage_to_curves(Mimage mimage, Curves curves);
Mimage mw_curves_to_mimage(Curves curves, Mimage mimage);
Fcurves mw_mimage_to_fcurves(Mimage mimage, Fcurves fcurves);
Mimage mw_fcurves_to_mimage(Fcurves fcurves, Mimage mimage);
Morpho_sets mw_morpho_set_to_morpho_sets(Morpho_set is,
                                         Morpho_sets morpho_sets);
Morpho_set mw_morpho_sets_to_morpho_set(Morpho_sets morpho_sets,
                                        Morpho_set is);
Mimage mw_morpho_sets_to_mimage(Morpho_sets iss, Mimage mimage);
Morpho_sets mw_mimage_to_morpho_sets(Mimage mimage, Morpho_sets iss);
Cmorpho_line mw_morpho_line_to_cmorpho_line(Morpho_line ll, Cmorpho_line cll);
Morpho_line mw_cmorpho_line_to_morpho_line(Cmorpho_line cll, Morpho_line ll);
Cfmorpho_line mw_fmorpho_line_to_cfmorpho_line(Fmorpho_line ll,
                                               Cfmorpho_line cll);
Fmorpho_line mw_cfmorpho_line_to_fmorpho_line(Cfmorpho_line cll,
                                              Fmorpho_line ll);
Cfmorpho_line mw_cmorpho_line_to_cfmorpho_line(Cmorpho_line ll,
                                               Cfmorpho_line fll);
Cmorpho_line mw_cfmorpho_line_to_cmorpho_line(Cfmorpho_line fll,
                                              Cmorpho_line ll);
Cmimage mw_cmorpho_line_to_cmimage(Cmorpho_line ll, Cmimage cmimage);
Cmorpho_line mw_cmimage_to_cmorpho_line(Cmimage cmimage, Cmorpho_line ll);
Cmorpho_line mw_curve_to_cmorpho_line(Curve curve, Cmorpho_line ll);
Curve mw_cmorpho_line_to_curve(Cmorpho_line ll, Curve cv);
Cfmorpho_line mw_fcurve_to_cfmorpho_line(Fcurve fcurve, Cfmorpho_line fll);
Fcurve mw_cfmorpho_line_to_fcurve(Cfmorpho_line fll, Fcurve cv);
Curves mw_cmimage_to_curves(Cmimage cmimage, Curves curves);
Cmimage mw_curves_to_cmimage(Curves curves, Cmimage cmimage);
Fcurves mw_cmimage_to_fcurves(Cmimage cmimage, Fcurves fcurves);
Cmimage mw_fcurves_to_cmimage(Fcurves fcurves, Cmimage cmimage);
Cmimage mw_cmorpho_sets_to_cmimage(Cmorpho_sets iss, Cmimage cmimage);
Cmorpho_sets mw_cmimage_to_cmorpho_sets(Cmimage cmimage, Cmorpho_sets iss);
Cmorpho_sets mw_cmorpho_set_to_cmorpho_sets(Cmorpho_set is,
                                            Cmorpho_sets cmorpho_sets);
Cmorpho_set mw_cmorpho_sets_to_cmorpho_set(Cmorpho_sets cmorpho_sets,
                                           Cmorpho_set is);
Cmimage mw_cfmorpho_line_to_cmimage(Cfmorpho_line fll, Cmimage cmimage);
Cfmorpho_line mw_cmimage_to_cfmorpho_line(Cmimage cmimage, Cfmorpho_line fll);
Cmimage mw_mimage_to_cmimage(Mimage mimage, Cmimage cmimage);
Mimage mw_cmimage_to_mimage(Cmimage cmimage, Mimage mimage);
Flist mw_flists_to_flist(Flists ls, Flist l);
Flists mw_flist_to_flists(Flist l, Flists ls);
Flist mw_fcurve_to_flist(Fcurve c, Flist l);
Fcurve mw_flist_to_fcurve(Flist l, Fcurve c);
Flists mw_fcurves_to_flists(Fcurves cs, Flists ls);
Fcurves mw_flists_to_fcurves(Flists ls, Fcurves cs);
Dlist mw_dlists_to_dlist(Dlists ls, Dlist l);
Dlists mw_dlist_to_dlists(Dlist l, Dlists ls);
Dlist mw_dcurve_to_dlist(Dcurve c, Dlist l);
Dcurve mw_dlist_to_dcurve(Dlist l, Dcurve c);
Dlists mw_dcurves_to_dlists(Dcurves cs, Dlists ls);
Dcurves mw_dlists_to_dcurves(Dlists ls, Dcurves cs);
Dlist mw_flist_to_dlist(Flist in, Dlist out);
Flist mw_dlist_to_flist(Dlist in, Flist out);
Dlists mw_flists_to_dlists(Flists in, Dlists out);
Flists mw_dlists_to_flists(Dlists in, Flists out);

/* src/file_type.c */
void _mw_lower_type(char type[]);
void _mw_put_range_type(char type[], char *mtype, int r);
void _mw_make_type(char type[], char type_in[], char *mtype);
void _mw_choose_type(char type[], char *type_force, char *mtype);
char *_mw_get_ftype_opt(char *ftype);
int _mw_is_of_ftype(char *in, char *type);
void _mw_make_comment(char comment[], char comment_in[]);
int _mw_get_file_type(char *fname, char *ftype, char *mtype, int *hsize,
                      float *version);

/* src/mwio.c */
long mw_fsize(FILE * fp);
void _mw_flip_image(register unsigned char *ptr, short size, short dx,
                    short dy, char flip);
FILE *_mw_write_header_file(char *fname, char *type, float IDvers);
int _search_filename(char *fname);
int _check_filename(const char *fname);
long _mw_find_pattern_in_file(FILE * fp, char *label);
int _mw_byte_ordering_is_little_endian(void);
short _mwload_cimage(char *name, char type[], char comment[], Cimage * im);
short _mwsave_cimage(char *name, char type[], char type_force[],
                     char comment[], Cimage im);
short _mwload_fimage(char *name, char type[], char comment[], Fimage * im);
short _mwsave_fimage(char *name, char type[], char type_force[],
                     char comment[], Fimage im);
short _mwload_cmovie(char *name, char type[], char comment[], Cmovie * movie);
short _mwsave_cmovie(char *name, char type[], char type_force[],
                     char comment[], Cmovie movie);
short _mwload_fmovie(char *name, char type[], char comment[], Fmovie * movie);
short _mwsave_fmovie(char *name, char type[], char type_force[],
                     char comment[], Fmovie movie);
short _mwload_ccmovie(char *name, char type[], char comment[],
                      Ccmovie * movie);
short _mwsave_ccmovie(char *name, char type[], char type_force[],
                      char comment[], Ccmovie movie);
short _mwload_cfmovie(char *name, char type[], char comment[],
                      Cfmovie * movie);
short _mwsave_cfmovie(char *name, char type[], char type_force[],
                      char comment[], Cfmovie movie);
short _mwload_curve(char *name, char type[], char comment[], Curve * cv);
short _mwsave_curve(char *name, char type[], char type_force[],
                    char comment[], Curve cv);
short _mwload_curves(char *name, char type[], char comment[], Curves * cv);
short _mwsave_curves(char *name, char type[], char type_force[],
                     char comment[], Curves cv);
short _mwload_polygon(char *name, char type[], char comment[],
                      Polygon * poly);
short _mwsave_polygon(char *name, char type[], char type_force[],
                      char comment[], Polygon poly);
short _mwload_polygons(char *name, char type[], char comment[],
                       Polygons * poly);
short _mwsave_polygons(char *name, char type[], char type_force[],
                       char comment[], Polygons poly);
short _mwload_fcurve(char *name, char type[], char comment[], Fcurve * cv);
short _mwsave_fcurve(char *name, char type[], char type_force[],
                     char comment[], Fcurve cv);
short _mwload_fcurves(char *name, char type[], char comment[], Fcurves * cv);
short _mwsave_fcurves(char *name, char type[], char type_force[],
                      char comment[], Fcurves cv);
short _mwload_fpolygon(char *name, char type[], char comment[],
                       Fpolygon * poly);
short _mwsave_fpolygon(char *name, char type[], char type_force[],
                       char comment[], Fpolygon poly);
short _mwload_fpolygons(char *name, char type[], char comment[],
                        Fpolygons * poly);
short _mwsave_fpolygons(char *name, char type[], char type_force[],
                        char comment[], Fpolygons poly);
short _mwload_fsignal(char *name, char type[], char comment[],
                      Fsignal * signal);
short _mwsave_fsignal(char *name, char type[], char type_force[],
                      char comment[], Fsignal signal);
short _mwload_wtrans1d(char *name, char type[], char comment[],
                       Wtrans1d * wtrans);
short _mwsave_wtrans1d(char *name, char type[], char type_force[],
                       char comment[], Wtrans1d wtrans);
short _mwload_wtrans2d(char *name, char type[], char comment[],
                       Wtrans2d * wtrans);
short _mwsave_wtrans2d(char *name, char type[], char type_force[],
                       char comment[], Wtrans2d wtrans);
short _mwload_vchain_wmax(char *name, char type[], char comment[],
                          Vchain_wmax * vchain);
short _mwsave_vchain_wmax(char *name, char type[], char type_force[],
                          char comment[], Vchain_wmax vchain);
short _mwload_vchains_wmax(char *name, char type[], char comment[],
                           Vchains_wmax * vchains);
short _mwsave_vchains_wmax(char *name, char type[], char type_force[],
                           char comment[], Vchains_wmax vchains);
short _mwload_ccimage(char *name, char type[], char comment[], Ccimage * im);
short _mwsave_ccimage(char *name, char type[], char type_force[],
                      char comment[], Ccimage im);
short _mwload_cfimage(char *name, char type[], char comment[], Cfimage * im);
short _mwsave_cfimage(char *name, char type[], char type_force[],
                      char comment[], Cfimage im);
short _mwload_morpho_line(char *name, char type[], char comment[],
                          Morpho_line * ll);
short _mwsave_morpho_line(char *name, char type[], char type_force[],
                          char comment[], Morpho_line ll);
short _mwload_fmorpho_line(char *name, char type[], char comment[],
                           Fmorpho_line * fll);
short _mwsave_fmorpho_line(char *name, char type[], char type_force[],
                           char comment[], Fmorpho_line fll);
short _mwload_morpho_set(char *name, char type[], char comment[],
                         Morpho_set * morpho_set);
short _mwsave_morpho_set(char *name, char type[], char type_force[],
                         char comment[], Morpho_set morpho_set);
short _mwload_morpho_sets(char *name, char type[], char comment[],
                          Morpho_sets * morpho_sets);
short _mwsave_morpho_sets(char *name, char type[], char type_force[],
                          char comment[], Morpho_sets morpho_sets);
short _mwload_mimage(char *name, char type[], char comment[],
                     Mimage * mimage);
short _mwsave_mimage(char *name, char type[], char type_force[],
                     char comment[], Mimage mimage);
short _mwload_cmorpho_line(char *name, char type[], char comment[],
                           Cmorpho_line * ll);
short _mwsave_cmorpho_line(char *name, char type[], char type_force[],
                           char comment[], Cmorpho_line ll);
short _mwload_cfmorpho_line(char *name, char type[], char comment[],
                            Cfmorpho_line * fll);
short _mwsave_cfmorpho_line(char *name, char type[], char type_force[],
                            char comment[], Cfmorpho_line fll);
short _mwload_cmorpho_set(char *name, char type[], char comment[],
                          Cmorpho_set * cmorpho_set);
short _mwsave_cmorpho_set(char *name, char type[], char type_force[],
                          char comment[], Cmorpho_set cmorpho_set);
short _mwload_cmorpho_sets(char *name, char type[], char comment[],
                           Cmorpho_sets * cmorpho_sets);
short _mwsave_cmorpho_sets(char *name, char type[], char type_force[],
                           char comment[], Cmorpho_sets cmorpho_sets);
short _mwload_cmimage(char *name, char type[], char comment[],
                      Cmimage * cmimage);
short _mwsave_cmimage(char *name, char type[], char type_force[],
                      char comment[], Cmimage cmimage);
short _mwload_shape(char *name, char type[], char comment[], Shape * shape);
short _mwsave_shape(char *name, char type[], char type_force[],
                    char comment[], Shape shape);
short _mwload_shapes(char *name, char type[], char comment[],
                     Shapes * shapes);
short _mwsave_shapes(char *name, char type[], char type_force[],
                     char comment[], Shapes shapes);
short _mwload_dcurve(char *name, char type[], char comment[], Dcurve * cv);
short _mwsave_dcurve(char *name, char type[], char type_force[],
                     char comment[], Dcurve cv);
short _mwload_dcurves(char *name, char type[], char comment[], Dcurves * cv);
short _mwsave_dcurves(char *name, char type[], char type_force[],
                      char comment[], Dcurves cv);
short _mwload_rawdata(char *name, char type[], char comment[], Rawdata * rd);
short _mwsave_rawdata(char *name, char type[], char type_force[],
                      char comment[], Rawdata rd);
short _mwload_flist(char *name, char type[], char comment[], Flist * lst);
short _mwsave_flist(char *name, char type[], char type_force[],
                    char comment[], Flist lst);
short _mwload_flists(char *name, char type[], char comment[], Flists * lsts);
short _mwsave_flists(char *name, char type[], char type_force[],
                     char comment[], Flists lsts);
short _mwload_dlist(char *name, char type[], char comment[], Dlist * lst);
short _mwsave_dlist(char *name, char type[], char type_force[],
                    char comment[], Dlist lst);
short _mwload_dlists(char *name, char type[], char comment[], Dlists * lsts);
short _mwsave_dlists(char *name, char type[], char type_force[],
                     char comment[], Dlists lsts);
short _mwload_wpack2d(char *name, char type[], char comment[],
                      Wpack2d * pack);
short _mwsave_wpack2d(char *name, char type[], char type_force[],
                      char comment[], Wpack2d pack);

/* src/type_conv.c */
void *mw_conv_internal_type(void *mwstruct, char *typein, char *typeout);
void *_mw_load_etype_to_itype(char *fname, char *typein, char *typeout,
                              char *Type);
short _mw_create_etype_from_itype(char *fname, void *mwstruct, char *typein,
                                  char *ftype);

/* src/mt19937ar.c */
void mw_srand_mt(unsigned long s);
void mw_srand_mt_array(unsigned long init_key[], int key_length);
unsigned long mw_rand32_mt(void);
double mw_drand53_mt(void);

#endif                          /* !_MW3_H_ */
