/**
 * @file definitions.h
 *
 * common definitions for libmw
 */

#ifndef _DEFINITIONS_H
#define _DEFINITIONS_H

/*
 * DEFINITIONS
 */

/* booleans */
/* TODO : remove?? */
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

/* M_PI is an Unix (XOPEN) specification */
/* TODO : handle in modules? MW_PI? */
#ifndef M_PI
#define M_PI           3.14159265358979323846   /* pi */
#endif
#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923   /* pi/2 */
#endif
#ifndef M_PI_4
#define M_PI_4         0.78539816339744830962   /* pi/4 */
#endif

/* error levels */
/* TODO : rename MW_XX */
/* TODO : use a better message lib */
#define WARNING  0
#define ERROR    1
#define FATAL    2
#define USAGE    3
#define INTERNAL 4

/* models of the color system */
#define MODEL_RGB 0
#define MODEL_YUV 1
/* Hue, Saturation, Intensity */
#define MODEL_HSI 2
/* Hue, Saturation, Value */
/* more or less same than above but without trigonometric transform */
#define MODEL_HSV 3

/* ASCII header */
#define _MW_DATA_ASCII_FILE_HEADER "MegaWave2 - DATA ASCII file -\n"

/* wmax */
/* Preset value for an unknown argument */
#define mw_not_an_argument 1.0e9

/* Preset value for an unknown magnitude */
#define mw_not_a_magnitude -1.0

/*
 * SIZES
 */

/* max size of
 * - the megawave memory types (such as "Cimage")
 * - the megawave file types (such as "IMG")
 * - some string fields
 */
/* TODO : capital */
#define mw_mtype_size 20
#define mw_ftype_size 20
#define mw_cmtsize 255
#define mw_namesize 255

#define SIZE_OF_MW2_BIN_TYPE_ID 6

/*
 * MACROS
 */

/* flip macros : invert bytes order */

/* get flipped value */
/* invert 2 bytes data : 01 -> 10 */
#define _mw_get_flip_b2(b2) (((((unsigned short) (b2)) & 0xff00)>>8) +  \
                             ((((unsigned short) (b2)) & 0x00ff)<<8) )

/* invert 4 bytes data : 0123 -> 3210 */
#define _mw_get_flip_b4(b4) (((((unsigned long) (b4)) & 0xff000000)>>24)+ \
                             ((((unsigned long) (b4)) & 0x00ff0000)>>8) + \
                             ((((unsigned long) (b4)) & 0x0000ff00)<<8) + \
                             ((((unsigned long) (b4)) & 0x000000ff)<<24) )

/* in-place flipping */
/* invert 2 bytes data : 01 -> 10 */
#define _mw_in_flip_b2(b2)  do                  \
     {                                          \
          (b2) = _mw_get_flip_b2(b2);           \
     }                                          \
     while (0)

/* invert 4 bytes data : 0123 -> 3210 */
#define _mw_in_flip_b4(b4)  do                  \
     {                                          \
          (b4) = _mw_get_flip_b4(b4);           \
     }                                          \
     while (0)

/*
 * for float, need to perform the flip by means of an unsigned long buffer
 */
#define _mw_in_flip_float(b4) do                \
     {                                          \
          unsigned long * flip_float;           \
          flip_float = (unsigned long *)(&b4);  \
          _mw_in_flip_b4(* flip_float);         \
     }                                          \
     while (0)

/*
 * for double=(u1, u2), performs (flip(u2),flip(u1))
 * where u1,u2 are long (b4 = 32 bits).
 */
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

/*
 * STRUCTURES
 */

/* TODO : unsigned, typedef, simplify */
/** unsigned char (byte) gray level image */
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

/** floating point gray level image */
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

/* FIXME : abstract the color model */
/** unsigned char (byte) color image */
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

/** floating point color image */
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

/** point in a curve */
typedef struct point_curve {
    int x, y;
               /**< coordinates of the point */

    /* for use in curve */
    struct point_curve *previous;   /**< pointer to the prev point */
    struct point_curve *next;       /**< pointer to the next point */
} *Point_curve;

/** curve(s) */
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

/** float point in a curve */
typedef struct point_fcurve {
    float x, y;  /**< coordinates of the point */

    /* for use in curve */
    struct point_fcurve *previous;   /**< pointer to the prev point */
    struct point_fcurve *next;       /**< pointer to the next point */
} *Point_fcurve;

/** float curve(s) */
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

/** double point in a curve */
typedef struct point_dcurve {
    double x, y;  /**< coordinates of the point */

    /* for use in curve */
    struct point_dcurve *previous;   /**< pointer to the prev point */
    struct point_dcurve *next;       /**< pointer to the next point */
} *Point_dcurve;

/** double curve(s) */
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

/** polygon(s) */
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

/** float polygon(s) */
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

/** type of a point in a curve */
typedef struct point_type {
    /* TODO : enum? */
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

/** morpho_line on the discrete grid */
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

/* TODO : morphosets? num? */
/** morpho_line on a pseudo-continuous plane */
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

/** horizontal segment on the discrete grid */
typedef struct hsegment {
    int xstart;
                 /**< left  x-coordinate of the segment */
    int xend;    /**< right x-coordinate of the segment */
    int y;       /**< y-coordinate of the segment       */
    struct hsegment *previous;   /**< pointer to the prev segment */
    struct hsegment *next;       /**< pointer to the next segment */
} *Hsegment;

/** morpho_set(s) on the discrete grid */
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

/** morphological image */
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

/** color */
/* TODO : abstract the model */
/* TODO : n-channels         */
typedef struct color {
    unsigned char model;  /**< model of the colorimetric system   */
    float red;            /**< the red   value if model=MODEL_RGB */
    float green;          /**< the green value if model=MODEL_RGB */
    float blue;           /**< the blue  value if model=MODEL_RGB */
} Color;

/** cmorpho_line on the discrete grid */
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

/** cmorpho_line on a pseudo-continuous plane */
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

/** cmorpho_set(s) on the discrete grid */
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

/** float list(s) */
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
/** double list(s) */
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

/** movie of unsigned char (byte) gray level images */
typedef struct cmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the image                     */
    Cimage first;            /**< pointer to the first image            */
} *Cmovie;

/** movie of unsigned char (byte) color images */
typedef struct ccmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the movie                     */
    Ccimage first;           /**< pointer to the first image            */
} *Ccmovie;

/** movie of floating point gray level images */
typedef struct fmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the image                     */
    Fimage first;            /**< pointer to the first image            */
} *Fmovie;

/** movie of floating point color images */
typedef struct cfmovie {
    float scale;             /**< time scale of the movie (should be 1) */
    char cmt[mw_cmtsize];    /**< comments                              */
    char name[mw_namesize];
                             /**< name of the image                     */
    Cfimage first;           /**< pointer to the first image            */
} *Cfmovie;

/** floating Point Signal */
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

/** maximum Number of levels (octaves) in a wavelet decomposition  */
/** # 0 is the original image                                      */
#define mw_max_nlevel 20
/** maximum Number of voices per octave in a wavelet decomposition */
#define mw_max_nvoice 50
/** maximum Number of filter file names (prefixs only)             */
#define mw_max_nfilter_1d 4
/** maximum Number of orientations in a wavelet decomposition      */
#define mw_max_norient 5
/** maximum Number of filter file names (prefixs only)             */
#define mw_max_nfilter 6

/** 1-dimensional wavelet representation */
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

/** 2-dimensional wavelet representation */
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

/** 2-dimensional wavelet packets */
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

/** virtual maxima point in a virtual chain */
typedef struct vpoint_wmax {
    int x, y;                  /**< location of the point in the image */
    float mag[mw_max_nlevel];  /**< magnitudes at scale 2^(n+1)        */
    float arg[mw_max_nlevel];  /**< arguments at scale 2^(n+1)         */
    float argp;                /**< prediction, if any done            */
    struct vpoint_wmax *previous, *next;  /**< link to others vpoints */
} *Vpoint_wmax;

/** virtual chain of maxima points */
typedef struct vchain_wmax {
    int size;           /**< number of vpoint in the vchain */
    Vpoint_wmax first;  /**< beginning of the vchain        */
    struct vchain_wmax *previous, *next; /**< link to others vchains */
} *Vchain_wmax;

/** set of virtual chains */
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

/** point in the plane */
typedef struct point_plane {
    short x, y;  /**< coordinates of the point */
} *Point_plane;

/** shape : a connected component of a level set, with filled holes */
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

/** set of shapes (complete representation of an image) */
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

/** raw data */
typedef struct rawdata {
    int size;                   /* number of samples */
    unsigned char *data;        /* data field        */
} *Rawdata;

/*
 * GLOBALS
 */

extern char *mw_native_ftypes[];
extern char *mw_ftypes_description[];
extern char *mw_type_conv_out[];
extern char *mw_type_conv_in[];

extern char *mwname;
extern char *mwgroup;

/* Types of the wtrans1d performed */

#define mw_orthogonal 1
#define mw_biorthogonal 2
#define mw_dyadic 3
#define mw_continuous 4

/* Types of the edges statments */

#define mw_edges_zeropad 1      /* image is Zero-padded (no statment) */
#define mw_edges_periodic 2     /* image is made periodic */
#define mw_edges_mirror 3       /* mirror effect */
#define mw_edges_wadapted 4     /* adapted wavelet on edges */

/* Maximum Number of levels in a Wavelet Decomposition */
/* # 0 is the original image */
#define mw_max_nlevel 20

/* Meaning of the orientation # */
#define mw_sample 0             /* Low-pass images */

/* Orthogonal, biorthogonal and dyadic cases */
#define mw_horizontal 1
#define mw_vertical 2

/* Orthogonal/biorthogonal cases only */
#define mw_diagonal 3

/* Dyadic case (polar coordinates representation) */
#define mw_magnitude 3
#define mw_argument 4

/* Types of the wtrans2d performed */
#define mw_orthogonal 1
#define mw_biorthogonal 2
#define mw_dyadic 3

/* Types of the edges statments */
#define mw_edges_zeropad 1      /* image is Zero-padded (no statment) */
#define mw_edges_periodic 2     /* image is made periodic */
#define mw_edges_mirror 3       /* mirror effect */
#define mw_edges_wadapted 4     /* adapted wavelet on edges */

/*
 * from error.h
 */

extern int mwdbg;

#endif                          /* !_DEFINITIONS_H */
