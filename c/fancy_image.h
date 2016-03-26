//
// FANCY IMAGE //
// -----------
//
// A data structure for dealing with large images


// The data structure.
//
// The fields of this struct are intended to be read by the users.
struct fancy_image {
	int w;  // width
	int h;  // height
	int pd; // pixel dimension (samples per pixel)
	int no; // number of octaves

	char pad[200000]; // padding for internal details
	// Note: the current padding is large (200Kb), because
	// it contains a potentially long list of filenames.
	// If this poses problems, It can be solved very easily, just ask me.
	//  --Enric
};


///////////////
// BASIC API //
///////////////

// open an image with the desired amount of cache
// (the cache size is honored only for tiled tiffs)
//
// The default options are the following:
// 	"r,megabytes=100,octaves=0,verbose=0"
//
// Options for reading and writing an existing file
// 	"rw"
//
// Options for creating a new file
// 	"c,width=*,height=*,pd=*,type=*[,tw=*,th=*]"
//
struct fancy_image *fancy_image_open(char *filename, char *options);


// close an image and free its associated cache
// (if the option "write" was given, it may write tiles upon this call)
void fancy_image_close(struct fancy_image *f);


// obtain a sample of the image at the given point and octave
// (if the coordinates or the octave are out of range, return NAN)
float fancy_image_getsample_oct(struct fancy_image *f, int octave, int i,int j, int l);


// set a sample of the given image
// return a boolean wether it failed or not
// the file will be actually written when calling "fancy_image_close"
int fancy_image_setsample(struct fancy_image *f, int i, int j, int l, float v);


// obtain a sample of the image at the given point
// (if the point is outside the image domain, return NAN)
float fancy_image_getsample(struct fancy_image *f, int i,int j, int l);











// UTILITY API: functions that have a more comfortable interface
// All these utility functions are defined in terms of the Basic Api above

struct fancy_image *fancy_image_create(char *filename, char *fmt, ...);

void fancy_image_fill_rectangle_float_vec(float *out, int w, int h,
		struct fancy_image *f, int octave, int x0, int y0);

void fancy_image_fill_rectangle_float_split(float *out, int w, int h,
		struct fancy_image *f, int octave, int x0, int y0);

// fill-in the array "out" with "f->pd" numbers
void fancy_image_getpixel(float *out, struct fancy_image *f, int i,int j);


// fill-in the array "out" with "f->pd" numbers
void fancy_image_getpixel_oct(float *out, struct fancy_image *f, int octave,
		int i,int j);

// getpixel with automatic scale selectoin, and trilinear interpolation
// and automatic choice of the octave
void fancy_image_trilinear(float *out, struct fancy_image *f,
		double x, double y, double dx, double dy);


// leaky abstraction
int fancy_image_leak_tiff_info(int *tw, int *th, int *fmt, int *bps,
		struct fancy_image *f);



// future API
// ----------
// default option string = "r,megabytes=100,octaves=0"
//
//
//
//struct fancy_image *fancy_image_open(char *filename, char *options);
//float fancy_image_getsample(struct fancy_image *f, int i, int, int);
//float fancy_image_getsample_oct(struct fancy_image *f, int o, int i, int, int);
//void fancy_image_setsample(struct fancy_image *f, int i, int j, int l, float v);
//int fancy_image_close(struct fancy_image *f);
//
//void f(void)
//{
//	struct fancy_image f = fancy_image_open("f.tiff", "r,megabytes=100");
//	struct fancy_image f = fancy_image_open("f.tiff", "w,megabytes=100");
//}
