/**
 * @file mw3-x11.h
 *
 * @author Nicolas Limare (2009)
 *
 * mw3-x11 API header
 */

#ifndef _MW_X11_H_
#define _MW_X11_H_

/* src/definitions.h */
#include <X11/Xlib.h>
#define  XK_MISCELLANY
#include <X11/keysymdef.h>

typedef struct {
    Window win;                 /*< X Window ID                             */
    int x, y;                   /*< Current Location of the Window          */
    int dx, dy;                 /*< Current Size of the Window              */
    unsigned char *pix;         /*< BitMap for the Window (client side)     */
    /*  Format of pixels is screen-dependent    */
    unsigned char *pic;         /*< BitMap - Format of pixels in 8 bits pp  */
    /*  Used only in case of not 8 bpp screens  */
    /*  (else pic = pix)                        */
    XImage *ximage;             /*< XImage structure of the Window          */
    /*  (client side)                           */
    int ix, iy;                 /*  Size of the window alloc.               */
    /*  for ximage and pix                      */
    Pixmap pixmap;              /*< Pixmap structure                        */
    /*  for buffering the graphics              */
    /*  (server side)                           */
    int px, py;                 /*< Size of the memory allocated for Pixmap */
    unsigned long event_mask;   /*< Event Mask for this Window              */
} Wframe;

#ifdef W_DEBUG_ON
#define WDEBUG(Function) (fprintf(stderr,"\n>>> Function <<<\n"))
#else
#define WDEBUG(Function)
#endif

#define WLIB_ERROR (fprintf(stderr,"Wdevice Library error: "))

extern int _W_XErrorOccured;
extern Display *_W_Display;
extern int _W_Screen;
extern int _W_Depth;
extern unsigned long _W_BlackPixel;
extern unsigned long _W_WhitePixel;
extern Colormap _W_Colormap;
extern GC _W_GC;
extern Visual *_W_Visual;
extern Cursor _W_Cursor;
extern XFontStruct *_W_Font;
extern int _W_nfcols;
extern unsigned long _W_freecols[256];
extern int _W_NumCols;
extern unsigned char _W_special_color;
extern unsigned char _W_Red[256], _W_Green[256], _W_Blue[256];
extern unsigned long _W_cols[256];
extern XColor _W_RGB[256];
extern int _W_KeyBuffer;
extern unsigned long _W_cols[256];

#define Wp_max_buttons 100      /* maximum number of buttons */
#define WP_STRSIZE 1000         /* maximum string size */

#define WP_BLACK    0
#define WP_GREY    40
#define WP_WHITE   63
#define WP_RED    164
#define WP_BLUE    68
#define WP_GREEN   84

#define WP_NULL    0
#define WP_TOGGLE  1
#define WP_INT     2
#define WP_FLOAT   3

typedef struct wp_toggle {
    char *text;                 /* text to display */
    int color;                  /* active color */
    short nbuttons;             /* number of buttons */
    short button;               /* current active button */
    char **button_text;         /* text for each button */
    int x, y;                   /* position on window (upleft corner) */
    int (*proc) (struct wp_toggle *, int);
    /* function to call when value changes */
    /* (may be NULL) */
} *Wp_toggle;

typedef struct wp_int {
    char *text;                 /* text to display */
    char *format;               /* format for int display (eg "%d") */
    int value;                  /* value */
    int strsize;                /* internal use (initialize to 0) */
    int scale;                  /* length of scale bar (0 means no bar) */
    int firstscale;             /* value of bar left border */
    int lastscale;              /* value of bar right border */
    int divscale;               /* number of bar scale divisions */
    int color;                  /* text color */
    short nbuttons;             /* number of buttons */
    char **button_text;         /* text for each button */
    int *button_inc;            /* increment for each button */
    int x, y;                   /* position on window (upleft corner) */
    int (*proc) (struct wp_int *, int);
    /* function to call when value changes */
    /* (may be NULL) */
} *Wp_int;

typedef struct wp_float {
    char *text;                 /* text to display */
    char *format;               /* format for int display (eg "%d") */
    float value;                /* value */
    int strsize;                /* internal use (initialize to 0) */
    int color;                  /* text color */
    short nbuttons;             /* number of buttons */
    char **button_text;         /* text for each button */
    float *button_inc;          /* increment for each button */
    int x, y;                   /* position on window (upleft corner) */
    int (*proc) (struct wp_float *, int);
    /* function to call when value changes */
    /* (may be NULL) */
} *Wp_float;

typedef struct wpanel {
    Wframe *window;             /* attached window */
    char state;                 /* -1 means that window should be closed */
    int nx, ny;                 /* size of bitmaps (initial window size) */
    char *type;                 /* bitmap (associated wp type) */
    void **action;              /* bitmap (pointer to wp structure) */
    short *button;              /* bitmap (associated button number) */
} *Wpanel;

extern int mwwindelay;
extern int mwrunmode;

/* src/wpanel.c */
int Wp_DrawButton(Wframe * window, int x, int y, char *str, int color);
void Wp_DrawScale(Wframe * window, int x, int y, int pos, int divisions,
                  int length, int color);
Wpanel Wp_Init(Wframe * window);
void Wp_SetButton(int type, Wpanel wp, void *b);
int Wp_handle(Wpanel wp, int event, int x, int y);
int Wp_notify(Wframe * window, void *wp);
void Wp_main_loop(Wpanel wp);

/* src/window.c */
Wframe *mw_get_window(Wframe * window, int dx, int dy, int x0, int y0,
                      char *title);
void mw_window_notify(Wframe * Win, void *param,
                      int (*proc) (Wframe *, void *));
void mw_window_main_loop(void);

/* src/wdevice.c */
int WIsAnActiveWindow(Wframe * window);
void WSetColorMap(void);
void WFlushWindow(Wframe * window);
void WFlushAreaWindow(Wframe * window, int x0, int y0, int x1, int y1);
int WColorsAvailable(void);
void WSetColorPencil(Wframe * window, int color);
void WSetForegroundColorPencil(Wframe * window);
void WSetBackgroundColorPencil(Wframe * window);
void WSetSpecialColorPencil(Wframe * window);
void WSetTypePencil(int opt);
void WDrawPoint(Wframe * window, int x, int y);
void WDrawLine(Wframe * window, int x0, int y0, int x1, int y1);
void WDrawString(Wframe * window, int x, int y, char *text);
void WDrawRectangle(Wframe * window, int x0, int y0, int x1, int y1);
void WFillRectangle(Wframe * window, int x0, int y0, int x1, int y1);
void WClearWindow(Wframe * window);
void WDestroyWdeviceWindow(Wframe * window);
void WDestroyWindow(Wframe * window);
void WMoveWindow(Wframe * window, int x, int y);
void WPutTitleWindow(Wframe * window, char *title);
void WSaveImageWindow(Wframe * window, int x, int y, int width, int height);
void WRestoreImageWindow(Wframe * window, int x, int y, int width,
                         int height);
int WLoadBitMapImage(Wframe * window, unsigned char *bitmap, int width,
                     int height);
int WLoadBitMapColorImage(Wframe * window, unsigned char *Red,
                          unsigned char *Green, unsigned char *Blue,
                          int width, int height);
void WSystemEvent(Wframe * window);
void WSetUserEvent(Wframe * window, long unsigned int user_event_mask);
int WUserEvent(Wframe * window);
int WGetStateMouse(Wframe * window, int *x, int *y, int *button_mask);
int WGetKeyboard(void);
Wframe *WNewImageWindow(void);
Wframe *WOpenImageWindow(int width, int height, int ltx, int lty,
                         char *label);
void WReOpenImageWindow(Wframe * window, int width, int height, int ltx,
                        int lty, char *label);

/* src/wdevice-misc.c */
int WX_ErrorHandler(Display * display, XErrorEvent * error);
int WX_Init(char *theDisplayName);
void WX_FreeColors(void);
void WX_AllocColors(void);
void WX_CreateXImage(Wframe * window, int dx, int dy);
void WX_AllocXImage(Wframe * window, int dx, int dy);
void WX_AllocXPixmap(Wframe * window, int dx, int dy);
void WX_Ditherize(Wframe * window, int dx, int dy);
int WX_KeyPress(XKeyEvent * event);

/* src/config.h */
#define WFONT1 "-misc-fixed-medium-r-normal-*-13-*"
#define WFONT2 "8x13"
#define WFONT3 "-*-courier-medium-r-*-*-12-*"
#define BORDER_WIDTH  2         /* width of the window border */

#define PLOT_RES_NAME   "Plot"  /* name of the plot window icon           */
#define PLOT_RES_CLASS  "Plot"  /* resource class of the plot window icon */
#define PLOT_MIN_WIDTH  50      /* minimum useful size of the plot window */
#define PLOT_MIN_HEIGHT 50
#define PLOT_MAX_WIDTH  2000    /* maximum useful size of the plot Window */
#define PLOT_MAX_HEIGHT 2000

#define IMAGE_RES_NAME   "View" /* name of the image window icon           */
#define IMAGE_RES_CLASS  "View" /* resource class of the image window icon */
#define IMAGE_MIN_WIDTH  50     /* minimum useful size of the image window */
#define IMAGE_MIN_HEIGHT 50
#define IMAGE_MAX_WIDTH  2000   /* maximum useful size of the image window */
#define IMAGE_MAX_HEIGHT 2000

#define SYSTEM_EVENT_MASK (ExposureMask                 \
                           | EnterWindowMask            \
                           | LeaveWindowMask            \
                           | StructureNotifyMask)

#define W_MS_LEFT    10         /* mouse buttons (not a mask) */
#define W_MS_RIGHT   11
#define W_MS_MIDDLE  12
#define W_MS_UP      13
#define W_MS_DOWN    14
#define W_MS_BUTTON  ButtonPressMask    /* mask for button scanning */

#define W_REPAINT  ExposureMask /* have to repaint the window    */
#define W_RESIZE   ResizeRedirectMask   /* have to resize the window     */
#define W_ENTER    EnterWindowMask      /* mouse enters the window       */
#define W_LEAVE    LeaveWindowMask      /* mouse leaves the window       */
#define W_KEYPRESS KeyPressMask /* a key has been pressed        */
#define W_DESTROY  StructureNotifyMask  /* the window has been destroyed */

#define W_COPY GXcopy
#define W_XOR  GXequiv

#endif                          /* !_MW3_X11_H_ */
