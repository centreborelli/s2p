cmake_minimum_required(VERSION 2.8)


include_directories("mw3")
SET(SRC   fgrain.c mw3/error.c   mw3/fimage.c  mw3/fsignal.c mw3/libmw.c   mw3/list.c    mw3/shape.c flst/flst.c flst/flst_reconstruct.c flst/flst_pixels.c   flst/flstb_dualchain.c  flst/fsaddles.c)

ADD_LIBRARY(libmw ${SRC})


