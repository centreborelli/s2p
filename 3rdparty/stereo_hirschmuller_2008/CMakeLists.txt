cmake_minimum_required(VERSION 2.6)

PROJECT( SGBM )

FIND_PACKAGE( OpenCV REQUIRED )
FIND_PACKAGE(PNG)
FIND_PACKAGE(JPEG)
FIND_PACKAGE(TIFF)


#### GF: HACK
add_library(IIOLIB iio.c iio.h) #Define that those files create the lib FooLib
target_link_libraries(IIOLIB ${PNG_LIBRARIES} ${JPEG_LIBRARIES} ${TIFF_LIBRARIES})
SET_TARGET_PROPERTIES(IIOLIB PROPERTIES
                             COMPILE_FLAGS "-O3 -std=c99 -funroll-loops -Wno-unused -DNDEBUG")



ADD_EXECUTABLE( ${PROJECT_NAME} sgbm.cpp )
set(CORELIBS ${OpenCV_LIBRARIES} ${OPENCV_LIB_COMPONENTS} m IIOLIB)
target_link_libraries(${PROJECT_NAME} ${CORELIBS})


