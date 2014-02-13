cmake_minimum_required(VERSION 2.8.3)

PROJECT( SGBM )

# DEFAULT BUILD: RELEASE
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   set(CMAKE_BUILD_TYPE "Release")
endif()


FIND_PACKAGE(OpenCV REQUIRED core highgui calib3d)
FIND_PACKAGE(PNG)
FIND_PACKAGE(JPEG)
FIND_PACKAGE(TIFF)

# OpenMP
FIND_PACKAGE( OpenMP )
IF(OPENMP_FOUND)
   MESSAGE("OPENMP FOUND")
   SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}" )
ENDIF()


#### GF: HACK
add_library(IIOLIB iio.c iio.h) #Define that those files create the lib FooLib
target_link_libraries(IIOLIB ${PNG_LIBRARIES} ${JPEG_LIBRARIES} ${TIFF_LIBRARIES})
SET_TARGET_PROPERTIES(IIOLIB PROPERTIES
                             COMPILE_FLAGS "-O3 -std=c99 -funroll-loops -Wno-unused -DNDEBUG")


set(CORELIBS ${OpenCV_LIBRARIES} ${OPENCV_LIB_COMPONENTS} m IIOLIB)

ADD_EXECUTABLE( ${PROJECT_NAME} sgbm.cpp )
target_link_libraries(${PROJECT_NAME} ${CORELIBS})


ADD_EXECUTABLE( sgbm2 sgbm2.cpp )
target_link_libraries(sgbm2 ${CORELIBS})
SET_TARGET_PROPERTIES(sgbm2 PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
