cmake_minimum_required(VERSION 2.8)
project(sift)

# The build type "Release" adds some optimizations
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release")
endif ()

# Are we using gcc?
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC on MacOs needs this option to use the clang assembler
    if (APPLE)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q")
    endif ()
    # Optimize to the current CPU and enable warnings
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -Wall -Wextra")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Wextra")
endif ()

# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
else ()
    set (CMAKE_C_STANDARD 99)
endif ()

# Enable OpenMP
find_package (OpenMP)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

# Link -lrt
find_library (RT_LIBRARY rt)
if (RT_LIBRARY)
    link_libraries (${RT_LIBRARY})
endif()

set_property(SOURCE iio.c PROPERTY COMPILE_FLAGS "-Wno-unused-variable -Wno-unused-parameter -Wno-pointer-sign -Wno-parentheses -Wno-deprecated-declarations -Wno-unused-function")

set(SOURCE_FILES_SIFT_ROI
    pickopt.c
    pickopt.h
    timing.c
    timing.h
    fancy_image.c
    fancy_image.h
    iio.h
    iio.c
    sift_anatomy_20141201/lib_description.c
    sift_anatomy_20141201/lib_description.h
    sift_anatomy_20141201/lib_discrete.c
    sift_anatomy_20141201/lib_discrete.h
    sift_anatomy_20141201/lib_keypoint.c
    sift_anatomy_20141201/lib_keypoint.h
    sift_anatomy_20141201/lib_scalespace.c
    sift_anatomy_20141201/lib_scalespace.h
    sift_anatomy_20141201/lib_sift.c
    sift_anatomy_20141201/lib_sift.h
    sift_anatomy_20141201/lib_sift_anatomy.c
    sift_anatomy_20141201/lib_sift_anatomy.h
    sift_anatomy_20141201/lib_util.c
    sift_anatomy_20141201/lib_util.h
    sift_roi.c)

set(SOURCE_FILES_MATCHING
    pickopt.c
    pickopt.h
    timing.c
    timing.h
    linalg.c
    linalg.h
    sift_anatomy_20141201/lib_util.c
    sift_anatomy_20141201/lib_util.h
    sift_anatomy_20141201/lib_keypoint.c
    sift_anatomy_20141201/lib_keypoint.h
    sift_anatomy_20141201/lib_matching.c
    sift_anatomy_20141201/lib_matching.h
    matching.c)

add_executable(sift_roi ${SOURCE_FILES_SIFT_ROI})
add_executable(matching ${SOURCE_FILES_MATCHING})

# Link LibTIFF, LibJPEG, LibPNG
find_package (TIFF REQUIRED)
find_package (JPEG REQUIRED)
find_package (PNG REQUIRED)
find_package (GDAL REQUIRED)
include_directories (PUBLIC ${TIFF_INCLUDE_DIR} PUBLIC ${JPEG_INCLUDE_DIR} PUBLIC ${PNG_INCLUDE_DIRS} PUBLIC ${GDAL_INCLUDE_DIRS}) 
target_link_libraries (sift_roi ${TIFF_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${GDAL_LIBRARIES})

# Link -lm
find_library (M_LIBRARY m)
if (M_LIBRARY)
    target_link_libraries (matching ${M_LIBRARY})
endif()
