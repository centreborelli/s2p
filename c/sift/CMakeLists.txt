cmake_minimum_required (VERSION 2.8)
project (sift)

# The build type "Release" adds some optimizations
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release")
endif ()

# GCC on MacOs needs this option to use the clang assembler
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (APPLE))
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q")
endif ()

# Optimize to the current CPU and enable warnings
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR
    (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -Wall -Wextra -ftree-vectorize -funroll-loops -msse2 -mavx")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Wextra -ftree-vectorize -funroll-loops -msse2 -mavx")
endif ()

# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
else ()
    set (CMAKE_C_STANDARD 99)
endif ()

# Enable C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()

# Enable OpenMP
find_package (OpenMP)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

# Link -lm
find_library (M_LIBRARY m)
if (M_LIBRARY)
    link_libraries (${M_LIBRARY})
endif()

# Link -lrt
find_library (RT_LIBRARY rt)
if (RT_LIBRARY)
    link_libraries (${RT_LIBRARY})
endif()

set(SOURCE_FILES_SIFT_ROI
    sift_roi.cpp
    pickopt.c
    LibSift/KeyPoint.cpp
    LibSift/LibSift.cpp
    LibSift/ScaleSpace.cpp
    LibSift/Octave.cpp
    LibImages/LibImages.cpp
    LibSSE/LibSSE.cpp
    Utilities/Memory.cpp
    Utilities/Parameters.cpp
    Utilities/Time.cpp
    Utilities/Utilities.cpp)

set(SOURCE_FILES_MATCHING
    pickopt.c
    timing.c
    linalg.c
    sift_anatomy_20141201/lib_util.c
    sift_anatomy_20141201/lib_keypoint.c
    sift_anatomy_20141201/lib_matching.c
    matching.c)

add_executable(sift_roi ${SOURCE_FILES_SIFT_ROI})
add_executable(matching ${SOURCE_FILES_MATCHING})

# Link GDAL
find_package (GDAL REQUIRED)
include_directories (PUBLIC ${GDAL_INCLUDE_DIRS})
target_link_libraries (sift_roi ${GDAL_LIBRARIES})
