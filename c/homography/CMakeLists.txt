cmake_minimum_required(VERSION 2.8)
project(homography)

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
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -march=native -ftree-vectorize -funroll-loops -msse2")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -march=native -ftree-vectorize -funroll-loops -msse2")
endif ()

# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else ()
  set (CMAKE_C_STANDARD 99)
endif ()

# Enable OpenMP
find_package (OpenMP)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

# Link GDAL
find_package (GDAL REQUIRED)
include_directories (PUBLIC ${GDAL_INCLUDE_DIRS} PUBLIC)
link_libraries (${GDAL_LIBRARIES})

set(SOURCE_FILES
    pickopt.c
    pickopt.h
    linalg.c
    linalg.h
    LibHomography/Homography.cpp
    LibHomography/Homography.h
    LibHomography/Splines.cpp
    LibHomography/Splines.h
    LibImages/LibImages.cpp
    LibImages/LibImages.h
    Utilities/Memory.cpp
    Utilities/Memory.h
    Utilities/Time.cpp
    Utilities/Time.h
    Utilities/Parameters.h
    Utilities/Parameters.cpp
    Utilities/Utilities.cpp
    Utilities/Utilities.h
    main.cpp)

add_executable(homography ${SOURCE_FILES})
