cmake_minimum_required(VERSION 2.8)
project(msmw)

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
endif ()

# Optimize to the current CPU and enable warnings
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -march=native")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -march=native")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -mavx")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -mavx")

# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else ()
    set (CMAKE_C_STANDARD 99)
endif ()

# Enable OpenMP
find_package (OpenMP)
if(OPENMP_FOUND)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Link LibTIFF, LibJPEG, LibPNG, LibFFTW
find_package (TIFF REQUIRED)
find_package (PNG REQUIRED)
find_package (JPEG REQUIRED)
find_path (FFTW_INCLUDE_DIR fftw3.h)
find_library (FFTWF_LIBRARIES NAMES fftw3f)
include_directories (PUBLIC ${TIFF_INCLUDE_DIR} PUBLIC ${JPEG_INCLUDE_DIR} PUBLIC ${PNG_INCLUDE_DIRS} PUBLIC ${FFTW_INCLUDE_DIR})
link_libraries (${TIFF_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${FFTWF_LIBRARIES})


set(SOURCE_FILES
    LibMSMW/ConnectedComponents.cpp
    LibMSMW/ConnectedComponents.h
    LibMSMW/LibMSMW.cpp
    LibMSMW/LibMSMW.h
    LibMSMW/UtilitiesMSMW.cpp
    LibMSMW/UtilitiesMSMW.h
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

add_executable(msmw ${SOURCE_FILES})
