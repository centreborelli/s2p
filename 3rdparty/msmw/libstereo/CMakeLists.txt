CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
   message("OPENMP FOUND")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


# debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 ")

INCLUDE_DIRECTORIES(..)




ADD_EXECUTABLE(iip_stereo_correlation_multi_win3 iip_stereo_correlation_multi_win2.cpp libstereo.h libstereoCV.cpp )
TARGET_LINK_LIBRARIES(iip_stereo_correlation_multi_win3 library2013 fftw3f libmw)

ADD_EXECUTABLE(iip_stereo_correlation_multi_win2 iip_stereo_correlation_multi_win2.cpp libstereo.h libstereo.cpp )
TARGET_LINK_LIBRARIES(iip_stereo_correlation_multi_win2 library2013 fftw3f libmw)

ADD_EXECUTABLE(iip_stereo_correlation iip_stereo_correlation.cpp libstereo.h libstereo.cpp )
TARGET_LINK_LIBRARIES(iip_stereo_correlation library2013 fftw3f libmw)


ADD_EXECUTABLE(iip_stereo_subpixel_computation iip_stereo_subpixel_computation.cpp libstereo.h libstereo.cpp )
TARGET_LINK_LIBRARIES(iip_stereo_subpixel_computation library2013 fftw3f libmw)

