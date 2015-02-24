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




ADD_EXECUTABLE(iip_stereo_correlation_multi_win2_newversion iip_stereo_correlation_multi_win2.cpp libstereo.h libstereo.cpp cubic.h ransacFIX.h vvector.h smartparameter.h)
TARGET_LINK_LIBRARIES(iip_stereo_correlation_multi_win2_newversion library2013 fftw3f libmw)


