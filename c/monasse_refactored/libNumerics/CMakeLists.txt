CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(libNumerics)

INCLUDE_DIRECTORIES(..)

SET(SRC homography.cpp homography.h
        matrix.cpp matrix.h
        numerics.cpp numerics.h
        vector.cpp)


ADD_LIBRARY(Numerics ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Numerics PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/libNumerics)
