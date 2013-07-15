CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(libTransform)

SET(SRC gauss_convol.cpp gauss_convol.h
        map_image.cpp map_image.h
        spline.h spline.cpp
        TransformSize.cpp TransformSize.h)

ADD_LIBRARY(Transform ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Transform PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

INCLUDE_DIRECTORIES(..)

TARGET_LINK_LIBRARIES(Transform Numerics)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/libTransform)
