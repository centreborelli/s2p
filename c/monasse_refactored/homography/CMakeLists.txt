CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(homography)

SET(SRC homography.cpp)

ADD_EXECUTABLE(homography ${SRC})

INCLUDE_DIRECTORIES(..)

IF(UNIX)
    SET_TARGET_PROPERTIES(homography PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

TARGET_LINK_LIBRARIES(homography IO Transform Numerics IIOLIB)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/homography)
