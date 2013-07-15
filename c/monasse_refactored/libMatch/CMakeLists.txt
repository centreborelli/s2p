CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(libMatch)

SET(SRC match.cpp match.h)
ADD_LIBRARY(Match ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Match PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/libMatch)
