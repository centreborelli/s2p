CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(rectify_mindistortion)

INCLUDE_DIRECTORIES(..)

SET(SRC rectify_mindistortion.cpp  utils.cpp utils.h rectify_mindistortion.hpp)

ADD_EXECUTABLE(rectify_mindistortion rectify_main.cpp ${SRC})
ADD_EXECUTABLE(show_epipoles show_epipoles.cpp ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(rectify_mindistortion PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
    SET_TARGET_PROPERTIES(show_epipoles PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

TARGET_LINK_LIBRARIES(rectify_mindistortion Numerics Match)
TARGET_LINK_LIBRARIES(show_epipoles IO Transform Match Numerics)

INSTALL(FILES CMakeLists.txt rectify_main.cpp show_epipoles.cpp ${SRC} DESTINATION src/minDistortionRectify)
