CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(sift)

SET(SRC library.cpp library.h
        domain.cpp domain.h
        splines.cpp splines.h
        filter.cpp filter.h
        demo_lib_sift.cpp demo_lib_sift.h
        demo_sift.cpp)

ADD_EXECUTABLE(sift ${SRC})
ADD_EXECUTABLE(sift_keypoints 
        library.cpp library.h
        domain.cpp domain.h
        splines.cpp splines.h
        filter.cpp filter.h
        demo_lib_sift.cpp demo_lib_sift.h
        main_sift_keypoints.cpp)

INCLUDE_DIRECTORIES(..)

IF(UNIX)
    SET_TARGET_PROPERTIES(sift PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
    SET_TARGET_PROPERTIES(sift_keypoints PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

TARGET_LINK_LIBRARIES(sift Numerics Match IO)
TARGET_LINK_LIBRARIES(sift_keypoints Numerics IO)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/sift)
INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/sift_keypoints)
