cmake_minimum_required (VERSION 2.6) #specify the required version of cmake
project (IIO)                        #give a name to the project/solution

# ADD LOCAL PATH TO SEARCH MODULES
list (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_PREFIX}")
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})



# SETUP A MAKE CLEAN TO CLEAN EVERYTHING
#
#
set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES 
   "iioConfig.cmake;CMakeCache.txt;CMakeFiles;Makefile;cmake_install.cmake;iio.o")


# LIBRARIES
#
#

set (LIBS "m")

# IF AVAILABLE USE OPENEXR
FIND_PACKAGE(OPENEXR) 
# OPENEXR_FOUND = OpenEXR found. 
# OPENEXR_INCLUDE_PATHS = OpenEXR include directories.
# OPENEXR_LIBRARIES = libraries that are needed to use OpenEXR.
IF (OPENEXR_FOUND) 
   include_directories(${OPENEXR_INCLUDE_PATHS})
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DI_CAN_HAS_LIBEXR")
   SET( LIBS ${LIBS} ${OPENEXR_LIBRARIES})
ENDIF()

# OTHER LIBRARIES
FIND_PACKAGE(TIFF)
IF(TIFF_FOUND)
   include_directories(${TIFF_INCLUDE_DIR})
   #  link_directories(${TIFF_LIBRARY})
   SET( LIBS ${LIBS} ${TIFF_LIBRARIES})
ENDIF()

FIND_PACKAGE(JPEG)
if(JPEG_FOUND)
   include_directories(${JPEG_INCLUDE_DIR})
   # link_directories(${JPEG_LIBRARY})
   SET( LIBS ${LIBS} ${JPEG_LIBRARIES})
ENDIF()

FIND_PACKAGE(PNG)
IF(PNG_FOUND)
   include_directories(${PNG_INCLUDE_DIR})
   #  link_directories(${PNG_LIBRARY})
   SET( LIBS ${LIBS} ${PNG_LIBRARIES})
ENDIF()


# DETECT LINUX
#
#
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   # Linux specific code flags
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DI_CAN_HAS_LINUX -D_GNU_SOURCE")
ENDIF()


# IIO
#
#
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=c99 -funroll-loops -Wno-unused -Wdeprecated-declarations -DNDEBUG")
message( "FLAGS: ${CMAKE_C_FLAGS}")
set(IIO_LIB_SRCS iio.h iio.c)       #Define a variable that contains files of the lib
add_library(IIOLIB ${IIO_LIB_SRCS}) #Define that those files create the lib FooLib
target_link_libraries(IIOLIB ${LIBS})



# generate custom iio.o
#
#
set(inc_list "")
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES) 
foreach(next_ITEM ${inc_dirs})
      set(inc_list  ${inc_list}  -I${next_ITEM})
endforeach(next_ITEM ${inc_dirs})


#HORRIBLE HACK TO CONCATENATE THE PARAMETERS
STRING(REGEX REPLACE " "  ";" FLGS "${CMAKE_C_FLAGS}")
foreach(next_ITEM ${FLGS})
      set(inc_list  ${inc_list}  ${next_ITEM})
endforeach(next_ITEM ${FLGS})

message("${inc_list}") 


ADD_CUSTOM_COMMAND(
   OUTPUT generate_iio.o 
   COMMAND ${CMAKE_C_COMPILER} ${inc_list} -c ${CMAKE_CURRENT_SOURCE_DIR}/iio.c -o ${CMAKE_CURRENT_BINARY_DIR}/iio.o
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iio.c ${CMAKE_CURRENT_SOURCE_DIR}/iio.h
   COMMENT "Generating iio.o"
   )
ADD_CUSTOM_TARGET(iio.o ALL DEPENDS generate_iio.o)




##setup Config.cmake
#
#
set(PARAMETROS_BASE_DIR "${PROJECT_SOURCE_DIR}")
set(PARAMETROS_CMAKE_DIR "${PROJECT_BINARY_DIR}")
set(IIO_LIBS ${LIBS} ${PROJECT_SOURCE_DIR}/libIIOLIB.a)
set(IIO_INCLUDE ${PROJECT_SOURCE_DIR})
message( "${PROJECT_BINARY_DIR}" )
message( "${PROJECT_SOURCE_DIR}" )
configure_file(iioConfig.cmake.in
   "${PROJECT_BINARY_DIR}/iioConfig.cmake" @ONLY)

