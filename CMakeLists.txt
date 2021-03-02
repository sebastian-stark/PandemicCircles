#######################################################################
#                                                                     #
# CMake script for configuration of PandemicCircles                   #
#                                                                     #
#######################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

PROJECT("pandemic_circles")

#name of executable and source files from which executable is compiled
SET(_target ${PROJECT_NAME})
SET(_target_src ${_target}.cc)

#locate GraphLib
MESSAGE(${GRAPH_LIB_DIR})
FIND_PACKAGE(graph_lib HINTS $ENV{GRAPH_LIB_DIR} ${GRAPH_LIB_DIR})
IF(NOT ${graph_lib_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of GraphLib. ***\n\n"
    "You may want to either pass a flag -DGRAPH_LIB_DIR=/path/to/GraphLib to cmake\n"
    "or set an environment variable \"GRAPH_LIB_DIR\" that contains this path."
    )
ENDIF()

#set up the target
ADD_EXECUTABLE(${_target} ${_target_src})

#link against GraphLib, include headers of GraphLib
#work around the cache of cmake
UNSET(GRAPH_LIB_LIBRARY CACHE)
FIND_LIBRARY(GRAPH_LIB_LIBRARY ${GRAPH_LIB_TARGET} HINTS ${GRAPH_LIB_INSTALL_PATH_LIBRARIES})
IF(${GRAPH_LIB_LIBRARY} MATCHES "GRAPH_LIB_LIBRARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate GraphLib ***\n")
ELSE()
  TARGET_LINK_LIBRARIES(${_target} ${GRAPH_LIB_LIBRARY})
  TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${GRAPH_LIB_INSTALL_PATH_HEADERS})
ENDIF()
