cmake_minimum_required(VERSION 3.15)

### Simple Matrix Operations Library (smol)
project (SMOL C)
set (CMAKE_BUILD_TYPE Debug)
set (LIB_NAME smol)

set (SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set (LIB_DIR "${CMAKE_SOURCE_DIR}/lib")
set (LIBS m)

set (SOURCE_FILES
  ${SRC_DIR}/matrix.h
  ${SRC_DIR}/matrix.c)

###
### The line of no edit.
###
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra -Wno-unused-variable")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O3")

set (LIBRARY_OUTPUT_PATH ${LIB_DIR})

add_library(${LIB_NAME} ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} ${LIBS})

enable_testing()
add_subdirectory(test)
