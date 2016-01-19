cmake_minimum_required(VERSION 2.8)
project(PVBARMS)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

enable_language (Fortran)

if (${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_BINARY_DIR})
  message(SEND_ERROR)
  message("-- ERROR: in-tree-build not allowed.")
  message("-- TRY: rm -f CMakeCache.txt; mkdir -p build; cd build; cmake ..; make")
  return()
endif()

set(PVBARMS_MAJOR_VERSION 1)
set(PVBARMS_MINOR_VERSION 0)
set(PVBARMS_BUILD_VERSION 0)
set(PVBARMS_VERSION "${PVBARMS_MAJOR_VERSION}.${PVBARMS_MINOR_VERSION}.${PVBARMS_BUILD_VERSION}")

#######################################################################################################################

set(PVBARMS_LIBRARY PVBARMS)
#message("PVBARMS ${PVBARMS_LIBRARY}")
set(BIN_BUILD_DIR ${CMAKE_BINARY_DIR}/bin)
set(INC_BUILD_DIR ${CMAKE_BINARY_DIR}/inc)
set(LIB_BUILD_DIR ${CMAKE_BINARY_DIR}/lib)
#set(MAT_DIR ${CMAKE_BINARY_DIR}/matrices)
#set(PLUGINS_BUILD_DIR ${CMAKE_BINARY_DIR}/inc/plug-ins)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_BUILD_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_BUILD_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_BUILD_DIR})

#######################################################################################################################
# Determine if library should be build shared or static
option(BUILD_SHARED "Build the PVBARMS library shared." OFF)
if(BUILD_SHARED)
  message("build shared lib")

  SET(LIBRARY_TYPE SHARED)
else()
  message("build static lib")

  SET(LIBRARY_TYPE STATIC)
endif()
#######################################################################################################################
# Find packages

# BLAS
find_package(BLAS)
if (NOT BLAS_FOUND)
  message("-- BLAS not found. Compiling WITHOUT BLAS support.")
else()
  option(SUPPORT_BLAS "Compile WITH BLAS support." ON)
  if (SUPPORT_BLAS)
    message("-- BLAS is supported")
  endif()
endif()

# LAPACK
find_package(LAPACK)
if (NOT LAPACK_FOUND)
  message("-- LAPACK not found. Compiling WITHOUT LAPACK support.")#TODO
else()
  option(SUPPORT_LAPACK "Compile WITH LAPACK support." ON)
  if (SUPPORT_LAPACK)
    message("-- LAPACK is supported")
  endif()
endif()

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
#######################################################################################################################


#######################################################################################################################
# Debug or RelWithDebInfo
set(PVARMS_BASIC_FLAGS " -DDBL -DHAS_BLAS -DVOID_POINTER_SIZE_8")

set(PVARMS_UNDERSCORE_FORTRAN_FLAGS " -DUSE_MPI -DREAL=double -DFORTRAN_UNDERSCORE")

#set(PVARMS_FORTRAN_FLAGS " -DVOID_POINTER_SIZE_8")


if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "debug")
  set(C_FLAGS "-pg -O0 -Wall -fPIC ${PVARMS_BASIC_FLAGS} ${PVARMS_UNDERSCORE_FORTRAN_FLAGS}")
  set(F_FLAGS "-pg -O0 ${PVARMS_BASIC_FLAGS}")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo" OR "${CMAKE_BUILD_TYPE}" STREQUAL "relwithdebinfo")
  set(C_FLAGS "-O3 -Wall -g -fPIC ${PVARMS_BASIC_FLAGS} ${PVARMS_UNDERSCORE_FORTRAN_FLAGS}")
  set(F_FLAGS "-O3 -g ${PVARMS_BASIC_FLAGS}")
else()
  set(C_FLAGS "-O3 -Wall -fPIC ${PVARMS_BASIC_FLAGS} ${PVARMS_UNDERSCORE_FORTRAN_FLAGS}")
  set(F_FLAGS "-O3 ${PVARMS_BASIC_FLAGS}")
endif()

#set(C_PVBARMS_FLAGS -DUSE_MPI -DREAL=double -DDBL -DHAS_BLAS)
#message("CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")

#######################################################################################################################
# Build examples
option(BUILD_EXAMPLES "Build examples." ON)

set(CMAKE_C_FLAGS "${C_FLAGS}")
set(CMAKE_Fortran_FLAGS "${F_FLAGS}")

add_subdirectory(src)
if (BUILD_EXAMPLES)
  add_subdirectory(examples/general)
endif()

#######################################################################################################################
# Testing



