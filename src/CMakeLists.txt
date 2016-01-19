# PVBARMS version 1.0.0 


set(PVBARMS_SOURCES)
set(PVBARMS_PUBLIC_HEADERS)
#  PVBARMS.h
#)

# Copy headers to binary include dir
# Fill PVBARMS lists with headers for each subdirectory
include(include/CMakeLists.txt)
foreach(i ${SRC_HEADERS})
  list(APPEND PVBARMS_PUBLIC_HEADERS include/${i})
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/${i}" "${INC_BUILD_DIR}/include/${i}" COPYONLY)
endforeach()

#include(../include/CMakeLists.txt)
#foreach(i ${HEADERS})
#  list(APPEND PVBARMS_PUBLIC_HEADERS ../include/${i})
#  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../include/${i}" "${INC_BUILD_DIR}/../include/${i}" COPYONLY)
#endforeach()

# Fill PVBARMS lists with sources for each subdirectory

set(SOURCES
 fgmres.c
 gmres.c
 parms_comm.c
 parms_complex.c
 parms_map.c
 parms_mat.c
 parms_mat_dvcsr.c
 parms_mat_vcsr.c
 parms_mem.c
 parms_operator.c
 parms_pc_bj.c
 parms_pc.c
 parms_pc_ras.c
 parms_pc_schur.c
 parms_qsplit.c
 parms_solver.c
 parms_table.c
 parms_timer.c
 parms_vbilu_vcsr.c
 parms_vec.c
 parms_viewer.c
)

foreach(i ${SOURCES})
  list(APPEND PVBARMS_SOURCES ${i})
endforeach()

include(DDPQ/CMakeLists.txt)
foreach(i ${DDPQ_SOURCES})
  list(APPEND PVBARMS_SOURCES DDPQ/${i})
endforeach()
foreach(i ${DDPQ_PUBLIC_HEADERS})
  list(APPEND PVBARMS_PUBLIC_HEADERS DDPQ/${i})
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/DDPQ/${i}" "${INC_BUILD_DIR}/DDPQ/${i}" COPYONLY)
endforeach()



#######################################################################################################################

# Create PVBARMS library

add_library(${PVBARMS_LIBRARY} ${LIBRARY_TYPE} ${PVBARMS_SOURCES})

#######################################################################################################################
#message("support_blas")

# Link platform specific libraries
#message("print out lib path ${CMAKE_LIBRARY_PATH}")

target_link_libraries(${PVBARMS_LIBRARY} dl m)#libmath.a and libdl.a

if (SUPPORT_BLAS)
#message("print out ${BLAS_LIBRARIES} and ${CUDA_LIBRARIES}")

  target_link_libraries(${PVBARMS_LIBRARY} ${BLAS_LIBRARIES})
#message("support_blas")
endif()

if (SUPPORT_LAPACK)
#message("print out LAPACK lib ${LAPACK_LIBRARIES}")
  target_link_libraries(${PVBARMS_LIBRARY} ${LAPACK_LIBRARIES})
endif()
