# PVBARMS version 1.0.0 


include_directories(${INC_BUILD_DIR})

#copy inputs and matrix path files into bin directory
configure_file(inputs ${BIN_BUILD_DIR}/inputs COPYONLY)
configure_file(matfileReal ${BIN_BUILD_DIR}/matfileReal COPYONLY)
configure_file(matfileCmplx ${BIN_BUILD_DIR}/matfileCmplx COPYONLY)
configure_file(../matrices/1138_bus.mtx ${CMAKE_BINARY_DIR}/matrices/1138_bus.mtx COPYONLY)

set(EXAMPLE_SOURCES
 aux.c
 mmio.c
 skit.f
)

set(EXAMPLE_HEARDERS
 aux.h
 mmio.h
)

set(EXAMPLES
  Bdd-bin-dse
  Bdd-mtx-dse
  Bdd-HB-dse
)

add_library(aux aux.c mmio.c skit.f)

set(AUX_LIB aux)

foreach(i ${EXAMPLE_HEARDERS})
  list(APPEND PVBARMS_PUBLIC_HEADERS ${i})
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${i}" "${INC_BUILD_DIR}/${i}" COPYONLY)
endforeach()

foreach(example ${EXAMPLES})
  add_executable(${example} ${example}.c)
  target_link_libraries(${example} ${PVBARMS_LIBRARY} ${AUX_LIB} ${MPI_C_LIBRARIES})
  install(PROGRAMS ${example} DESTINATION bin)
endforeach()
