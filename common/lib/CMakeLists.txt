file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../" projectdir)

#
# General options for all targets
#
set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${CMAKE_CURRENT_SOURCE_DIR}/../include -DRELEASE="'${RELEASE}'")

set(sources-f90
  accuracy.F90
  anglib.F90
  constants.F90
  coordtrans.F90
  eigensolver.F90
  fdiff.F90
  fifo.F90
  fifobase.F90
  fifo_real1.F90
  fifo_real2.F90
  gridgenerator.F90
  interpolation.F90
  lapack.F90
  lebedev_laikov.F90
  lebedev_laikov_f77.f
  message.F90
  partition.F90
  poisson.F90
  quadrature.F90
  sphericalharmonics.F90
  splines.F90
  taggedout.F90
  utils.F90)

add_library(skprogs-common ${sources-f90})

set(moddir ${CMAKE_CURRENT_BINARY_DIR}/modfiles)
set_target_properties(skprogs-common PROPERTIES Fortran_MODULE_DIRECTORY ${moddir})
target_include_directories(skprogs-common PUBLIC
  $<BUILD_INTERFACE:${moddir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

if(BUILD_SHARED_LIBS)
  install(TARGETS skprogs-common EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
