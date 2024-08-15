set(projectdir ${PROJECT_SOURCE_DIR})

#
# General options for all targets
#
set(fypp_flags ${FYPP_BUILD_FLAGS} ${FYPP_CONFIG_FLAGS})
list(APPEND fypp_flags -I${projectdir}/common/include -DRELEASE="'${RELEASE}'")

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

set(sources-fpp
  assert.F90
  environment.F90
  globalenv.F90
  mpifx.F90
  schedule.F90)

if(WITH_MPI)
  list(APPEND sources-fpp
    mpienv.F90)
endif()

skprogs_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_library(skprogs-common ${sources-f90} ${sources-f90-preproc})

if(WITH_MPI)
  target_link_libraries(skprogs-common MPI::MPI_Fortran)
  target_link_libraries(skprogs-common MpiFx::MpiFx)
endif()

set(moddir ${CMAKE_CURRENT_BINARY_DIR}/modfiles)
set_target_properties(skprogs-common PROPERTIES Fortran_MODULE_DIRECTORY ${moddir})
target_include_directories(skprogs-common PUBLIC
  $<BUILD_INTERFACE:${moddir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

if(BUILD_SHARED_LIBS)
  install(TARGETS skprogs-common EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
