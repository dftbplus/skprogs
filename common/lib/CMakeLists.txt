set(sources-f90
  accuracy.F90
  constants.F90
  fifo.F90
  fifo_real1.F90
  fifo_real2.F90
  fifobase.F90
  taggedout.F90)

add_library(skprogs-common ${sources-f90})

set(moddir ${CMAKE_CURRENT_BINARY_DIR}/modfiles)
set_target_properties(skprogs-common PROPERTIES Fortran_MODULE_DIRECTORY ${moddir})
target_include_directories(skprogs-common PUBLIC
  $<BUILD_INTERFACE:${moddir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

if(BUILD_SHARED_LIBS)
  install(TARGETS skprogs-common EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
