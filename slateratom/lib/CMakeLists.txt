set(sources-f90
  broyden.f90
  core_overlap.f90
  coulomb_hfex.f90
  coulomb_potential.f90
  density.f90
  densitymatrix.f90
  dft.f90
  diagonalizations.f90
  globals.f90
  hamiltonian.f90
  input.f90
  integration.f90
  numerical_differentiation.f90
  output.f90
  total_energy.f90
  utilities.f90
  xcfunctionals.f90
  zora_routines.f90)

add_library(skprogs-slateratom ${sources-f90})

target_link_libraries(skprogs-slateratom skprogs-common Libxc::xcf90 Libxc::xc)
target_link_libraries(skprogs-slateratom skprogs-common Libxc::xcf03 Libxc::xc)
target_link_libraries(skprogs-slateratom skprogs-common LAPACK::LAPACK)

set(moddir ${CMAKE_CURRENT_BINARY_DIR}/modfiles)
set_target_properties(skprogs-slateratom PROPERTIES Fortran_MODULE_DIRECTORY ${moddir})
target_include_directories(skprogs-slateratom PUBLIC
  $<BUILD_INTERFACE:${moddir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

if(BUILD_SHARED_LIBS)
  install(TARGETS skprogs-slateratom EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
#install(DIRECTORY ${moddir}/ DESTINATION ${CMAKE_INSTALL_MODULEDIR})
