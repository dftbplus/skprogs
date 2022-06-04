set(sources-f90
  bisection.f90
  coordtrans.f90
  dftbxc.f90
  gridgenerator.f90
  gridorbital.f90
  interpolation.f90
  partition.f90
  quadrature.f90
  sphericalharmonics.f90
  twocnt.f90)

add_library(skprogs-sktwocnt ${sources-f90})

target_link_libraries(skprogs-sktwocnt skprogs-common)

# for a potential libxc integration:
# target_link_libraries(skprogs-sktwocnt skprogs-common Libxc::xcf90 Libxc::xc)

set(moddir ${CMAKE_CURRENT_BINARY_DIR}/modfiles)
set_target_properties(skprogs-sktwocnt PROPERTIES Fortran_MODULE_DIRECTORY ${moddir})
target_include_directories(skprogs-sktwocnt PUBLIC
  $<BUILD_INTERFACE:${moddir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

if(BUILD_SHARED_LIBS)
  install(TARGETS skprogs-sktwocnt EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
#install(DIRECTORY ${moddir}/ DESTINATION ${CMAKE_INSTALL_MODULEDIR})
