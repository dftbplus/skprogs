set(sources-f90
  bisection.f90
  gridorbital.f90
  quadrature.f90
  twocnt.f90
  xcfunctionals.f90)

add_library(skprogs-sktwocnt ${sources-f90})

target_link_libraries(skprogs-sktwocnt skprogs-common)

target_link_libraries(skprogs-sktwocnt skprogs-common Libxc::xcf03 Libxc::xc)
target_link_libraries(skprogs-sktwocnt skprogs-common LAPACK::LAPACK)

set(moddir ${CMAKE_CURRENT_BINARY_DIR}/modfiles)
set_target_properties(skprogs-sktwocnt PROPERTIES Fortran_MODULE_DIRECTORY ${moddir})
target_include_directories(skprogs-sktwocnt PUBLIC
  $<BUILD_INTERFACE:${moddir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

if(BUILD_SHARED_LIBS)
  install(TARGETS skprogs-sktwocnt EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
#install(DIRECTORY ${moddir}/ DESTINATION ${CMAKE_INSTALL_MODULEDIR})
