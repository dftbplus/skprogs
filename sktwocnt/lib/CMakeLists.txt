set(projectdir ${PROJECT_SOURCE_DIR})

#
# General options for all targets
#
set(fypp_flags ${FYPP_BUILD_FLAGS} ${FYPP_CONFIG_FLAGS})
list(APPEND fypp_flags -I${projectdir}/common/include -DRELEASE="'${RELEASE}'")

set(sources-f90
  bisection.f90
  gridorbital.f90
  quadrature.f90
  xcfunctionals.f90)

set(sources-fpp
  twocnt.F90)

skprogs_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_library(skprogs-sktwocnt ${sources-f90} ${sources-f90-preproc})

target_link_libraries(skprogs-sktwocnt skprogs-common)

target_link_libraries(skprogs-sktwocnt Libxc::xcf03 Libxc::xc)
target_link_libraries(skprogs-sktwocnt LAPACK::LAPACK)

if(WITH_MPI)
  target_link_libraries(skprogs-sktwocnt MPI::MPI_Fortran)
  target_link_libraries(skprogs-sktwocnt MpiFx::MpiFx)
endif()

set(moddir ${CMAKE_CURRENT_BINARY_DIR}/modfiles)
set_target_properties(skprogs-sktwocnt PROPERTIES Fortran_MODULE_DIRECTORY ${moddir})
target_include_directories(skprogs-sktwocnt PUBLIC
  $<BUILD_INTERFACE:${moddir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>)

if(BUILD_SHARED_LIBS)
  install(TARGETS skprogs-sktwocnt EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
#install(DIRECTORY ${moddir}/ DESTINATION ${CMAKE_INSTALL_MODULEDIR})
