set(projectdir ${PROJECT_SOURCE_DIR})

#
# General options for all targets
#
set(fypp_flags ${FYPP_BUILD_FLAGS} ${FYPP_CONFIG_FLAGS})
list(APPEND fypp_flags -I${projectdir}/common/include -DRELEASE="'${RELEASE}'")

set(sources-f90
  cmdargs.f90
  input.f90
  output.f90)

set(sources-fpp
  main.F90)

skprogs_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(sktwocnt ${sources-f90} ${sources-f90-preproc})

target_link_libraries(skprogs-sktwocnt skprogs-common Libxc::xcf03 Libxc::xc)
target_link_libraries(sktwocnt skprogs-sktwocnt)

install(TARGETS sktwocnt EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_BINDIR})
