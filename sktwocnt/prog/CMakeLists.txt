set(sources-f90
  cmdargs.f90
  input.f90
  main.f90
  output.f90)

add_executable(sktwocnt ${sources-f90})

target_link_libraries(sktwocnt skprogs-sktwocnt)

install(TARGETS sktwocnt EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_BINDIR})
