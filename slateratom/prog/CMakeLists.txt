set(sources-f90
  cmdargs.f90
  main.f90)

add_executable(slateratom ${sources-f90})

target_link_libraries(slateratom skprogs-slateratom)
  
install(TARGETS slateratom EXPORT skprogs-targets DESTINATION ${CMAKE_INSTALL_BINDIR})
