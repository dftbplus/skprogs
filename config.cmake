#
# Global architecture independent build settings
#

#set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (Release|RelWithDebInfo|Debug|MinSizeRel)")
# CMAKE_BUILD_TYPE is commented out in order to allow for multi-configuration builds. It will
# automatically default to RelWithDebInfo if used in a single configuration build. Uncomment or
# override it only if you want a non-default single configuration build.

option(WITH_MPI "Whether SkProgs should support MPI-parallelism" FALSE)

#
# Test environment settings
#
set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of MPI processes used for testing")

# Command line used to launch the test code.
# The escaped variables (\${VARIABLE}) will be substituted by the corresponding CMake variables.
if(WITH_MPI)
  set(TEST_RUNNER_TEMPLATE "mpiexec -n \${TEST_MPI_PROCS}" CACHE STRING "How to run the tests")
else()
  set(TEST_RUNNER_TEMPLATE " " CACHE STRING "How to run the tests")
endif()

#
# Installation options
#
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_install" CACHE STRING
  "Directory to install the compiled code into")

set(INSTALL_INCLUDEDIR "skprogs" CACHE PATH
  "Name of the project specific sub-folder within the install folder for include files")

set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within the install folder for include files)")


#
# Advanced options (e.g. for developers and packagers)
#

#set(TOOLCHAIN "gnu" CACHE STRING "Prefix of the toolchain file to be read from the sys/ folder")
# Uncomment and set it if you want to override the automatic, compiler based toolchain file
# selection.

set(HYBRID_CONFIG_METHODS "Submodule;Find;Fetch" CACHE STRING
  "Configuration methods to try in order to satisfy hybrid dependencies")
