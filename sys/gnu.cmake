#
# Toolchain file for
#
# GNU compiler
#
# Notes:
#
#  * Settings here should work out of the box on Ubuntu (tested on 18.4). Other build environments
#    may need some fine tuning.
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually
#


#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O2 -funroll-all-loops"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG
  "-g -Wall -pedantic -std=f2018 -fcheck=all -ffpe-trap=invalid,zero -finit-real=nan -finit-integer='huge(1)'"
  CACHE STRING "Fortran compiler flags for Debug build")

set(Fortran_FLAGS_COVERAGE "-O0 -g --coverage")

# Use intrinsic Fortran 2008 erf/erfc functions
set(INTERNAL_ERFC CACHE BOOL 0)


#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O2 -funroll-all-loops"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_RELWITDEBINFO "-g ${C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for RelWithDebInfo build")

set(C_FLAGS_DEBUG "-g -Wall -pedantic -fbounds-check"
  CACHE STRING "C compiler flags for Debug build")

set(C_FLAGS_COVERAGE "-O0 -g --coverage")


#
# External libraries
#

# NOTE: Libraries with CMake export files (e.g. ELSI and if the HYBRID_CONFIG_METHODS variable
# contains the "Find" method also libNEGF, libMBD, ScalapackFx and MpiFx) are included by searching
# for the export file in the paths defined in the CMAKE_PREFIX_PATH **environment** variable. Make
# sure your CMAKE_PREFIX_PATH variable is set up accordingly.

# LAPACK and BLAS
# (if the BLAS library contains the LAPACK functions, set LAPACK_LIBRARY to "NONE")
#set(BLAS_LIBRARY "openblas" CACHE STRING "BLAS libraries to link")
#set(BLAS_LIBRARY_DIR "" CACHE STRING "Directories where BLAS libraries can be found")
#set(LAPACK_LIBRARY "NONE" CACHE STRING "LAPACK libraries to link")
#set(LAPACK_LIBRARY_DIR "" CACHE STRING "Directories where LAPACK libraries can be found")
