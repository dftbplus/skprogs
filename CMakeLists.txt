cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

project(SkProgs VERSION 0.1 LANGUAGES Fortran)

include(GNUInstallDirs)

set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Build type: ${default_build_type} (default single-config)")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
elseif(CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS
    "Build type: ${CMAKE_BUILD_TYPE} (manually selected single-config)")
else()
  message(STATUS "Build type: Multi-Config (build type selected at the build step)")
endif()


#################### External dependencies ####################

find_package(Libxc QUIET)
if (NOT Libxc_FOUND)
  message(STATUS "Libxc: No CMake export file found, trying to find with pkg-config")
  find_package(PkgConfig QUIET)
  pkg_check_modules(pc_libxc REQUIRED libxc)
  pkg_check_modules(pc_libxcf90 REQUIRED libxcf90)
  add_library(Libxc::xc INTERFACE IMPORTED)
  target_link_libraries(Libxc::xc INTERFACE ${pc_libxc_LINK_LIBRARIES})
  target_include_directories(Libxc::xc INTERFACE ${pc_libxc_INCLUDE_DIRS})
  add_library(Libxc::xcf90 INTERFACE IMPORTED)
  target_link_libraries(Libxc::xcf90 INTERFACE ${pc_libxcf90_LINK_LIBRARIES})
  target_include_directories(Libxc::xc INTERFACE ${pc_libxcf90_INCLUDE_DIRS})
elseif(NOT TARGET Libxc::xcf90)
  message(FATAL_ERROR "Libxc CMake export file found, but target Libxc::xcf90 is missing "
    "(maybe Libxc was built without the -DENABLE_FORTRAN=True switch?")
endif()

find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(PYTHON_INTERPRETER "${Python3_EXECUTABLE}")
set(PYTHON_VERSION_MAJOR_MINOR "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")

find_package(CustomBlas REQUIRED)
find_package(CustomLapack REQUIRED)

#################### Add source components ####################

add_subdirectory(common/lib)
add_subdirectory(slateratom)
add_subdirectory(sktwocnt)
add_subdirectory(sktools)

#################### Extra install ####################

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/export/skprogs-activate.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/skprogs-activate.sh
  @ONLY)

install(
  PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/skprogs-activate.sh"
  DESTINATION "${CMAKE_INSTALL_BINDIR}/")

