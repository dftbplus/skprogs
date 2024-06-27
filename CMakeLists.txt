cmake_minimum_required(VERSION 3.16)

project(SkProgs VERSION 0.3 LANGUAGES Fortran)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(SkProgsUtils)

skprogs_ensure_out_of_source_build()
skprogs_load_build_settings()
skprogs_setup_build_type()
skprogs_load_toolchain_settings()
skprogs_setup_global_compiler_flags()

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

if (Libxc_FOUND AND Libxc_VERSION VERSION_LESS "6.0.0")
  message(FATAL_ERROR "SkProgs requires libXC version >=6.0.0 to work properly.")
endif()

if (NOT Libxc_FOUND)
  message(STATUS "Libxc: No CMake export file found, trying to find with pkg-config")
  find_package(PkgConfig QUIET)
  pkg_check_modules(pc_libxc REQUIRED libxc>=6.0.0)
  pkg_check_modules(pc_libxcf03 REQUIRED libxcf03>=6.0.0)
  add_library(Libxc::xc INTERFACE IMPORTED)
  target_link_libraries(Libxc::xc INTERFACE ${pc_libxc_LINK_LIBRARIES})
  target_include_directories(Libxc::xc INTERFACE ${pc_libxc_INCLUDE_DIRS})
  add_library(Libxc::xcf03 INTERFACE IMPORTED)
  target_link_libraries(Libxc::xcf03 INTERFACE ${pc_libxcf03_LINK_LIBRARIES})
  target_include_directories(Libxc::xc INTERFACE ${pc_libxcf03_INCLUDE_DIRS})
elseif(NOT TARGET Libxc::xcf03)
  message(FATAL_ERROR "Libxc CMake export file found, but target Libxc::xcf03 is missing "
    "(maybe Libxc was built without the -DENABLE_FORTRAN=True switch?")
endif()

find_package(CustomBlas REQUIRED)
find_package(CustomLapack REQUIRED)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(PYTHON_INTERPRETER "${Python3_EXECUTABLE}")
set(PYTHON_VERSION_MAJOR_MINOR "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")

#################### Preprocessor details #####################

set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")
skprogs_add_fypp_defines(FYPP_FLAGS)
set(FYPP_FLAGS "${FYPP_FLAGS}")

set(FYPP_CONFIG_FLAGS "${FYPP_FLAGS}")
# Make sure, the line-marker option is not set
list(REMOVE_ITEM FYPP_CONFIG_FLAGS "-n")
set(FYPP_BUILD_FLAGS "${FYPP_FLAGS}" "--file-var-root=${CMAKE_SOURCE_DIR}"
  "$<IF:$<CONFIG:Debug>,-DDEBUG=1,-DDEBUG=0>")

set(PYTHON_INTERPRETER "python3" CACHE STRING
  "Python interpreter to use for installing and test python components")

#################### Add source components ####################

add_subdirectory(common)
add_subdirectory(slateratom)
add_subdirectory(sktwocnt)
add_subdirectory(sktools)

########################### Testing ###########################

string(CONFIGURE "${TEST_RUNNER_TEMPLATE}" TEST_RUNNER)
enable_testing()
add_subdirectory(test)

######################## Extra install ########################

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/export/skprogs-activate.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/skprogs-activate.sh
  @ONLY)

install(
  PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/skprogs-activate.sh"
  DESTINATION "${CMAKE_INSTALL_BINDIR}/")
