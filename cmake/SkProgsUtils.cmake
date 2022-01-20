include(FetchContent)

# Stops the code if the source and the build folders are identical.
#
function(skprogs_ensure_out_of_source_build)

  get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build SkProgs from its source folder. Please, create a \
separate build directory and invoke CMake from that directory. See the INSTALL.rst file for \
detailed build instructions.")
  endif()

endfunction()


# Makes sure, that a compiler has been already defined for a given language
#
# Args:
#     language [in]: The language to look at.
#
function(skprogs_ensure_compiler_def language)

  if(NOT DEFINED CMAKE_${language}_COMPILER)
    message(FATAL_ERROR "Undefined ${language} compiler. The automatic detection of compilers, \
flags and libraries is disabled. You must provide configuration parameters explicitely (e.g. in a \
toolchain file). See the INSTALL.rst file for detailed instructions.")
  endif()

endfunction()


# Loads global build settings (either from config.cmake or from user defined file)
#
macro (skprogs_load_build_settings)

  if(NOT DEFINED BUILD_CONFIG_FILE)
    if(DEFINED ENV{SKPROGS_BUILD_CONFIG_FILE}
        AND NOT "$ENV{SKPROGS_BUILD_CONFIG_FILE}" STREQUAL "")
      set(BUILD_CONFIG_FILE "$ENV{SKPROGS_BUILD_CONFIG_FILE}")
    else()
      set(BUILD_CONFIG_FILE "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake")
    endif()
  endif()
  message(STATUS "Reading global build config file: ${BUILD_CONFIG_FILE}")
  include(${BUILD_CONFIG_FILE})

endmacro()


# Sets up the build type.
function (skprogs_setup_build_type)

  set(default_build_type "RelWithDebInfo")
  get_property(_multiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if(_multiConfig)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;Coverage")
    message(STATUS "Build type: Multi-Config (build type selected at the build step)")
  else()
    if(NOT CMAKE_BUILD_TYPE)
      message(STATUS "Build type: ${default_build_type} (default single-config)")
      set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Build type" FORCE)
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY HELPSTRING "Choose the type of build")
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo" "Coverage")
    else()
      message(STATUS "Build type: ${CMAKE_BUILD_TYPE} (manually selected single-config)")
    endif()
  endif()

endfunction()


# Tries to guess which toolchain to load based on the environment.
#
# Args:
#     toolchain [out]: Name of the selected toolchain or undefined if it could not be selected
#
function(skprogs_guess_toolchain toolchain)

  if("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "GNU|GNU")
    set(_toolchain "gnu")
  elseif("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "Intel|Intel")
    set(_toolchain "intel")
  elseif("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "NAG|GNU")
    set(_toolchain "nag")
  else()
    set(_toolchain "generic")
  endif()

  set(${toolchain} "${_toolchain}" PARENT_SCOPE)

endfunction()


# Loads toolchain settings.
#
macro(skprogs_load_toolchain_settings)

  if(NOT DEFINED TOOLCHAIN_FILE AND NOT "$ENV{SKPROGS_TOOLCHAIN_FILE}" STREQUAL "")
    set(TOOLCHAIN_FILE "$ENV{SKPROGS_TOOLCHAIN_FILE}")
  endif()
  if(NOT DEFINED TOOLCHAIN AND NOT "$ENV{SKPROGS_TOOLCHAIN}" STREQUAL "")
    set(TOOLCHAIN "$ENV{SKPROGS_TOOLCHAIN}")
  endif()
  if(NOT DEFINED TOOLCHAIN_FILE OR TOOLCHAIN_FILE STREQUAL "")
    if(NOT DEFINED TOOLCHAIN OR TOOLCHAIN STREQUAL "")
      skprogs_guess_toolchain(TOOLCHAIN)
    endif()
    set(TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/sys/${TOOLCHAIN}.cmake)
  endif()
  message(STATUS "Reading build environment specific toolchain file: ${TOOLCHAIN_FILE}")
  include(${TOOLCHAIN_FILE})
endmacro()


# Sets up the global compiler flags
#
macro(skprogs_setup_global_compiler_flags)

  if(CMAKE_BUILD_TYPE)
    set(_buildtypes ${CMAKE_BUILD_TYPE})
  else()
    set(_buildtypes ${CMAKE_CONFIGURATION_TYPES})
  endif()
  foreach(_buildtype IN LISTS _buildtypes)
    foreach (lang IN ITEMS Fortran C)
      string(TOUPPER "${_buildtype}" _buildtype_upper)
      set(CMAKE_${lang}_FLAGS " ${${lang}_FLAGS}")
      set(CMAKE_${lang}_FLAGS_${_buildtype_upper} " ${${lang}_FLAGS_${_buildtype_upper}}")
      message(STATUS "Flags for ${lang}-compiler (build type: ${_buildtype}): "
        "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${_buildtype_upper}}")
    endforeach()
  endforeach()
  unset(_buildtypes)
  unset(_buildtype)
  unset(_buildtype_upper)
endmacro()


# Checks whether the current compiler fullfills minimal version requirements.
#
#
# Arguments:
#   lang [in]: Language for which the compiler should be checked (e.g. Fortran, C, CXX)
#   compiler_versions [in]: List with alternating compiler ids and minimal version numbers, e.g.
#       "Intel;19.0;GNU;9.0". If the compiler is amoung the listed ones and its version number is
#       less than the specified one, a fatal error message will be issued. Otherwise the function
#       returns silently.
#
function (skprogs_check_minimal_compiler_version lang compiler_versions)
  while(compiler_versions)
    list(POP_FRONT compiler_versions compiler version)
    if("${CMAKE_${lang}_COMPILER_ID}" STREQUAL "${compiler}"
        AND CMAKE_${lang}_COMPILER_VERSION VERSION_LESS "${version}")
      message(FATAL_ERROR "${compiler} ${lang} compiler is too old "
          "(found \"${CMAKE_${lang}_COMPILER_VERSION}\", required >= \"${version}\")")
    endif()
  endwhile()
endfunction()
