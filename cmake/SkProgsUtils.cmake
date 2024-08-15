include(FetchContent)

# Replaces the extension of a given file
#
# Args:
#     oldext [in]: Old extension
#     newext [in]: New extension
#     fname [in]: File name in which extension should be replaced.
#     newfname [out]: File name after extension replacement.
#
function(skprogs_replace_extension oldext newext fname newfname)

  string(REGEX REPLACE "\\.${oldext}$" ".${newext}" _newfname ${fname})
  set(${newfname} ${_newfname} PARENT_SCOPE)

endfunction()


# Registers files for preprocessing
#
# Args:
#     preproc [in]: Preprocessor to use
#     preprocopts [in]:  Preprocessor command line arguments (but not in/out file)
#     oldext [in]: Extension of the unpreprocessed files.
#     newext [in]: Extension of the preprocessed files.
#     oldfiles [in]: List of unpreprocessed file names.
#     newfiles [out]: List of preprocessed file names.
#
function(skprogs_preprocess preproc preprocopts oldext newext oldfiles newfiles)

  set(_newfiles)
  foreach(oldfile IN LISTS oldfiles)
    skprogs_replace_extension(${oldext} ${newext} ${oldfile} newfile)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      COMMAND ${preproc} ${preprocopts} ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile} ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile})
    list(APPEND _newfiles ${CMAKE_CURRENT_BINARY_DIR}/${newfile})
  endforeach()
  set(${newfiles} ${_newfiles} PARENT_SCOPE)

endfunction()


# Build -D command line arguments for Fypp preprocessor based on current configuration
#
# Args:
#     fyppflags [inout]: Current Fypp flags on enter, with -D options extended flags on exit.
#
function (skprogs_add_fypp_defines fyppflags)

  set(_fyppflags "${${fyppflags}}")

  if(WITH_MPI)
    list(APPEND _fyppflags -DWITH_MPI)
  endif()

  set(${fyppflags} ${_fyppflags} PARENT_SCOPE)

endfunction()

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

  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(_toolchain "gnu")
  elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel|IntelLLVM")
    set(_toolchain "intel")
  elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
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
    foreach (lang IN ITEMS Fortran)
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


# Handles a hybrid dependency.
#
# Depending on the list items in the config_methods variable, it will try to:
#
# - checkout the source as a submodule within the origin sub-folder ("Submodule")
# - find the package as external dependency ("Find")
# - fetch the source from a git repository ("Fetch") into the build folder
#
# The methods are tried in the order of their appearance until success the first eligible one.
#
# The methods "Submodule" and "Fetch" would call the passed sub-directory with add_subdirectory()
# passing two variables with the source and binary directory.
#
# Args:
#     package [in]: Name of the dependency to look for.
#     config_methods [in]: Config methods to try
#     target [in]: Name of the target, which must be exported after the configuration.
#     findpkgopts [in]: Options to pass to find_package()
#     subdir [in]: Subdirectory with CMakeFiles.txt for integrating package source.
#     subdiropts [in]: Options to pass to the add_subdir() command.
#     git_repository [in]: Git repository to fetch the package from.
#     git_tag [in]: Git tag to use when fetching the source.
#
# Variables:
#     <UPPER_PACKAGE_NAME>_SOURCE_DIR, <UPPER_PACKAGE_NAME>_BINARY_DIR:
#         Source and binary directories for the build (to pass to add_subdirectory())
#
macro(skprogs_config_hybrid_dependency package target config_methods findpkgopts subdir subdiropts
    git_repository git_tag)

  set(_allowed_methods "submodule;find;fetch;pkgconf")
  string(TOLOWER "${package}" _package_lower)
  string(TOUPPER "${package}" _package_upper)

  foreach(_config_method IN ITEMS ${config_methods})

    string(TOLOWER "${_config_method}" _config_lower)
    if(NOT ${_config_lower} IN_LIST _allowed_methods)
      message(FATAL_ERROR "${package}: Unknown configuration method '${_config_method}'")
    endif()

    if("${_config_lower}" STREQUAL "find")

      message(STATUS "${package}: Trying to find installed package")
      find_package(${package} ${findpkgopts})
      if(${package}_FOUND)
        message(STATUS "${package}: Installed package found")
        break()
      else()
        message(STATUS "${package}: Installed package could not be found")
      endif()

    elseif("${_config_lower}" STREQUAL "pkgconf")
      message(STATUS "${package}: Trying to find installed package (pkg-config)")

      find_package(PkgConfig QUIET)
      if(PkgConfig_FOUND)
        pkg_check_modules("${_package_upper}" QUIET "${package}")
        if("${${_package_upper}_FOUND}")
          message(STATUS "${package}: Installed package found (pkg-config)")
          add_library("${target}" INTERFACE IMPORTED)
          target_link_libraries(
            "${target}"
            INTERFACE
            "${${_package_upper}_LINK_LIBRARIES}"
          )
          target_include_directories(
            "${target}"
            INTERFACE
            "${${_package_upper}_INCLUDE_DIRS}"
          )
          break()
        else()
          message(STATUS "${package}: Installed package could not be found (pkg-config)")
        endif()
      endif()

    elseif("${_config_lower}" STREQUAL "submodule")

      if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/origin/CMakeLists.txt
          AND GIT_WORKING_COPY)
        message(STATUS "${package}: Downloading via git submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ${subdir}/origin
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      endif()

      if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/origin/CMakeLists.txt)
        message(STATUS "${package}: Using source in ${subdir}/origin")
        set(${_package_upper}_SOURCE_DIR "origin")
        set(${_package_upper}_BINARY_DIR)
        add_subdirectory(${subdir} ${subdiropts})
        break()
      endif()

    elseif("${_config_lower}" STREQUAL "fetch")

      message(STATUS "${package}: Fetching from repository ${git_repository}@${git_tag}")
      FetchContent_Declare(${_package_lower} GIT_REPOSITORY ${git_repository} GIT_TAG ${git_tag})
      FetchContent_GetProperties(${_package_lower})
      if(NOT ${_package_lower}_POPULATED)
        FetchContent_Populate(${_package_lower})
      endif()
      set(${_package_upper}_SOURCE_DIR "${${_package_lower}_SOURCE_DIR}")
      set(${_package_upper}_BINARY_DIR "${${_package_lower}_BINARY_DIR}")
      add_subdirectory(${subdir} ${subdiropts})
      break()

    endif()

  endforeach()

  if(NOT TARGET ${target})
    message(FATAL_ERROR "Could not configure ${package} to export target '${target}'")
  endif()

  unset(_allowed_methods)
  unset(_package_lower)
  unset(_package_upper)
  unset(_config_method)
  unset(_config_lower)

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
