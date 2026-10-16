# -*- cmake -*-

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(FATAL_ERROR "PyROL can only be built with shared libraries. Building of shared libraries is currently set to OFF. To enable shared libraries please set the cache variable \"BUILD_SHARED_LIBS\" to ON")
ENDIF()

# Set the package version number
SET(PyROL_VERSION ${Trilinos_VERSION})

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_PIP_INSTALL
  PYROL_SCIKIT
  "Logic specific to the scikit-build-core backend."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_ENABLE_BINDER
  PYROL_GENERATE_SRC
  "Enable regeneration of source files with Binder."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_BINDER_SUPPRESS_ERRORS
  PYROL_SUPPRESS_ERRORS
  "Enable the suppress errors option of Binder."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_BINDER_USE_ONE_FILE
  PYROL_USE_ONE_FILE
  "Enable the use of one file by Binder."
  OFF )

SET(PYROL_BINDER_NUM_FILES "100" CACHE STRING "Number of files to use when PYROL_USE_ONE_FILE=OFF.")

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_BINDER_CMAKE_ERROR
  PYROL_CMAKE_ERROR
  "Stop the configuration if Binder fails."
  ON )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_BINDER_VERBOSE
  PYROL_B_VERBOSE
  "Increase the verbosity of Binder"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_ENABLE_BINDER_UPDATE
  PYROL_UPDATE_GENERATED_SRC
  "Enable the update of the generated source files with Binder."
  OFF )

##################################################
# Helpers                                        #
##################################################

# Append the include directories of all targets that LIBRARY_NAME links to
# (recursively, with $<BUILD_INTERFACE:...> unwrapped) to all_include_dirs.
function(get_all_include_dirs  LIBRARY_NAME all_include_dirs all_visited_libs)
  if (TARGET ${LIBRARY_NAME})
    get_property(depend_libs TARGET ${LIBRARY_NAME} PROPERTY INTERFACE_LINK_LIBRARIES)
    foreach(depend_lib IN LISTS depend_libs)
      if (TARGET ${depend_lib} AND (NOT ${depend_lib} IN_LIST all_visited_libs))
        list(APPEND all_visited_libs "${depend_lib}") # Update list in the current scope only
        get_property(current_includes TARGET ${depend_lib} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
        foreach(include IN LISTS current_includes)
          STRING(REPLACE "$<BUILD_INTERFACE:" "" new_tmp_include ${include})
          STRING(REPLACE ">" "" new_include ${new_tmp_include})
          list(APPEND all_include_dirs "${new_include}") # Update list in the current scope only
        endforeach()
        get_all_include_dirs(${depend_lib} "${all_include_dirs}" "${all_visited_libs}")
      endif()
    endforeach()
    set(all_include_dirs ${all_include_dirs} PARENT_SCOPE)
    set(all_visited_libs ${all_visited_libs} PARENT_SCOPE)
  endif()
endfunction()

# Set result to the list of subdirectories of curdir.
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

# Write the items of the list variable list_var to filename, one per line.
function(write_lines filename list_var)
  set(contents "")
  foreach(line IN LISTS ${list_var})
    string(APPEND contents "${line}\n")
  endforeach()
  file(WRITE ${filename} "${contents}")
endfunction()

##################################################
# Python and pybind11                            #
##################################################

find_program(Python3_EXECUTABLE NAMES python3 python)
MESSAGE("-- Python3_EXECUTABLE: ${Python3_EXECUTABLE}")

# Interpreter lets FindPython match the interpreter's ABI; without it, CMake 4.4+
# cannot find Development.Module for free-threaded Python (e.g. 3.14t).
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
SET(PYTHON_VERSION ${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR})

EXECUTE_PROCESS(COMMAND
  ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE pybind11_DIR
  ERROR_VARIABLE  pybind11_CMAKE_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
MESSAGE(STATUS "pybind11 CMake path: ${pybind11_DIR}")
find_package(pybind11 REQUIRED)

EXECUTE_PROCESS(COMMAND
  ${Python3_EXECUTABLE} -c "import mpi4py; print(mpi4py.get_include())"
  OUTPUT_VARIABLE Mpi4Py_INCLUDE_DIR
  ERROR_VARIABLE  Mpi4Py_INCLUDE_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# The pure-Python part of the package: top-level modules and subpackages.
FILE(GLOB PyROLPyFiles ${CMAKE_CURRENT_SOURCE_DIR}/python/*.py)
SUBDIRLIST(PyROL_SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/python)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pyrol)

##################################################
# Include directories                            #
##################################################

GET_PROPERTY( PTR_IMPL             GLOBAL PROPERTY PTR_IMPL  )
GET_PROPERTY( PTR_DIR              GLOBAL PROPERTY PTR_DIR  )
GET_PROPERTY( STACKTRACE_DIR       GLOBAL PROPERTY STACKTRACE_DIR  )
GET_PROPERTY( PARAMETERLIST_DIR    GLOBAL PROPERTY PARAMETERLIST_DIR  )
GET_PROPERTY( LA_DIR               GLOBAL PROPERTY LA_DIR  )
GET_PROPERTY( LAPACK_DIR           GLOBAL PROPERTY LAPACK_DIR  )
GET_PROPERTY( BLAS_DIR             GLOBAL PROPERTY BLAS_DIR  )

# ROL's include directories, in search order. They are used to compile the
# module and, when Binder runs, to collect the headers it may need.
SET(ROL_SRC_DIR "${${PACKAGE_NAME}_SOURCE_DIR}/src")
SET(ROL_SRC_SUBDIRS
  algorithm
  algorithm/TypeB
  algorithm/TypeE
  algorithm/TypeG
  algorithm/TypeP
  algorithm/TypeU
  algorithm/TypeB/pqn
  algorithm/TypeG/augmentedlagrangian/
  algorithm/TypeG/fletcher
  algorithm/TypeG/interiorpoint
  algorithm/TypeG/moreauyosida
  algorithm/TypeG/stabilizedlcl
  algorithm/TypeU/bundle
  algorithm/TypeU/linesearch
  algorithm/TypeU/linesearch/descent
  algorithm/TypeU/trustregion
  elementwise
  function
  function/dynamic
  function/boundconstraint
  function/constraint
  function/nlls
  function/projection
  function/objective
  function/operator
  function/simopt
  function/sketching
  function/std
  sol/function
  sol/function/distribution
  sol/function/expectationquad
  sol/function/randvarfunctional
  sol/function/randvarfunctional/risk
  sol/function/randvarfunctional/risk/fdivergence
  sol/function/randvarfunctional/risk/spectral
  sol/function/randvarfunctional/error
  sol/function/randvarfunctional/regret
  sol/function/randvarfunctional/deviation
  sol/function/randvarfunctional/probability
  sol/function/progressivehedging
  sol/algorithm
  sol/algorithm/PrimalDual
  sol/vector
  sol/sampler
  sol/sampler/SROM
  sol/status
  status
  step
  step/augmentedlagrangian
  step/bundle
  step/fletcher
  step/interiorpoint
  step/krylov
  step/linesearch
  step/moreauyosidapenalty
  step/nonlinearcg
  step/secant
  step/trustregion
  utils
  utils/function_bindings
  vector
  zoo
  )
SET(ROL_OED_SUBDIRS
  oed
  oed/algorithm
  oed/constraint
  oed/factors
  oed/objective
  oed/objective/het
  oed/objective/hom
  oed/utilities
  oed/vector
  )
list(TRANSFORM ROL_SRC_SUBDIRS PREPEND "${ROL_SRC_DIR}/")
list(TRANSFORM ROL_OED_SUBDIRS PREPEND "${ROL_SRC_DIR}/")

SET(ROL_all_include_dirs
  "${CMAKE_CURRENT_BINARY_DIR}/../src"
  "${ROL_SRC_DIR}"
  ${ROL_SRC_SUBDIRS}
  "${PTR_DIR}"
  "${STACKTRACE_DIR}"
  "${PARAMETERLIST_DIR}"
  "${LA_DIR}"
  "${LAPACK_DIR}"
  "${BLAS_DIR}"
  ${ROL_OED_SUBDIRS}
  )

##################################################
# Binder-generated sources                       #
##################################################

IF (PYROL_GENERATE_SRC)
  MESSAGE("-- PyROL_BINDER_EXECUTABLE:")
  IF(NOT DEFINED PyROL_BINDER_EXECUTABLE)
    find_program(PyROL_BINDER_EXECUTABLE
        NAMES binder
        )
    MESSAGE("  -- CMake has set: PyROL_BINDER_EXECUTABLE = ${PyROL_BINDER_EXECUTABLE}")
  ELSE()
    MESSAGE("  -- User has set: PyROL_BINDER_EXECUTABLE = ${PyROL_BINDER_EXECUTABLE}")
  ENDIF()
  file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
  file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/binder)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/binder)

  MESSAGE("PTR_IMPL = ${PTR_IMPL}")
  IF( PTR_IMPL STREQUAL "std::shared_ptr")
    SET( BINDER_CFG "${CMAKE_CURRENT_SOURCE_DIR}/scripts/PyROL_shared_ptr.cfg" )
  ELSE()
    SET( BINDER_CFG "${CMAKE_CURRENT_SOURCE_DIR}/scripts/PyROL_RCP.cfg" )
  ENDIF()

  set(binder_include_name "${CMAKE_CURRENT_BINARY_DIR}/trilinos_includes.hpp")
  set(all_header_with_dir_list "${CMAKE_CURRENT_BINARY_DIR}/list_with_dir.txt")
  set(all_header_without_dir_list "${CMAKE_CURRENT_BINARY_DIR}/list_without_dir.txt")
  set(all_include_list "${CMAKE_CURRENT_BINARY_DIR}/list_include.txt")

  # Include directories of ROL's dependencies, followed by ROL's own.
  set(all_include_dirs "")
  set(all_visited_libs "")
  foreach(depPkg IN LISTS ROL_LIB_ENABLED_DEPENDENCIES)
    get_all_include_dirs(${depPkg}::all_libs "${all_include_dirs}" "${all_visited_libs}")
  endforeach()

  list(REMOVE_DUPLICATES all_include_dirs)
  list(REMOVE_ITEM all_include_dirs "")

  FOREACH(include_dir IN LISTS ROL_all_include_dirs)
    list(APPEND all_include_dirs "${include_dir}")
  ENDFOREACH(include_dir)

  # The headers in those directories, with and without the directory prefix.
  set(PyROL_all_include_files_with_dir "")
  set(PyROL_all_include_files_without_dir "")
  foreach(include_dir IN LISTS all_include_dirs)
    file(GLOB include_files
        "${include_dir}/*.hpp"
        "${include_dir}/*.h"
        "${include_dir}/*/*.hpp"
        "${include_dir}/*/*.h"
        "${include_dir}/*/*/*.hpp"
        "${include_dir}/*/*/*.h"
        "${include_dir}/*/*/*/*.hpp"
        "${include_dir}/*/*/*/*.h"
        "${include_dir}/*/*/*/*.inc"
        "${include_dir}/*/*/*/*.inc_predicate"
    )
    foreach(include_file IN LISTS include_files)
      list(APPEND PyROL_all_include_files_with_dir "${include_file}")
      string(REPLACE "${include_dir}/" "" include_file_without_dir "${include_file}")
      list(APPEND PyROL_all_include_files_without_dir "${include_file_without_dir}")
    endforeach()
  endforeach()

  write_lines(${all_include_list} all_include_dirs)
  write_lines(${all_header_with_dir_list} PyROL_all_include_files_with_dir)
  write_lines(${all_header_without_dir_list} PyROL_all_include_files_without_dir)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)

  # Copy the headers to include_tmp with #include "..." rewritten to <...>,
  # and write Binder's input file.
  EXECUTE_PROCESS(COMMAND
  ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/gather_includes.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${all_header_with_dir_list} ${all_header_without_dir_list} ${binder_include_name}
  )

  set(BINDER_OPTIONS "")
  list(APPEND BINDER_OPTIONS ${binder_include_name})
  list(APPEND BINDER_OPTIONS --root-module pyrol)
  list(APPEND BINDER_OPTIONS --prefix ${CMAKE_CURRENT_BINARY_DIR}/binder)
  IF(PYROL_USE_ONE_FILE)
    list(APPEND BINDER_OPTIONS -single-file)
  ELSE()
    list(APPEND BINDER_OPTIONS -max-file-size=1000000)
    list(APPEND BINDER_OPTIONS -flat)
  ENDIF()
  list(APPEND BINDER_OPTIONS --bind Teuchos)
  list(APPEND BINDER_OPTIONS --bind ROL)
  list(APPEND BINDER_OPTIONS --bind pyrol)
  IF(PYROL_B_VERBOSE)
    list(APPEND BINDER_OPTIONS -v)
  ENDIF()
  list(APPEND BINDER_OPTIONS --config ${BINDER_CFG})
  IF(PYROL_SUPPRESS_ERRORS)
    list(APPEND BINDER_OPTIONS --suppress-errors)
  ENDIF()
  list(APPEND BINDER_OPTIONS --)
  IF(TPL_ENABLE_CUDA)
    list(APPEND BINDER_OPTIONS -x cuda --cuda-host-only)
  ENDIF()
  list(APPEND BINDER_OPTIONS ${PyROL_BINDER_FLAGS})
  list(APPEND BINDER_OPTIONS -std=c++17)
  list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
  list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/src)
  list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/src)

  list(APPEND BINDER_OPTIONS -I${PyROL_BINDER_clang_include_dirs})
  list(APPEND BINDER_OPTIONS -iwithsysroot${PyROL_BINDER_LibClang_include_dir})
  IF(PyROL_BINDER_GCC_TOOLCHAIN)
    list(APPEND BINDER_OPTIONS --gcc-toolchain=${PyROL_BINDER_GCC_TOOLCHAIN})
  ENDIF()
  list(APPEND BINDER_OPTIONS -DNDEBUG)

  message("BINDER_OPTIONS='${BINDER_OPTIONS}'")

  # Binder writes as many pyrol_<i>.cpp files as it needs; pre-create the
  # fixed set that src/CMakeLists.txt compiles, so unused ones stay empty.
  IF(NOT PYROL_USE_ONE_FILE)
    MATH(EXPR NUMBER_FILE "${PYROL_BINDER_NUM_FILES}")

    foreach(index RANGE 0 ${NUMBER_FILE})
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/binder/pyrol_${index}.cpp "")
    endforeach()
  ENDIF()

  EXECUTE_PROCESS(COMMAND
    ${PyROL_BINDER_EXECUTABLE} ${BINDER_OPTIONS}
    RESULT_VARIABLE STATUS
    OUTPUT_VARIABLE OUTPUT_BINDER
  )

  if(STATUS AND NOT STATUS EQUAL 0)
    message("${OUTPUT_BINDER}")
    if(PYROL_CMAKE_ERROR)
      message(FATAL_ERROR "BINDER FAILED: ${STATUS}")
    else()
      message("BINDER FAILED: ${STATUS}")
    endif()
  else()
    message(STATUS "BINDER SUCCESS:")
    message("${OUTPUT_BINDER}")
  endif()

  IF(PYROL_UPDATE_GENERATED_SRC)
    FILE(REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/binder)
    FILE(COPY ${CMAKE_CURRENT_BINARY_DIR}/binder DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
  ENDIF()
ELSE()
  # Copy binder files to ${CMAKE_CURRENT_BINARY_DIR}
  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/binder DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()

##################################################
# Module, installation and tests                 #
##################################################

# Determine the install directory
IF(PYROL_SCIKIT)
  SET(PyROL_INSTALL_DIR
    "."
    CACHE STRING "The path where PyROL will be installed."
  )
ELSE()
  SET(PyROL_INSTALL_DIR
    lib/python${PYTHON_VERSION}/site-packages/pyrol
    CACHE STRING "The path where PyROL will be installed"
    )
ENDIF()

cmake_path(IS_RELATIVE PyROL_INSTALL_DIR is_pyrol_install_relative)

IF(NOT is_pyrol_install_relative)
  MESSAGE(FATAL_ERROR "PyROL install dir is not relative.")
ENDIF()

cmake_path(GET PyROL_INSTALL_DIR RELATIVE_PART PyROL_INSTALL_DIR_RELATIVE)

MESSAGE(STATUS "PyROL installation path: ${PyROL_INSTALL_DIR_RELATIVE}")

INSTALL(FILES
  ${PyROLPyFiles}
  DESTINATION ${PyROL_INSTALL_DIR})

FOREACH(PyROL_subdir ${PyROL_SUBDIRS})
  INSTALL(DIRECTORY
    ${PyROL_subdir}
    DESTINATION ${PyROL_INSTALL_DIR})
ENDFOREACH()

ADD_SUBDIRECTORY( src )

# Lay out the package in the build tree too, so that tests can import it.
FILE(COPY ${PyROLPyFiles} ${PyROL_SUBDIRS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pyrol)

SET(PyROL_PYTHONPATH "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

TRIBITS_ADD_TEST_DIRECTORIES(test)
