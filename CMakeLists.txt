#####################################################################
# asyncde - Asynchronous Differential Evolution Package
# Top-Level CMake file
#####################################################################

cmake_minimum_required(VERSION 3.13)

project(asyncde)

set(CMAKE_VERBOSE_MAKEFILE on)

find_package(OpenMP)
# find_package(OpenMP QUIET)
find_package(MPI)

# Adjust paths for installation
if(UNIX)
  # Compiler options
#  set(CMAKE_CXX_COMPILER clang++)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED on)
  set(CMAKE_CXX_FLAGS "-g")
  add_compile_options(-O3) # optimization
  add_compile_options(-Wall -Wshadow -Woverloaded-virtual) # warnings
#  add_compile_options(-Werror) # no -Werror with openmpi
  add_compile_options(-Wextra) # extra warnings

  find_library(FIND_PROFILER_LIB NAMES profiler)
  if(FIND_PROFILER_LIB)
    message(STATUS "Found profiler library: " ${FIND_PROFILER_LIB})
    #    link_libraries("-lprofiler")
  elseif()
    message(STATUS "profiler library: not found")
  endif()
endif(UNIX)

message(STATUS "C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "Build directory: " ${PROJECT_BINARY_DIR})
message(STATUS "Install into: " ${CMAKE_INSTALL_PREFIX})

# paths for libraries and binaries
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# add sources
include_directories(include)
add_subdirectory(include)

# OMP multi-threading
if (OPENMP_FOUND)
  add_subdirectory(omp)
endif(OPENMP_FOUND)

# MPI part
if (MPI_CXX_FOUND)
  add_subdirectory(mpi)
endif(MPI_CXX_FOUND)

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(multithreads)

# uninstall according to http://www.cmake.org/Wiki/CMake_FAQ
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
