cmake_minimum_required(VERSION 3.16.0)

project(arpack VERSION 3.9.1
               DESCRIPTION "ARPACK-NG: solve large scale eigenvalue problems"
               LANGUAGES C)

set(arpack_ng_MAJOR_VERSION 3)
set(arpack_ng_MINOR_VERSION 9)
set(arpack_ng_PATCH_VERSION 1)
set(arpack_ng_VERSION ${arpack_ng_MAJOR_VERSION}.${arpack_ng_MINOR_VERSION}.${arpack_ng_PATCH_VERSION})

option(ENABLE_TESTS "Compile ARPACK tests" ON)
option(ENABLE_TIMING "Enable ARPACK timing feature" ON)
option(ENABLE_TRACING "Enable ARPACK tracing feature" ON)

if(MSVC)
  add_compile_definitions(NOMINMAX _COMPLEX_DEFINED)
endif()

if(NOT ENABLE_TIMING)
  add_compile_definitions(NO_TIMER)
endif()

if(NOT ENABLE_TRACING)
  add_compile_definitions(NO_TRACE)
endif()

# Find BLAS

if (NOT TARGET BLAS::BLAS) # Search only if not already found by upper CMakeLists.txt
  find_package(BLAS REQUIRED)

  # BLAS::BLAS target was already created at this point by FindBLAS.cmake if cmake version >= 3.18
  if (NOT TARGET BLAS::BLAS) # Create target "at hand" to ensure compatibility if cmake version < 3.18
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
  endif()
endif()

# Find LAPACK

if (NOT TARGET LAPACK::LAPACK) # Search only if not already found by upper CMakeLists.txt
  find_package(LAPACK REQUIRED)

  # LAPACK::LAPACK target was already created at this point by FindLAPACK.cmake if cmake version >= 3.18
  if (NOT TARGET LAPACK::LAPACK) # Create target "at hand" to ensure compatibility if cmake version < 3.18
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
  endif()
endif()

if(ENABLE_TESTS)
  enable_testing()
endif()

add_subdirectory(src)
