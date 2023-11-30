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

include(GNUInstallDirs)

set(ARPACK_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/arpack")
set(ARPACK_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/arpackng")

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

# Install

set(PACKAGE_NAME ${PROJECT_NAME})
set(PACKAGE_VERSION ${arpack_ng_VERSION})
set(PACKAGE_URL "https://github.com/wo80/arpack-ds/")

set(installTargets arpack)

install(TARGETS ${installTargets}
  EXPORT arpackngTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES
  src/arpack-ng/arpack.h
  src/arpack-ng/arpack_types.h
  DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})

install(EXPORT arpackngTargets
  DESTINATION "${ARPACK_INSTALL_CMAKEDIR}"
)

# Provide find_package for arpack-ng to users.
configure_file(cmake/arpackng-config.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config.cmake" @ONLY)
configure_file(cmake/arpackng-config-version.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake" @ONLY)
install(
  FILES
  "${PROJECT_BINARY_DIR}/arpackng-config.cmake"
  "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arpackng)
