add_library(arpack)
add_library(ARPACK::ARPACK ALIAS arpack)

set_target_properties(arpack PROPERTIES DEBUG_POSTFIX _d)

target_link_libraries(arpack
  PUBLIC
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,LAPACK::LAPACK>>
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,BLAS::BLAS>>
  $<BUILD_INTERFACE:LAPACK::LAPACK>
  $<BUILD_INTERFACE:BLAS::BLAS>
)

add_subdirectory(arpack-ng)
