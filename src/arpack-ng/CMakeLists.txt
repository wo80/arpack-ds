target_include_directories(arpack
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

file(GLOB src_SRC SRC/*.c)
file(GLOB src_UTIL UTIL/*.c)

target_sources(arpack PRIVATE f2c.c arpack_types.c ${src_SRC} ${src_UTIL})

if(ENABLE_TESTS)
  add_subdirectory(EXAMPLES)
endif()
