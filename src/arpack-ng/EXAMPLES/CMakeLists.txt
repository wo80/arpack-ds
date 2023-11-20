function(examples list_name)
    foreach(name ${${list_name}})
        add_executable(${name} "${arpackexample_DIR}/${name}.c")
        target_link_libraries(${name} arpack band)
        add_test(NAME ${name}_test COMMAND ${name})
    endforeach()
endfunction(examples)

set(band_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BAND/)

add_library(band
    "${band_DIR}/cnband.c"
    "${band_DIR}/dnband.c" "${band_DIR}/dsband.c"
    "${band_DIR}/snband.c" "${band_DIR}/ssband.c"
    "${band_DIR}/znband.c")

target_link_libraries(band arpack)

# examples/band
set(arpackexample_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BAND/)
set(examples_SRCS
    cnbdr1 cnbdr2 cnbdr3 cnbdr4
    dnbdr1 dnbdr2 dnbdr3 dnbdr4 dnbdr5 dnbdr6
    dsbdr1 dsbdr2 dsbdr3 dsbdr4 dsbdr5 dsbdr6
    snbdr1 snbdr2 snbdr3 snbdr4 snbdr5 snbdr6
    ssbdr1 ssbdr2 ssbdr3 ssbdr4 ssbdr5 ssbdr6
    znbdr1 znbdr2 znbdr3 znbdr4)
examples(examples_SRCS)

# examples/complex
set(arpackexample_DIR ${CMAKE_CURRENT_SOURCE_DIR}/COMPLEX/)
set(examples_SRCS
    cndrv1 cndrv2 cndrv3 cndrv4
    zndrv1 zndrv2 zndrv3 zndrv4)
examples(examples_SRCS)

# examples/nonsym
set(arpackexample_DIR ${CMAKE_CURRENT_SOURCE_DIR}/NONSYM/)
set(examples_SRCS
    dndrv1 dndrv2 dndrv3 dndrv4 dndrv5 dndrv6
    sndrv1 sndrv2 sndrv3 sndrv4 sndrv5 sndrv6)
examples(examples_SRCS)

# examples/simple
set(arpackexample_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SIMPLE/)
set(examples_SRCS
    cnsimp dnsimp dssimp
    snsimp sssimp znsimp)
examples(examples_SRCS)

# examples/svd
set(arpackexample_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SVD/)
set(examples_SRCS dsvd ssvd)
examples(examples_SRCS)

# examples/sym
set(arpackexample_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SYM/)
set(examples_SRCS
    dsdrv1 dsdrv2 dsdrv3 dsdrv4 dsdrv5 dsdrv6
    ssdrv1 ssdrv2 ssdrv3 ssdrv4 ssdrv5 ssdrv6)
examples(examples_SRCS)
