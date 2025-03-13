
file(GLOB SRCS "*.cpp")

add_library(asyncdethreadsStatic STATIC ${SRCS})
set_target_properties(asyncdethreadsStatic PROPERTIES OUTPUT_NAME asyncdethreads)
add_library(asyncdethreadsShared SHARED ${SRCS})
set_target_properties(asyncdethreadsShared PROPERTIES OUTPUT_NAME asyncdethreads)

install(TARGETS asyncdethreadsStatic asyncdethreadsShared DESTINATION lib)
