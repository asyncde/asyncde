
file(GLOB SRCS "*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_library(asyncdeStatic STATIC ${SRCS})
set_target_properties(asyncdeStatic PROPERTIES OUTPUT_NAME asyncde)

add_library(asyncdeShared SHARED ${SRCS})
set_target_properties(asyncdeShared PROPERTIES OUTPUT_NAME asyncde)

install(TARGETS asyncdeStatic asyncdeShared DESTINATION lib)
