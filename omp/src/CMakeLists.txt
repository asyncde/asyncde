
file(GLOB SRCS "*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_library(asyncdeompStatic STATIC ${SRCS})
set_target_properties(asyncdeompStatic PROPERTIES OUTPUT_NAME asyncdeomp)
add_library(asyncdeompShared SHARED ${SRCS})
set_target_properties(asyncdeompShared PROPERTIES OUTPUT_NAME asyncdeomp)

install(TARGETS asyncdeompStatic asyncdeompShared DESTINATION lib)
