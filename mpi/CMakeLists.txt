
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})

message(STATUS "MPI C++ compiler: " ${CMAKE_CXX_COMPILER})

# add sources
include_directories(include)

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(examples)
