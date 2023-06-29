
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# add sources
include_directories(include)

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(examples)
