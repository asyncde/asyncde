
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(omp_master omp_master.cpp)
target_include_directories(omp_master PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tests/include/>)
target_link_libraries(omp_master asyncdeompShared asyncdeShared)
