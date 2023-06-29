
add_executable(mpi_master mpi_master.cpp)
target_include_directories(mpi_master PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tests/include/>)
target_link_libraries(mpi_master asyncdempiShared asyncdeShared ${MPI_CXX_LIBRARIES})
