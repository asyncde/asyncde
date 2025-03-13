
add_executable(threads_master threads_master.cpp)
target_include_directories(threads_master PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tests/include/>)
target_link_libraries(threads_master asyncdethreadsShared asyncdeShared)
