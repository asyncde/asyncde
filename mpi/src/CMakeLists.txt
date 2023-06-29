
file(GLOB SRCS "*.cpp")

add_library(asyncdempiStatic STATIC ${SRCS})
set_target_properties(asyncdempiStatic PROPERTIES OUTPUT_NAME asyncdempi)
add_library(asyncdempiShared SHARED ${SRCS})
set_target_properties(asyncdempiShared PROPERTIES OUTPUT_NAME asyncdempi)

install(TARGETS asyncdempiStatic asyncdempiShared DESTINATION lib)
