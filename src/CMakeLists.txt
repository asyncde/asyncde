
file(GLOB SRCS "*.cpp")

add_library(asyncdeStatic STATIC ${SRCS})
set_target_properties(asyncdeStatic PROPERTIES OUTPUT_NAME asyncde)

add_library(asyncdeShared SHARED ${SRCS})
set_target_properties(asyncdeShared PROPERTIES OUTPUT_NAME asyncde)

install(TARGETS asyncdeStatic asyncdeShared DESTINATION lib)
