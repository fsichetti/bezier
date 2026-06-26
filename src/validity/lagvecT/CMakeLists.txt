file(GLOB SOURCES "*.cpp")
if(NOT INSTANTIATE_3_3_4)
    list(FILTER SOURCES EXCLUDE REGEX "lagrangeVectorT_3_3_4.*\.cpp$")
endif()

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})