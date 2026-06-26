include(catch2)

add_executable(unit_tests tests.cpp)
target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain)
target_include_directories(unit_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_property(TARGET unit_tests PROPERTY COMPILE_WARNING_AS_ERROR ON)
target_link_libraries(unit_tests PRIVATE bezier::bezier)
target_link_libraries(unit_tests PRIVATE bezier_interface)

if(EIGEN_INTERFACE)
	add_executable(eigen_tests eigen_tests.cpp)
	target_link_libraries(eigen_tests PRIVATE Catch2::Catch2WithMain)
	target_include_directories(eigen_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)
	set_property(TARGET eigen_tests PROPERTY COMPILE_WARNING_AS_ERROR ON)
	target_link_libraries(eigen_tests PRIVATE bezier::bezier)
	target_link_libraries(eigen_tests PRIVATE bezier_interface)
endif()

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(unit_tests)

if(EIGEN_INTERFACE)
catch_discover_tests(eigen_tests)
endif()
