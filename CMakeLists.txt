# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(BEZIER_TOPLEVEL_PROJECT OFF)
else()
    set(BEZIER_TOPLEVEL_PROJECT ON)
endif()

cmake_minimum_required(VERSION 3.18)
project(bezier)


file(GLOB BEZIER_SOURCES
    "${BEZIER_SOURCE_DIR}/*.cpp"
    "${BEZIER_SOURCE_DIR}/numbers/*.cpp"
    "${BEZIER_SOURCE_DIR}/utils/*.cpp"
    "${BEZIER_SOURCE_DIR}/validity/*.cpp"
    "${BEZIER_SOURCE_DIR}/validity/*/*.cpp"
)
message("${BEZIER_SOURCES}")
add_library(bezier ${GLOB_BEZIER_SOURCES})
add_library(bezier::bezier ALIAS bezier)
add_library(bezier_interface INTERFACE)
target_link_libraries(bezier PRIVATE bezier_interface)

# Set the C++ standard
# target_compile_features(bezier_interface INTERFACE cxx_std_20)

# Options
option(IPRED_ARITHMETIC "Use the efficient Indirect Predicates library" ON)
option(UNIT_TESTS "Run unit tests" ON)
option(HDF5_INTERFACE "Process HDF5 datasets" ON)
option(EIGEN_INTERFACE "Read data from Eigen classes" ON)

# Experimental flags and stuff for benchmarks (leave these OFF)
option(EXPERIM_NAIVE_SUBDIV "Use a naive subdivision scheme" OFF)
option(EXPERIM_FLOAT_ARITHMETIC "Use non-robust floating point arithmetic" OFF)

# Add definitions
if(EXPERIM_NAIVE_SUBDIV)
    add_compile_definitions(EXPERIM_NAIVE_SUBDIV)
endif()
if(EXPERIM_FLOAT_ARITHMETIC)
    add_compile_definitions(EXPERIM_FLOAT_ARITHMETIC)
endif()

# Cmake recipes & CPM Cache 
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
include(bezier_cpm_cache)

set_target_properties(bezier PROPERTIES COMPILE_WARNING_AS_ERROR ON)
target_compile_options(bezier INTERFACE
    -Wall
    -Wno-unknown-pragmas
    -fno-strict-aliasing
)

# Use fast interval arithmetic
if(DEFINED IPRED_ARITHMETIC)
    target_compile_options(bezier_interface INTERFACE
        -frounding-math
        -msse2
        -mavx2
        -mfma
    )
    set(USE_SIMD_INSTRUCTIONS ON)
endif()

# Ensure the build type is set
if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "CMAKE_BUILD_TYPE not specified. Defaulting to Release.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()


set(BEZIER_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(BEZIER_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

add_subdirectory(${BEZIER_SOURCE_DIR})

target_include_directories(bezier PUBLIC ${BEZIER_INCLUDE_DIR})

################################################################################
# Required libraries
################################################################################

##### OpenMP ###################################################################
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(bezier_interface INTERFACE OpenMP::OpenMP_CXX)
else()
endif()

##### GMP ######################################################################
find_package(GMP REQUIRED)
target_link_libraries(bezier PRIVATE gmp::gmp)

##### Eigen ####################################################################
if (EIGEN_INTERFACE)
    target_compile_definitions(bezier_interface INTERFACE EIGEN_INTERFACE)
    include(eigen)
    target_link_libraries(bezier_interface INTERFACE Eigen3::Eigen)
endif()

##### INDIRECT_PREDICATES ######################################################
if (IPRED_ARITHMETIC)
    target_compile_definitions(bezier_interface INTERFACE IPRED_ARITHMETIC)

    include(indirect_predicates)
endif()


##### HDF5 #####################################################################
if (HDF5_INTERFACE AND BEZIER_TOPLEVEL_PROJECT)
    target_compile_definitions(bezier_interface INTERFACE HDF5_INTERFACE)
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    target_include_directories(bezier_interface INTERFACE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(bezier_interface INTERFACE ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
endif()


##### EXECUTABLES ##############################################################
if (BEZIER_TOPLEVEL_PROJECT)
    add_executable(bin ${PROJECT_SOURCE_DIR}/src/main.cpp)
    target_link_libraries(bin PRIVATE bezier::bezier)
    target_link_libraries(bin PRIVATE bezier_interface)
endif()

################################################################################
# Tests
################################################################################

# Compile extras only if this is a top-level project
if(BEZIER_TOPLEVEL_PROJECT AND UNIT_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()