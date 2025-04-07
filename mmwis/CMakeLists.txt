cmake_minimum_required(VERSION 3.10)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(MMWIS C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()


add_definitions(-w)
# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
  add_definitions(-fno-stack-limit)
endif()
#CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
#if(COMPILER_SUPPORTS_WALL)
  #add_definitions(-Wall)
#endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  add_definitions(-march=native)
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
  add_definitions(-fpermissive)
endif()

CHECK_CXX_COMPILER_FLAG(-Wno-unused-value COMPILER_SUPPORTS_NOUNUSED)
if(COMPILER_SUPPORTS_NOUNUSED)
add_definitions(-Wno-unused-value)
endif()
CHECK_CXX_COMPILER_FLAG(-Wno-unused-value COMPILER_SUPPORTS_NOUNUSEDRES)
if(COMPILER_SUPPORTS_NOUNUSEDRES)
add_definitions(-Wno-unused-result)
endif()

# check dependencies
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

# 64 Bit option
option(64BITMODE "64 bit mode" OFF)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()

add_subdirectory(lib/modKaHIPFiles)


set(MMWIS ${CMAKE_CURRENT_SOURCE_DIR})
set(MMWIS_INCLUDES
    ${MMWIS}/app
    ${MMWIS}/lib
    ${MMWIS}/lib/data_structure
    ${MMWIS}/lib/data_structure/priority_queues
    ${MMWIS}/lib/mis
    ${MMWIS}/lib/mis/ils
    ${MMWIS}/lib/mis/hils
    ${MMWIS}/lib/mis/initial_mis
    ${MMWIS}/lib/mis/kernel
    ${MMWIS}/lib/mis/evolutionary
    ${MMWIS}/lib/mis/evolutionary/combine
    ${MMWIS}/lib/mis/max_flow
    ${MMWIS}/lib/tools
    ${MMWIS}/extern/
    ${MMWIS}/extern/argtable3-3.2.1
    ${MMWIS}/extern/struction/lib/mis/kernel
    ${MMWIS}/extern/struction/lib/data_structure
    ${KaHIP_INCLUDES}
)

set(LIB_MMWIS
    ${MMWIS}/lib/tools/mmwis_log.cpp
    ${MMWIS}/lib/mis/ils/ils.cpp
    ${MMWIS}/lib/mis/ils/local_search.cpp
    ${MMWIS}/lib/mis/hils/Solution.cpp
    ${MMWIS}/lib/mis/initial_mis/greedy_mis.cpp
    ${MMWIS}/lib/mis/initial_mis/greedy_vertex.cpp
    ${MMWIS}/lib/mis/initial_mis/greedy_mwis.cpp
    ${MMWIS}/lib/mis/initial_mis/greedy_weighted_vertex.cpp
    ${MMWIS}/lib/mis/initial_mis/random_mis.cpp
    ${MMWIS}/lib/mis/initial_mis/initial_mis.cpp
    ${MMWIS}/lib/mis/initial_mis/cyclicFast.cpp
    ${MMWIS}/lib/data_structure/mis_permutation.cpp
    ${MMWIS}/lib/data_structure/candidate_list.cpp
    ${MMWIS}/lib/data_structure/operation_log.cpp
    ${MMWIS}/lib/data_structure/priority_queues/bucket_array.cpp
    ${MMWIS}/lib/mis/kernel/branch_and_reduce_algorithm.cpp
    ${MMWIS}/lib/mis/kernel/cout_handler.cpp
    ${MMWIS}/lib/mis/kernel/reductions.cpp
    ${MMWIS}/extern/struction/lib/mis/kernel/struction_reductions.cpp
    ${MMWIS}/extern/struction/lib/mis/kernel/struction_branch_and_reduce_algorithm.cpp
    ${MMWIS}/extern/struction/lib/mis/kernel/extended_struction.cpp
    ${MMWIS}/extern/struction/lib/mis/kernel/key_functions.cpp
    ${MMWIS}/extern/struction/lib/mis/kernel/mwis_finder.cpp
    ${MMWIS}/extern/struction/lib/mis/kernel/original_struction.cpp
    ${MMWIS}/lib/mis/max_flow/bipartite_weighted_vertex_cover.cpp
    ${MMWIS}/lib/mis/evolutionary/population_mis.cpp
    ${MMWIS}/lib/mis/evolutionary/reduction_evolution.cpp
    ${MMWIS}/lib/mis/evolutionary/separator_pool.cpp
    ${MMWIS}/lib/mis/evolutionary/combine/combine.cpp
    ${MMWIS}/lib/mis/evolutionary/combine/cover_combine.cpp
    ${MMWIS}/lib/mis/evolutionary/combine/separator_combine.cpp
    ${MMWIS}/lib/mis/evolutionary/combine/multiway_combine.cpp
)   
add_library(lib_mmwis OBJECT ${LIB_MMWIS})
target_include_directories(lib_mmwis PRIVATE ${KaHIP_INCLUDES})
target_include_directories(lib_mmwis PRIVATE ${MMWIS_INCLUDES})

add_dependencies(lib_mmwis libkaffpa)

option(STRUCTION_ENABLE_APPS "enable/disable struction executables")
add_subdirectory(extern/struction EXCLUDE_FROM_ALL)

add_executable(mmwis "app/mmwis.cpp" $<TARGET_OBJECTS:libkaffpa> 
                                     $<TARGET_OBJECTS:libmapping> 
                                     $<TARGET_OBJECTS:libnodeordering>
                                     $<TARGET_OBJECTS:lib_mmwis> )
target_include_directories(mmwis PRIVATE ${MMWIS_INCLUDES})
target_link_libraries(mmwis ${OpenMP_CXX_LIBRARIES})

install(TARGETS mmwis DESTINATION bin)

