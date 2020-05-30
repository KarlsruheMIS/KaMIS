cmake_minimum_required(VERSION 2.8.9)
project(components)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release) 
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fno-omit-frame-pointer -march=native -DNDEBUG -funroll-loops -fno-stack-limit -fPIC -O3 ")

CHECK_CXX_COMPILER_FLAG(-Wno-unused-value COMPILER_SUPPORTS_NOUNUSED)
if(COMPILER_SUPPORTS_NOUNUSED)
add_definitions(-Wno-unused-value)
endif()
CHECK_CXX_COMPILER_FLAG(-Wno-unused-value COMPILER_SUPPORTS_NOUNUSEDRES)
if(COMPILER_SUPPORTS_NOUNUSEDRES)
add_definitions(-Wno-unused-result)
endif()

add_subdirectory(extern/KaHIP)

include_directories(
        app
        lib
        lib/data_structure
        lib/data_structure/priority_queues
        lib/mis
        lib/mis/ils
        lib/mis/initial_mis
        lib/mis/kernel
        lib/tools
        extern/KaHIP
        extern/KaHIP/app
        extern/KaHIP/interface
        extern/KaHIP/lib
        extern/KaHIP/lib/algorithms
        extern/KaHIP/lib/data_structure
        extern/KaHIP/lib/data_structure/matrix
        extern/KaHIP/lib/data_structure/priority_queues
        extern/KaHIP/lib/io
        extern/KaHIP/lib/parallel_mh
        extern/KaHIP/lib/parallel_mh/exchange
        extern/KaHIP/lib/parallel_mh/galinier_combine
        extern/KaHIP/lib/partition
        extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement
        extern/KaHIP/lib/tools
)


set(SOURCES
        lib/tools/mis_log.cpp
        lib/mis/ils/ils.cpp
        lib/mis/ils/local_search.cpp
        lib/mis/initial_mis/greedy_mis.cpp
        lib/mis/initial_mis/greedy_vertex.cpp
        lib/mis/initial_mis/random_mis.cpp
        lib/mis/initial_mis/initial_mis.cpp
        lib/data_structure/mis_permutation.cpp
        lib/data_structure/candidate_list.cpp
        lib/data_structure/operation_log.cpp
        lib/data_structure/priority_queues/bucket_array.cpp
        lib/mis/kernel/branch_and_reduce_algorithm.cpp
        lib/mis/kernel/reductions.cpp
        ../extern/argtable3-3.0.3/argtable3.c)

add_library(libsources OBJECT ${SOURCES})

#add_executable(graph_checker "app/graphchecker.cpp" $<TARGET_OBJECTS:libsources>)
add_executable(branch_reduce "app/branch_reduce.cpp" $<TARGET_OBJECTS:libsources>)
add_executable(weighted_ls "app/weighted_ls.cpp" $<TARGET_OBJECTS:libsources>)
add_executable(merge_graph_weights "app/merge_graph_weights.cpp" $<TARGET_OBJECTS:libsources>)

set(LIBS libkaffpa2)

#target_link_libraries(graph_checker LINK_PUBLIC ${LIBS})
target_link_libraries(branch_reduce LINK_PUBLIC $<TARGET_OBJECTS:libkaffpa2>)
target_link_libraries(weighted_ls LINK_PUBLIC $<TARGET_OBJECTS:libkaffpa2>)
target_link_libraries(merge_graph_weights LINK_PUBLIC $<TARGET_OBJECTS:libkaffpa2>)
