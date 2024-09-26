

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fno-omit-frame-pointer -march=native -DNDEBUG -fPIC -O3 ")

include_directories(
        app
        lib
        lib/data_structure
        ../../lib
        ../../lib/data_structure
        ../../lib/data_structure/priority_queues
        ../../lib/mis
        ../../lib/mis/ils
        ../../lib/mis/hils
        ../../lib/mis/initial_mis
        lib/mis/kernel
        ../KaHIP
        ../KaHIP/app
        ../KaHIP/interface
        ../KaHIP/lib
        ../KaHIP/lib/algorithms
        ../KaHIP/lib/data_structure
        ../KaHIP/lib/data_structure/matrix
        ../KaHIP/lib/data_structure/priority_queues
        ../KaHIP/lib/io
        ../KaHIP/lib/parallel_mh
        ../KaHIP/lib/parallel_mh/exchange
        ../KaHIP/lib/parallel_mh/galinier_combine
        ../KaHIP/lib/partition
        ../KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement
        ../KaHIP/lib/tools
        ../argtable-2.10
)

set(LIBSTRUCTION_SOURCE_FILES
# lib/tools/struction_mis_log.cpp
# lib/data_structure/candidate_list.cpp
# lib/data_structure/mis_permutation.cpp 
# lib/data_structure/operation_log.cpp
# lib/mis/initial_mis/greedy_mis.cpp
# lib/mis/initial_mis/greedy_vertex.cpp
# lib/mis/initial_mis/initial_mis.cpp
# lib/mis/initial_mis/random_mis.cpp
# lib/mis/kernel/struction_branch_and_reduce_algorithm.cpp
# lib/mis/kernel/struction_reductions.cpp
# lib/mis/kernel/extended_struction.cpp
# lib/mis/kernel/key_functions.cpp 
# lib/mis/kernel/mwis_finder.cpp 
# lib/mis/kernel/original_struction.cpp 
# lib/mis/ils/local_search.cpp
# ../../wmis/lib/mis/hils/Solution.cpp
)
add_library(libstruction OBJECT ${LIBSTRUCTION_SOURCE_FILES})

