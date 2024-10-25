cmake_minimum_required(VERSION 3.10)

set(KaHIP ${CMAKE_CURRENT_SOURCE_DIR}/../../extern/KaHIP)
set(modKaHIP ${CMAKE_CURRENT_SOURCE_DIR})
set(KaHIP_INCLUDES
    ${modKaHIP}/interface
    ${modKaHIP}/lib
    ${modKaHIP}/lib/data_structure
    ${modKaHIP}/lib/algorithms
    ${modKaHIP}/lib/tools
    ${modKaHIP}/lib/io
    ${modKaHIP}/../data_structure
    ${modKaHIP}/../data_structure/priority_queues
    ${KaHIP}/app
    ${KaHIP}/lib
    ${KaHIP}/lib/io
    ${KaHIP}/lib/tools
    ${KaHIP}/lib/algorithms
    ${KaHIP}/lib/partition
    ${KaHIP}/lib/data_structure/priority_queues
    ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement
)
set(KaHIP_INCLUDES ${KaHIP_INCLUDES} PARENT_SCOPE)

set(LIBKAFFPA_SOURCE_FILES
  # ${modKaHIP}/lib/io/graph_io.cpp
  # ${modKaHIP}/lib/tools/graph_extractor.cpp
  ${modKaHIP}/lib/algorithms/strongly_connected_components.cpp
  ${modKaHIP}/interface/kaHIP_interface.cpp
  ${KaHIP}/lib/data_structure/graph_hierarchy.cpp
  ${KaHIP}/lib/algorithms/topological_sort.cpp
  ${KaHIP}/lib/algorithms/push_relabel.cpp
  ${KaHIP}/lib/tools/quality_metrics.cpp
  ${KaHIP}/lib/tools/misc.cpp
  ${KaHIP}/lib/tools/partition_snapshooter.cpp
  ${KaHIP}/lib/tools/random_functions.cpp
  ${KaHIP}/lib/partition/graph_partitioner.cpp
  ${KaHIP}/lib/partition/w_cycles/wcycle_partitioner.cpp
  ${KaHIP}/lib/partition/coarsening/coarsening.cpp
  ${KaHIP}/lib/partition/coarsening/contraction.cpp
  ${KaHIP}/lib/partition/coarsening/edge_rating/edge_ratings.cpp
  ${KaHIP}/lib/partition/coarsening/matching/matching.cpp
  ${KaHIP}/lib/partition/coarsening/matching/random_matching.cpp
  ${KaHIP}/lib/partition/coarsening/matching/gpa/path.cpp
  ${KaHIP}/lib/partition/coarsening/matching/gpa/gpa_matching.cpp
  ${KaHIP}/lib/partition/coarsening/matching/gpa/path_set.cpp
  ${KaHIP}/lib/partition/coarsening/clustering/node_ordering.cpp
  ${KaHIP}/lib/partition/coarsening/clustering/size_constraint_label_propagation.cpp
  ${KaHIP}/lib/partition/initial_partitioning/initial_partitioning.cpp
  ${KaHIP}/lib/partition/initial_partitioning/initial_partitioner.cpp
  ${KaHIP}/lib/partition/initial_partitioning/initial_partition_bipartition.cpp
  ${KaHIP}/lib/partition/initial_partitioning/initial_refinement/initial_refinement.cpp
  ${KaHIP}/lib/partition/initial_partitioning/bipartition.cpp
  ${KaHIP}/lib/partition/initial_partitioning/initial_node_separator.cpp
  ${KaHIP}/lib/partition/uncoarsening/uncoarsening.cpp
  ${KaHIP}/lib/partition/uncoarsening/separator/area_bfs.cpp
  ${KaHIP}/lib/partition/uncoarsening/separator/vertex_separator_algorithm.cpp
  ${KaHIP}/lib/partition/uncoarsening/separator/vertex_separator_flow_solver.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/cycle_improvements/greedy_neg_cycle.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/cycle_improvements/problem_factory.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/mixed_refinement.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/refinement.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/two_way_flow_refinement.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/boundary_bfs.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/cut_flow_problem_solver.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/most_balanced_minimum_cuts/most_balanced_minimum_cuts.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/quotient_graph_scheduling.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/active_block_quotient_graph_scheduler.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph_fabric.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/cycle_improvements/advanced_models.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/node_separators/greedy_ns_local_search.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/node_separators/fm_ns_local_search.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/node_separators/localized_fm_ns_local_search.cpp
  ${KaHIP}/lib/algorithms/cycle_search.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.cpp
  ${KaHIP}/lib/partition/uncoarsening/refinement/tabu_search/tabu_search.cpp
  ${KaHIP}/extern/argtable3-3.2.2/argtable3.c)
add_library(libkaffpa OBJECT ${LIBKAFFPA_SOURCE_FILES})
target_include_directories(libkaffpa PRIVATE ${KaHIP_INCLUDES})

set(LIBMAPPING_SOURCE_FILES
  ${KaHIP}/lib/mapping/local_search_mapping.cpp
  ${KaHIP}/lib/mapping/full_search_space.cpp
  ${KaHIP}/lib/mapping/full_search_space_pruned.cpp
  ${KaHIP}/lib/mapping/communication_graph_search_space.cpp
  ${KaHIP}/lib/mapping/fast_construct_mapping.cpp
  ${KaHIP}/lib/mapping/construct_distance_matrix.cpp
  ${KaHIP}/lib/mapping/mapping_algorithms.cpp
  ${KaHIP}/lib/mapping/construct_mapping.cpp)
add_library(libmapping OBJECT ${LIBMAPPING_SOURCE_FILES})
target_include_directories(libmapping PRIVATE ${KaHIP_INCLUDES})

set(NODE_ORDERING_SOURCE_FILES
  ${KaHIP}/lib/node_ordering/min_degree_ordering.cpp
  ${KaHIP}/lib/node_ordering/nested_dissection.cpp
  ${KaHIP}/lib/node_ordering/ordering_tools.cpp
  ${KaHIP}/lib/node_ordering/reductions.cpp)
add_library(libnodeordering OBJECT ${NODE_ORDERING_SOURCE_FILES})
target_include_directories(libnodeordering PRIVATE ${KaHIP_INCLUDES})