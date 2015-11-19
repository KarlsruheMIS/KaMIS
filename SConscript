import platform
import sys

# Get the current platform.
SYSTEM = platform.uname()[0]

Import('env')

# Library files needed for EvoMIS
libfiles = [    
                'lib/tools/mis_log.cpp',
                'lib/mis/ils/ils.cpp',
                'lib/mis/ils/local_search.cpp',
                'lib/mis/initial_mis/greedy_mis.cpp',
                'lib/mis/initial_mis/greedy_vertex.cpp',
                'lib/mis/initial_mis/random_mis.cpp',
                'lib/mis/initial_mis/initial_mis.cpp',
                'lib/data_structure/mis_permutation.cpp',
                'lib/data_structure/candidate_list.cpp',
                'lib/data_structure/operation_log.cpp',
                'lib/data_structure/priority_queues/bucket_array.cpp',
                'lib/mis/evolutionary/population_mis.cpp',
                'lib/mis/evolutionary/reduction_evolution.cpp',
                'lib/mis/hopcroft/bipartite_vertex_cover.cpp',
                'lib/mis/kernel/branch_and_reduce_algorithm.cpp',
                'lib/mis/kernel/modified.cpp',
                'lib/mis/evolutionary/separator_pool.cpp',
                'lib/mis/evolutionary/combine/combine.cpp',
                'lib/mis/evolutionary/combine/cover_combine.cpp',
                'lib/mis/evolutionary/combine/separator_combine.cpp',
                'lib/mis/evolutionary/combine/multiway_combine.cpp',
                ]

libkahip_files = [   'extern/KaHIP//lib/data_structure/graph_hierarchy.cpp',
                      'extern/KaHIP//lib/algorithms/strongly_connected_components.cpp',
                      'extern/KaHIP//lib/algorithms/topological_sort.cpp',
                      'extern/KaHIP//lib/algorithms/push_relabel.cpp',
                      'extern/KaHIP//lib/io/graph_io.cpp',
                      'extern/KaHIP//lib/tools/quality_metrics.cpp',
                      'extern/KaHIP//lib/tools/random_functions.cpp',
                      'extern/KaHIP//lib/tools/graph_extractor.cpp',
                      'extern/KaHIP//lib/tools/misc.cpp',
                      'extern/KaHIP//lib/tools/partition_snapshooter.cpp',
                      'extern/KaHIP//lib/partition/graph_partitioner.cpp',
                      'extern/KaHIP//lib/partition/w_cycles/wcycle_partitioner.cpp',
                      'extern/KaHIP//lib/partition/coarsening/coarsening.cpp',
                      'extern/KaHIP//lib/partition/coarsening/contraction.cpp',
                      'extern/KaHIP//lib/partition/coarsening/edge_rating/edge_ratings.cpp',
                      'extern/KaHIP//lib/partition/coarsening/clustering/node_ordering.cpp',
                      'extern/KaHIP//lib/partition/coarsening/clustering/size_constraint_label_propagation.cpp',
                      'extern/KaHIP//lib/partition/coarsening/matching/matching.cpp',
                      'extern/KaHIP//lib/partition/coarsening/matching/random_matching.cpp',
                      'extern/KaHIP//lib/partition/coarsening/matching/gpa/path.cpp',
                      'extern/KaHIP//lib/partition/coarsening/matching/gpa/gpa_matching.cpp',
                      'extern/KaHIP//lib/partition/coarsening/matching/gpa/path_set.cpp',
                      'extern/KaHIP//lib/partition/initial_partitioning/initial_partitioning.cpp',
                      'extern/KaHIP//lib/partition/initial_partitioning/initial_partitioner.cpp',
                      'extern/KaHIP//lib/partition/initial_partitioning/initial_partition_bipartition.cpp',
                      'extern/KaHIP//lib/partition/initial_partitioning/initial_refinement/initial_refinement.cpp',
                      'extern/KaHIP//lib/partition/initial_partitioning/bipartition.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/uncoarsening.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/separator/vertex_separator_algorithm.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/separator/vertex_separator_flow_solver.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/cycle_improvements/greedy_neg_cycle.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/cycle_improvements/problem_factory.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/refinement.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/mixed_refinement.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/two_way_flow_refinement.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/boundary_bfs.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/most_balanced_minimum_cuts/most_balanced_minimum_cuts.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/cut_flow_problem_solver.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/quotient_graph_scheduling.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/active_block_quotient_graph_scheduler.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph_fabric.cpp', 
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/cycle_improvements/advanced_models.cpp', 
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.cpp', 
                      'extern/KaHIP//lib/algorithms/cycle_search.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.cpp',
                      'extern/KaHIP//lib/partition/uncoarsening/refinement/tabu_search/tabu_search.cpp',
                      'extern/KaHIP/interface/kaHIP_interface.cpp'
                      ]

if SYSTEM == 'Darwin':
    env['CXX'] = 'g++-4.8'

libkahip = env.Library('libkahip', libkahip_files)

if env['program'] == 'redumis':
        env.Program('redumis', ['app/reduction_evomis.cpp']+libfiles, LIBS=[libkahip, 'libargtable2', 'gomp'])
if env['program'] == 'graph_checker':
        env.Program('graphchecker', ['app/graphchecker.cpp']+libfiles, LIBS=[libkahip, 'libargtable2', 'gomp'])
