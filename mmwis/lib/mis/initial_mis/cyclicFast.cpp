/**
 * cyclicFast.cpp
 * Purpose: Compute an_and_reduce_algorithm reducer(G, mis_config);initial solution (maximum weight independent set)
 *          by using a greedy algorithm, that always picks the node 
 *          with the largest residual weight.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#include "cyclicFast.h"
#include "struction_branch_and_reduce_algorithm.h"
#include "struction/app/configuration_struction.h"
#include "random_functions.h"

#include "data_structure/priority_queues/bucket_array.h"

using namespace mmwis;

cyclicFast::cyclicFast() {

}

cyclicFast::~cyclicFast() {

}

void cyclicFast::initial_partition(const unsigned int seed, graph_access & G) {
    std::cerr << "for cyclicFast initial solutions use initial_partition_struction" << std::endl;
    exit(1);
}


bool cyclicFast::initial_partition_struction(MISConfig & config, graph_access & G, double remaining_time) {
    random_functions::setSeed(config.seed);

    //set config for cyclicFast:
    MISConfig struction_misconfig;
    configuration_struction cfg;
    cfg.cyclicFast(struction_misconfig);

    struction_misconfig.seed=config.seed;
    if (remaining_time > 60)
        struction_misconfig.time_limit=60;
    else 
        struction_misconfig.time_limit=remaining_time;
    struction_misconfig.ils_time_limit=struction_misconfig.time_limit;;

    struction::cout_handler::disable_cout();


    struction::branch_and_reduce_algorithm reducer(G, struction_misconfig);
    reducer.run_branch_reduce();
    reducer.apply_branch_reduce_solution(G);

    struction::cout_handler::enable_cout();
    return !reducer.timeout;
}


void cyclicFast::generate_permutation(graph_access & G, NodePermutationMap & permutation) {
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);
}

