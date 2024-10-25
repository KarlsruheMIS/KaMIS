/**
 * random_mis.cpp
 * Purpose: Compute an initial solution (maximum independent set)
 *          by using a simple random selection of vertices. 
 *
 *****************************************************************************/

#include "random_mis.h"

#include "random_functions.h"

random_mis::random_mis() {

}

random_mis::~random_mis() {

}

void random_mis::initial_partition(const unsigned int seed, graph_access & G) {
    random_functions::setSeed(seed);
    NodePermutationMap permutation;
    generate_permutation(G, permutation);
    G.set_partition_count(2);

    // Initialize the partition index
    forall_nodes (G, n) {
        G.setPartitionIndex(n, 0);
    } endfor

    // Random selection
    forall_nodes (G, n) {
        NodeID node = permutation[n];
        bool insert = true;
        forall_out_edges(G, edge, node) {
            NodeID target = G.getEdgeTarget(edge);
            if (G.getPartitionIndex(target) == 1) insert = false;
        } endfor
        if (insert) G.setPartitionIndex(node, 1);
        else G.setPartitionIndex(node, 0);
    } endfor
}

void random_mis::generate_permutation(graph_access & G, NodePermutationMap & permutation) {
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);
}

