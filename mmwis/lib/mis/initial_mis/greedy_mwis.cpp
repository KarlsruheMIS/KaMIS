/**
 * greedy_mwis.cpp
 * Purpose: Compute an initial solution (maximum weight independent set)
 *          by using a greedy algorithm, that always picks the node 
 *          with the largest residual weight.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#include "greedy_mwis.h"

#include "random_functions.h"

#include "data_structure/priority_queues/bucket_array.h"
#include "data_structure/priority_queues/MaxHeap.h"

greedy_mwis::greedy_mwis() {

}

greedy_mwis::~greedy_mwis() {

}

void greedy_mwis::initial_partition(const unsigned int seed, graph_access & G) {
    random_functions::setSeed(seed);
    NodePermutationMap permutation;
    generate_permutation(G, permutation);

    MaxHeap<NodeWeight> *pq = new MaxHeap<NodeWeight>;
    G.set_partition_count(2);

    // Initialize the priority queue
    forall_nodes (G, n) {
        NodeID node = permutation[n];
        NodeWeight node_weight = G.getNodeWeight(node);
        std::pair<NodeID, NodeWeight> node_pair(node, node_weight);
        pq->insert(node, node_weight);
        G.setPartitionIndex(node, 0);
    } endfor

    // While there are still free vertices left
    while (pq->size()) {
        // Select the vertex with the largest weight
        NodeID start_node = pq->maxElement();
        NodeWeight w = pq->maxValue();
        pq->deleteMax();
        G.setPartitionIndex(start_node, 1);
        
        forall_out_edges (G, edge, start_node) {
           NodeID target = G.getEdgeTarget(edge);
           // Check if the neighbor is still free
           if (pq->contains(target)) {
               // Remove neighbor since its no longer free
               pq->deleteNode(target);
           }
        } endfor
    }
    
    delete pq;

}

void greedy_mwis::generate_permutation(graph_access & G, NodePermutationMap & permutation) {
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);
}

