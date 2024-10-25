/**
 * greedy_weighted_vertex.cpp
 * Purpose: Compute an initial solution (maximum independent set) by 
 *          building a minimal vertex cover.
 *
 *****************************************************************************/

#include "greedy_weighted_vertex.h"

#include "random_functions.h"
#include "data_structure/priority_queues/MaxHeap.h"

greedy_weighted_vertex::greedy_weighted_vertex() {

}

greedy_weighted_vertex::~greedy_weighted_vertex() {

}

void greedy_weighted_vertex::initial_partition(const unsigned int seed, graph_access & G) {
    random_functions::setSeed(seed);
    NodePermutationMap permutation;
    generate_permutation(G, permutation);

    MaxHeap<NodeWeight> *pq = new MaxHeap<NodeWeight>;
    int *uncovered_edges = new int[G.number_of_nodes()];
    G.set_partition_count(2);

    // Initialize the priority queue
    NodeWeight max_weight = G.getMaxWeight();
    forall_nodes (G, n) {
        NodeID node = permutation[n];
        pq->insert(node, max_weight - G.getNodeWeight(node));
        uncovered_edges[node] = G.getNodeDegree(node);
        G.setPartitionIndex(node, 0);
    } endfor

    // While there are still vertices left
    while (pq->size()) {
        // Pick the node with the smallest weight
        NodeID start_node = pq->maxElement();
        NodeWeight w = pq->maxValue();
        pq->deleteMax();
        // Insert it into the solution
        if (uncovered_edges[start_node] > 0) {
                G.setPartitionIndex(start_node, 1);
                uncovered_edges[start_node] = 0;
        } 
        
        forall_out_edges(G, edge, start_node) {
            NodeID target = G.getEdgeTarget(edge);
            // Remove any neighbors that are now completely covered
            if (pq->contains(target)) {
                uncovered_edges[target]--;
                if (uncovered_edges[target] <= 0) {
                    pq->deleteNode(target);
                }
            }
        } endfor
    }

    // Build complement
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 0) G.setPartitionIndex(node, 1); 
        else G.setPartitionIndex(node, 0);
    } endfor

    delete pq;
    delete[] uncovered_edges;
}

void greedy_weighted_vertex::generate_permutation(graph_access & G, NodePermutationMap & permutation) {
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);
}

