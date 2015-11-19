/**
 * greedy_vertex.cpp
 * Purpose: Compute an initial solution (maximum independent set) by 
 *          building a minimal vertex cover.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "greedy_vertex.h"

#include "random_functions.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/bucket_array.h"

greedy_vertex::greedy_vertex() {

}

greedy_vertex::~greedy_vertex() {

}

void greedy_vertex::initial_partition(const unsigned int seed, graph_access & G) {
    random_functions::setSeed(seed);
    NodePermutationMap permutation;
    generate_permutation(G, permutation);

    bucket_array *buckets = new bucket_array(G.number_of_nodes());
    int uncovered_edges[G.number_of_nodes()];
    G.set_partition_count(2);

    // Initialize the priority queue
    forall_nodes (G, n) {
        NodeID node = permutation[n];
        buckets->increment(node, G.getMaxDegree() - G.getNodeDegree(node));
        uncovered_edges[node] = G.getNodeDegree(node);
        G.setPartitionIndex(node, 0);
    } endfor

    // While there are still vertices left
    while (1) {
        // Pick the node with the most uncovered edges
        int start_node = buckets->pickSmallest();
        if (start_node == -1) break;
        buckets->remove(start_node);
        uncovered_edges[start_node] = 0;
        // Insert it into the solution
        G.setPartitionIndex(start_node, 1);
        
        forall_out_edges(G, edge, start_node) {
            NodeID target = G.getEdgeTarget(edge);
            // Remove any neighbors that are now completely covered
            if (buckets->contains(target)) {
                buckets->increment(target);
                uncovered_edges[target]--;
                if (uncovered_edges[target] <= 0) {
                    buckets->remove(target);
                }
            }
        } endfor
    }

    // Build complement
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 0) G.setPartitionIndex(node, 1); 
        else G.setPartitionIndex(node, 0);
    } endfor

    delete buckets;
}

void greedy_vertex::generate_permutation(graph_access & G, NodePermutationMap & permutation) {
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);
}

