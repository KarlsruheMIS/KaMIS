/**
 * greedy_mis.cpp
 * Purpose: Compute an initial solution (maximum independent set)
 *          by using a greedy algorithm, that always picks the node 
 *          with the smallest residual degree.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
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

#include "greedy_mis.h"

#include "random_functions.h"

#include "data_structure/priority_queues/bucket_array.h"

greedy_mis::greedy_mis() {

}

greedy_mis::~greedy_mis() {

}

void greedy_mis::initial_partition(const unsigned int seed, graph_access & G) {
    random_functions::setSeed(seed);
    NodePermutationMap permutation;
    generate_permutation(G, permutation);

    bucket_array *buckets = new bucket_array(G.number_of_nodes());
    G.set_partition_count(2);

    // Initialize the priority queue
    forall_nodes (G, n) {
        NodeID node = permutation[n];
        EdgeWeight node_degree = G.getNodeDegree(node);
        buckets->increment(node, node_degree);
        G.setPartitionIndex(node, 0);
    } endfor

    // While there are still free vertices left
    while (1) {
        // Select the vertex with the least number of free neighbors 
        int start_node = buckets->pickSmallest();
        if (start_node == -1) break;
        G.setPartitionIndex(start_node, 1);
        buckets->remove(start_node);
        
        forall_out_edges (G, edge, start_node) {
           NodeID target = G.getEdgeTarget(edge);
           // Check if the neighbor is still free
           if (buckets->contains(target)) {
               // Remove neighbor since its no longer free
               buckets->remove(target);
               // Update the surrounding vertices by decreasing
               // its count of free neighbors
               forall_out_edges (G, target_edge, target) {
                    NodeID target_neighbor = G.getEdgeTarget(target_edge);
                    if (buckets->contains(target_neighbor)) buckets->decrement(target_neighbor);
               } endfor
           }
        } endfor
    }
    
    delete buckets;
}

void greedy_mis::generate_permutation(graph_access & G, NodePermutationMap & permutation) {
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);
}

