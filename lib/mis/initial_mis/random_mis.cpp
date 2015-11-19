/**
 * random_mis.cpp
 * Purpose: Compute an initial solution (maximum independent set)
 *          by using a simple random selection of vertices. 
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

