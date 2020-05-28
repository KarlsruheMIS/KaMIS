/**
 * greedy_mis_online.cpp
 * Purpose: Compute an initial solution (maximum independent set)
 *          by using a greedy algorithm, that always picks the node
 *          with the smallest residual degree.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
  *****************************************************************************/

#include "greedy_mis_online.h"

#include "random_functions.h"

#include "data_structure/priority_queues/bucket_array.h"

greedy_mis_online::greedy_mis_online() {}

greedy_mis_online::~greedy_mis_online() {}

void greedy_mis_online::initial_partition(const unsigned int seed, graph_access &G) {
  random_functions::setSeed(seed);
  NodePermutationMap permutation;
  generate_permutation(G, permutation);

  bucket_array *buckets = new bucket_array(G.number_of_nodes());
  G.set_partition_count(2);

  // Initialize the priority queue
  forall_nodes(G, n) {
    NodeID node = permutation[n];
    EdgeWeight node_degree = G.getNodeDegree(node);
    buckets->increment(node, node_degree);
    G.setPartitionIndex(node, 0);
  }
  endfor

      // While there are still free vertices left
      while (1) {
    // Select the vertex with the least number of free neighbors
    int start_node = buckets->pickSmallest();
    if (start_node == -1) break;
    G.setPartitionIndex(start_node, 1);
    buckets->remove(start_node);

    forall_out_edges(G, edge, start_node) {
      NodeID target = G.getEdgeTarget(edge);
      // Check if the neighbor is still free
      if (buckets->contains(target)) {
        // Remove neighbor since its no longer free
        buckets->remove(target);
        // Update the surrounding vertices by decreasing
        // its count of free neighbors
        forall_out_edges(G, target_edge, target) {
          NodeID target_neighbor = G.getEdgeTarget(target_edge);
          if (buckets->contains(target_neighbor))
            buckets->decrement(target_neighbor);
        }
        endfor
      }
    }
    endfor
  }

  delete buckets;
}

void greedy_mis_online::initial_partition(const unsigned int seed, graph_access &G,
                                   std::vector<long long> &marked_degree,
                                   long long degree_limit) {
  random_functions::setSeed(seed);
  NodePermutationMap permutation;
  generate_permutation(G, permutation);

  bucket_array *buckets = new bucket_array(G.number_of_nodes());
  G.set_partition_count(2);

  // Initialize the priority queue
  // Modification: "Skip" reduced and high degree vertices
  forall_nodes(G, n) {
    NodeID node = permutation[n];
    EdgeWeight node_degree = marked_degree[node];
    if (node_degree == -1) G.setPartitionIndex(node, 1);
    // else if (node_degree == -2 || node_degree > degree_limit)
    // G.setPartitionIndex(node, 0);
    else if (node_degree == -2)
      G.setPartitionIndex(node, 0);
    else {
      buckets->increment(node, node_degree);
      G.setPartitionIndex(node, 0);
    }
  }
  endfor

      // While there are still free vertices left
      while (1) {
    // Select the vertex with the least number of free neighbors
    int start_node = buckets->pickSmallest();
    if (start_node == -1) break;
    G.setPartitionIndex(start_node, 1);
    buckets->remove(start_node);

    forall_out_edges(G, edge, start_node) {
      NodeID target = G.getEdgeTarget(edge);
      // Check if the neighbor is still free
      if (buckets->contains(target)) {
        // Remove neighbor since its no longer free
        buckets->remove(target);
        // Update the surrounding vertices by decreasing
        // its count of free neighbors
        forall_out_edges(G, target_edge, target) {
          NodeID target_neighbor = G.getEdgeTarget(target_edge);
          if (buckets->contains(target_neighbor))
            buckets->decrement(target_neighbor);
        }
        endfor
      }
    }
    endfor
  }

  delete buckets;
}

void greedy_mis_online::generate_permutation(graph_access &G,
                                      NodePermutationMap &permutation) {
  permutation.resize(G.number_of_nodes());
  random_functions::permutate_vector_good(permutation, true);
}
