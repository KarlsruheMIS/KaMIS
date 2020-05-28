/**
 * mis_permutation.cpp
 * Purpose: Data structure used for the local search algorithm.
 *          The structure itself is a permutation of all nodes divided into
 *three blocks. The first block contains the solution nodes. The second block
 *contains free nodes. The third block contains all remaining
 *(non-free, non-solution) nodes.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
  *****************************************************************************/

#include "mis_permutation_online.h"

#include <stdio.h>
#include <algorithm>

#include "data_structure/operation_log.h"
#include "macros_assertions.h"

mis_permutation_online::mis_permutation_online() {}

mis_permutation_online::~mis_permutation_online() {}

void mis_permutation_online::construct(graph_access &G,
                                std::vector<long long> &marked_degree,
                                long long degree_limit) {
  solution_size = 0;
  total_size = G.number_of_nodes();
  nodes.clear();
  tightness.clear();
  position.clear();
  nodes.resize(total_size);
  tightness.resize(total_size);
  position.resize(total_size);
  onetight_all.init(G.number_of_nodes());

  // Insert solution nodes
  // Modification: Unnecessary to iterate over next neighbors
  folded_vertices = 0;
  free_size = total_size;
  forall_nodes(G, n) {
    nodes[n] = n;
    position[n] = n;
    tightness[n] = 0;
  } endfor
}

int mis_permutation_online::calculate_tightness(NodeID node, graph_access &G,
                                         std::vector<long long> &marked_degree,
                                         long long degree_limit) {
  int tightness = 0;
  forall_out_edges(G, edge, node) {
    NodeID target = G.getEdgeTarget(edge);

    // Modification: Skip
    if (marked_degree[target] < 0) continue;
    bool target_index = is_solution_node(target);
    if (target_index == 1) {
      tightness++;
    }
  } endfor 

  return tightness;
}

NodeID mis_permutation_online::get_solution_node(unsigned int i) {
  ASSERT_LT(i, solution_size);
  return nodes[i];
}

NodeID mis_permutation_online::get_free_node(unsigned int i) {
  ASSERT_LT(i, free_size);
  return nodes[solution_size + i];
}

NodeID mis_permutation_online::get_non_free_node(unsigned int i) {
  ASSERT_LT(i, total_size - (free_size + solution_size));
  return nodes[solution_size + free_size + i];
}

NodeID mis_permutation_online::get_non_solution_node(unsigned int i) {
  ASSERT_LT(i, total_size - solution_size);
  return nodes[solution_size + i];
}

unsigned int mis_permutation_online::get_solution_size() {
  return solution_size;
}

unsigned int mis_permutation_online::get_free_size() { return free_size; }

bool mis_permutation_online::is_maximal() { return (free_size == 0); }

bool mis_permutation_online::search_adj(NodeID source, NodeID target,
                                 graph_access &G) {
  // Binary search on boundaries
  EdgeID low = G.get_first_edge(source);
  EdgeID high = G.get_first_invalid_edge(source) - 1;
  EdgeID mid = 0;
  while (low <= high) {
    mid = low + (high - low) / 2;
    if (target == G.getEdgeTarget(mid))
      return true;
    else if (target < G.getEdgeTarget(mid))
      high = mid - 1;
    else
      low = mid + 1;
  }
  return false;
}

bool mis_permutation_online::add_and_reduce(NodeID node,
                                     std::vector<long long> &marked_degree,
                                     long long degree_limit, graph_access &G) {
  if (marked_degree[node] < 0) return true;

  // Modification: Pendant vertex
  if (marked_degree[node] <= 2 && marked_degree[node] >= 0) {
    // Degree zero
    if (marked_degree[node] == 0) {
      marked_degree[node] = -1;
    }
    // Degree one
    else if (marked_degree[node] == 1) {
#ifndef NDEBUG
      NodeID u = 0;
      forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);
        if (marked_degree[target] >= 0) {
          u = target;
          break;
        }
      } endfor
      if (is_solution_node(u)) {
        fprintf(stderr, "Degree one - This shouldn't happen\n");
        fprintf(stderr, "x=%d m=%d.\n", node, marked_degree[node]);
        fprintf(stderr, "u=%d m=%d.\n", u, marked_degree[u]);
        exit(1);
      }
#endif
      marked_degree[node] = -1;
    }
    // Degree two
    else {
      NodeID u = 0;
      EdgeID u_edge = 0;
      NodeID w = 0;
      EdgeID w_edge = 0;
      forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);
        if (u > 0 && marked_degree[target] >= 0) {
          w = target;
          break;
        }
        if (marked_degree[target] >= 0) {
          u = target;
        }
      } endfor
      if (is_solution_node(u) || is_solution_node(w)) {
#ifndef NDEBUG
        fprintf(stderr, "Degree two - This shouldn't happen\n");
        fprintf(stderr, "x=%d m=%d s=%d.\n", node, marked_degree[node], is_solution_node(node));
        fprintf(stderr, "u=%d m=%d s=%d.\n", u, marked_degree[u], is_solution_node(u));
        fprintf(stderr, "w=%d m=%d s=%d.\n", w, marked_degree[w], is_solution_node(w));
        exit(1);
#endif
      }
      else {
        // Search w in adj[u], then u - w
        // u - v and w - v by definition of edges
        if (search_adj(w, u, G)) {
          marked_degree[node] = -1;
        }
        // Path reduction
        // else if (!search_adj(w, u, G) && marked_degree[w] == 2 &&
        //          marked_degree[u] == 2) {
        //   // Gather other neighbors
        //   NodeID u_prime = G.number_of_nodes();
        //   EdgeID u_prime_edge = 0;
        //   NodeID w_prime = G.number_of_nodes();
        //   EdgeID w_prime_edge = 0;

        //   // Find u_prime
        //   forall_out_edges(G, edge, u) {
        //     NodeID target = G.getEdgeTarget(edge);
        //     if (target != node && marked_degree[target] >= 0) {
        //       u_prime = target;
        //       break;
        //     }
        //   } endfor

        //   // Find w_prime
        //   forall_out_edges(G, edge, w) {
        //     NodeID target = G.getEdgeTarget(edge);
        //     if (target != node && marked_degree[target] >= 0) {
        //       w_prime = target;
        //       break;
        //     }
        //   } endfor

        //   // Stop if no valid next neighbors
        //   if (u_prime == G.number_of_nodes() || w_prime == G.number_of_nodes()) return false;

        //   // Find (node, u)
        //   forall_out_edges(G, edge, node) {
        //     NodeID target = G.getEdgeTarget(edge);
        //     if (target == u) {
        //       u_edge = edge;
        //       break;
        //     }
        //   } endfor

        //   // Find (u_prime, u)
        //   forall_out_edges(G, edge, u_prime) {
        //     NodeID target = G.getEdgeTarget(edge);
        //     if (target == u) {
        //       u_prime_edge = edge;
        //       break;
        //     }
        //   } endfor

        //   // Relink (u_prime, u) -> (u_prime, node);
        //   G.setEdgeTarget(u_prime_edge, u_prime, node);
        //   // Relink (node, u) -> (node, u_prime);
        //   G.setEdgeTarget(u_edge, node, u_prime);

        //   // Find (node, w)
        //   forall_out_edges(G, edge, node) {
        //     NodeID target = G.getEdgeTarget(edge);
        //     if (target == w) {
        //       w_edge = edge;
        //       break;
        //     }
        //   } endfor

        //   // Find (w_prime, w)
        //   forall_out_edges(G, edge, w_prime) {
        //     NodeID target = G.getEdgeTarget(edge);
        //     if (target == w) {
        //       w_prime_edge = edge;
        //       break;
        //     }
        //   } endfor

        //   if (w_prime != u_prime) {
        //     // Relink (w_prime, w) -> (w_prime, node);
        //     G.setEdgeTarget(w_prime_edge, w_prime, node);
        //     // Relink (node, w) -> (node, w_prime);
        //     G.setEdgeTarget(w_edge, node, w_prime);
        //   }
        //   else {
        //     // "Diamond" reduction
        //     marked_degree[node]--;
        //     marked_degree[u_prime]--;
        //   }

        //   // Remove u and w and store fold
        //   marked_degree[u] = -2;
        //   marked_degree[w] = -2;
        //   move_to_non_free(u, G);
        //   move_to_non_free(w, G);
        //   tightness[u] = 0;
        //   tightness[w] = 0;
        //   if (tightness[u] == 1) onetight_all.remove(u);
        //   if (tightness[w] == 1) onetight_all.remove(w);
        //   folded_vertices++;
        //   folds.emplace(node, u, w);

        //   tightness[node] = 0;
        //   if (is_solution_node(u_prime) || is_solution_node(w_prime)) {
        //     move_to_non_free(node, G);
        //     if (is_solution_node(u_prime)) tightness[node]++;
        //     if (is_solution_node(w_prime)) tightness[node]++;
        //     if (tightness[node] == 1 && !onetight_all.contains(node))
        //       onetight_all.insert(node);
        //     else if (onetight_all.contains(node))
        //       onetight_all.remove(node);
        //   } else move_to_free(node, G);

        //   return true;
        // }
      }
    }
  }

  // Skip if vertex could not be reduced
  if (marked_degree[node] != -1) return false;

  move_to_solution(node, G);
  if (tightness[node] == 1) onetight_all.remove(node);
  tightness[node] = 0;

  // Update information about v's neighbors
  forall_out_edges(G, edge, node) {
    NodeID target = G.getEdgeTarget(edge);

    // Modification: Set marked degree
    if (marked_degree[target] > 0) {
      marked_degree[target] = -2;
      bool target_in_solution = is_solution_node(target);
      // Modification: Make sure neighbor is not free anymore
      move_to_non_free(target, G);
      if (tightness[target] == 1) onetight_all.remove(target);
      tightness[target] = 0;
      // Modification: Update degree of next neighbors
      forall_out_edges(G, target_edge, target) {
        NodeID ntarget = G.getEdgeTarget(target_edge);
        if (marked_degree[ntarget] > 0) {
          marked_degree[ntarget]--;
          // Update tightness of next neighbor
          if (target_in_solution) {
            if (tightness[ntarget] == 1) onetight_all.remove(ntarget);
            tightness[ntarget]--;
            if (tightness[ntarget] == 1) onetight_all.insert(ntarget);
            if (tightness[ntarget] == 0) move_to_free(ntarget, G);
          }
        }
      } endfor
    }
  } endfor

  return true;
}

void mis_permutation_online::add_to_solution(NodeID node,
                                      std::vector<long long> &marked_degree,
                                      long long degree_limit, graph_access &G) {
  // Node is already in the solution
  if (is_solution_node(node)) return;

  // Modification: Reduction
  if (add_and_reduce(node, marked_degree, degree_limit, G)) return;

  // Add node to the solution
  move_to_solution(node, G);

  // Update neighbors
  forall_out_edges(G, edge, node) {
    NodeID target = G.getEdgeTarget(edge);

    // Modification: Skip
    if (marked_degree[target] < 0) continue;

    ASSERT_TRUE(!is_solution_node(target));
    if (tightness[target] == 0) move_to_non_free(target, G);

    if (tightness[target] == 1) onetight_all.remove(target);
    tightness[target]++;
    if (tightness[target] == 1) onetight_all.insert(target);
  } endfor

      ASSERT_LT(position[node], solution_size);
  operation_log::instance()->report_insert(node);
  return;
}

void mis_permutation_online::remove_from_solution(
    NodeID node, std::vector<long long> &marked_degree, long long degree_limit,
    graph_access &G) {
  // Node isn't in the solution
  if (!is_solution_node(node)) return;

  // Modification: Skip
  if (marked_degree[node] < 0) return;

  // Remove node from the solution
  bool is_free = (tightness[node] == 0);
  if (is_free)
    move_to_free(node, G);
  else
    move_to_non_free(node, G);

  // Update neighbors
  forall_out_edges(G, edge, node) {
    NodeID target = G.getEdgeTarget(edge);

    // Modification: Skip
    if (marked_degree[target] < 0) continue;

    if (tightness[target] == 1) onetight_all.remove(target);
    tightness[target]--;
    if (tightness[target] == 1) onetight_all.insert(target);

    ASSERT_TRUE(!is_solution_node(node));
    if (tightness[target] == 0) move_to_free(target, G);
  } endfor

  ASSERT_GEQ(position[node], solution_size);
  operation_log::instance()->report_remove(node);
}

void mis_permutation_online::move_to_solution(NodeID node, graph_access &G) {
  unsigned int node_position = position[node];

  // Node already in the solution then skip
  if (node_position < solution_size) return;
  // Is the node free?
  if (node_position < solution_size + free_size) {
    swap_nodes(solution_size, node_position);
    free_size--;
  }
  // Non free node
  else {
    swap_nodes(solution_size + free_size, position[node]);
    swap_nodes(solution_size, solution_size + free_size);
  }
  solution_size++;
}

void mis_permutation_online::move_to_free(NodeID node, graph_access &G) {
  unsigned int node_position = position[node];

  // Is the node in the solution?
  if (node_position < solution_size) {
    swap_nodes(solution_size - 1, node_position);
    solution_size--;
  }
  // Node already free then skip
  else if (node_position < solution_size + free_size)
    return;
  // Non free node
  else {
    swap_nodes(solution_size + free_size, position[node]);
  }
  free_size++;
}

void mis_permutation_online::move_to_non_free(NodeID node, graph_access &G) {
  unsigned int node_position = position[node];

  // printf("Move %d to non free\n", node);
  // Is the node already non free
  if (node_position >= solution_size + free_size) return;
  // Is the node free?
  if (node_position >= solution_size) {
    swap_nodes(solution_size + free_size - 1, node_position);
    free_size--;
  }
  // Solution node
  else {
    swap_nodes(solution_size - 1, position[node]);
    swap_nodes(solution_size + free_size - 1, solution_size - 1);
    solution_size--;
  }
}

void mis_permutation_online::swap_nodes(unsigned int first_pos,
                                 unsigned int second_pos) {
  if (first_pos == second_pos) return;

  NodeID first_node = nodes[first_pos];
  NodeID second_node = nodes[second_pos];
  nodes[first_pos] = second_node;
  nodes[second_pos] = first_node;
  position[first_node] = second_pos;
  position[second_node] = first_pos;
}

bool mis_permutation_online::is_solution_node(NodeID node) {
  return (position[node] < solution_size);
}

bool mis_permutation_online::is_free_node(NodeID node) {
  return (position[node] >= solution_size &&
          position[node] < free_size + solution_size);
}

bool mis_permutation_online::is_non_solution_node(NodeID node) {
  return (position[node] >= solution_size + free_size &&
          position[node] < total_size);
}

void mis_permutation_online::print(bool details) {
  printf("\n");
  printf("********************\n");
  printf("Solution size: %d\n", solution_size);
  printf("Free size: %d\n", free_size);
  printf("Maximal: %d\n", is_maximal());
  if (details) {
    printf("\n");
    printf("Solution nodes: \n");
    for (unsigned int i = 0; i < solution_size; ++i) printf(" %d ", nodes[i]);
    printf("\n");

    printf("Free nodes: \n");
    for (unsigned int i = solution_size; i < solution_size + free_size; ++i)
      printf(" %d ", nodes[i]);
    printf("\n");

    printf("Remaining nodes: \n");
    for (unsigned int i = solution_size + free_size; i < total_size; ++i)
      printf(" %d ", nodes[i]);
    printf("\n");
  }
  printf("********************\n");
  printf("\n");
}

void mis_permutation_online::print_position() {
  printf("\n");
  for (unsigned int i = 0; i < solution_size; ++i)
    printf("Position: %d\n", position[i]);
  for (unsigned int i = solution_size; i < solution_size + free_size; ++i)
    printf("Position: %d\n", position[i]);
  for (unsigned int i = solution_size + free_size; i < total_size; ++i)
    printf("Position: %d\n", position[i]);
  printf("\n");
}

void mis_permutation_online::print_tightness() {
  printf("\n");
  for (unsigned int i = 0; i < solution_size; ++i)
    printf("Tightness: %d\n", tightness[i]);
  for (unsigned int i = solution_size; i < solution_size + free_size; ++i)
    printf("Tightness: %d\n", tightness[i]);
  for (unsigned int i = solution_size + free_size; i < total_size; ++i)
    printf("Tightness: %d\n", tightness[i]);
  printf("\n");
}

bool mis_permutation_online::check_permutation() {
  bool solution_correct = true;

  for (unsigned int i = 0; i < solution_size; ++i) {
    NodeID node = nodes[i];
    if (tightness[node] != 0) {
      printf("Tightness error in solution\n");
      solution_correct = false;
    }
  }

  for (unsigned int i = solution_size; i < solution_size + free_size; ++i) {
    NodeID node = nodes[i];
    if (tightness[node] != 0) {
      printf("Tightness error in free\n");
      solution_correct = false;
    }
  }

  for (unsigned int i = solution_size + free_size; i < total_size; ++i) {
    NodeID node = nodes[i];
    if (tightness[node] < 1) {
      printf("Tightness error in non free\n");
      solution_correct = false;
    }
  }

  return solution_correct;
}

bool mis_permutation_online::check_consistency(graph_access &G) {
  forall_nodes(G, node) {
    if (is_solution_node(node) != 1) return false;
    if (!is_solution_node(node) != 0) return false;
  } endfor 
  return true;
}
