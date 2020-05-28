/**
 * local_search_online.cpp
 * Purpose: Apply the local search algorithm to a maximum independent set.
 *
 * The original code from Andrade et. al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#include "ils/local_search_online.h"

#include <algorithm>

#include "random_functions.h"
#include "timer.h"

local_search_online::local_search_online() { long_runs = 0; }

local_search_online::~local_search_online() {}

void local_search_online::direct_improvement(graph_access &G,
                                      std::vector<long long> &marked_degree,
                                      long long degree_limit, bool forced,
                                      long long forced_node) {
  if (candidates.is_empty()) build_candidates(G, marked_degree, degree_limit);

  while (1) {
    NodeID x;
    if (candidates.is_empty()) {
      if (forced && forced_node > -1) {
        x = forced_node;
        forced = false;
      } else
        break;
    } else {
      x = candidates.remove_random();
      // Modification: Check
      if (marked_degree[x] < 0 || marked_degree[x] > degree_limit) {
#ifndef NDEBUG
        // This might trigger for path reduced vertices
        fprintf(stderr, "This shouldn't happen.\n");
        fprintf(stderr, "v=%d m=%d.\n", x, marked_degree[x]);
#endif
        continue;
      }
      if (!perm.is_solution_node(x)) continue;
      if (forced && x == forced_node) continue;
    }

    // Build L(x)
    build_onetight(x, G, marked_degree, degree_limit);

    if (onetight_size < 2) continue;
    int solutions_found = 0;
    int max_found = 100;
    NodeID improv_v = 0;
    NodeID improv_w = 0;

    // run_t.restart();
    for (int i = onetight_size - 1; i >= 0; i--) {
      NodeID v = onetight[i];
      // No need to last neighbor if no swap was found so far
      if (i == 0) {
        // No swap can be found
        if (solutions_found == 0) break;
        // Swap would have been found already
        if (onetight_size == 2) break;
      }

      // Fix from Renato
      if (solutions_found >= max_found) break;

      // Build A(v)
      build_neighbors(v, G, marked_degree, degree_limit);
      unsigned int v_pos = 0;
      unsigned int x_pos = 0;

      while (1) {
        // No more neighbors to check
        if (x_pos >= onetight_size) break;

        NodeID x_neighbor = onetight[x_pos];
        // Skip if the neighbor is v itself
        if (x_neighbor == v) {
          x_pos++;
          continue;
        }

        // Reached the end of A(v) so all remaining nodes are possible
        // improvements
        if (v_pos >= neighbors_size) {
          solutions_found++;
          if (solutions_found == 1 ||
              random_functions::nextInt(1, solutions_found) == 1) {
            improv_v = v;
            improv_w = x_neighbor;
          }
          x_pos++;
          continue;
        }

        NodeID v_neighbor = neighbors[v_pos];
        // Skip if the neighbor is x itself
        if (v_neighbor == x) {
          v_pos++;
          continue;
        }

        // Skip if non candidate neighbor
        if (v_neighbor < x_neighbor) {
          v_pos++;
          continue;
        }

        // Skip if both neighbors are the same
        if (v_neighbor == x_neighbor) {
          v_pos++;
          x_pos++;
          continue;
        }

        // Candidate found but still continue looking
        if (v_neighbor > x_neighbor) {
          solutions_found++;
          if (random_functions::nextInt(1, solutions_found) == 1) {
            improv_v = v;
            improv_w = x_neighbor;
          }
          x_pos++;
        }
      }
    }
    // fprintf(stderr, "candidates t=%f\n", run_t.elapsed());
    // run_t.restart();

    if (solutions_found > 0) {
      perm.remove_from_solution(x, marked_degree, degree_limit, G);
      perm.add_to_solution(improv_v, marked_degree, degree_limit, G);
      perm.add_to_solution(improv_w, marked_degree, degree_limit, G);

      // Modification: Update
      if (marked_degree[improv_v] < 0 ||
          marked_degree[improv_v] > degree_limit) {
        if (marked_degree[improv_v] < 0) {
          forall_out_edges(G, edge, improv_v) {
            NodeID w = G.getEdgeTarget(edge);
            if (marked_degree[w] < 0 && candidates.contains(w)) candidates.remove(w);
          } endfor
        }
      } else {
        if (!candidates.contains(improv_v)) candidates.insert(improv_v);
      }

      // Modification: Update
      if (marked_degree[improv_w] < 0 ||
          marked_degree[improv_w] > degree_limit) {
        if (marked_degree[improv_w] < 0) {
          forall_out_edges(G, edge, improv_w) {
            NodeID w = G.getEdgeTarget(edge);
            if (marked_degree[w] < 0 && candidates.contains(w)) candidates.remove(w);
          } endfor
        }
      } else {
        if (!candidates.contains(improv_w)) candidates.insert(improv_w);
      }

      if (!perm.is_maximal())
        make_maximal_and_update(G, marked_degree, degree_limit);

      // Incremental
      update_candidates(x, G, marked_degree, degree_limit);
    }
    // fprintf(stderr, "insert t=%f\n", run_t.elapsed());
    // ASSERT_TRUE(perm.check_permutation());
  }
}

void local_search_online::make_maximal(graph_access &G,
                                std::vector<long long> &marked_degree,
                                long long degree_limit) {
  while (perm.get_free_size() > 0) {
    int random = random_functions::nextInt(0, perm.get_free_size() - 1);
    NodeID free_node = perm.get_free_node(random);
    perm.add_to_solution(free_node, marked_degree, degree_limit, G);
  }
}

void local_search_online::make_maximal_and_update(
    graph_access &G, std::vector<long long> &marked_degree,
    long long degree_limit) {
  while (perm.get_free_size() > 0) {
    int random = random_functions::nextInt(0, perm.get_free_size() - 1);
    NodeID free_node = perm.get_free_node(random);
    perm.add_to_solution(free_node, marked_degree, degree_limit, G);

    // Modification: Update
    if (marked_degree[free_node] < 0) {
      if (candidates.contains(free_node)) candidates.remove(free_node);
      forall_out_edges(G, edge, free_node) {
        NodeID w = G.getEdgeTarget(edge);
        if (marked_degree[w] < 0 && candidates.contains(w)) candidates.remove(w);
      } endfor
    }

    // Modification: Update
    if (marked_degree[free_node] >= 0 &&
        marked_degree[free_node] <= degree_limit)
      if (!candidates.contains(free_node)) {
        candidates.insert(free_node);
      }
  }
}

void local_search_online::build_onetight(NodeID node, graph_access &G,
                                  std::vector<long long> &marked_degree,
                                  long long degree_limit) {
  onetight_size = 0;
  forall_out_edges(G, edge, node) {
    NodeID target = G.getEdgeTarget(edge);

    // Modification: Skip
    if (marked_degree[target] < 0) continue;

    if (perm.get_tightness(target) == 1) onetight[onetight_size++] = target;
  } endfor
}

void local_search_online::print_onetight() {
  for (unsigned int i = 0; i < onetight_size; ++i) {
    printf("Node: %d\n", onetight[i]);
  }
}

void local_search_online::build_neighbors(NodeID node, graph_access &G,
                                   std::vector<long long> &marked_degree,
                                   long long degree_limit) {
  neighbors_size = 0;
  forall_out_edges(G, edge, node) {
    NodeID target = G.getEdgeTarget(edge);

    // Modification: Skip
    if (marked_degree[target] < 0) continue;

    neighbors[neighbors_size++] = target;
  } endfor
}

void local_search_online::build_candidates(graph_access &G,
                                    std::vector<long long> &marked_degree,
                                    long long degree_limit) {
  candidates.init(G.number_of_nodes());
  unsigned int solution_size = perm.get_solution_size();
  for (unsigned int i = 0; i < solution_size; ++i) {
    NodeID v = perm.get_solution_node(i);

    // Modification: Skip marked nodes and high degree nodes
    if (marked_degree[v] < 0 || marked_degree[v] > degree_limit) {
      if (candidates.contains(v)) candidates.remove(v);
    } else {
      if (!candidates.contains(v)) {
        candidates.insert(v);
      }
    }
  }
}

void local_search_online::update_candidates(NodeID node, graph_access &G,
                                     std::vector<long long> &marked_degree,
                                     long long degree_limit) {
  forall_out_edges(G, edge, node) {
    NodeID target = G.getEdgeTarget(edge);
    // Skip if neighbor is not 1-tight
    if (perm.get_tightness(target) != 1) continue;

    // Modification: Skip
    if (marked_degree[target] < 0 || marked_degree[target] > degree_limit)
      continue;

    forall_out_edges(G, target_edge, target) {
      NodeID candidate = G.getEdgeTarget(target_edge);

      // Modification: Skip
      if (marked_degree[candidate] < 0 ||
          marked_degree[candidate] > degree_limit)
        continue;

      if (perm.get_position(candidate) < perm.get_solution_size()) {
        if (!candidates.contains(candidate)) {
          candidates.insert(candidate);
        }
        // There can only be one valid candidate
        break;
      }
    } endfor
  } endfor
}

void local_search_online::print_permutation() {
  perm.print(0);
  perm.check_permutation();
}
