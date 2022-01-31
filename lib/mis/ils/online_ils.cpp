/**
 * online_ils.cpp
 * Purpose: Perform the iterated local search (ILS) as described by Andrade et
 *al.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
  *****************************************************************************/

#include "online_ils.h"

#include <algorithm>
#include <fstream>

#include "data_structure/operation_log.h"
#include "greedy_mis_online.h"
#include "mis_log.h"
#include "random_functions.h"

online_ils::online_ils() {
  // Set config parameters
  plateau_down = 1;
  plateau_up = 0;
  plateau_best = 0;
  plateau_best_again = 0;
  pden_floor = 1;
  delta_penalty = 0;
  limit_plateau = false;
  swap_on_failure = true;

  force_list = NULL;
  best_solution = 0;
}

online_ils::~online_ils() { reset(); }

void online_ils::perform_ils(MISConfig &config, graph_access &G,
                      unsigned int iteration_limit) {
  reset();
  // Init operation log
  operation_log::instance()->init(G.number_of_nodes());
  operation_log::instance()->reset();
  operation_log::instance()->deactivate();

  // Init last forced
  last_forced.clear();
  last_forced.resize(G.number_of_nodes(), 0);

  // Init RNG
  srand(config.seed);
  random_functions::setSeed(config.seed);

  // Initialize "friends"
  cand = &local.candidates;
  perm = &local.perm;
  one = &local.perm.onetight_all;

  // Start timer
  t.restart();
  mis_log::instance()->restart_online_timer();

  // Modification: Init marks
  degree_limit = G.getMaxDegree();
  marked_degree.resize(G.number_of_nodes(), -3);
  local.init(G, marked_degree, degree_limit);

  // Modification: Insert eagerly
  forall_nodes(G, node) { 
    marked_degree[node] = G.getNodeDegree(node); 
  } endfor

  // Modification: Calculate 1% percentile degree limit
  // Use adaptive approach suggested by Prof. Sanders
  // Immediately add vertices with degree < 2 to the solution
  // ASSERT_TRUE(perm->check_consistency(G));
  double percentage = 0.01;
  long long cut_size = ceil((double)G.number_of_nodes() * percentage);
  std::vector<long long> percentile(marked_degree);
  std::nth_element(percentile.begin(), percentile.begin() + cut_size,
                   percentile.end(), std::greater<long long>());
  degree_limit = percentile[cut_size];

  forall_nodes(G, node) {
    if (marked_degree[node] < 0 || marked_degree[node] > 2) continue;
    perm->add_and_reduce(node, marked_degree, degree_limit, G);
  } endfor

  // Init local search data structures
  if (config.start_greedy_adaptive) {
          greedy_mis_online init_mis;
          init_mis.initial_partition(config.seed, G, marked_degree, degree_limit);
          forall_nodes(G, node) {
                  if (marked_degree[node] >= 0 && G.getPartitionIndex(node))
                          perm->add_to_solution(node, marked_degree, degree_limit, G);
          } endfor
  } else {
          local.make_maximal(G, marked_degree, degree_limit);
  }

  local.direct_improvement(G, marked_degree, degree_limit);

  best_solution = perm->get_solution_size() + perm->get_folded_vertices();
  mis_log::instance()->set_best_size_online(config, best_solution);

  unsigned int denominator = 4;
  unsigned int plateau = plateau_best * (perm->get_solution_size() + perm->get_folded_vertices());
  force_list = new candidate_list();
  force_list->init(G.number_of_nodes());

  timer run_t;
  run_t.restart();
  unsigned int iterations = 0;
  while (1) {
    // fprintf(stderr, "total t=%f\n", run_t.elapsed());
    // std::cin.get();
    run_t.restart();
    iterations++;
    // Stop if the time limit was passed
    if (t.elapsed() > config.time_limit) break;
    if (iterations % 10000 == 0) {
      mis_log::instance()->set_best_size_online(config, best_solution);
    }

    plateau--;
    if (plateau < 0) plateau = 0;

    // Activate the operation log
    operation_log::instance()->reset();
    operation_log::instance()->activate();

    int non_solution = G.number_of_nodes() - (perm->get_solution_size() + perm->get_folded_vertices());
    if (non_solution <= 0) break;

    long long should_skip = 0;
    forall_nodes(G, node) {
      if (marked_degree[node] < 0 || marked_degree[node] > degree_limit) {
        should_skip++;
      }
    } endfor
    if (non_solution - should_skip <= 0) break;

    denominator = 2 * (perm->get_solution_size() + perm->get_folded_vertices());
    unsigned int solution_before = perm->get_solution_size() + perm->get_folded_vertices();
    unsigned int forced = 1;

    // Check if more than one vertex should be forced and if so
    // determine how many
    unsigned int outer_denominator = denominator;
    if (random_functions::nextInt(1, outer_denominator) == 1) {
      unsigned int inner_denominator = 2;
      forced++;
      while (random_functions::nextInt(1, inner_denominator) == 1) forced++;
    }

    // Look at a certain number of candidates for the forceful insertion
    // and pick the one that hasn't been in the solution the longest time.
    long long v;
    long long best_v = -1;
    non_solution = G.number_of_nodes() - (perm->get_solution_size() + perm->get_folded_vertices());
    for (int i = 0; i < config.force_cand; i++) {
      unsigned int position = random_functions::nextInt(0, non_solution - 1);
      v = perm->get_non_solution_node(position);

      // Modification: Skip nodes with degree higher than the limit
      // Modification: Don't touch neighbors of solution nodes
      if (marked_degree[v] < 0 || marked_degree[v] > degree_limit) {
        i--;
        continue;
      }

      if (best_v == -1)
        best_v = v;
      else if (last_forced[v] < last_forced[best_v])
        best_v = v;
    }
    v = best_v;

    // Actual insertion
    if (forced == 1)
      force(config, G, v, NULL);
    else {
      force_list->reset();
      for (int i = 0; i < forced; i++) {
        // No vertex left then stop
        if (v == -1) {
          forced = i;
          break;
        }
        // Force the node
        force(config, G, v, force_list);
        v = -1;
        // Modification: Remove high degree and marked vertices
        unsigned int num_cand = force_list->get_size();
        for (int j = 0; j < num_cand; j++) {
          NodeID w = force_list->pick(j);
          if (marked_degree[w] < 0 || marked_degree[w] > degree_limit) {
            force_list->remove(w);
            continue;
          }
        }
        // Pick a node thats close to the ones we already removed
        num_cand = force_list->get_size();
        if (num_cand == 0) break;
        force_list->random_permute();
        unsigned int valid_count = 0;
        for (int j = 0; j < num_cand; j++) {
          NodeID w = force_list->pick(j);
          forall_out_edges(G, edge, w) {
            NodeID u = G.getEdgeTarget(edge);

            // Modification: Skip nodes with degree higher than the limit
            // Modification: Don't touch neighbors of marked nodes
            if (marked_degree[u] < 0 || marked_degree[u] > degree_limit)
              continue;

            if (perm->is_solution_node(u)) continue;
            if (force_list->contains(u)) continue;

            valid_count++;
            if (random_functions::nextInt(1, valid_count) == 1) v = u;
          } endfor
        }
      }
    }

#ifndef NDEBUG
    if (forced == 0) printf("Should have at least inserted one vertex\n");
#endif
    if (forced != 1) v = -1;

    // See if the new solution can be further improved
    local.make_maximal_and_update(G, marked_degree, degree_limit);

    // Modification: Only process unreduced nodes
    if (marked_degree[v] >= 0)
      local.direct_improvement(G, marked_degree, degree_limit, true, v);

    // Check the difference
    unsigned int solution_after = perm->get_solution_size() + perm->get_folded_vertices();
    int delta = solution_before - solution_after;
    int delta_best = best_solution - solution_after;

    // New best?
    if (solution_after > best_solution) {
      best_solution = perm->get_solution_size() + perm->get_folded_vertices();
      plateau = plateau_best * (perm->get_solution_size() + perm->get_folded_vertices());
      // mis_log::instance()->set_best_size(config, best_solution);
    } else if (solution_after > solution_before) {
      if (solution_after == best_solution) {
        plateau = plateau_best_again * (perm->get_solution_size() + perm->get_folded_vertices());
      } else {
        plateau = plateau_up * (perm->get_solution_size() + perm->get_folded_vertices());
      }
    }

    // If the new solution solution isn't better
    // it's possible to stay put
    if (delta > 0 || (delta == 0 && limit_plateau)) {
      int c = delta_penalty;
      int pden = pden_floor + (delta_best + c) * (delta + c);

      // Keep the solution with a small probability
      if (plateau > 0 || random_functions::nextInt(1, pden) != 1) {
        // Undo operations that led to the new solution
        unwind(G);
        if (swap_on_failure) {
          // Perform a random 1-swap
          int x = -1;
          if (!one->is_empty()) {
            for (int i = 0; i < config.force_cand; i++) {
              int y = one->pick_random();

              // Modification: Skip nodes with degree higher than the limit
              // Modification: Don't touch neighbors of marked nodes
              if (marked_degree[y] < 0 || marked_degree[y] > degree_limit) {
                i--;
                continue;
              }

              if (x == -1 || last_forced[y] < last_forced[x]) x = y;
            }
            force(config, G, x, NULL);
            // Try to improve the solution
            // fprintf(stderr, "delta v=%d s=%d m=%d\n", x,
            // perm->is_solution_node(x), marked_degree[x]); Modification: Only
            // process unreduced nodes
            if (marked_degree[x] >= 0)
              local.direct_improvement(G, marked_degree, degree_limit, true, x);
            // Improvement found?
            if (perm->get_solution_size() + perm->get_folded_vertices() > best_solution) {
              best_solution = perm->get_solution_size() + perm->get_folded_vertices();
              plateau = plateau_best * (perm->get_solution_size() + perm->get_folded_vertices());
            }
          }
        }
      } else {
        plateau = plateau_down * (perm->get_solution_size() + perm->get_folded_vertices());
      }
    }

    // Update the "Forced list"
    for (unsigned int pos = operation_log::instance()->get_size(); pos > 0; pos--) {
      int v = operation_log::instance()->peek(pos - 1);
      if (v < 0) v = -v;
      if (!perm->is_solution_node(v)) last_forced[v] = iterations;
    }

    // Candidates should be empty anyway
    if (cand && !cand->is_empty()) {
      // Modification: Candidates can be empty if node was reduced
      // printf("Candidates should be empty\n");
      cand->reset();
    }
    // ASSERT_TRUE(perm->check_consistency(G));
  }

  forall_nodes(G, node) {
    G.setPartitionIndex(node, perm->is_solution_node(node));
  } endfor

  //std::cout << "folded " << perm->get_folded_vertices() << std::endl;
}

void online_ils::force(MISConfig &config, graph_access &G, NodeID v,
                candidate_list *force_list) {
  // Modification: Skip
  if (marked_degree[v] < 0 || marked_degree[v] > degree_limit) return;

#ifndef NDEBUG
  if (perm->is_solution_node(v)) printf("Vertex already in the solution.\n");
#endif
  if (force_list) {
    if (!force_list->contains(v)) force_list->insert(v);
#ifndef NDEBUG
    else
      printf("Vertex already in the force list.\n");
#endif
  }

  forall_out_edges(G, edge, v) {
    NodeID w = G.getEdgeTarget(edge);

    // Modification: Skip
    if (marked_degree[w] < 0) continue;

    if (force_list) {
      if (!force_list->contains(w)) force_list->insert(w);
    }
    if (perm->is_solution_node(w)) {
      perm->remove_from_solution(w, marked_degree, degree_limit, G);
      if (cand) local.update_candidates(w, G, marked_degree, degree_limit);
    }
  } endfor

  perm->add_to_solution(v, marked_degree, degree_limit, G);
  if (marked_degree[v] >= 0 && marked_degree[v] <= degree_limit){
    if (!cand->contains(v)) cand->insert(v);
  } else {
    if (cand->contains(v)) cand->remove(v);
  }
}

void online_ils::unwind(graph_access &G) {
  operation_log::instance()->deactivate();
  while (!operation_log::instance()->is_empty()) {
    int v = operation_log::instance()->unwind();
    // Element was removed
    // if (marked_degree[std::abs(v)] < 0) continue;
    if (v < 0) perm->add_to_solution(-v, marked_degree, degree_limit, G);
    // Element was added
    else
      perm->remove_from_solution(v, marked_degree, degree_limit, G);
  }
  operation_log::instance()->activate();
}

void online_ils::reset() {
  if (force_list != NULL) {
    delete force_list;
    force_list = NULL;
  }
}
