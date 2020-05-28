/**
 * online_ils.h
 * Purpose: Perform the iterated local search (ILS) as described by Andrade et
 *al.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _ONLINEILS_H_
#define _ONLINEILS_H_

#include <stack>
#include <tuple>
#include <vector>

#include "data_structure/graph_access.h"
#include "local_search_online.h"
#include "mis_config.h"
#include "data_structure/mis_permutation_online.h"
#include "timer.h"

struct is_solution {
  NodeID *solution;
  unsigned int solution_size;
};

class online_ils {
 public:
  /**
   * Default Constructor.
   */
  online_ils();

  /**
   * Default Destructor.
   */
  virtual ~online_ils();

  /**
   * Main algorithm of the ILS.
   * More detailed information can be found in the original paper
   * of Andrade et al.
   *
   * @param config Config used for the ILS.
   * @param G Graph representation.
   * @param iteration_limit Maximum number of iterations.
   */
  void perform_ils(MISConfig &config, graph_access &G,
                   unsigned int iteration_limit = 0);

  // Stack for folds.
  std::stack<std::tuple<NodeID, NodeID, NodeID>> &get_folds() {
    return perm->get_folds();
  };

  /**
   * Reset the ILS.
   * Clears the force-list and best solution.
   */
  void reset();

 private:
  // Array for storing the last time a node was forced in the solution.
  std::vector<NodeID> last_forced;
  // Main solution data structure.
  mis_permutation_online *perm;
  // List of candidates for the local search.
  candidate_list *cand;
  // List of nodes that were forced in the solution.
  candidate_list *force_list;
  // List of onetight nodes.
  candidate_list *one;
  // Best solution found so far
  NodeID best_solution;
  // Local search algorithm
  local_search_online local;
  // Timer for measuring the time taken for the ILS.
  timer t;
  // Degree marks
  std::vector<long long> marked_degree;
  // Degree limit
  long long degree_limit;

  // Fold stack
  std::stack<std::tuple<NodeID, NodeID, NodeID>> folds;

  // ILS config
  unsigned int plateau_down;
  unsigned int plateau_up;
  unsigned int plateau_best;
  unsigned int plateau_best_again;
  unsigned int pden_floor;
  unsigned int delta_penalty;
  bool limit_plateau;
  bool swap_on_failure;

  /**
   * Force the given node in the solution.
   *
   * @param config Config used for the ILS.
   * @param G Graph representation.
   * @param v Node to force into the solution.
   * @param force_list List storing the forced nodes.
   */
  void force(MISConfig &config, graph_access &G, NodeID v,
             candidate_list *force_list = NULL);

  /**
   * Undo all operations currently stored in the operation log.
   *
   * @param G Graph representation
   */
  void unwind(graph_access &G);
};

#endif
