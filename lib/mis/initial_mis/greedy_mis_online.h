/**
 * greedy_mis.h
 * Purpose: Compute an initial solution (maximum independent set)
 *          by using a greedy algorithm, that always picks the node
 *          with the smallest residual degree.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _GREEDY_MIS_ONLINE_H_
#define _GREEDY_MIS_ONLINE_H_

#include <vector>

#include "data_structure/graph_access.h"
#include "initial_mis.h"

class greedy_mis_online : public initial_mis {
 public:
  /**
   * Default Constructor.
   */
  greedy_mis_online();

  /**
   * Default Destructor.
   */
  virtual ~greedy_mis_online();

  /**
   * Generate an initial solution.
   * Use a greedy algorithm, that always picks the node
   * with the smalles residual degree.
   *
   * @param seed Seed for the RNG.
   * @param G Graph representation.
   */
  void initial_partition(const unsigned int seed, graph_access &G);
  void initial_partition(const unsigned int seed, graph_access &G,
                         std::vector<long long> &marked_degree,
                         long long degree_limit);

 private:
  /**
   * Generate a permutation of the nodes.
   *
   * @param G Graph representation
   * @param permutation Permutation that is created
   */
  void generate_permutation(graph_access &G, NodePermutationMap &permutation);
};

#endif
