/**
 * local_search_online.h
 * Purpose: Apply the local search algorithm to a maximum independent set.
 *
 * The original code from Andrade et. al. was kindly provided by Renato Werneck.
 *
  *****************************************************************************/

#ifndef _DIRECT_LOCAL_ONLINE_H_
#define _DIRECT_LOCAL_ONLINE_H_

#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/mis_permutation_online.h"

class local_search_online {
  friend class online_ils;

 public:
  /**
   * Constructor.
   */
  local_search_online();

  /**
   * Destructor.
   */
  virtual ~local_search_online();

  /**
   * Initialize the underlying permutation
   * data structure.
   *
   * @param G Graph representation.
   */
  void init(graph_access &G, std::vector<long long> &marked_degree,
            long long degree_limit) {
    perm.construct(G, marked_degree, degree_limit);
    onetight.clear();
    neighbors.clear();
    onetight.resize(G.getMaxDegree());
    neighbors.resize(G.getMaxDegree());
  }

  /**
   * Forcefully insert a given amount of random nodes into the solution.
   *
   * @param G Graph representation.
   * @param k Number of nodes.
   */
  void force(graph_access &G, unsigned int k);

  /**
   * Forcefully insert the given node.
   *
   * @param G Graph representation.
   * @param node Node to insert.
   */
  void force_node(graph_access &G, NodeID node);

  /**
   * Run the direct version of the two-improvement algorithm proposed by Andrade
   * et al.
   *
   * @param G Graph representation.
   * @param marked_degree Degree marks for skipping.
   * @param degree_limit Degree limit for skipping.
   * @param forced Whether or not the candidates contain a forced node.
   * @param forced_node The forced node.
   */
  void direct_improvement(graph_access &G,
                          std::vector<long long> &marked_degree,
                          long long degree_limit, bool forced = false,
                          long long forced_node = 0);

  /**
   * Make the solution maximal by inserting free vertices.
   *
   * @param G Graph representation.
   * @param marked_degree Degree marks for skipping.
   * @param degree_limit Degree limit for skipping.
   */
  void make_maximal(graph_access &G, std::vector<long long> &marked_degree,
                    long long degree_limit);

  /**
   * Make the solution maximal by inserting free vertices.
   * Additionaly update candidate list.
   *
   * @param G Graph representation.
   * @param marked_degree Degree marks for skipping.
   * @param degree_limit Degree limit for skipping.
   */
  void make_maximal_and_update(graph_access &G,
                               std::vector<long long> &marked_degree,
                               long long degree_limit);

  /**
   * Print the current state of the permutation without details
   * and check if it's valid.
   */
  void print_permutation();

 private:
  // List of 1-tight neighbors of a node.
  std::vector<NodeID> onetight;
  // Adjacency list of a node.
  std::vector<NodeID> neighbors;
  // Candidate list
  candidate_list candidates;
  // Permutation data-structure.
  mis_permutation_online perm;
  // Amount of 1-tight neighbors.
  unsigned int onetight_size;
  // Amount of neighbors.
  unsigned int neighbors_size;

  unsigned int long_runs;

  /**
   * Build a candidate list based on the solution nodes of the graph.
   *
   * @param G Graph representation.
   * @param marked_degree Degree marks for skipping.
   * @param degree_limit Degree limit for skipping.
   */
  void build_candidates(graph_access &G, std::vector<long long> &marked_degree,
                        long long degree_limit);

  /**
   * Updates the candidate list after removal of a node based
   * on the incremental two-improvement algorithm proposed by
   * Andrade et al.
   *
   * @param node Node that was removed from the candidates.
   * @param G Graph representation.
   */
  void update_candidates(NodeID node, graph_access &G,
                         std::vector<long long> &marked_degree,
                         long long degree_limit);

  /**
   * Build a list of neighbors for a given node
   *
   * @param node Node for which the list of neighbors should be created.
   * @param G Graph representation.
   */
  void build_neighbors(NodeID node, graph_access &G,
                       std::vector<long long> &marked_degree,
                       long long degree_limit);

  /**
   * Build a list of 1-tight neighbors for a given node
   *
   * @param node Node for which the list of neighbors should be created.
   * @param G Graph representation.
   */
  void build_onetight(NodeID node, graph_access &G,
                      std::vector<long long> &marked_degree,
                      long long degree_limit);

  /**
   * Print the current list of 1-tight nodes.
   */
  void print_onetight();
};

#endif
