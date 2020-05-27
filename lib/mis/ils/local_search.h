/**
 * local_search.h
 * Purpose: Apply the local search algorithm to a maximum independent set.
 *
 * The original code from Andrade et. al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _DIRECT_LOCAL_H_
#define _DIRECT_LOCAL_H_

#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/mis_permutation.h"

class local_search {
    friend class ils;
    public:
        /**
         * Constructor.
         */
        local_search();

        /**
         * Destructor.
         */
        virtual ~local_search();

        /**
         * Preprocess the given graph.
         * This involves creating the permutation and
         * building a candidate list.
         * The graph has to contain an mis.
         *
         * @param G Graph representation.
         */
        void preprocess_graph(graph_access & G);

        /**
         * Preprocess the given graph with a certain list of candidates.
         * The graph has to contain an mis.
         *
         * @param G Graph representation.
         * @param cand Candidate list.
         * @param cand_size Size of the candidate list.
         */
        void preprocess_graph_with_candidates(graph_access & G, std::vector<NodeID> cand, unsigned int cand_size);

        /**
         * Forcefully insert a given amount of random nodes into the solution.
         *
         * @param G Graph representation.
         * @param k Number of nodes.
         */
        void force(graph_access & G, unsigned int k);

        /**
         * Forcefully insert the given node.
         *
         * @param G Graph representation.
         * @param node Node to insert.
         */
        void force_node(graph_access & G, NodeID node);

        /**
         * Run the two-improvement algorithm proposed by Andrade et al.
         *
         * @param G Graph representation.
         * @param forced Whether or not the candidates contain a forced node.
         * @param forced_node The forced node.
         */
        void simple_improvement(graph_access & G, bool forced = false, NodeID forced_node = 0);

        /**
         * Run the direct version of the two-improvement algorithm proposed by Andrade et al.
         *
         * @param G Graph representation.
         * @param forced Whether or not the candidates contain a forced node.
         * @param forced_node The forced node.
         */
        void direct_improvement(graph_access & G, bool forced = false, NodeID forced_node = 0);

        /**
         * Make the solution maximal by inserting free vertices.
         *
         * @param G Graph representation.
         */
        void make_maximal(graph_access & G);

        /**
         * Add a given amount of candidates to the list.
         *
         * @param G Graph representation
         * @param cand List of candidates to be added.
         * @param num_cand Number of candidates to be added.
         */
        void insert_candidates(graph_access & G, std::vector<NodeID> cand, unsigned int num_cand);

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
        mis_permutation perm;
        // Amount of 1-tight neighbors.
        unsigned int onetight_size;
        // Amount of neighbors.
        unsigned int neighbors_size;

        /**
         * Build a candidate list based on the solution nodes of the graph.
         *
         * @param G Graph representation.
         */
        void build_candidates(graph_access & G);

        /**
         * Updates the candidate list after removal of a node based
         * on the incremental two-improvement algorithm proposed by
         * Andrade et al.
         *
         * @param node Node that was removed from the candidates.
         * @param G Graph representation.
         */
        void update_candidates(NodeID node, graph_access & G);

        /**
         * Build a list of neighbors for a given node
         *
         * @param node Node for which the list of neighbors should be created.
         * @param G Graph representation.
         */
        void build_neighbors(NodeID node, graph_access & G);

        /**
         * Build a list of 1-tight neighbors for a given node
         *
         * @param node Node for which the list of neighbors should be created.
         * @param G Graph representation.
         */
        void build_onetight(NodeID node, graph_access & G);

        /**
         * Print the current list of 1-tight nodes.
         */
        void print_onetight();
};

#endif

