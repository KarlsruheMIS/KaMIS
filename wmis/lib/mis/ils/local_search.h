/*
 * local_search.h
 * Purpose: Apply the local search algorithm to a maximum independent set.
 *
 * The original code from Andrade et. al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _DIRECT_LOCAL_H_
#define _DIRECT_LOCAL_H_

#include <vector>

#include "mis_config.h"
#include "data_structure/graph_access.h"
#include "data_structure/mis_permutation.h"
#include "data_structure/priority_queues/maxNodeHeap.h"

class local_search {
    friend class ils;
    public:
        /**
         * Constructor.
         *
         * @param sort_freenodes Wheter the free nodes should be sorted
         */
        local_search(const MISConfig& config);

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
		// insertion node of a (1,2)-swap
		using Swap_1_2 = std::pair<NodeID, NodeID>;
		static constexpr NodeID INVALID_NODE = std::numeric_limits<NodeID>::max();

        // List of 1-tight neighbors of a node.
        std::vector<NodeID> onetight;
		// (1,2)-swaps for each of the candidates
		std::unordered_map <NodeID, Swap_1_2> candidate_swaps;
        // Adjacency list of a node.
        std::vector<NodeID> neighbors;
        // Candidate list as a pq with gains as key
		// TODO: add template parameter to maxNodeHeap
		maxNodeHeap candidates;
        // Permutation data-structure.
        mis_permutation perm;
        // Amount of 1-tight neighbors.
        unsigned int onetight_size;
        // Amount of neighbors.
        unsigned int neighbors_size;
		
		bool sort_freenodes;

        /**
         * Build a candidate list based on the solution nodes of the graph.
         *
         * @param G Graph representation.
         */
        void build_candidates(graph_access & G);

		/**
		* Add a node to the candidate list.
		*
		* @param node The new candidate node.
        * @param G Graph representation.
		*/
		void add_candidate(NodeID node, graph_access & G);

		/**
		* Add a (1,2)-swap to the candidate_swap list.
		*
		* @param node The candidate node.
		* @param swap The new swap.
		* @param G Graph representation.
		* @return The gain of the swap nodes.
		*/
		int insert_swap(NodeID node, Swap_1_2 swap, graph_access & G);

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
		* Updates a candidate node because of a newly available
		* 1-tight neighboring node.
		*
		* @param node Candidate node which is to be updated.
		* @param one_tight 1-tight node which caused the update.
		* @param G Graph representation.
		*/
		void update_candidate(NodeID node, NodeID one_tight, graph_access & G);

		/**
		* Updates the (1,2)-swaps for candidates because of conflicts
		* while inserting the new solution node.
		*
		* @param node New solution node.
		* @param G Graph representation.
		*/
		void update_swaps(NodeID node, graph_access & G);


		/**
		* Finds the best (1,2)-swap for a candidate node.
		*
		* @param node The candidate node.
		* @param G Graph representation.
		* @return The best (1,2)-swap.
		*/
		Swap_1_2 find_best_swap(NodeID node, graph_access & G);

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


        /**
        * Sorts the given range by node weights in descending order using counting sort.
        *
        * @param begin Beginning of the range in vector
        * @param end End of the range in vector
        */
        void sort_by_weight(graph_access & G, std::vector<NodeID>::iterator begin, std::vector<NodeID>::iterator end);
};

#endif

