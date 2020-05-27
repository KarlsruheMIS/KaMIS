/**
 * multiway_combine.h
 * Purpose: Combine a list of individuals using k-separators/partitions.
 *
 *****************************************************************************/

#ifndef _MULTIWAY_COMBINE_H_
#define _MULTIWAY_COMBINE_H_

#include "combine.h"
#include "mis_config.h"
#include "ils/ils.h"
#include "ils/local_search.h"
#include "definitions.h"
#include "population_mis.h"
#include "separator_pool.h"
#include "graph_extractor.h"
#include "data_structure/graph_access.h"

class multiway_combine : public combine {
    public:
        /**
         * Default Constructor.
         */
        multiway_combine();

        /**
         * Default Destructor.
         */
        virtual ~multiway_combine();

        /**
         * Combine a given list of parents.
         * A new offsrping is gained by selecting a k-separator/k-partition and
         * then calculating the scores for each partition for each parent.
         * A combined solution is build using the highest scoring parents.
         *
         * @param config Config for the combine operator.
         * @param G Graph representation.
         * @param pool Separator pool containing the k-separators.
         * @param pop List of parents (whole population).
         * @param out Resulting offspring.
         */
        void combine(MISConfig & config, graph_access & G, separator_pool *pool, population_mis & pop, individuum_mis & out);

    private:
        // Local search algorithms
        ils iterate;
        local_search local;

        /**
         * Build a combined solution using a vertex cover approach.
         * Uncovered nodes get added using a greedy algorithm, that
         * always picks the node with the most uncovered edges.
         *
         * @param config Config for the combine operator.
         * @param G Graph representation.
         * @param best_parents List of highest scoring parents.
         * @param candidates Nodes added from boundary.
         */
        void build_combined_solution_vc(MISConfig & config, graph_access & G, std::vector<individuum_mis> & best_parents, std::vector<NodeID> & candidates);

        /**
         * Build a combined solution using a node separator approach.
         * The separator nodes get added using a greedy algorithm, that
         * always picks the node with the least residual degree.
         *
         * @param config Config for the combine operator.
         * @param G Graph representation.
         * @param best_parents List of highest scoring parents.
         * @param candidates Nodes added from separator.
         */
        void build_combined_solution_ns(MISConfig & config, graph_access & G, std::vector<individuum_mis> & best_parents, std::vector<NodeID> & candidates);

        /**
         * Creates the complement of the given individuum.
         *
         * @param G Graph representation.
         * @param in Individuum to build the complement for.
         * @param out Individuum containing the complement.
         */
        void create_complement(graph_access & G, individuum_mis & in, individuum_mis & out);

        /**
         * Greedily add uncovered nodes to the solutions.
         * Uses a greedy algorithm, which always chooses the 
         * node with the maximal residual degree.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param candidates Nodes added from boundary.
         */
        void build_vertex_cover_candidates(MISConfig & config, graph_access & G, std::vector<NodeID> & candidates);

        /**
         * Add the nodes of the separator to the solutions.
         * Uses a greedy algorithm, which always chooses the 
         * node with the least residual degree.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param candidates Nodes added from separator.
         */
        void build_separator_candidates(MISConfig & config, graph_access & G, std::vector<NodeID> & candidates);

        /**
         * Get a random k-partition from the pool build with the KaHIP-interface.
         * Apply that partition to the graph.
         *
         * @param config Config used for the combinator.
         * @param G Graph representation.
         * @param pool Pool of k-partitions.
         * @param part The partition that was applied.
         */
        void apply_k_partition_kahip(MISConfig & config, graph_access & G, separator_pool *pool, partition & part);

        /**
         * Get a random k-separator from the pool build with the KaHIP-interface.
         * Apply that separator to the graph.
         *
         * @param config Config used for the combinator.
         * @param G Graph representation.
         * @param pool Pool of k-separators.
         * @param sep The separator that was applied.
         */
        void apply_k_separator_kahip(MISConfig & config, graph_access & G, separator_pool *pool, separator & sep);
};

#endif
