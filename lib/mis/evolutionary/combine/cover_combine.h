/**
 * cover_combine.h
 * Purpose: Combine two individuals with a vertex cover based mechanism.
 *
 *****************************************************************************/

#ifndef _COVER_COMBINE_H_
#define _COVER_COMBINE_H_

#include "combine.h"
#include "mis_config.h"
#include "ils/ils.h"
#include "ils/local_search.h"
#include "definitions.h"
#include "population_mis.h"
#include "separator_pool.h"
#include "data_structure/graph_access.h"

class cover_combine : public combine {
    public:
        /**
         * Default Constructor.
         */
        cover_combine();

        /**
         * Default Destructor.
         */
        virtual ~cover_combine();

        /**
         * Combine two individuals.
         * Combination works by first building the complements of the 
         * input individuals (which are vertex covers).
         * Then a partition for the graph is build and for each partition 
         * the individuum with lesser solution nodes is chosen.
         * The vertex covers are then combined based on this process.
         * To ensure the final vertex cover is valid the boundary nodes
         * get added greedily (or using the Hopcroft-Karp algorithm).
         * The complement of the final vertex cover represents a new
         * maximum independent set.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param pool Pool of partitions.
         * @param first First individuum.
         * @param second Second individuum.
         * @param out_first Individuum resulting from the combination
         * @param out_second Second individuum resulting from the combination
         */
        void combine(MISConfig & config, graph_access & G, separator_pool *pool, individuum_mis & first, individuum_mis & second, individuum_mis & out_first, individuum_mis & out_second);

    private:
        // Local search algorithms
        ils iterate;
        local_search local;
    
        /**
         * Creates the complement of the given individuum.
         *
         * @param G Graph representation.
         * @param in Individuum to build the complement for.
         * @param out Individuum containing the complement.
         */
        void create_complement(graph_access & G, individuum_mis & in, individuum_mis & out);

        /**
         * Fix the partition cut boundaries.
         * First extracts all boundaries and then creates
         * corresponding vertex cover using the Hopcroft-Karp algorithm.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param cover Vertex cover representation.
         */
        void build_hopcroft_cover(MISConfig & config, graph_access & G, graph_access & cover);

        /**
         * Turn a given vertex cover to a maximum independent set.
         * Also improve the individual by using a local search algorithm.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param cover Individual storing a vertex cover.
         * @param ind Individual containing a maximum indepent set.
         */
        void vertex_cover_to_mis(MISConfig & config, graph_access & G, individuum_mis & cover, individuum_mis & ind);

        /**
         * Extract all possible boundaries found in the given graph.
         * A boundary is represented by two arrays. Each array indicates if a node is
         * on the left or right side of the boundary.
         *
         * @param G Graph representation.
         * @param boundaries Resulting array of boundaries.
         */
        void extract_boundaries(graph_access & G, std::vector<std::pair<std::vector<NodeID>, std::vector<NodeID>>> & boundaries);

        /**
         * Greedily add uncovered nodes to the solutions.
         * Uses a greedy algorithm, which always chooses the 
         * node with the maximal residual degree.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param partition_index Use first (0) or second (1) partition.
         */
        void build_vertex_cover_candidates(MISConfig & config, graph_access & G, unsigned int partition_index);

        /**
         * Checks which individual has a smaller vertex cover for the given graph.
         *
         * @param G Graph representation.
         * @param mapping Node mapping for the graph.
         * @param first First individuum.
         * @param second Second individuum.
         * @return 1 if the first cover is smaller, 2 if the second one.
         */
        unsigned int get_smaller_individuum(graph_access & G, std::vector<NodeID> mapping, individuum_mis & first, individuum_mis & second);

        /**
         * Get a random partition from the pool build with the KaHIP-interface.
         *
         * @param config Config used for the combinator.
         * @param G Graph representation.
         * @param pool Pool of partitions.
         */
        void apply_partition_kahip(MISConfig & config, graph_access & G, separator_pool *pool);
};

#endif
