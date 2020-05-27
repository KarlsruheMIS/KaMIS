/**
 * separator_combine.h
 * Purpose: Combine two individuals based on node separators.
 *
 *****************************************************************************/

#ifndef _SEPARATOR_COMBINE_H_
#define _SEPARATOR_COMBINE_H_

#include "combine.h"
#include "mis_config.h"
#include "ils/ils.h"
#include "ils/local_search.h"
#include "definitions.h"
#include "population_mis.h"
#include "separator_pool.h"
#include "graph_extractor.h"
#include "data_structure/graph_access.h"

class separator_combine : public combine {
    public:
        /**
         * Default Constructor.
         */
        separator_combine();

        /**
         * Default Destructor.
         */
        virtual ~separator_combine();

        /**
         * Combine two individuals.
         * Combination works by first building a partition of the graph. 
         * After that a node separator S (G = S + A + B) is created from this partition.
         * Then the nodes from the first individuum that are in A and
         * the nodes from the second individuum that are in B are taken and combined.
         * Finally local search is applied to this combination using 
         * the node separator as a candidate list.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param pool Pool of node separators.
         * @param first First individuum.
         * @param second Second individuum.
         * @param out_first Individuum resulting from the combination
         * @param out_second Individuum resulting from the combination
         */
        void combine(MISConfig & config, graph_access & G, separator_pool *pool, individuum_mis & first, individuum_mis & second, individuum_mis & out_first, individuum_mis & out_second);

    private:
        // Local search algorithms
        ils iterate;
        local_search local;

        /**
         * Add the nodes of the separator to the solutions.
         *
         * @param config Config for the combinator.
         * @param G Graph representation.
         * @param candidates Nodes that got added from the separator.
         */
        void build_separator_candidates(MISConfig & config, graph_access & G, std::vector<NodeID> & candidates);

        /**
         * Get a random separator from the pool build with the KaHIP-interface.
         *
         * @param config Config used for the combinator.
         * @param G Graph representation.
         * @param pool Pool of separators.
         */
        void apply_separator_kahip(MISConfig & config, graph_access & G, separator_pool *pool);

};

#endif
