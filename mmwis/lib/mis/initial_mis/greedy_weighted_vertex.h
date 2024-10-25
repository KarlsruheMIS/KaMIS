/**
 * greedy_weighted_vertex.h
 * Purpose: Compute an initial solution (maximum independent set) by 
 *          building a minimal vertex cover.
 *
 *****************************************************************************/

#ifndef _GREEDY_WEIGHTED_VERTEX_H
#define _GREEDY_WEIGHTED_VERTEX_H


#include "initial_mis.h"
#include "mmwis_graph_access.h"

class greedy_weighted_vertex : public initial_mis {
    public:
        /**
         * Default Constructor.
         */
        greedy_weighted_vertex();

        /**
         * Default Destructor.
         */
        virtual ~greedy_weighted_vertex();

        /**
         * Generate an initial solution.
         * Build a minimal vertex cover of the given graph
         * and use the complement as a maximum independent set.
         *
         * @param seed Seed for the RNG.
         * @param G Graph representation.
         */
        void initial_partition( const unsigned int seed,
                                mmwis::graph_access & G );

    private:
        /**
         * Generate a permutation of the nodes.
         *
         * @param G Graph representation
         * @param permutation Permutation that is created
         */
        void generate_permutation ( mmwis::graph_access & G,
                                    NodePermutationMap & permutation);
};

#endif

