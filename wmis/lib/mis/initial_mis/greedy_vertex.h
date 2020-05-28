/**
 * greedy_vertex.h
 * Purpose: Compute an initial solution (maximum independent set) by 
 *          building a minimal vertex cover.
 *
 *****************************************************************************/

#ifndef _GREEDY_VERTEX_H_
#define _GREEDY_VERTEX_H_

#include "initial_mis.h"
#include "data_structure/graph_access.h"

class greedy_vertex : public initial_mis {
    public:
        /**
         * Default Constructor.
         */
        greedy_vertex();

        /**
         * Default Destructor.
         */
        virtual ~greedy_vertex();

        /**
         * Generate an initial solution.
         * Build a minimal vertex cover of the given graph
         * and use the complement as a maximum independent set.
         *
         * @param seed Seed for the RNG.
         * @param G Graph representation.
         */
        void initial_partition( const unsigned int seed,
                                graph_access & G );

    private:
        /**
         * Generate a permutation of the nodes.
         *
         * @param G Graph representation
         * @param permutation Permutation that is created
         */
        void generate_permutation ( graph_access & G,
                                    NodePermutationMap & permutation);
};

#endif

