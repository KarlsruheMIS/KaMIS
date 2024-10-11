/**
 * greedy_mwis.h
 * Purpose: Compute an initial solution (maximum weight independent set)
 *          by using a greedy algorithm, that always picks the node 
 *          with the largest residual weight.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _GREEDY_MWIS_H_
#define _GREEDY_MWIS_H_

#include "initial_mis.h"
#include "data_structure/graph_access.h"

class greedy_mwis : public initial_mis {
    public:
        /**
         * Default Constructor.
         */
        greedy_mwis();

        /**
         * Default Destructor.
         */
        virtual ~greedy_mwis();

        /**
         * Generate an initial solution.
         * Use a greedy algorithm, that always picks the node
         * with the smalles residual degree.
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
        void generate_permutation( graph_access & G,
                                   NodePermutationMap & permutation);
};

#endif

