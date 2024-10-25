/**
 * greedy_mis.h
 * Purpose: Compute an initial solution (maximum independent set)
 *          by using a greedy algorithm, that always picks the node 
 *          with the smallest residual degree.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _GREEDY_MIS_H_
#define _GREEDY_MIS_H_

#include "initial_mis.h"
#include "mmwis_graph_access.h"

class greedy_mis : public initial_mis {
    public:
        /**
         * Default Constructor.
         */
        greedy_mis();

        /**
         * Default Destructor.
         */
        virtual ~greedy_mis();

        /**
         * Generate an initial solution.
         * Use a greedy algorithm, that always picks the node
         * with the smalles residual degree.
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
        void generate_permutation( mmwis::graph_access & G,
                                   NodePermutationMap & permutation);
};

#endif

