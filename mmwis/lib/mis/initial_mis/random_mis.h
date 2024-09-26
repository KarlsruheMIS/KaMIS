/**
 * random_mis.h
 * Purpose: Compute an initial solution (maximum independent set)
 *          by using a simple random selection of vertices. 
 *
 *****************************************************************************/

#ifndef _RANDOM_MIS_H_
#define _RANDOM_MIS_H_

#include "initial_mis.h"
#include "data_structure/graph_access.h"

class random_mis : public initial_mis {
    public:
        /**
         * Default Constructor.
         */
        random_mis();

        /**
         * Default Destructor.
         */
        virtual ~random_mis();

        /**
         * Generate an initial solution.
         * Use random selection of nodes.
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

