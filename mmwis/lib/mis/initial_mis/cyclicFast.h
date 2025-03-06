/**
 * greedy_mwis.h
 * Purpose: Compute an initial solution (maximum weight independent set)
 *          by using a greedy algorithm, that always picks the node 
 *          with the largest residual weight.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _CYCLIC_FAST_H_
#define _CYCLIC_FAST_H_

#include "initial_mis.h"
#include "mmwis_config.h"
#include "graph_access.h"

namespace mmwis {
class cyclicFast: public initial_mis {
    public:
        /**
         * Default Constructor.
         */
        cyclicFast();

        /**
         * Default Destructor.
         */
        virtual ~cyclicFast();

        /**
         * Generate an initial solution.
         * Use a greedy algorithm, that always picks the node
         * with the smalles residual degree.
         *
         * @param seed Seed for the RNG.
         * @param G Graph representation.
         */
        bool initial_partition_struction(MISConfig & config, 
                                graph_access & G,
                                double remaining_time); 

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

}
#endif

