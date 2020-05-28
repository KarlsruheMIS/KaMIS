/**
 * initial_mis.h
 * Purpose: Interface for different initial solution algorithms.
 *
 *****************************************************************************/

#ifndef _INITIAL_MIS_H_
#define _INITIAL_MIS_H_

#include "data_structure/graph_access.h"

class initial_mis {
    public:
        /**
         * Default Constructor.
         */
        initial_mis();

        /**
         * Default Destructor.
         */
        virtual ~initial_mis();

        /**
         * Interface for performing the initial partitioning.
         *
         * @param seed Seed for the RNG.
         * @param G Graph representation.
         */
        virtual void initial_partition( const unsigned int seed,
                                        graph_access & G ) = 0;
};

#endif

