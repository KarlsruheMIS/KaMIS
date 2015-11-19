/**
 * greedy_vertex.h
 * Purpose: Compute an initial solution (maximum independent set) by 
 *          building a minimal vertex cover.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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

