/**
 * initial_mis.h
 * Purpose: Interface for different initial solution algorithms.
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

