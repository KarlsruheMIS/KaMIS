/**
 * configuration.h
 * Purpose: Contains preset configurations for the evolutionary algorithms.
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

#ifndef _CONFIGURATION_MIS_H_
#define _CONFIGURATION_MIS_H_

#include "definitions.h"
#include "mis_config.h"
#include "interface/kaHIP_interface.h"

class configuration_mis {
    public:
        /**
         * Default Constructor.
         */
        configuration_mis() {} ;

        /**
         * Default Destructor.
         */
        virtual ~configuration_mis() {};

        /**
         * Set the standard configuration.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void standard( MISConfig & config );

        /**
         * Set the configuration for social network graphs.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void social( MISConfig & config );

        /**
         * Set the configuration for the experimental evaluation.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void full_standard( MISConfig & config );

        /**
         * Set the configuration for the experimental evaluation 
         * for social network graphs.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void full_social( MISConfig & config );
};

inline void configuration_mis::standard( MISConfig & mis_config ) {
    // Reductions
    mis_config.all_reductions                         = true;
    mis_config.num_reps                               = 5;
}

inline void configuration_mis::social( MISConfig & mis_config ) {
    standard(mis_config);
}

inline void configuration_mis::full_standard( MISConfig & mis_config ) {
    standard(mis_config);
}

inline void configuration_mis::full_social( MISConfig & mis_config ) {
    full_standard(mis_config);
}

#endif 
