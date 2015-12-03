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
    // Basic 
    mis_config.population_size                        = 50;
    mis_config.repetitions                            = 50;
    mis_config.time_limit                             = 1000.0;
    // KaHIP
    mis_config.kahip_mode                             = FAST;
    // Randomization
    mis_config.seed                                   = 0;
    mis_config.diversify                              = true;
    mis_config.imbalance                              = 0.03;
    mis_config.randomize_imbalance                    = true;
    // Selection
    mis_config.enable_tournament_selection            = true;
    mis_config.tournament_size                        = 2;
    // Mutation
    mis_config.flip_coin                              = 1;
    // Combination
    mis_config.use_hopcroft                           = false;
    mis_config.optimize_candidates                    = true;
    // Multiway
    mis_config.use_multiway_vc                        = false;
    mis_config.multiway_blocks                        = 64;
    // Thresholds
    mis_config.insert_threshold                       = 150;
    mis_config.pool_threshold                         = 250;
    mis_config.pool_renewal_factor                    = 10.0;
    // Separator pool
    mis_config.number_of_separators                   = 10;
    mis_config.number_of_partitions                   = 10;
    mis_config.number_of_k_separators                 = 10;
    mis_config.number_of_k_partitions                 = 10;
    // Output
    mis_config.print_repetition                       = true;
    mis_config.print_population                       = false;
    mis_config.console_log                            = false;
    mis_config.check_sorted                           = true;
    // ILS
    mis_config.ils_iterations                         = 15000;
    mis_config.force_k                                = 1;
    mis_config.force_cand                             = 4;
    // Reductions
    mis_config.all_reductions                         = true;
    // Convergence
    mis_config.reduction_threshold                    = 350;
    mis_config.remove_fraction                        = 0.10;
    mis_config.extract_best_nodes                     = true;
}

inline void configuration_mis::social( MISConfig & mis_config ) {
    standard(mis_config);
    mis_config.kahip_mode                             = FASTSOCIAL;
}

inline void configuration_mis::full_standard( MISConfig & mis_config ) {
    standard(mis_config);
    mis_config.population_size                        = 250;
    mis_config.time_limit                             = 36000.0;
    mis_config.number_of_separators                   = 30;
    mis_config.number_of_partitions                   = 30;
    mis_config.number_of_k_separators                 = 30;
    mis_config.number_of_k_partitions                 = 30;
    mis_config.flip_coin                              = 10;
    mis_config.pool_threshold                         = 200;
    mis_config.reduction_threshold                    = 1000;
}

inline void configuration_mis::full_social( MISConfig & mis_config ) {
    full_standard(mis_config);
    mis_config.kahip_mode                             = FASTSOCIAL;
}

#endif 
