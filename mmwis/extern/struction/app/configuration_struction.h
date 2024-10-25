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

#ifndef _STRUCTION_CONFIGURATION_MIS_H_
#define _STRUCTION_CONFIGURATION_MIS_H_

#include "definitions.h"
#include "mmwis_config.h"

class configuration_struction {
    public:
        /**
         * Default Constructor.
         */
        configuration_struction() {} ;

        /**
         * Default Destructor.
         */
        virtual ~configuration_struction() {};

        /**
         * Set the standard configuration.
         *
         * @param config Config to be initialized.
         */
        void standard( ::mmwis::MISConfig & config );
        void cyclicFast( ::mmwis::MISConfig & config );
        void cyclicStrong( ::mmwis::MISConfig & config );
};

inline void configuration_struction::standard( ::mmwis::MISConfig & mis_config ) {
    // Basic
    mis_config.time_limit                             = 1000.0;
    // Randomization
    mis_config.seed                                   = 0;
    // Output
    mis_config.console_log                            = false;
    mis_config.check_sorted                           = true;
    // ILS
    mis_config.ils_iterations                         = 15000;
    mis_config.ils_time_limit                         = 1000.0;
    mis_config.force_cand                             = 4;
	mis_config.sort_freenodes                         = true;
    // Reductions
	mis_config.perform_reductions                     = true;
    mis_config.struction_reduction_style                        = ::mmwis::MISConfig::StructionReduction_Style::NORMAL;
    // Weights
    mis_config.weight_source                          = ::mmwis::MISConfig::Weight_Source::FILE;

    mis_config.set_limit                              = 1024;
    mis_config.struction_degree                       = 256;
    mis_config.struction_type                         = ::mmwis::MISConfig::Struction_Type::EXTENDED;
    mis_config.key_type                               = ::mmwis::MISConfig::Key_Type::APPROXIMATE_INCREASE;
    mis_config.key_reinsert_factor                    = 2;
    mis_config.global_blow_up_factor                  = 9999;
    mis_config.phase_blow_up_factor                   = 2;
    mis_config.phase_blow_ups                         = 1;
    mis_config.max_unimproving_phases                 = 100;
    mis_config.backtrack_style                        = ::mmwis::MISConfig::Backtrack_Type::IMMEDIATE_EXCLUDE;
    mis_config.reduce_and_peel                        = false;
    mis_config.disable_clique_neighborhood            = false;
    mis_config.disable_generalized_fold               = false;
    mis_config.disable_critical_set                   = false;
    mis_config.disable_clique                         = false;
    mis_config.disable_blow_up                        = false;
    mis_config.plain_struction                        = false;
    mis_config.perform_hils                           = true;
}

inline void configuration_struction::cyclicFast( ::mmwis::MISConfig & mis_config ) {
    standard(mis_config);
    mis_config.struction_degree                       = 64;
    mis_config.max_unimproving_phases                 = 25;
    mis_config.set_limit                              = 512;
    mis_config.disable_generalized_fold               = true;
    mis_config.disable_clique_neighborhood            = true;

}

inline void configuration_struction::cyclicStrong( ::mmwis::MISConfig & mis_config ) {
    standard(mis_config);
    mis_config.struction_degree                       = 512;
    mis_config.max_unimproving_phases                 = 64;
    mis_config.set_limit                              = 2048;
    mis_config.disable_generalized_fold               = true;
    mis_config.disable_clique_neighborhood            = true;
}
#endif
