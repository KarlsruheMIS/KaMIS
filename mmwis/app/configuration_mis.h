/**
 * configuration.h
 * Purpose: Contains preset configurations for the evolutionary algorithms.
 *
 *****************************************************************************/

#pragma once

#include "definitions.h"
#include "mmwis_config.h"
#include "kaHIP_interface.h"

#ifdef OMIS
const int FAST           = 0;
const int ECO            = 1;
const int STRONG         = 2;
const int FASTSOCIAL     = 3;
const int ECOSOCIAL      = 4;
const int STRONGSOCIAL   = 5;
#endif

namespace mmwis {

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
         * Set the configuration for the algorithm mmwis
         *
         * @param config Config to be initialized.
         */
        void mmwis( MISConfig & config );

        /**
         * Set the configuration for the algorithm mmwis+s
         *
         * @param config Config to be initialized.
         */
        void mmwiss( MISConfig & config );

    private:
        void standard( MISConfig & config );

};

inline void configuration_mis::standard( MISConfig & mis_config ) {
    // Basic
    mis_config.population_size                        = 250;
    mis_config.repetitions                            = 50;
    mis_config.time_limit                             = 1000.0;
    mis_config.evo_time_limit                         = 1000.0;
    // KaHIP
    mis_config.kahip_mode                             = FAST;
    // Randomization
    mis_config.seed                                   = 0;
    mis_config.diversify                              = true;
    mis_config.imbalance                              = 0.03;
    mis_config.randomize_imbalance                    = true;

    //Struction CyclicFast
    mis_config.set_limit                              = 1024;
    mis_config.struction_degree                       = 256;
    mis_config.struction_type                         = MISConfig::Struction_Type::EXTENDED;
    mis_config.key_type                               = MISConfig::Key_Type::APPROXIMATE_INCREASE;
    mis_config.key_reinsert_factor                    = 2;
    mis_config.global_blow_up_factor                  = 2;
    mis_config.phase_blow_up_factor                   = 2;
    mis_config.phase_blow_ups                         = 1;
    mis_config.max_unimproving_phases                 = 100;
    mis_config.backtrack_style                        = MISConfig::Backtrack_Type::IMMEDIATE_EXCLUDE;
    mis_config.reduce_and_peel                        = false;
    mis_config.disable_generalized_fold               = false;
    mis_config.disable_clique_neighborhood            = false;
    mis_config.disable_critical_set                   = false;
    mis_config.disable_clique                         = false;
    mis_config.disable_blow_up                        = false;
    mis_config.plain_struction                        = false;
    mis_config.perform_hils                           = true;
    mis_config.use_struction_initial_sol              = 0;


    // Weights
    mis_config.weight_source                          = MISConfig::Weight_Source::FILE;

    // Selection
    mis_config.enable_tournament_selection            = true;
    mis_config.tournament_size                        = 2;
    // Mutation
    mis_config.flip_coin                              = 1;
    // Combination
    mis_config.use_max_flow                           = true;
    mis_config.optimize_candidates                    = true;
    // Multiway
    mis_config.use_multiway_vc                        = true;
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
    mis_config.local_search_threshold                 = mis_config.ils_iterations;
    mis_config.force_k                                = 1;
    mis_config.force_cand                             = 4;
    mis_config.sort_freenodes                         = true;
    // Reductions
    mis_config.all_reductions                         = true;
	mis_config.perform_reductions                     = true; // run on kernel or on full graph
    mis_config.struction_reduction_style              = MISConfig::StructionReduction_Style::NORMAL;
    mis_config.reduction_style                        = MISConfig::Reduction_Style::initial;
    // Convergence
    mis_config.reduction_threshold                    = 350;
    mis_config.fraction                               = 0;
    mis_config.extract_best_nodes                     = true;
    // Initial solution
    mis_config.start_greedy_adaptive = false;

}

inline void configuration_mis::mmwis( MISConfig & mis_config ) {
    standard(mis_config);
    mis_config.evo_time_limit            = 3600;
    mis_config.time_limit                = 36000;
    mis_config.use_struction_initial_sol = 0;
    mis_config.time_solve_exact          = 100;
    mis_config.V_solve_exact             = 15000;
}

inline void configuration_mis::mmwiss( MISConfig & mis_config ) {
    standard(mis_config);
    mis_config.evo_time_limit            = 3600;
    mis_config.time_limit                = 36000;
    mis_config.use_struction_initial_sol = 250;
    mis_config.time_solve_exact          = 1000;
    mis_config.V_solve_exact             = 100;
}
}

