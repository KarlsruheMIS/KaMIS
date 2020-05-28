/**
 * configuration.h
 * Purpose: Contains preset configurations for the evolutionary algorithms.
 *
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
         *
         * @param config Config to be initialized.
         */
        void standard( MISConfig & config );
};

inline void configuration_mis::standard( MISConfig & mis_config ) {
    // Basic
    mis_config.time_limit                             = 1000.0;
    // Randomization
    mis_config.seed                                   = 0;
    // Output
    mis_config.console_log                            = false;
    mis_config.check_sorted                           = true;
    // ILS
    mis_config.ils_iterations                         = 15000;
    mis_config.force_cand                             = 4;
	mis_config.sort_freenodes                         = true;
    // Reductions
	mis_config.perform_reductions                     = true;
    mis_config.reduction_style                         = MISConfig::Reduction_Style::NORMAL;
    // Weights
    mis_config.weight_source                          = MISConfig::Weight_Source::FILE;
}

#endif
